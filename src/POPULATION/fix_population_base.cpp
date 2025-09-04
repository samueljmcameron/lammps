/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Sam Cameron
------------------------------------------------------------------------- */

#include "fix_population_base.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "random_mars.h"
#include "pair.h"
#include "update.h"
#include "modify.h"

#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstring>

#include<unistd.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPopulationBase::FixPopulationBase(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), alive_indices(nullptr),
  dead_indices(nullptr),nalive_per_proc(nullptr),
  ndead_per_proc(nullptr),divdeath_array(nullptr),shift(0.0)
{


  vector_flag = 1;   // fix calculates a vector having global nalive and ndead in it
  size_vector = 6;
  global_freq = 1;
  extvector = 0;
  
  peratom_flag = 1; // fix calculates a peratom array having division and death rates
  size_peratom_cols = 2;
  peratom_freq = 1;

  cleanevery = 0;
  
  nspecified_args = 3;
  alivetype = utils::inumeric(FLERR, arg[nspecified_args++], false, lmp);
  if (alivetype <= 0 || alivetype > atom->ntypes)
    error->all(FLERR, "Invalid atom type in fix population command");
  
  deadtype = utils::inumeric(FLERR, arg[nspecified_args++], false, lmp);
  if (deadtype <= 0 || deadtype > atom->ntypes)
    error->all(FLERR, "Invalid atom type in fix population command");

  seed = utils::inumeric(FLERR, arg[nspecified_args++], false, lmp);

  if (seed <= 0)
    error->all(FLERR, "Seed must be positive in fix population command");
  
  if (strcmp(arg[nspecified_args], "delete") == 0) {
    nspecified_args += 1;
    if (strcmp(arg[nspecified_args], "efficient") == 0) {
      delete_flag = EFFICIENT;
    } else if (strcmp(arg[nspecified_args], "every") == 0) {
      delete_flag = EVERY;
      nspecified_args += 1;
      cleanevery = utils::inumeric(FLERR, arg[nspecified_args], false, lmp);
    } else if (strcmp(arg[nspecified_args], "never") == 0) {
      delete_flag = NEVER;
    } else
      error->all(FLERR, "delete argument must be either efficient, every, or never in "
		 "fix population command");

    nspecified_args += 1;
  } else delete_flag = EFFICIENT;

  if (strcmp(arg[nspecified_args], "recycle") == 0) {
    nspecified_args += 1;
    if (strcmp(arg[nspecified_args], "yes") == 0) {
      recycle_flag = 1;
    } else if (strcmp(arg[nspecified_args], "no") == 0) {
      recycle_flag = 0;
    } else
      error->all(FLERR, "recycle argument must be either yes or no in "
		 "fix population command");
    nspecified_args += 1;
  } else recycle_flag = 0;

  if (delete_flag == NEVER && !recycle_flag) {
    error->warning(FLERR, "Never deleting and never recycling atoms in fix population command "
		   "will lead to unbounded accumulation of dead atoms.");
  }
  
  rng = new RanMars(lmp, seed + comm->me*100);


  // communicate atom types to neighboring processors (since types can be switched by this fix)
  comm_forward = 1;

  // must reneighbor if new atoms are created or dead atoms are deleted (see next_reneighbor flag in relevant functions below)
  force_reneighbor = 1;
  
  
  next_reneighbor = 0;

  // add arrays which store nalive for local number of processors (superfluous to store this but good debug check)
  memory->create(nalive_per_proc,comm->nprocs, "fix_population:nalive_per_proc");
  memory->create(ndead_per_proc,comm->nprocs, "fix_population:ndead_per_proc");


  // allocate other arrays
  FixPopulationBase::grow_arrays(atom->nmax);
  // resize arrays in this fix by adding a callback to atom (so atom will handle resizing)
  atom->add_callback(Atom::GROW);
  
  for (int i = 0; i < atom->nlocal; i++) {
    divdeath_array[i][0] = 0.0;
    divdeath_array[i][1] = 0.0;
  }
 

  
}

/* ---------------------------------------------------------------------- */

FixPopulationBase::~FixPopulationBase()
{

  memory->destroy(alive_indices);
  memory->destroy(dead_indices);
  memory->destroy(nalive_per_proc);
  memory->destroy(ndead_per_proc);

  memory->destroy(divdeath_array);
  if (modify->get_fix_by_id(id)) atom->delete_callback(id, Atom::GROW);
  
  
  delete rng;

}

void FixPopulationBase::reset_dt()
{
  
  double dt = update->dt;

}

 
/* ---------------------------------------------------------------------- */

int FixPopulationBase::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

void FixPopulationBase::init()
{
  count_vitals(); 

}


/* ----------------------------------------------------------------------
   allocate atom-based arrays, and point array_atom to divdeath_array
------------------------------------------------------------------------- */

void FixPopulationBase::grow_arrays(int nmax)
{
  memory->grow(divdeath_array, nmax,2, "fix_population:divdeath_array");
  memory->grow(alive_indices ,nmax, "fix_population:alive_indices");
  memory->grow(dead_indices, nmax, "fix_population:dead_indices");
  array_atom = divdeath_array;

}
void FixPopulationBase::copy_arrays(int i, int j, int /*delflag*/)
{
  divdeath_array[j][0] = divdeath_array[i][0];
  divdeath_array[j][1] = divdeath_array[i][1];
}



/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixPopulationBase::pack_exchange(int i, double *buf)
{
  buf[0] = divdeath_array[i][0];
  buf[1] = divdeath_array[i][1];
  return 2;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixPopulationBase::unpack_exchange(int nlocal, double *buf)
{

  divdeath_array[nlocal][0] = buf[0];
  divdeath_array[nlocal][1] = buf[1];
  return 2;
}


/* ----------------------------------------------------------------------
   actually perform the birth and death of atoms.

   this will force a reneighboring if any of the following are true:
   1) new atoms are created
   2) deadtype atoms are deleted from the system

   it will NOT force a reneighboring if a deadtype atom is swapped
   to an alivetype atom to become the daughter of another atom.
   This means that if there is at least one (local) deadtype atom
   remaining and a (local) alivetype atom is ready to divide, the
   deadtype atom's position will be moved according to the
   divide(i,j) function. This will likely mean the deadtype atom
   moves larger than half the skin distance, and will therefore
   be considered a "dangerous build" by neighbor.

   additionally, there might be issues with using neigh_modify exclude
   deadtype atoms?

------------------------------------------------------------------------- */

void FixPopulationBase::post_integrate()
{

  int *mask = atom->mask;
  double dt = update->dt;

  int i;
  double ran;
  
  count_vitals();  // create alive_indices and compute nalive, ndead

  // array to store parents which have no nearby free atoms
  std::vector<int> dividing_from_scratch_atoms; 



  int number_of_recycled_atoms = 0;
  int number_of_dieing_atoms = 0;
  // update divdeath_array


  compute_division_and_death_rates();


  // iterate over alive atoms only
  for (int ia = 0; ia < localalive; ia++) {
    i = alive_indices[ia];
      
    ran = rng->uniform();

    
    if (ran  <= divdeath_array[i][0]*dt) { // if true then division event will occur

      if (recycle_flag) {
	// recycle local dead atom into an alive atom if possible
	bool recycled = recycle_from_dead(i);

	if (! recycled) { // alive atom must be created from scratch 
	  dividing_from_scratch_atoms.push_back(i);
	  
	} else number_of_recycled_atoms += 1;
      } else dividing_from_scratch_atoms.push_back(i);
      
      
    } else if (ran <= (divdeath_array[i][0] + divdeath_array[i][1])*dt) { // if true then death event occurs
      number_of_dieing_atoms += 1;      
      atom->type[i] = deadtype;
      
    }
    
  }


  // sum over processors for events that don't necessarily create/delete atoms

  MPI_Allreduce(&number_of_recycled_atoms,&total_atoms_recycled,1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&number_of_dieing_atoms,&total_atoms_killed,1, MPI_INT, MPI_SUM, world);


  // now move on to events that create/delete atoms
  
  create_new_atoms(dividing_from_scratch_atoms);
  nalive += total_atoms_from_scratch+total_atoms_recycled-total_atoms_killed;

  total_atoms_deleted = 0;
  if ((delete_flag == EFFICIENT && total_atoms_from_scratch > 0) ||
      (delete_flag == EVERY && update->ntimestep % cleanevery == 0)) {
    delete_dead_atoms(); // computes total_atoms_deleted;
    ndead = 0;
  } else // if atoms aren't deleted 
    ndead += total_atoms_killed-total_atoms_recycled;

  
  if (total_atoms_from_scratch > 0 || total_atoms_deleted > 0 || total_atoms_recycled > 0) {
    // an atom has either been created, deleted, or moved an arbitrary amount
    //  (within a processor), respectively, meaning a reneighboring must be done.

    next_reneighbor = update->ntimestep;

    if (atom->map_style != Atom::MAP_NONE && (total_atoms_from_scratch > 0 || total_atoms_deleted > 0)) {
      atom->nghost = 0;
      atom->map_init();
      atom->map_set();
    }
  } else  // if either atoms are labelled dead but not deleted, or nothing at all has happened.
    comm->forward_comm(this);

}


/* ---------------------------------------------------------------------- */
/* Store atoms which are alive in alive_indices and count total
   number of alive atoms and dead atoms across all processors. */
/* ---------------------------------------------------------------------- */

void FixPopulationBase::count_vitals()
{
  
  localalive = 0;
  localdead = 0;

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      if (atom->type[i] == alivetype) {
	alive_indices[localalive++] = i;
      } else if (atom->type[i] == deadtype) {
	dead_indices[localdead++] = i;
      }
    }
  }


  MPI_Allgather(&localalive,1,MPI_INT,&(nalive_per_proc[0]),1,MPI_INT,world);
  MPI_Allgather(&localdead,1,MPI_INT,&(ndead_per_proc[0]),1,MPI_INT,world);
  
  
  MPI_Allreduce(&localalive, &nalive, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&localdead, &ndead, 1, MPI_INT, MPI_SUM, world);

  int sum = 0;

  for (int i = 0; i < comm->nprocs; i++)
    sum += nalive_per_proc[i];

  if (sum != nalive)
    printf("ERROR sum = %d but nalive = %d !!!!\n\n\n\n\n\n",sum,nalive);

  sum = 0;

  for (int i = 0; i < comm->nprocs; i++)
    sum += ndead_per_proc[i];

  if (sum != ndead)
    printf("ERROR sum = %d but ndead = %d !!!!\n\n\n\n\n\n",sum,ndead);
  

  return;
}


/* ---------------------------------------------------------------------- */
/* For atom i which is alivetype try and divide by converting a
   (local) deadtype atom into one of its daughters.  */
/* ---------------------------------------------------------------------- */

bool FixPopulationBase::recycle_from_dead(int i) {
  
  bool recycled = false;

  for (int j = 0;  j < atom->nlocal; j++) {
    
    
    if (atom->type[j] == deadtype) {
      divide(i,j);
      recycled = true;
      break;
    }

  }

  return recycled;
}



/* ----------------------------------------------------------------------
   Create a set of new daughter atoms from the list of parent atoms.
   Also, tally up the total_atoms_from_scratch.
---------------------------------------------------------------------- */
void FixPopulationBase::create_new_atoms(const std::vector<int> &new_atoms)
{
  // store number of local atoms before depleted parents generate new gametes
  bigint nlocal = atom->nlocal;
  double x[3];


  // create new atoms (overwrites ghost atoms so need to rebuild neighbor list next step).
  // total atoms created will = (new_atoms.size() * atom->numgametes)
  int i,n;
  for (auto i: new_atoms) {
    x[0] = atom->x[i][0];
    x[1] = atom->x[i][1];
    x[2] = atom->x[i][2];
    atom->avec->create_atom(atom->type[i],x);
    n = atom->nlocal - 1;
    atom->avec->copy(i,n,0);
    atom->type[n] = deadtype;
    atom->tag[n] = 0;
    atom->f[n][0] = 0.0;
    atom->f[n][1] = 0.0;
    atom->f[n][2] = 0.0;
    modify->create_attribute(n);
    divide(i,n);
  }


  
  bigint newlocal = atom->nlocal;
  bigint natoms_previous = atom->natoms;
  
  MPI_Allreduce(&newlocal, &atom->natoms, 1, MPI_INT, MPI_SUM, world);
  
  total_atoms_from_scratch = atom->natoms - natoms_previous;

  if (atom->natoms < 0 || atom->natoms >= MAXBIGINT)
    error->all(FLERR, "Too many total atoms");

  if (total_atoms_from_scratch > 0) {
    
    // add IDs for newly created atoms
    // check that atom IDs are valid
    
    if (atom->tag_enable) atom->tag_extend();
    atom->tag_check();

  }

  return;

}

/* ----------------------------------------------------------------------
   Delete dead atoms from the simulation. Also, tally up the
   total_atoms_deleted.
---------------------------------------------------------------------- */
void FixPopulationBase::delete_dead_atoms()
{

  bigint natoms_previous = atom->natoms;
  int nlocal = atom->nlocal;


  std::vector<int> dlist(nlocal);   // vector of atoms to be deleted (1 if should delete, 0 if shouldn't delete)
  
  for (int i = 0; i < nlocal; i++) {
    if (atom->type[i] == deadtype)
      dlist[i] = 1;
    else
      dlist[i] = 0;
  }

  int i = 0;

  while (i < nlocal) {  // proceed to delete dead atoms from simulation
    if (dlist[i]) {
      atom->avec->copy(nlocal - 1, i, 1);
      dlist[i] = dlist[nlocal - 1];
      nlocal--;
    } else
      i++;
  }
  atom->nlocal = nlocal;
  bigint nblocal = atom->nlocal;

  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  total_atoms_deleted = natoms_previous - atom->natoms;

  return;
  
}

/* ---------------------------------------------------------------------- 
   divide by using atom[i] (which should be alivetype) and atom[j]
   (which should be deadtype)
---------------------------------------------------------------------- */

void FixPopulationBase::divide(int i, int j)
{

  double **x = atom->x;
  double **v = atom->v;

  
  double polar,azim;
  double cx,cy,cz,dx,dy,dz;

  cx = x[i][0];
  cy = x[i][1];
  cz = x[i][2];
  
  if (domain->dimension == 3) {
    polar = M_PI*rng->uniform();
    azim = 2*M_PI*rng->uniform();
    dx = shift*cos(azim)*sin(polar)/2.;
    dy = shift*sin(azim)*sin(polar)/2.;
    dz = shift*cos(polar)/2;
    
  } else {
    azim = 2*M_PI*rng->uniform();
    dx = shift*cos(azim)/2.;
    dy = shift*sin(azim)/2.;
    dz = 0.0;
  }

  atom->type[j] = alivetype;
	    
  x[j][0] = cx+dx;
  x[j][1] = cy+dy;
  x[j][2] = cz+dz;
  
  x[i][0] = cx-dx;
  x[i][1] = cy-dy;
  x[i][2] = cz-dz;

  v[j][0] = v[i][0]/2.;
  v[j][1] = v[i][1]/2.;
  v[j][2] = v[i][2]/2.;

  v[i][0] /= 2;
  v[i][1] /= 2;
  v[i][2] /= 2;

  divdeath_array[j][0] = divdeath_array[i][0];
  divdeath_array[j][1] = divdeath_array[i][1];
  return;
}


/* ----------------------------------------------------------------------
   birth and death rates for each atom
------------------------------------------------------------------------- */

double FixPopulationBase::compute_vector(int i)
{

  //if (count_vitals_flag)
  //  count_vitals(); 

  if (i == 0) return nalive;
  else if (i == 1) return ndead;
  else if (i == 2) return total_atoms_from_scratch;
  else if (i == 3) return total_atoms_deleted;
  else if (i == 4) return total_atoms_recycled;
  else if (i == 5) return total_atoms_killed;


  return -1;
}



/* ----------------------------------------------------------------------
   Pack atom->type (since this
   fix may change alivetype atoms to deadtype atoms and vice versa)
---------------------------------------------------------------------- */

int FixPopulationBase::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  tagint *tag = atom->tag;
  int *type = atom->type;

  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = ubuf(type[j]).d;
  }

  return m;
}



/* ----------------------------------------------------------------------
   Unpack atom->type
---------------------------------------------------------------------- */
void FixPopulationBase::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  tagint *tag = atom->tag;
  int *type = atom->type;

  m = 0;
  last = first + n;

  for (i = first; i < last; i++) {
    type[i] = (int) ubuf(buf[m++]).i; 
  }

}
