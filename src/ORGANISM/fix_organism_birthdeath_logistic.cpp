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

#include "fix_organism_birthdeath_logistic.h"

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

FixOrganismBirthDeathLogistic::FixOrganismBirthDeathLogistic(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), local_alive_list(nullptr),
  local_dead_list(nullptr),proc_dist(0,comm->nprocs-1),atom_dist(0,1)
{

  // required args
  vector_flag = 1;
  extvector = 0;
  size_vector = 2;
  global_freq = 1;
  
  alivetype = utils::inumeric(FLERR, arg[3], false, lmp);
  if (alivetype <= 0 || alivetype > atom->ntypes)
    error->all(FLERR, "Invalid atom type in fix organism/birthdeath/logistic command");
  
  deadtype = utils::inumeric(FLERR, arg[4], false, lmp);
  if (deadtype <= 0 || deadtype > atom->ntypes)
    error->all(FLERR, "Invalid atom type in fix organism/birthdeath/logistic command");

  seed = utils::inumeric(FLERR, arg[5], false, lmp);

  if (seed <= 0)
    error->all(FLERR, "Seed must be positive in fix organism/birthdeath/logistic command");
  
  birthrate = utils::numeric(FLERR, arg[6], false, lmp);
  deathrate = utils::numeric(FLERR, arg[7], false, lmp);

  if (birthrate <= 0)
    error->all(FLERR, "Birth rate must be positive in fix organism/birthdeath/logistic command");

  if (deathrate <= 0)
    error->all(FLERR, "Death rate must be positive in fix organism/birthdeath/logistic command");
  
  shift = utils::numeric(FLERR, arg[8], false, lmp);
  
  if (shift < 0)
    error->all(FLERR, "Shift must be positive in fix organism/birthdeath/logistic command");

  cleanevery = utils::inumeric(FLERR, arg[9], false, lmp);
  
  rng = new RanMars(lmp, seed + comm->me*100);

  comm_forward = 1;  
  force_reneighbor = 1;
  
  atom_swap_nmax = 0;
  next_reneighbor = 0;

  gen.seed(seed + comm->me*100);

  global_alive_list =
    (int *) memory->smalloc(comm->nprocs * sizeof(int), "MCSWAP:global_alive_list");
  global_dead_list =
    (int *) memory->smalloc(comm->nprocs * sizeof(int), "MCSWAP:global_dead_list");


  if (birthrate/deathrate < 1.0)
    error->all(FLERR, "Steady-state number of particles < 1 in fix organism/birthdeath/logistic command");  
  

  
}

/* ---------------------------------------------------------------------- */

FixOrganismBirthDeathLogistic::~FixOrganismBirthDeathLogistic()
{

  memory->destroy(local_alive_list);
  memory->destroy(local_dead_list);
  memory->destroy(global_alive_list);
  memory->destroy(global_dead_list);
  
  delete rng;

}

void FixOrganismBirthDeathLogistic::reset_dt()
{
  
  double Nav = birthrate/deathrate;

  double dt = update->dt;
  
  if (Nav*birthrate*dt > 1.0)
    error->all(FLERR,"Birth rate in fix organism/birthdeath/logistic command is too large, multiple events will occur for given dt.");    
  else if (Nav*dt*(deathrate*(Nav-1) + birthrate) > 1.0) 
    error->all(FLERR,"Death rate in fix organism/birthdeath/logistic command is too large, multiple events will occur for given dt.");

}

 
/* ---------------------------------------------------------------------- */

int FixOrganismBirthDeathLogistic::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

void FixOrganismBirthDeathLogistic::init()
{

}

/* ----------------------------------------------------------------------
   birth death
------------------------------------------------------------------------- */

void FixOrganismBirthDeathLogistic::post_integrate()
{

  
  int *mask = atom->mask;
  double dt = update->dt;



  int i;
  double ran;
  
  count_vitals();  // create local_alive_list and compute nalive, ndead

  std::vector<int> new_atoms; // array to store parents which have no nearby free atoms
  //                                     but need to procreate

  int proc = -1;
  int birth = 0;

  // decide whether birth or death happens only on one processor
  if (comm->me == 0) {
    ran = rng->uniform();

    // choose processor and which atom  
    if (ran  <= nalive*birthrate*dt) {
      birth = 1;
      proc = proc_dist(gen);
      while (global_alive_list[proc] == 0) 
	proc = proc_dist(gen);

    } else if (ran <= nalive*dt*(deathrate*(nalive-1) + birthrate)) {
      proc = proc_dist(gen);
      while (global_alive_list[proc] == 0) 
	proc = proc_dist(gen);

      
    }

    
  }

  MPI_Bcast(&proc,1,MPI_INT,0,world);
  MPI_Bcast(&birth,1,MPI_INT,0,world);

  if (proc == comm->me) {
    atom_dist.param(decltype(atom_dist)::param_type(0,localalive-1));

    i = local_alive_list[atom_dist(gen)];

    if (birth) {
      bool procreated = birth_from_dead(i);

      if (! procreated) {
	new_atoms.push_back(i);
      }
      
    } else {

      atom->type[i] = deadtype;
    }
  
  }

  /*
  for (int ia = 0; ia < localalive; ia++) {
    i = local_alive_list[ia];
    
    ran = rng->uniform();

    if (ran  <= birthrate*dt) {

      bool procreated = birth_from_dead(i);

      if (! procreated) {
	new_atoms.push_back(i);
      }
      
    } else if (ran <= deathrate*(nalive-1)*dt + birthrate*dt) {

      atom->type[i] = deadtype;
      
    }

  }
  */

  create_new_atoms(new_atoms);
  
  comm->forward_comm(this);
  delete_dead_atoms();

  

}
      

/* ---------------------------------------------------------------------- */

int FixOrganismBirthDeathLogistic::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
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

/* ---------------------------------------------------------------------- */

void FixOrganismBirthDeathLogistic::unpack_forward_comm(int n, int first, double *buf)
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


void FixOrganismBirthDeathLogistic::delete_dead_atoms()
{

  if (update->ntimestep % cleanevery != 0)

    return;

  next_reneighbor = update->ntimestep;
  
  bigint natoms_previous = atom->natoms;
  int nlocal = atom->nlocal;

  std::vector<int> dlist(nlocal);
  
  for (int i = 0; i < nlocal; i++) {
    if (atom->type[i] == deadtype)
      dlist[i] = 1;
    else
      dlist[i] = 0;
  }

  int i = 0;

  while (i < nlocal) {
    if (dlist[i]) {
      atom->avec->copy(nlocal - 1, i, 1);
      dlist[i] = dlist[nlocal - 1];
      nlocal--;
    } else
      i++;
  }
  
  atom->nlocal = nlocal;

  
  // reset atom->natoms and also topology counts
  
  bigint nblocal = atom->nlocal;
  MPI_Allreduce(&nblocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  
  // reset atom->map if it exists
  // set nghost to 0 so old ghosts of deleted atoms won't be mapped
  
  if (atom->map_style != Atom::MAP_NONE) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }
  
  
  bigint ndelete = natoms_previous - atom->natoms;

}


void FixOrganismBirthDeathLogistic::create_new_atoms(const std::vector<int> &new_atoms)
{
  // store number of local atoms before depleted parents generate new gametes
  bigint nlocal = atom->nlocal;


  // create new atoms (overwrites ghost atoms so need to rebuild neighbor list next step).
  // total atoms created will = (new_atoms.size() * atom->numgametes)
  int i,n;
  for (auto i: new_atoms) {
    atom->avec->create_atom(atom->type[i], atom->x[i]);
    n = atom->nlocal - 1;
    atom->avec->copy(i,n,0);
    atom->type[n] = deadtype;
    atom->tag[n] = 0;
    atom->f[n][0] = 0.0;
    atom->f[n][1] = 0.0;
    atom->f[n][2] = 0.0;
    modify->create_attribute(n);
    procreate(i,n);
  }


  
  int reneigh =  new_atoms.size();
  int globalreneigh;
  MPI_Allreduce(&reneigh, &globalreneigh, 1, MPI_INT, MPI_SUM, world);
  if (globalreneigh > 0) {
    next_reneighbor = update->ntimestep;
  

    bigint newlocal = atom->nlocal;
    MPI_Allreduce(&newlocal, &atom->natoms, 1, MPI_LMP_BIGINT, MPI_SUM, world);
    if (atom->natoms < 0 || atom->natoms >= MAXBIGINT)
      error->all(FLERR, "Too many total atoms");
    
    // add IDs for newly created atoms
    // check that atom IDs are valid
    
    if (atom->tag_enable) atom->tag_extend();
    atom->tag_check();
    
    // if global map exists, reset it
    // invoke map_init() b/c atom count has grown
    
    if (atom->map_style != Atom::MAP_NONE) {
      atom->map_init();
      atom->map_set();
    }
  }

  return;
}


/* ---------------------------------------------------------------------- */
/* For atom i which is ALIVE but has no gametes, try and absorb a (local)
   FREE neighbor, turn it into a gamete, and then procreate from it.  */
/* ---------------------------------------------------------------------- */

bool FixOrganismBirthDeathLogistic::birth_from_dead(int i) {
  
  bool procreated = false;

  for (int j = 0;  j < atom->nlocal; j++) {
    
    
    if (atom->type[j] == deadtype) {
      procreate(i,j);
      procreated = true;
      break;
    }

  }

  return procreated;
}

/* ---------------------------------------------------------------------- */
/* For atom i which is ALIVE and atom j which is GAMETE (stored in atom
   i's gametes at gindex), convert atom j into ALIVE and split remaining
   gametes between atoms i and j. ASSUMES gindex IS LOCATION OF FIRST
   GAMETE (i.e. atom->gametes[i][g] == -1 for g < gindex) AND
   atom->gametes[i][gindex] == atom->tag[j].  */
/* ---------------------------------------------------------------------- */

void FixOrganismBirthDeathLogistic::procreate(int i, int j)
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
  
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  return;
}

/* ---------------------------------------------------------------------- */
/* Store atoms which are alive in local_alive_list and count total
   number of alive atoms and dead atoms across all processors. */
/* ---------------------------------------------------------------------- */

void FixOrganismBirthDeathLogistic::count_vitals()
{
  
  localalive = 0;
  localdead = 0;

  if (atom->nmax > atom_swap_nmax) {
    memory->sfree(local_alive_list);
    memory->sfree(local_dead_list);
    atom_swap_nmax = atom->nmax;
    local_alive_list =
        (int *) memory->smalloc(atom_swap_nmax * sizeof(int), "MCSWAP:local_alive_list");
    local_dead_list =
        (int *) memory->smalloc(atom_swap_nmax * sizeof(int), "MCSWAP:local_dead_list");
  }

  for (int i = 0; i < atom->nlocal; i++) {
    if (atom->mask[i] & groupbit) {
      if (atom->type[i] == alivetype) {
	local_alive_list[localalive++] = i;
      } else if (atom->type[i] == deadtype) {
	local_dead_list[localdead++] = i;
      }
    }
  }


  MPI_Allgather(&localalive,1,MPI_INT,&(global_alive_list[0]),1,MPI_INT,world);
  MPI_Allgather(&localdead,1,MPI_INT,&(global_dead_list[0]),1,MPI_INT,world);
  
  
  MPI_Allreduce(&localalive, &nalive, 1, MPI_INT, MPI_SUM, world);
  MPI_Allreduce(&localdead, &ndead, 1, MPI_INT, MPI_SUM, world);

  int sum = 0;

  for (int i = 0; i < comm->nprocs; i++)
    sum += global_alive_list[i];

  if (sum != nalive)
    printf("ERROR sum = %d but nalive = %d !!!!\n\n\n\n\n\n",sum,nalive);

  sum = 0;

  for (int i = 0; i < comm->nprocs; i++)
    sum += global_dead_list[i];

  if (sum != ndead)
    printf("ERROR sum = %d but ndead = %d !!!!\n\n\n\n\n\n",sum,ndead);
  

  return;
}


/* ----------------------------------------------------------------------
   number of alive particles
------------------------------------------------------------------------- */

double FixOrganismBirthDeathLogistic::compute_vector(int i)
{

  count_vitals();  // create local_alive_list and compute nalive, ndead

  if (i == 0) return nalive;
  if (i == 1) return ndead;

  return -1;
}
