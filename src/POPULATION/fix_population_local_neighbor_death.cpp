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

#include "fix_population_local_neighbor_death.h"

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
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstring>

#include<unistd.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPopulationLocalNeighborDeath::FixPopulationLocalNeighborDeath(LAMMPS *lmp, int narg, char **arg) :
  FixPopulationBase(lmp, narg, arg),ncoeff(7)
{


  cutflag = 0;
  comm_reverse = 2;


  shift = utils::numeric(FLERR, arg[nspecified_args++], false, lmp);
  
  if (shift < 0)
    error->all(FLERR, "shift must be positive in fix population/sensing command");

  int iarg = nspecified_args;
  
  if (strcmp(arg[nspecified_args], "coeff") == 0)  {
    cutflag = 0;
  } else if (strcmp(arg[nspecified_args], "cutoff") == 0) {
    cutflag = 1;
    cutoff_user = utils::numeric(FLERR, arg[nspecified_args+1], false, lmp);
    iarg = nspecified_args+2;
  } else 
    error->all(FLERR, "Invalid keyword in fix population/sensing command (must be coeff or cutoff).");
  

  if ((narg - iarg  ) % ncoeff != 0)
    error->all(FLERR, "Invalid coefficients in fix population/sensing command (need multiple of 8).");


  allocate();
  const int rounds = (narg-iarg )/ncoeff;

  for (int i = 0; i < rounds; i++) {
      
    if (strcmp(arg[iarg + i*ncoeff], "coeff") == 0) 
      coeff(&arg[iarg + 1 + i*ncoeff]);
    else {
      error->all(FLERR,
		 "Invalid coefficients in fix population/sensing command.");
    }
  }

  

}


FixPopulationLocalNeighborDeath::~FixPopulationLocalNeighborDeath()
{


  memory->destroy(setflag);

  memory->destroy(b0);
  memory->destroy(d0);
  memory->destroy(sigma);
  memory->destroy(width);
  memory->destroy(cg_volume);    

}




/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void FixPopulationLocalNeighborDeath::allocate()
{

  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "fix:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;

  memory->create(d0, n, n, "fix:d0");
  memory->create(b0, n, n, "fix:b0");
  memory->create(sigma, n, n, "fix:sigma");
  memory->create(width, n, n, "fix:width");
  memory->create(cg_volume, n, n, "fix:cg_volume");



  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = 1; j <= atom->ntypes; j++) {
      setflag[i][j] = 0;
      b0[i][j] = 0.0;
      d0[i][j] = 0.0;
      sigma[i][j] = 0.0;
      width[i][j] = 1.0; // 1.0 to avoid divide by zero
      cg_volume[i][j] = 1.0;  // 1.0 to avoid divide by zero
    }
  }
  

}



void FixPopulationLocalNeighborDeath::init()
{

  FixPopulationBase::init();
  
  if (force->pair && cutflag)
    error->all(FLERR,"Fix population/sensing requires a pair style be defined "
               "exclusively or cutoff specified exclusively, not both");

  if (cutflag) {
    double skin = neighbor->skin;
    mycutneigh = cutoff_user + skin;

    double cutghost;            // as computed by Neighbor and Comm
    cutghost = comm->cutghostuser;

    if (mycutneigh > cutghost)
      error->all(FLERR,"Fix population/sensing cutoff exceeds ghost atom range - "
                 "use comm_modify cutoff command");

  }

  // need a neighbor list
  // if user specified, request a cutoff = cutoff_user + skin
  // skin is included b/c Neighbor uses this value similar
  //   to its cutneighmax = force cutoff + skin
  auto req = neighbor->add_request(this, NeighConst::REQ_DEFAULT);

  if (cutflag) {
    req->set_cutoff(mycutneigh);
  }

  full_neigh_list = false;
  
}



void FixPopulationLocalNeighborDeath::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

void FixPopulationLocalNeighborDeath::coeff(char **arg)
{

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double b0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double d0_one = utils::numeric(FLERR, arg[3], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[4], false, lmp);
  double width_one = utils::numeric(FLERR, arg[5], false, lmp);
  double cg_volume_one;

  if (domain->dimension == 2)
    cg_volume_one = M_PI*sigma_one*sigma_one;
    //cg_volume_one = 2*M_PI*sigma_one*sigma_one;
    //cg_volume_one = 1.0;
  else
    cg_volume_one = 4*M_PI/3.*sigma_one*sigma_one*sigma_one;
    //cg_volume_one = sqrt(2*M_PI*sigma_one*sigma_one)*sqrt(2*M_PI*sigma_one*sigma_one)*sqrt(2*M_PI*sigma_one*sigma_one);
    //cg_volume_one = 1.0;
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      b0[i][j] = b0_one;
      d0[i][j] = d0_one;
      sigma[i][j] = sigma_one;
      width[i][j] = width_one;
      cg_volume[i][j] = cg_volume_one;
      setflag[i][j] = 1;
      count++;
      
    }
  }


  if (count == 0) error->all(FLERR, "Incorrect args for fix population/sensing coefficients");
}


/* ---------------------------------------------------------------------- */


void FixPopulationLocalNeighborDeath::compute_division_and_death_rates()
{

  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double r_dist,rsq;
  int *ilist, *jlist, *numneigh, **firstneigh;


  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton = force->newton;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;


  // set divdeath_array components to zero (very heavily plagiarised from force_clear() method of verlet.cpp)
  // note that ghost atoms must be included if force->newton = 1
  size_t nbytes;
  if (neighbor->includegroup == 0) {
    nbytes = sizeof(double) * nlocal;
    if (true) nbytes += sizeof(double) * atom->nghost;

    if (nbytes) {
      memset(&divdeath_array[0][0],0,2*nbytes);
    }

  // neighbor includegroup flag is set
  // clear force only on initial nfirst particles
  // if either newton flag is set, also include ghosts

  } else {
    nbytes = sizeof(double) * atom->nfirst;

    if (nbytes) {
      memset(&divdeath_array[0][0],0,2*nbytes);
    }

    if (true) {
      nbytes = sizeof(double) * atom->nghost;

      if (nbytes) {
        memset(&divdeath_array[nlocal][0],0,2*nbytes);
      }
    }
  }

  
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];


    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      j &= NEIGHMASK;


      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];



      jtype = type[j];
      rsq = delx * delx + dely * dely + delz * delz;
      r_dist = sqrt(rsq);
      if (rsq < sigma[itype][jtype]) {//force->pair->cutsq[itype][jtype]) {


	//divdeath_array[i][1] += d0[itype][jtype]/(1+exp((r_dist-sigma[itype][jtype])/width[itype][jtype]))/cg_volume[itype][jtype];
	//divdeath_array[i][1] += d0[itype][jtype]/cg_volume[itype][jtype]*exp(-rsq/(2*sigma[itype][jtype]*sigma[itype][jtype]));
	divdeath_array[i][1] += d0[itype][jtype]/cg_volume[itype][jtype];
	
	// only do this bit if half neighbor list is being used (this is the default,
	//  but child classes of sensing might use a full neighbor list)
	if (!full_neigh_list) {
	  if (newton || j < nlocal) {
	    //divdeath_array[j][1] += d0[jtype][itype]/(1+exp((r_dist-sigma[jtype][itype])/width[jtype][itype]))/cg_volume[jtype][itype];
	    //divdeath_array[j][1] += d0[jtype][itype]/cg_volume[jtype][itype]*exp(-rsq/(2*sigma[jtype][itype]*sigma[jtype][itype]));
	    divdeath_array[j][1] += d0[jtype][itype]/cg_volume[jtype][itype];
	  }
	}
      }

    }

  }





  // communicate division and death contributions from ghost atoms to other processors

  if (!full_neigh_list) {
    comm->reverse_comm(this);
  }


  for (int i = 0; i < atom->nlocal; i++) {
    itype = type[i];
    divdeath_array[i][0] += b0[itype][itype];
    divdeath_array[i][1] += d0[itype][itype]/cg_volume[itype][itype];
  }

}


int FixPopulationLocalNeighborDeath::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m, last;

  
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = divdeath_array[i][0];
    buf[m++] = divdeath_array[i][1];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixPopulationLocalNeighborDeath::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    divdeath_array[j][0] += buf[m++];
    divdeath_array[j][1] += buf[m++];
  }
}


