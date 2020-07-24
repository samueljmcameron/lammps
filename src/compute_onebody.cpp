/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Sam Cameron
------------------------------------------------------------------------- */

#include "compute_onebody.h"
#include <mpi.h>
#include <stdio.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "domain.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeOneBody::ComputeOneBody(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  hist(NULL), histall(NULL), typecount(NULL), icount(NULL)
{
  if (narg != 5)
    error->all(FLERR,"Illegal compute three body command");

  array_flag = 1;
  extarray = 0;

  nbin = force->inumeric(FLERR,arg[3]);
  if (nbin < 1) error->all(FLERR,"Illegal compute three body command");

  radius = force->numeric(FLERR,arg[4]);

  delr = radius/nbin;
  delrinv = 1.0/delr;
  size_array_rows = nbin;

  
  int ntypes = atom->ntypes;

  // store densities for each type of atom, plus the density of all.
  if (ntypes > 1) {
    hist_rows = 1 + ntypes;
    size_array_cols = 2 + ntypes;
  } else {
    hist_rows = 1;
    size_array_cols = 2 ;
  }

  memory->create(hist,hist_rows,nbin,"rdf:hist");
  memory->create(histall,hist_rows,nbin,"rdf:histall");
  memory->create(array,size_array_rows,size_array_cols,"rdf:array");

  // type count will count the number of atoms of each type found
  // locally, indexed using atom->type
  typecount = new int[ntypes + 1];
  // icount counts the number of atoms found (first locally, then globally)
  // of each type, indexed from 0 to ntypes + 1. icount[0] is the number
  // of total atoms, and icount[i>0] is number of atoms of type i. If
  // ntypes ==1 , then this array will be one too big.
  icount = new int[ntypes+1];
  
  dynamic = 0;
  natoms_old = 0;
}

/* ---------------------------------------------------------------------- */

ComputeOneBody::~ComputeOneBody()
{
  memory->destroy(hist);
  memory->destroy(histall);
  memory->destroy(array);
  delete [] typecount;
  delete [] icount;
}

/* ---------------------------------------------------------------------- */

void ComputeOneBody::init()
{

  for (int i = 0; i < nbin; i++)
    array[i][0] = (i+0.5) * delr;


  natoms_old = atom->natoms;
  dynamic = group->dynamic[igroup];
  if (dynamic_user) dynamic = 1;
  init_norm();


  // NO NEEED FOR NEIGHBOR LIST!!
  // need an occasional half neighbor list
  // if user specified, request a cutoff = cutoff_user + skin
  // skin is included b/c Neighbor uses this value similar
  //   to its cutneighmax = force cutoff + skin
  // also, this NeighList may be used by this compute for multiple steps
  //   (until next reneighbor), so it needs to contain atoms further
  //   than cutoff_user apart, just like a normal neighbor list does


  //int irequest = neighbor->request(this,instance_me);
  //neighbor->requests[irequest]->pair = 0;
  //neighbor->requests[irequest]->compute = 1;
  //neighbor->requests[irequest]->occasional = 1;
  //if (cutflag) {
  //  neighbor->requests[irequest]->cut = 1;
  //  neighbor->requests[irequest]->cutoff = mycutneigh;
  //}

  
}

/* ---------------------------------------------------------------------- */

void ComputeOneBody::init_norm()
{
  int i;

  // count atoms of each type that are also in group

  const int nlocal = atom->nlocal;
  const int ntypes = atom->ntypes;
  const int * const mask = atom->mask;
  const int * const type = atom->type;


  // typecount[i] = # number of atoms of each type found locally
  for (i = 1; i <= ntypes; i++) typecount[i] = 0;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) typecount[type[i]]++;
  }
    
  // icount[m] = # of atoms of type m found locally

  for (i = 1; i <= ntypes; i++) icount[i] = typecount[i];

  // icount[0] = total # found locally
  icount[0] = 0;
  for (i = 1; i <= ntypes; i++) icount[0] += icount[i];

  int *scratch = new int[ntypes+1];
  MPI_Allreduce(icount,scratch,ntypes+1,MPI_INT,MPI_SUM,world);
  for (i = 0; i <= ntypes; i++) {
    icount[i] = scratch[i];
    //printf("total number of atoms of type %d = %d\n",i,icount[i]);
  }

  delete [] scratch;
}

/* ---------------------------------------------------------------------- */

void ComputeOneBody::compute_array()
// ======================================================================
// Only counting triplets where all three distances are less than cutoff
// away.
//
// ======================================================================
{
  
  if (natoms_old != atom->natoms) {
    dynamic = 1;
    natoms_old = atom->natoms;
  }

  // if the number of atoms has changed or we have a dynamic group
  // or dynamic updates are requested (e.g. when changing atom types)
  // we need to recompute some normalization parameters

  if (dynamic) init_norm();

  invoked_array = update->ntimestep;

  // compute 
  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;

  // zero the histogram counts

  for (int i = 0; i < hist_rows; i++) {
    for (int j = 0; j < nbin; j++) {
      hist[i][j] = 0;
    }
  }
  
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;

  double xp,yp,zp, r;
  int ibin;
  int itype;

  // loop over full neighbor list of my atoms
  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    itype = type[i];
    
    xp = x[i][0];
    yp = x[i][1];
    zp = x[i][2];

    r = sqrt(xp*xp+yp*yp+zp*zp);

    ibin = static_cast<int> (r*delrinv);
	
    if (ibin >= nbin) continue;

    hist[0][ibin] += 1.0;
    if (ntypes > 1)
      hist[itype][ibin] += 1.0;

  }

  // sum histograms across procs

  MPI_Allreduce(hist[0],histall[0],hist_rows*nbin,MPI_DOUBLE,MPI_SUM,world);

  // convert counts to g(r) and coord(r) and copy into output array
  // vfrac = fraction of volume in shell m

  
  double prefac, normfac;

  if (domain->dimension==2) {
    prefac = MY_PI*radius*radius;
    prefac = MY_PI/prefac;

    //normfac = (icount[0] > 0) ? static_cast<double>(jcount[0])
    // - static_cast<double>(duplicates[0])/icount[0] : 0.0;
    normfac = static_cast<double>(icount[0]);
    set_array(prefac, normfac,0);

    if (ntypes > 1) {
      for (int i = 1; i < ntypes; i++) {
	normfac = static_cast<double>(icount[i]);
	set_array(prefac, normfac,i);
      }
    }
  }//  else (domain->dimension == 3) {

}

void ComputeOneBody::set_array(double prefac, double normfac,int a_type)
{
  double ulower,uupper;
  double vfrac;

  int ibin;
  double rhor;
  int col = a_type + 1;
      
  for (ibin = 0; ibin < nbin; ibin++) {
      ulower = ibin*delr;
      uupper = (ibin+1)*delr;
      vfrac = (prefac * (uupper*uupper - ulower*ulower)/2.0);
      if (vfrac * normfac != 0.0) {
	rhor = histall[a_type][ibin]/(vfrac *normfac);
      } else {
	rhor = 0.0;
      }
      array[ibin][col] = rhor;
  }

  return;
}
