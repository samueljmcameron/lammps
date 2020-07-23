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
  if (ntypes != 1) {
    error->all(FLERR,"Cannot compute onebody for system with multiple "
	       "atom types.");
  }

  // store densities for each type of atom, plus the density of all.
  size_array_cols = 2;

  memory->create(hist,1,nbin,"rdf:hist");
  memory->create(histall,1,nbin,"rdf:histall");
  memory->create(array,size_array_rows,size_array_cols,"rdf:array");  
  
  dynamic = 0;
  natoms_old = 0;
}

/* ---------------------------------------------------------------------- */

ComputeOneBody::~ComputeOneBody()
{
  memory->destroy(hist);
  memory->destroy(histall);
  memory->destroy(array);
  memory->destroy(typecount);
  memory->destroy(icount);
}

/* ---------------------------------------------------------------------- */

void ComputeOneBody::init()
{

  natoms_old = atom->natoms;
  dynamic = group->dynamic[igroup];
  if (dynamic_user) dynamic = 1;
  init_norm();

  // need an occasional half neighbor list
  // if user specified, request a cutoff = cutoff_user + skin
  // skin is included b/c Neighbor uses this value similar
  //   to its cutneighmax = force cutoff + skin
  // also, this NeighList may be used by this compute for multiple steps
  //   (until next reneighbor), so it needs to contain atoms further
  //   than cutoff_user apart, just like a normal neighbor list does

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->occasional = 1;
  if (cutflag) {
    neighbor->requests[irequest]->cut = 1;
    neighbor->requests[irequest]->cutoff = mycutneigh;
  }
  // need full neighbor list
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  
}

/* ---------------------------------------------------------------------- */

void ComputeOneBody::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeOneBody::init_norm()
{
  int i,j,m;

  // count atoms of each type that are also in group

  const int nlocal = atom->nlocal;
  const int ntypes = atom->ntypes;
  const int * const mask = atom->mask;
  const int * const type = atom->type;

  //here
  typecount[1] = 0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) typecount[type[i]]++;

  // icount = # of I atoms participating in density

  //here
  icount[0] = 0;
  icount[0] += typecount[1];

  //here
  int *scratch = new int[1];
  MPI_Allreduce(icount,scratch,1,MPI_INT,MPI_SUM,world);
  icount[0] = scratch[0];

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
  int i,j,m,ii,jj,inum,jnum,itype,jtype,ipair,jpair;
  int i_bin;
  int k,kk,knum,ktype,kpair;
  double xtmp,ytmp,ztmp;
  double xij,yij,zij,xik,yik,zik,rij,rik,rjk,theta;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double dumcostheta;
  double factor_lj,factor_coul;
  double r;

  double denom;
  
  if (natoms_old != atom->natoms) {
    dynamic = 1;
    natoms_old = atom->natoms;
  }

  // if the number of atoms has changed or we have a dynamic group
  // or dynamic updates are requested (e.g. when changing atom types)
  // we need to recompute some normalization parameters

  if (dynamic) init_norm();

  invoked_array = update->ntimestep;

  // invoke full neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero the histogram counts

  for (i = 0; i < 1; i++) {
    for (j = 0; j < nbin; j++) {
      hist[i][j] = 0;
    }
  }
  // tally the three body
  // both atom i and j must be in fix group
  // itype,jtype must have been specified by user
  // consider I,J as one interaction even if neighbor pair is stored on 2 procs
  // tally I,J pair each time I is central atom, and each time J is central

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;

  // loop over full neighbor list of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];



    r = sqrt(xtmp*xtmp+ytmp*ytmp+ztmp*ztmp);

    i_bin = static_cast<int> (r*delrinv);


	
    if (i_bin >= nbin) continue;

    hist[0][i_bin] += 1.0;

  }

  // sum histograms across procs

  MPI_Allreduce(hist[0],histall[0],1*nbin,MPI_DOUBLE,MPI_SUM,world);

  // convert counts to g(r) and coord(r) and copy into output array
  // vfrac = fraction of volume in shell m

  
  double constant,vfrac,gr,ulower,uupper,vlower,vupper,normfac;

  if (domain->dimension==2) {
    constant = (domain->xprd*domain->yprd);
    constant = MY_PI/constant;

    //normfac = (icount[0] > 0) ? static_cast<double>(jcount[0])
    // - static_cast<double>(duplicates[0])/icount[0] : 0.0;
    normfac = static_cast<double>(icount[0]);
    normfac = normfac;

    set_array(constant, normfac);
	
  }
  /*
  else (domain->dimension == 3) {
    constant = 4.0*MY_PI / (3.0*domain->xprd*domain->yprd*domain->zprd);

    for (m = 0; m < npairs; m++) {
      normfac = (icount[m] > 0) ? static_cast<double>(jcount[m])
                - static_cast<double>(duplicates[m])/icount[m] : 0.0;
      ncoord = 0.0;
      for (ibin = 0; ibin < nbin; ibin++) {
        rlower = ibin*delr;
        rupper = (ibin+1)*delr;
        vfrac = constant * (rupper*rupper*rupper - rlower*rlower*rlower);
        if (vfrac * normfac != 0.0)
          gr = histall[m][ibin] / (vfrac * normfac * icount[m]);
        else gr = 0.0;
        if (icount[m] != 0)
          ncoord += gr * vfrac * normfac;
        array[ibin][1+2*m] = gr;
        array[ibin][2+2*m] = ncoord;
      }
    }

    }
  */
}

void ComputeOneBody::set_array(double constant, double normfac)
{
  double ulower,uupper;
  double vfrac;

  int i_bin;
  double gr;
      
  for (i_bin = 0; i_bin < nbin; i_bin++) {
      ulower = i_bin*delr;
      uupper = (i_bin+1)*delr;
      vfrac = (constant * (uupper*uupper - ulower*ulower)/2.0);
      if (vfrac * normfac != 0.0) {
	gr = histall[0][i_bin]/(vfrac *normfac);
      } else {
	gr = 0.0;
      }
      array[i_bin][0] = (i_bin + 0.5)*delr;
      array[i_bin][1] = gr;
  }

  return;
}

double ComputeOneBody::compute_rjk(double rij, double rik, double theta)
{
  return sqrt(rij*rij + rik*rik - 2*rij*rik*cos(theta));
}
