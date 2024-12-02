// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Paul Crozier (SNL), Jeff Greathouse (SNL)
------------------------------------------------------------------------- */

#include "compute_structurefactor_base.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeStructureFactorBase::ComputeStructureFactorBase(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  sqpair(nullptr), nsqpair(nullptr), ilo(nullptr), ihi(nullptr), jlo(nullptr), jhi(nullptr),
  hist(nullptr), histall(nullptr), typecount(nullptr), icount(nullptr), jcount(nullptr),
  duplicates(nullptr)
{
  if (narg < 4) error->all(FLERR,"Illegal compute sq command");

  array_flag = 1;
  extarray = 0;

  nbin = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nbin < 1) error->all(FLERR,"Illegal compute sq command");

  // optional args
  // nargpair = # of pairwise args, starting at iarg = 4

  int nargpair = narg - 4;
  // pairwise args

  if (nargpair == 0) npairs = 1;
  else {
    if (nargpair % 2) error->all(FLERR,"Illegal compute sq command");
    npairs = nargpair/2;
  }

  hist_bins = 1;
  for (int idim = 0; idim < domain->dimension-1; idim ++ ) hist_bins *= nbin;

  hist_bins *= (nbin/2+1);

  int ntypes = atom->ntypes;
  memory->create(sqpair,npairs,ntypes+1,ntypes+1,"sq:sqpair");
  memory->create(nsqpair,ntypes+1,ntypes+1,"sq:nsqpair");
  ilo = new int[npairs];
  ihi = new int[npairs];
  jlo = new int[npairs];
  jhi = new int[npairs];

  if (nargpair == 0) {
    ilo[0] = 1; ihi[0] = ntypes;
    jlo[0] = 1; jhi[0] = ntypes;
  } else {
    int iarg = 4;
    for (int ipair = 0; ipair < npairs; ipair++) {
      utils::bounds(FLERR,arg[iarg],1,atom->ntypes,ilo[ipair],ihi[ipair],error);
      utils::bounds(FLERR,arg[iarg+1],1,atom->ntypes,jlo[ipair],jhi[ipair],error);
      if (ilo[ipair] > ihi[ipair] || jlo[ipair] > jhi[ipair])
        error->all(FLERR,"Illegal compute sq command");
      iarg += 2;
    }
  }

  int i,j;
  for (i = 1; i <= ntypes; i++)
    for (j = 1; j <= ntypes; j++)
      nsqpair[i][j] = 0;

  int ihisto;
  for (int m = 0; m < npairs; m++)
    for (i = ilo[m]; i <= ihi[m]; i++)
      for (j = jlo[m]; j <= jhi[m]; j++) {
        ihisto = nsqpair[i][j]++;
        sqpair[ihisto][i][j] = m;
      }

  
  
  memory->create(hist,2*npairs,hist_bins,"sq:hist");
  memory->create(histall,2*npairs,hist_bins,"sq:histall");

  typecount = new int[ntypes+1];
  icount = new int[npairs];
  jcount = new int[npairs];
  duplicates = new int[npairs];
  delta_q = new double[domain->dimension];

  dynamic = 0;
  natoms_old = 0;
}

/* ---------------------------------------------------------------------- */

ComputeStructureFactorBase::~ComputeStructureFactorBase()
{
  memory->destroy(sqpair);
  memory->destroy(nsqpair);
  delete [] ilo;
  delete [] ihi;
  delete [] jlo;
  delete [] jhi;
  memory->destroy(hist);
  memory->destroy(histall);
  delete [] typecount;
  delete [] icount;
  delete [] jcount;
  delete [] duplicates;
  delete [] delta_q;
}

/* ---------------------------------------------------------------------- */

void ComputeStructureFactorBase::init()
{

  // set 1st column of output array to bin coords

  

  for (int idim = 0; idim < domain->dimension; idim ++ )
    delta_q[idim] = 2*M_PI/domain->prd[idim];
      
  // initialize normalization, finite size correction, and changing atom counts

  write_qs();

  natoms_old = atom->natoms;
  dynamic = group->dynamic[igroup];
  if (dynamic_user) dynamic = 1;
  init_norm();

}


/* ---------------------------------------------------------------------- */

void ComputeStructureFactorBase::init_norm()
{

  
  int i,j,m;

  // count atoms of each type that are also in group

  const int nlocal = atom->nlocal;
  const int ntypes = atom->ntypes;
  const int * const mask = atom->mask;
  const int * const type = atom->type;

  for (i = 1; i <= ntypes; i++) typecount[i] = 0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) typecount[type[i]]++;

  // icount = # of I atoms participating in I,J pairs for each histogram
  // jcount = # of J atoms participating in I,J pairs for each histogram
  // duplicates = # of atoms in both groups I and J for each histogram

  for (m = 0; m < npairs; m++) {
    icount[m] = 0;
    for (i = ilo[m]; i <= ihi[m]; i++) icount[m] += typecount[i];
    jcount[m] = 0;
    for (i = jlo[m]; i <= jhi[m]; i++) jcount[m] += typecount[i];
    duplicates[m] = 0;
    for (i = ilo[m]; i <= ihi[m]; i++)
      for (j = jlo[m]; j <= jhi[m]; j++)
        if (i == j) duplicates[m] += typecount[i];
  }

  int *scratch = new int[npairs];
  MPI_Allreduce(icount,scratch,npairs,MPI_INT,MPI_SUM,world);
  for (i = 0; i < npairs; i++) icount[i] = scratch[i];
  MPI_Allreduce(jcount,scratch,npairs,MPI_INT,MPI_SUM,world);
  for (i = 0; i < npairs; i++) jcount[i] = scratch[i];
  MPI_Allreduce(duplicates,scratch,npairs,MPI_INT,MPI_SUM,world);
  for (i = 0; i < npairs; i++) duplicates[i] = scratch[i];
  delete [] scratch;
}

/* ---------------------------------------------------------------------- */

void ComputeStructureFactorBase::compute_histogram()
{

  int *mask = atom->mask;
  int *type = atom->type;
  double **x = atom->x;

  if (natoms_old != atom->natoms) {
    dynamic = 1;
    natoms_old = atom->natoms;
  }

  // if the number of atoms has changed or we have a dynamic group
  // or dynamic updates are requested (e.g. when changing atom types)
  // we need to recompute some normalization parameters

  if (dynamic) init_norm();

  invoked_array = update->ntimestep;


  for (int i = 0; i < 2*npairs; i++)
    for (int j = 0; j < hist_bins; j++)
      hist[i][j] = 0;


  double *list;  // list has entries for x,y,z,mask, and type
  memory->create(list,atom->nlocal*5,"sq:list");
  
  // need to store
  int n = 0;
  for (int i = 0; i < atom->nlocal; i++) {
    list[n++] = ubuf(mask[i]).d;
    list[n++] = ubuf(type[i]).d;
    list[n++] = x[i][0];
    list[n++] = x[i][1];
    list[n++] = x[i][2];

  }

    

  comm->ring(n,sizeof(double),list,1,callback,nullptr,(void *) this);



  memory->destroy(list);


  MPI_Allreduce(hist[0],histall[0],2*npairs*hist_bins,MPI_DOUBLE,MPI_SUM,world);

}


void ComputeStructureFactorBase::callback(int n, char *cbuf, void *ptr)
{
  auto sqptr = (ComputeStructureFactorBase *) ptr;
  auto list = (double *) cbuf;

  int groupbit = sqptr->groupbit;
  int nlocal = sqptr->atom->nlocal;
  
  double ** hist = sqptr->hist;
  double **x = sqptr->atom->x;
  int *mask = sqptr->atom->mask;
  int *type = sqptr->atom->type;
  int ***sqpair = sqptr->sqpair;
  int **nsqpair = sqptr->nsqpair;
  int nbin = sqptr->nbin;
  double *delta_q = sqptr->delta_q;
  int dimension = sqptr->domain->dimension;
  bigint tstep = sqptr->update->ntimestep;
  int npairs = sqptr->npairs;
  int hist_bins = sqptr->hist_bins;

  double xtmp,ytmp,ztmp,fac,delx,dely,delz;
  int itype,jtype,ipair,jpair;


  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    
    int j = 0;
    while (j < n) {
      if (!((int) ubuf(list[j++]).i & groupbit)) {
	j += 4;
	continue;
      }
      jtype = (int) ubuf(list[j++]).i;

      ipair = nsqpair[itype][jtype];
      jpair = nsqpair[jtype][itype];
      if (!ipair && !jpair) {
	j += 3;
	continue;
      }
      
      delx = xtmp - list[j++];
      dely = ytmp - list[j++];
      delz = ztmp - list[j++];

      int ibin,m;
      if (dimension == 2) {
	for (int nx = -nbin/2; nx < nbin/2; nx++ )
	  for (int ny = 0; ny < nbin/2+1; ny++ ) {
	    fac = delta_q[0]*nx*delx + delta_q[1]*ny*dely;
	    ibin = ny  + (nbin/2+1)*(nx+nbin/2);
	    for (int ihisto = 0; ihisto < ipair; ihisto++) {
	      m = sqpair[ihisto][itype][jtype];
	      hist[2*m][ibin] += cos(fac);
	      hist[2*m+1][ibin] += sin(fac);
	    }
	  }

	
	
      } else {
	for (int nx = -nbin/2; nx < nbin/2; nx++ )
	  for (int ny = -nbin/2; ny < nbin/2; ny++ )
	    for (int nz = 0; nz < nbin/2+1; nz++ ) {
	      fac = delta_q[0]*nx*delx + delta_q[1]*ny*dely + delta_q[2]*nz*delz;
	      ibin = nz  + (nbin/2+1)*((ny+nbin/2) + nbin*(nx+nbin/2));
	      for (int ihisto = 0; ihisto < ipair; ihisto++) {
		m = sqpair[ihisto][itype][jtype];
		hist[2*m][ibin] += cos(fac);
		hist[2*m+1][ibin] += sin(fac);
	      }
	    }
      }
      
      
    }

    
  }
  
}
