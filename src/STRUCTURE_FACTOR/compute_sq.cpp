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

#include "compute_sq.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

ComputeSQ::ComputeSQ(LAMMPS *lmp, int narg, char **arg) :
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

  cutflag = 0;

  int iarg;
  for (iarg = 4; iarg < narg; iarg++)
    if (strcmp(arg[iarg],"cutoff") == 0) break;

  int nargpair = iarg - 4;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"cutoff") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute sq command");
      cutoff_user = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (cutoff_user <= 0.0) cutflag = 0;
      else cutflag = 1;
      iarg += 2;
    } else error->all(FLERR,"Illegal compute sq command");
  }

  // pairwise args

  if (nargpair == 0) npairs = 1;
  else {
    if (nargpair % 2) error->all(FLERR,"Illegal compute sq command");
    npairs = nargpair/2;
  }

  size_array_rows = nbin;
  size_array_cols = 1 + 2*npairs;

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
    iarg = 4;
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

  
  memory->create(hist,2*npairs,nbin,"sq:hist");
  memory->create(histall,2*npairs,nbin,"sq:histall");
  memory->create(array,nbin,1+2*npairs,"sq:array");
  typecount = new int[ntypes+1];
  icount = new int[npairs];
  jcount = new int[npairs];
  duplicates = new int[npairs];
  ibins = new int[nbin*(nbin-1)];

  dynamic = 0;
  natoms_old = 0;
}

/* ---------------------------------------------------------------------- */

ComputeSQ::~ComputeSQ()
{
  memory->destroy(sqpair);
  memory->destroy(nsqpair);
  delete [] ilo;
  delete [] ihi;
  delete [] jlo;
  delete [] jhi;
  memory->destroy(hist);
  memory->destroy(histall);
  memory->destroy(array);
  delete [] typecount;
  delete [] icount;
  delete [] jcount;
  delete [] duplicates;
  delete [] ibins;
}

/* ---------------------------------------------------------------------- */

void ComputeSQ::init()
{

  if (!force->pair && !cutflag)
    error->all(FLERR,"Compute sq requires a pair style be defined "
               "or cutoff specified");

  if (cutflag) {
    double skin = neighbor->skin;
    mycutneigh = cutoff_user + skin;

    double cutghost;            // as computed by Neighbor and Comm
    if (force->pair)
      cutghost = MAX(force->pair->cutforce+skin,comm->cutghostuser);
    else
      cutghost = comm->cutghostuser;

    if (mycutneigh > cutghost)
      error->all(FLERR,"Compute sq cutoff exceeds ghost atom range - "
                 "use comm_modify cutoff command");
    if (force->pair && mycutneigh < force->pair->cutforce + skin)
      if (comm->me == 0)
        error->warning(FLERR,"Compute sq cutoff less than neighbor cutoff - "
                       "forcing a needless neighbor list build");

    delr = cutoff_user / nbin;
  } else delr = force->pair->cutforce / nbin;

  delrinv = 1.0/delr;

  // set 1st column of output array to bin coords

  double min_L = domain->xprd;
  min_L =  min_L > domain->yprd ? domain->yprd : min_L;
  if (domain->dimension != 2)
    min_L =  min_L > domain->zprd ? domain->zprd : min_L;

  delta_q = 2*M_PI/min_L;

  // only count first quadrant excluding the x = 0 line
  // the q = 0 point is uninteresting at is a delta function
  for (int nx = 0; nx < nbin; nx++) 
    for (int ny = 1; ny < nbin; ny++)
      ibins[nx+nbin*(ny-1)] = static_cast<int>(sqrt(nx*nx + ny*ny)/sqrt(2)-0.5);

    
  
  for (int i = 0; i < nbin; i++)
    array[i][0] = (i+1) * delta_q*sqrt(2);

  // initialize normalization, finite size correction, and changing atom counts

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

  auto req = neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
  if (cutflag) req->set_cutoff(mycutneigh);
}

/* ---------------------------------------------------------------------- */

void ComputeSQ::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeSQ::init_norm()
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

void ComputeSQ::compute_array()
{
  int i,j,m,ii,jj,inum,jnum,itype,jtype,ipair,jpair,ihisto,ibin;
  double xtmp,ytmp,ztmp,delx,dely,delz,r,qx,qy,Stmp_real;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double factor_lj,factor_coul;

  if (natoms_old != atom->natoms) {
    dynamic = 1;
    natoms_old = atom->natoms;
  }

  // if the number of atoms has changed or we have a dynamic group
  // or dynamic updates are requested (e.g. when changing atom types)
  // we need to recompute some normalization parameters

  if (dynamic) init_norm();

  invoked_array = update->ntimestep;

  // invoke half neighbor list (will copy or build if necessary)

  neighbor->build_one(list);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero the histogram counts

  for (i = 0; i < 2*npairs; i++)
    for (j = 0; j < nbin; j++)
      hist[i][j] = 0;

  // tally the SQ
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
  int newton_pair = force->newton_pair;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      // if both weighting factors are 0, skip this pair
      // could be 0 and still be in neigh list for long-range Coulombics
      // want consistency with non-charged pairs which wouldn't be in list

      if (factor_lj == 0.0 && factor_coul == 0.0) continue;

      if (!(mask[j] & groupbit)) continue;
      jtype = type[j];
      ipair = nsqpair[itype][jtype];
      jpair = nsqpair[jtype][itype];
      if (!ipair && !jpair) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      for (int nx = 0; nx < nbin; nx ++ ) {
	qx = delta_q*nx;
	for (int ny = 1; ny < nbin; ny ++ ) {

	  ibin = ibins[nx + nbin*(ny-1)];
	  if (ibin > nbin) continue;
	  qy = delta_q*ny;

	  Stmp_real = 2*(cos(qx*delx+qy*dely) + cos(qy*delx - qx*dely));
	  //Stmp_imag = sin(qx*delx+qy*dely);
	  
	  for (ihisto = 0; ihisto < ipair; ihisto++) {
	    m = sqpair[ihisto][itype][jtype];
	    hist[m][ibin] += 1.0;
	    hist[m+npairs][ibin] += Stmp_real;
	  }
	  if (newton_pair || j < nlocal) {
	    for (ihisto = 0; ihisto < jpair; ihisto++) {
	      m = sqpair[ihisto][jtype][itype];
	      hist[m][ibin] += 1.0;
	      hist[m+npairs][ibin] += Stmp_real;
	    }
	  }
	}
      }
    }
  }

  // sum histograms across procs

  MPI_Allreduce(hist[0],histall[0],npairs*nbin,MPI_DOUBLE,MPI_SUM,world);

  for (m = 0; m < npairs; m++) {
    for (ibin = 0; ibin < nbin; ibin++) {
      array[ibin][1+2*m] = hist[m][ibin];
      array[ibin][2+2*m] = hist[m+npairs][ibin];
    }
  }

}
