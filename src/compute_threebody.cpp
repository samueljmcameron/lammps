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

#include "compute_threebody.h"
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

ComputeThreeBody::ComputeThreeBody(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  rdfpair(NULL), ilo(NULL), ihi(NULL), jlo(NULL), jhi(NULL),
  hist(NULL), histall(NULL), typecount(NULL), icount(NULL), jcount(NULL),
  duplicates(NULL)
{
  if (narg != 5 && narg != 7 && narg != 9 && narg != 11)
    error->all(FLERR,"Illegal compute three body command");

  array_flag = 1;
  extarray = 0;

  nbin_dist = force->inumeric(FLERR,arg[3]);
  if (nbin_dist < 1) error->all(FLERR,"Illegal compute three body command");

  nbin_theta = force->inumeric(FLERR,arg[4]);
  if (nbin_theta < 1) error->all(FLERR,"Illegal compute three body command");

  nbin_total = nbin_dist*nbin_dist*nbin_theta;
  // optional args
  // nargpair = # of pairwise args, starting at iarg = 5

  cutflag = 0;
  lower_cut = 0;
  condenseflag = 0;

  if (narg > 5) {
    if (strcmp(arg[5],"cutoff") == 0) {
      cutoff_user = force->numeric(FLERR,arg[6]);
      if (cutoff_user <= 0.0) cutflag = 0;
      else cutflag = 1;
      if (narg > 7) {
	if (strcmp(arg[7],"lower_cutoff") == 0) {
	  lower_cut = force->numeric(FLERR,arg[8]);
	  if (lower_cut <= 0.0) lower_cut = 0;
	  if (narg == 11) {
	    if (strcmp(arg[9],"condense")==0) {
	      if (strcmp(arg[10],"yes")==0) { 
		condenseflag = 1;
	      }
	    } else error->all(FLERR,"Illegal compute threebody command");
	  }
	} else if (strcmp(arg[7],"condense")==0) {
	  if (strcmp(arg[8],"yes")==0) { 
	    condenseflag = 1;
	  }
	  if (narg == 11) {
	    if (strcmp(arg[9],"lower_cutoff") == 0) {
	      lower_cut = force->numeric(FLERR,arg[10]);
	      if (lower_cut <= 0.0) lower_cut = 0;
	    } else error->all(FLERR,"Illegal compute threebody command");
	  }
	} else error->all(FLERR,"Illegal compute threebody command");
      }
    } else if (strcmp(arg[5],"lower_cutoff") == 0) {
      lower_cut = force->numeric(FLERR,arg[6]);
      if (lower_cut <= 0.0) lower_cut = 0;
      if (narg > 7) {
	if (strcmp(arg[7],"cutoff") == 0) {
	  cutoff_user = force->numeric(FLERR,arg[8]);
	  if (cutoff_user <= 0.0) cutflag = 0;
	  else cutflag = 1;
	  if (narg == 11) {
	    if (strcmp(arg[9],"condense")==0) {
	      if (strcmp(arg[10],"yes")==0) { 
		condenseflag = 1;
	      }
	    } else error->all(FLERR,"Illegal compute threebody command");
	  }
	} else if (strcmp(arg[7],"condense")==0) {
	  if (strcmp(arg[8],"yes")==0) { 
	    condenseflag = 1;
	  }
	  if (narg == 11) {
	    if (strcmp(arg[9],"cutoff") == 0) {
	      cutoff_user = force->numeric(FLERR,arg[10]);
	      if (cutoff_user <= 0.0) cutflag = 0;
	      else cutflag = 1;
	    } else error->all(FLERR,"Illegal compute threebody command");
	  }
	} else error->all(FLERR,"Illegal compute threebody command");
      }
    } else if (strcmp(arg[5],"condense") == 0) {
      if (strcmp(arg[6],"yes")==0) { 
	condenseflag = 1;
      }
      if (narg > 7) {
	if (strcmp(arg[7],"cutoff") == 0) {
	  cutoff_user = force->numeric(FLERR,arg[8]);
	  if (cutoff_user <= 0.0) cutflag = 0;
	  else cutflag = 1;
	  if (narg == 11) {
	    if (strcmp(arg[9],"lower_cutoff") == 0) {
	      lower_cut = force->numeric(FLERR,arg[10]);
	      if (lower_cut <= 0.0) lower_cut = 0;
	    } else error->all(FLERR,"Illegal compute threebody command");
	  }
	} else if (strcmp(arg[7],"lower_cutoff") == 0) {
	  lower_cut = force->numeric(FLERR,arg[8]);
	  if (lower_cut <= 0.0) lower_cut = 0;
	  if (narg == 11) {
	    if (strcmp(arg[9],"cutoff") == 0) {
	      cutoff_user = force->numeric(FLERR,arg[10]);
	      if (cutoff_user <= 0.0) cutflag = 0;
	      else cutflag = 1;
	    } else error->all(FLERR,"Illegal compute threebody command");
	  }
	} else error->all(FLERR,"Illegal compute threebody command");
      }
    } else error->all(FLERR,"Illegal compute threebody command");
  }

  printf("lower_cut = %f, cutoff_user = %f, condenseflag = %d\n",
	 lower_cut,cutoff_user,condenseflag);
  if (force->newton_pair) {
    error->all(FLERR,"compute threebody command is incompatible with "
	       "newton pair being on.");
  }

  // pairwise args, fix it to be always one set of pairs for now,
  // but might change this later to deal with density differences

  npairs = 1;


  int length;

  deltheta = MY_PI/nbin_theta;

  delthetainv = 1.0/deltheta;


  
  if (condenseflag && cutflag) {
    deldist = (cutoff_user - lower_cut) / nbin_dist;
    length = length_of_array();
    size_array_rows = length;
  } else {
    if (condenseflag) {
      error->all(FLERR,"compute threebody command requires a "
		 "user-specified cutoff in order to condense the "
		 "output.");
    }
    size_array_rows = nbin_total;
  }

  size_array_cols = 4;

  
  int ntypes = atom->ntypes;
  if (ntypes != 1) {
    error->all(FLERR,"Cannot compute threebody for system with multiple "
	       "atom types.");
  }
    
  memory->create(rdfpair,npairs,ntypes+1,ntypes+1,"rdf:rdfpair");

  ilo = new int[npairs];
  ihi = new int[npairs];
  jlo = new int[npairs];
  jhi = new int[npairs];

  ilo[0] = 1; ihi[0] = 1;
  jlo[0] = 1; jhi[0] = 1;
  

  nrdfpair = 1;
  

  rdfpair[0][1][1] = 1;


  memory->create(hist,nbin_theta,nbin_dist,nbin_dist,"rdf:hist");
  memory->create(histall,nbin_theta,nbin_dist,nbin_dist,"rdf:histall");
  memory->create(array,size_array_rows,size_array_cols,"rdf:array");  
  typecount = new int[ntypes+1];
  icount = new int[npairs];
  jcount = new int[npairs];
  duplicates = new int[npairs];
  
  dynamic = 0;
  natoms_old = 0;
}

/* ---------------------------------------------------------------------- */

ComputeThreeBody::~ComputeThreeBody()
{
  memory->destroy(rdfpair);
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
}

/* ---------------------------------------------------------------------- */

void ComputeThreeBody::init()
{

  if (!force->pair && !cutflag)
    error->all(FLERR,"Compute threebody requires a pair style be defined "
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
      error->all(FLERR,"Compute threebody cutoff exceeds ghost atom range - "
                 "use comm_modify cutoff command");
    if (force->pair && mycutneigh < force->pair->cutforce + skin)
      if (comm->me == 0)
        error->warning(FLERR,"Compute threebody cutoff less than neighbor "
		       "cutoff - "
                       "forcing a needless neighbor list build");

    deldist = (cutoff_user - lower_cut) / nbin_dist;
  } else deldist = (force->pair->cutforce - lower_cut) / nbin_dist;

  deldistinv = 1.0/deldist;


  // since this file is going to be huge, don't include coordinates in it.

  //  for (int i = 0; i < nbin; i++)
  //array[i][0] = (i+0.5) * delr;

  // initialize normalization, finite size correction, and changing atom counts
  // ...not sure what this does

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

void ComputeThreeBody::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeThreeBody::init_norm()
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

  // icount = # of I atoms participating in I,J pairs for each histogram
  // jcount = # of J atoms participating in I,J pairs for each histogram
  // duplicates = # of atoms in both groups I and J for each histogram

  //here
  icount[0] = 0;
  for (i = ilo[0]; i <= ihi[0]; i++) icount[0] += typecount[i];
  jcount[0] = 0;
  for (i = jlo[0]; i <= jhi[0]; i++) jcount[0] += typecount[i];
  duplicates[0] = 0;
  for (i = ilo[0]; i <= ihi[0]; i++)
    for (j = jlo[0]; j <= jhi[0]; j++)
      if (i == j) duplicates[0] += typecount[i];


  //here
  int *scratch = new int[1];
  MPI_Allreduce(icount,scratch,npairs,MPI_INT,MPI_SUM,world);
  for (i = 0; i < npairs; i++) icount[i] = scratch[i];
  MPI_Allreduce(jcount,scratch,npairs,MPI_INT,MPI_SUM,world);
  for (i = 0; i < npairs; i++) jcount[i] = scratch[i];
  MPI_Allreduce(duplicates,scratch,npairs,MPI_INT,MPI_SUM,world);
  for (i = 0; i < npairs; i++) duplicates[i] = scratch[i];

  delete [] scratch;
}

/* ---------------------------------------------------------------------- */

void ComputeThreeBody::compute_array()
// ======================================================================
// Only counting triplets where all three distances are less than cutoff
// away.
//
// ======================================================================
{
  int i,j,m,ii,jj,inum,jnum,itype,jtype,ipair,jpair;
  int ij_bin,ik_bin,theta_bin,dum_jk_bin;
  int k,kk,knum,ktype,kpair;
  double xtmp,ytmp,ztmp;
  double xij,yij,zij,xik,yik,zik,rij,rik,rjk,theta;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double factor_lj,factor_coul;

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

  for (i = 0; i < nbin_theta; i++)
    for (j = 0; j < nbin_dist; j++)
      for (k = 0; k < nbin_dist; k++)
	hist[i][j][k] = 0;

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

      for (kk = 0; kk < jnum; kk++) {
	if (jj == kk) continue;
	k = jlist[kk];
	factor_lj = special_lj[sbmask(k)];
	factor_coul = special_coul[sbmask(k)];
	
	k &= NEIGHMASK;


	if (factor_lj == 0.0 && factor_coul == 0.0) continue;

	if (!(mask[k] & groupbit)) continue;

	
	xij = x[j][0]-xtmp;
	yij = x[j][1]-ytmp;
	zij = x[j][2]-ztmp;

	xik = x[k][0]-xtmp;
	yik = x[k][1]-ytmp;
	zik = x[k][2]-ztmp;

	rij = sqrt(xij*xij+yij*yij+zij*zij);
	rik = sqrt(xik*xik+yik*yik+zik*zik);

	rjk = sqrt((xik-xij)*(xik-xij)
		   +(yik-yij)*(yik-yij)+(zik-zij)*(zik-zij));
	
	theta = acos((xij*xik + yij*yik + zij*zik)/(rij*rik));

	ij_bin = static_cast<int> ((rij-lower_cut)*deldistinv);
	ik_bin = static_cast<int> ((rik-lower_cut)*deldistinv);
	dum_jk_bin = static_cast<int> ((rjk-lower_cut)*deldistinv);
	theta_bin = static_cast<int> (theta*delthetainv);
	
	if (ij_bin >= nbin_dist || ik_bin >=nbin_dist
	    || dum_jk_bin >= nbin_dist) continue;
	if (theta_bin >= nbin_theta) {
	  printf("theta = %f, theta_bin = %d\n",theta,theta_bin);
	  error->all(FLERR,"theta > 3.1415 somehow? "
		     "Error in compute threebody.");
	}


	hist[theta_bin][ij_bin][ik_bin] += 1.0;

      }
    }
  }

  // sum histograms across procs

  MPI_Allreduce(hist[0][0],histall[0][0],nbin_total,MPI_DOUBLE,MPI_SUM,world);

  // convert counts to g(r) and coord(r) and copy into output array
  // vfrac = fraction of volume in shell m
  // npairs = number of pairs, corrected for duplicates
  // duplicates = pairs in which both atoms are the same

  
  double constant,vfrac,gr,ulower,uupper,vlower,vupper,normfac;
  int flat_index;

  if (domain->dimension==2) {
    constant = (domain->xprd*domain->yprd);
    constant = constant*constant;
    // the 2 in front of the deltheta below is necessary because
    // deltatheta = pi/nbins over the [0,pi] region, but really
    // the g_3bod is defined from 0 to 2pi. So we are double-counting
    // in theta since we are only letting theta [0,pi).
    constant = 2*MY_PI*2*deltheta/constant;

    //normfac = (icount[0] > 0) ? static_cast<double>(jcount[0])
    // - static_cast<double>(duplicates[0])/icount[0] : 0.0;
    normfac = static_cast<double>(icount[0]);
    normfac = normfac*normfac*normfac;

    if (condenseflag) {
      set_condensed_array(constant, normfac);
    } else {
      set_array(constant, normfac);
    }
	
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

int ComputeThreeBody::length_of_array()
{
  double cutoff = cutoff_user;

  double rij, rik, theta;
  double rjk;
  
  int count = 0;
  for (int i = 0; i < nbin_theta; i++) {
    for (int j = 0; j < nbin_dist; j++) {
      for (int k = 0; k < nbin_dist; k++) {
	theta = (i + 0.5)*deltheta;
	rij = (j + 0.5)*deldist+lower_cut;
	rik = (k + 0.5)*deldist+lower_cut;
	rjk = compute_rjk(rij,rik,theta);
	if (rjk < cutoff) {
	  count += 1;
	}
      }
    }
  }
  return count;
}

void ComputeThreeBody::set_array(double constant, double normfac)
{
  double ulower,uupper,vlower,vupper;
  double vfrac;

  int theta_bin,ij_bin,ik_bin;
  double gr;

  int flat_index;
  
  for (theta_bin = 0; theta_bin < nbin_theta; theta_bin ++ ) {
    
    for (ij_bin = 0; ij_bin < nbin_dist; ij_bin++) {
      ulower = ij_bin*deldist;
      uupper = (ij_bin+1)*deldist;
      for (ik_bin = 0; ik_bin < nbin_dist; ik_bin++) {
	vlower = ik_bin*deldist;
	vupper = (ik_bin+1)*deldist;	  
	vfrac = (constant * (uupper*uupper - ulower*ulower)/2.0
		 * (vupper*vupper - vlower*vlower)/2.0);
	//if (histall[theta_bin][ij_bin][ik_bin] < -1e-9) {
	//  printf("histogram  is less than zero!\n");
	//}
	if (vfrac * normfac != 0.0) {
	  gr = histall[theta_bin][ij_bin][ik_bin]/(vfrac *normfac);
	} else {
	  gr = 0.0;
	}
	flat_index = ik_bin + (ij_bin+ theta_bin * nbin_dist) *nbin_dist;
	array[flat_index][0] = (theta_bin + 0.5)*deltheta;
	array[flat_index][1] = (ij_bin + 0.5)*deldist+lower_cut;
	array[flat_index][2] = (ik_bin + 0.5)*deldist+lower_cut;
	array[flat_index][3] = gr;
      }
    }
  }

  return;
}

void ComputeThreeBody::set_condensed_array(double constant, double normfac)
{
  double ulower,uupper,vlower,vupper;
  double vfrac;
  double cutoff = cutoff_user;

  int theta_bin,ij_bin,ik_bin;
  double gr;
  double rij,rik,rjk,theta;

  int flat_index = 0;
  for (theta_bin = 0; theta_bin < nbin_theta; theta_bin ++ ) {
    
    for (ij_bin = 0; ij_bin < nbin_dist; ij_bin++) {
      ulower = ij_bin*deldist;
      uupper = (ij_bin+1)*deldist;
      for (ik_bin = 0; ik_bin < nbin_dist; ik_bin++) {
	theta = (theta_bin + 0.5)*deltheta;
	rij = (ij_bin + 0.5)*deldist+lower_cut;
	rik = (ik_bin + 0.5)*deldist+lower_cut;
	rjk = compute_rjk(rij,rik,theta);
	if (rjk < cutoff) {
	  vlower = ik_bin*deldist;
	  vupper = (ik_bin+1)*deldist;	  
	  vfrac = (constant * (uupper*uupper - ulower*ulower)/2.0
		   * (vupper*vupper - vlower*vlower)/2.0);

	  if (vfrac * normfac != 0.0) {
	    gr = histall[theta_bin][ij_bin][ik_bin]/(vfrac *normfac);
	  } else {
	    gr = 0.0;
	  }
	  array[flat_index][0] = theta;
	  array[flat_index][1] = rij;
	  array[flat_index][2] = rik;
	  array[flat_index][3] = gr;
	  flat_index += 1;
	}
      }
    }
  }

  return;
}


double ComputeThreeBody::compute_rjk(double rij, double rik, double theta)
{
  return sqrt(rij*rij + rik*rik - 2*rij*rik*cos(theta));
}
