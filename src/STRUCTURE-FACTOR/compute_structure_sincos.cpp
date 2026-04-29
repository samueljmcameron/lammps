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


#include "compute_structure_sincos.h"

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

ComputeStructureSinCos::ComputeStructureSinCos(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg), hist(nullptr), q_counts(nullptr),histall(nullptr)
{
  if (narg != 4) error->all(FLERR,"Illegal compute structure/sincos command");

  int triclinic = domain->triclinic;
  if (triclinic == 1)
     error->all(FLERR,"Compute structure/sincos does not work with triclinic structures");

  array_flag = 1;
  extarray = 0;

  nq_bin = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nq_bin < 1) error->all(FLERR,"Illegal compute intermediatescattering command");

  dynamic = 0;
  
  size_array_rows = nq_bin;
  size_array_cols = 4;

  memory->create(q_counts,nq_bin,"sq:q_counts");
  memory->create(hist,2,nq_bin,"sq:hist");
  memory->create(histall,2,nq_bin,"sq:histall");
  memory->create(array,size_array_rows,size_array_cols,"sq:array");
  


    
}

/* ---------------------------------------------------------------------- */

ComputeStructureSinCos::~ComputeStructureSinCos()
{
  
  memory->destroy(q_counts);
  memory->destroy(hist);
  memory->destroy(histall);
  memory->destroy(array);
}


/* ---------------------------------------------------------------------- */

void ComputeStructureSinCos::init()
{
  init_norm();
}

/* ---------------------------------------------------------------------- */

void ComputeStructureSinCos::init_norm()
{

  // set 1st column of output array to bin coords

  
  delta_q = 2*M_PI/domain->prd[0];

  
  for (int idim = 1; idim < domain->dimension; idim ++ ) {
    
    double tmp = 2*M_PI/domain->prd[idim];
    
    if (tmp > delta_q)
      delta_q = 2*M_PI/domain->prd[idim];
    
  }


  double q,qx,qy,qz;
  int qbin;

  // initialize normalization, finite size correction, and changing atom counts

  // this will eventually be the qvalues, but calculating them in a precise way below.
  for (int i = 0; i < nq_bin; i++)
    array[i][0] = 0.0;

  for (int j = 0; j < nq_bin; j++)
    q_counts[j] = 0;

  
  if (domain->dimension == 2) {
    for (int iq = 0; iq < nq_bin; iq++) {
      qx = iq *delta_q;
      for (int jq = 0; jq < nq_bin; jq++) {
	qy = jq*delta_q;
	  
	q = sqrt(qx*qx+qy*qy);
	qbin = static_cast<int> (q/delta_q);
	  
	if (qbin < nq_bin) {
	  q_counts[qbin] += 1;
	  array[qbin][0] += q;
	}
      }
    }

  } else {
    for (int iq = 0; iq < nq_bin; iq++) {
      qx = iq *delta_q;
      for (int jq = 0; jq < nq_bin; jq++) {
	qy = jq*delta_q;
	for (int kq = 0; kq < nq_bin; kq++) {
	  qz = kq*delta_q;
	  q = sqrt(qx*qx+qy*qy+qz*qz);
	  qbin = static_cast<int> (q/delta_q);
	  
	  if (qbin < nq_bin) { 
	    q_counts[qbin] += 1;
	    array[qbin][0] += q;
	  }
	}
      }
    }
  }

  for (int i = 0; i < nq_bin; i++)
    array[i][0] /= q_counts[i];


}


/* ---------------------------------------------------------------------- */

void ComputeStructureSinCos::compute_array()
{

  int *mask = atom->mask;
  int *type = atom->type;
  double **x = atom->x;

  double q,qx,qy,qz,qdotr;
  int qbin;


  if (domain->box_change)
    init_norm();

  
  invoked_array = update->ntimestep;


  for (int i = 0; i < 2; i++)
    for (int j = 0; j < nq_bin; j++)
      hist[i][j] = 0;



  if (domain->dimension == 2) {
    for (int iq = 0; iq < nq_bin; iq++) {
      qx = iq *delta_q;
      for (int jq = 0; jq < nq_bin; jq++) {
	qy = jq*delta_q;
	  
	q = sqrt(qx*qx+qy*qy);
	qbin = static_cast<int> (q/delta_q);
	  
	if (qbin < nq_bin) {
	    
	  for (int i = 0; i < atom->nlocal; i++) {
	    qdotr = qx*x[i][0]+qy*x[i][1];
	    
	    hist[0][qbin] += cos(qdotr);
	    hist[1][qbin] += sin(qdotr);
	      
	  }
	}
      }
    }

  } else {
    for (int iq = 0; iq < nq_bin; iq++) {
      qx = iq *delta_q;
      for (int jq = 0; jq < nq_bin; jq++) {
	qy = jq*delta_q;
	for (int kq = 0; kq < nq_bin; kq++) {
	  qz = kq*delta_q;
	  
	  q = sqrt(qx*qx+qy*qy+qz*qz);
	  qbin = static_cast<int> (q/delta_q);
	  
	  if (qbin < nq_bin) {
	    
	    for (int i = 0; i < atom->nlocal; i++) {
	      qdotr = qx*x[i][0]+qy*x[i][1]+qz*x[i][2];
	      
	      hist[0][qbin] += cos(qdotr);
	      hist[1][qbin] += sin(qdotr);
	      
	    }
	  }
	}
      }
    }
  }
  MPI_Allreduce(hist[0],histall[0],2*nq_bin,MPI_DOUBLE,MPI_SUM,world);

  for (int i = 0 ; i < nq_bin; i++) {
    array[i][1] = q_counts[i];
    array[i][2] = histall[0][i]/sqrt(q_counts[i]);
    array[i][3] = histall[1][i]/sqrt(q_counts[i]);
  }

}
