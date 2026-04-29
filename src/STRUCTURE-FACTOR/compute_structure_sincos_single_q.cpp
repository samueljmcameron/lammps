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


#include "compute_structure_sincos_single_q.h"

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

ComputeStructureSinCosSingleQ::ComputeStructureSinCosSingleQ(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 4) error->all(FLERR,"Illegal compute structure/sincos/single command");

  int triclinic = domain->triclinic;
  if (triclinic == 1)
     error->all(FLERR,"Compute structure/sincos/single does not work with triclinic structures");

  vector_flag = 1;
  size_vector = 4;
  extvector = 0;
  single_bin = utils::inumeric(FLERR,arg[3],false,lmp);
  if (single_bin < 0) error->all(FLERR,"Illegal compute intermediatescattering/single command");

  

  vector = new double[size_vector];
    
}

/* ---------------------------------------------------------------------- */

ComputeStructureSinCosSingleQ::~ComputeStructureSinCosSingleQ()
{
  delete[] vector;  
}



/* ---------------------------------------------------------------------- */

void ComputeStructureSinCosSingleQ::init()
{
  init_norm();
}


/* ---------------------------------------------------------------------- */

void ComputeStructureSinCosSingleQ::init_norm()
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

  q_counts = 0;
  binned_q = 0.0;

  
  if (domain->dimension == 2) {
    for (int iq = 0; iq <= single_bin; iq++) {
      qx = iq *delta_q;
      for (int jq = 0; jq <= single_bin; jq++) {
	qy = jq*delta_q;
	  
	q = sqrt(qx*qx+qy*qy);
	qbin = static_cast<int> (q/delta_q);
	  
	if (qbin == single_bin) {
	  q_counts += 1;
	  binned_q += q;
	}
      }
    }

  } else {
    for (int iq = 0; iq <= single_bin; iq++) {
      qx = iq *delta_q;
      for (int jq = 0; jq <= single_bin; jq++) {
	qy = jq*delta_q;
	for (int kq = 0; kq <= single_bin; kq++) {
	  qz = kq*delta_q;
	  q = sqrt(qx*qx+qy*qy+qz*qz);
	  qbin = static_cast<int> (q/delta_q);
	  
	  if (qbin == single_bin)  {
	      q_counts += 1;
	      binned_q += q;
	  }
	}
      }
    }
  }

  binned_q /= q_counts;

}


/* ---------------------------------------------------------------------- */

void ComputeStructureSinCosSingleQ::compute_vector()
{

  int *mask = atom->mask;
  int *type = atom->type;
  double **x = atom->x;

  double q,qx,qy,qz,qdotr;
  int qbin;


  if (domain->box_change)
    init_norm();

  
  invoked_array = update->ntimestep;


  cosval = 0;
  sinval = 0;



  if (domain->dimension == 2) {
    for (int iq = 0; iq <= single_bin; iq++) {
      qx = iq *delta_q;
      for (int jq = 0; jq <= single_bin; jq++) {
	qy = jq*delta_q;
	  
	q = sqrt(qx*qx+qy*qy);
	qbin = static_cast<int> (q/delta_q);
	  
	if (qbin == single_bin) {
	    
	  for (int i = 0; i < atom->nlocal; i++) {
	    qdotr = qx*x[i][0]+qy*x[i][1];
	    
	    cosval += cos(qdotr);
	    sinval += sin(qdotr);
	      
	  }
	}
      }
    }

  } else {
    for (int iq = 0; iq <= single_bin; iq++) {
      qx = iq *delta_q;
      for (int jq = 0; jq <= single_bin; jq++) {
	qy = jq*delta_q;
	for (int kq = 0; kq <= single_bin; kq++) {
	  qz = kq*delta_q;
	  
	  q = sqrt(qx*qx+qy*qy+qz*qz);
	  qbin = static_cast<int> (q/delta_q);
	  
	  if (qbin == single_bin) {
	    
	    for (int i = 0; i < atom->nlocal; i++) {
	      qdotr = qx*x[i][0]+qy*x[i][1]+qz*x[i][2];
	      
	      cosval += cos(qdotr);
	      sinval += sin(qdotr);
	      
	    }
	  }
	}
      }
    }
  }


  vector[0] = binned_q;
  vector[1] = q_counts;
  MPI_Allreduce(&cosval,&vector[2],1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&sinval,&vector[3],1,MPI_DOUBLE,MPI_SUM,world);
  vector[2] /= sqrt(q_counts);
  vector[3] /= sqrt(q_counts);


}
