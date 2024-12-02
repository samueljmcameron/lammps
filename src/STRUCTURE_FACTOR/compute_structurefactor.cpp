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

#include "compute_structurefactor.h"

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

ComputeStructureFactor::ComputeStructureFactor(LAMMPS *lmp, int narg, char **arg) :
  ComputeStructureFactorBase(lmp, narg, arg)
{

  size_array_rows = hist_bins;
  size_array_cols = domain->dimension + 2*npairs;

  memory->create(array,size_array_rows,size_array_cols,"sq:array");
}

/* ---------------------------------------------------------------------- */

ComputeStructureFactor::~ComputeStructureFactor()
{

  memory->destroy(array);
}


void ComputeStructureFactor::write_qs()
{

  double qx,qy;
  if (domain->dimension == 2)
    for (int i = 0; i < nbin; i++) {
      qx = (i-nbin/2) * delta_q[0];
      for (int j = 0; j < nbin/2+1; j++) {
	array[j + (nbin/2+1)*i][0] = qx;
	array[j + (nbin/2+1)*i][1] = j*delta_q[1];
      }
    }
  
  else
    
    for (int i = 0; i < nbin; i++) {
      qx = (i-nbin/2) * delta_q[0];
      for (int j = 0; i < nbin; i++) {
	qy = (j-nbin/2) * delta_q[1];
	for (int k = 0; k < nbin/2+1; k++) {
	  array[k + (nbin/2+1)*(j+nbin*i)][0] = qx;
	  array[k + (nbin/2+1)*(j+nbin*i)][1] = qy;
	  array[k + (nbin/2+1)*(j+nbin*i)][2] = k*delta_q[2];
	  
	}
      }
    }
}


/* ---------------------------------------------------------------------- */

void ComputeStructureFactor::compute_array()
{
  ComputeStructureFactorBase::compute_array();
  int m;
  for (m = 0; m < npairs; m++) {
    for (int ibin = 0; ibin < size_array_rows; ibin++) {
      array[ibin][domain->dimension+2*m] = hist[2*m][ibin];
      array[ibin][domain->dimension+1+2*m] = hist[2*m+1][ibin];
    }
  }

}
