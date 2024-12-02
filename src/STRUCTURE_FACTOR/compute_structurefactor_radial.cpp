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

#include "compute_structurefactor_radial.h"

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

ComputeStructureFactorRadial::ComputeStructureFactorRadial(LAMMPS *lmp, int narg, char **arg) :
  ComputeStructureFactorBase(lmp, narg, arg)
{
  
  size_array_rows = nbin/2+1;


  // first col of array stores |\bf{q}|, second col stores counts
  size_array_cols = 2 + npairs;
  memory->create(array,size_array_rows,size_array_cols,"sq:array");
  counts = new int[nbin/2+1];
}

/* ---------------------------------------------------------------------- */

ComputeStructureFactorRadial::~ComputeStructureFactorRadial()
{

  memory->destroy(array);
  delete [] counts;
}


void ComputeStructureFactorRadial::write_qs()
{


  sqrt_dq2 = 0;

  for (int i = 0; i < domain->dimension; i++)
    sqrt_dq2 += delta_q[i]*delta_q[i];

  sqrt_dq2 = sqrt(sqrt_dq2/domain->dimension);

  for (int i = 0; i < nbin/2+1; i++) {
    array[i][0] = sqrt_dq2*i;
    counts[i] = 0;
  }
  
  double qx,qy,qz;

  int ibin;
  if (domain->dimension == 2)
    for (int i = 0; i < nbin; i++) {
      qx = (i-nbin/2) * delta_q[0];
      for (int j = 0; j < nbin/2+1; j++) {
	qy = j*delta_q[1];
	ibin = q_to_bin(sqrt(qx*qx+qy*qy));
	if (ibin <= nbin/2) counts[ibin] += 1;
      }
    }
  
  else
    
    for (int i = 0; i < nbin; i++) {
      qx = (i-nbin/2) * delta_q[0];
      for (int j = 0; i < nbin; i++) {
	qy = (j-nbin/2) * delta_q[1];
	for (int k = 0; k < nbin/2+1; k++) {
	  qz = k*delta_q[2];
	  ibin = q_to_bin(sqrt(qx*qx+qy*qy+qz*qz));
	  if (ibin <= nbin/2) counts[ibin] += 1; 
	}
      }
    }

}
  
  
/* ---------------------------------------------------------------------- */

void ComputeStructureFactorRadial::compute_array()
{

  compute_histogram();

  
  int histbin,ibin,m;
  double sqr,sqi;
  double qx,qy,qz;


  for (m = 0; m < npairs; m++) 
    for (ibin = 0; ibin < size_array_rows; ibin++) 
      array[ibin][2+m] = 0;
  
  if (domain->dimension == 2) {

    
    for (m = 0; m < npairs; m++) {

    
      for (int nx = -nbin/2; nx < nbin/2; nx++ )
	for (int ny = 0; ny < nbin/2+1; ny++ ) {
	  qx = delta_q[0]*nx;
	  qy = delta_q[1]*ny;
	  ibin = q_to_bin(sqrt(qx*qx+qy*qy));
	  if (ibin > nbin/2) continue;
	  histbin = ny  + (nbin/2+1)*(nx+nbin/2);
	  sqr = (hist[2*m][histbin]);
	  sqi = (hist[1+2*m][histbin]);	
	  array[ibin][2+m] += sqr*sqr+sqi*sqi;
	}
    }
    
    
  } else {
    
    for (m = 0; m < npairs; m++) {
      
      for (int nx = -nbin/2; nx < nbin/2; nx++ )
	for (int ny = -nbin/2; ny < nbin/2; ny++ )
	  for (int nz = 0; nz < nbin/2+1; nz++ ) {
	    qx = delta_q[0]*nx;
	    qy = delta_q[1]*ny;
	    qz = delta_q[2]*nz;
	    ibin = q_to_bin(sqrt(qx*qx+qy*qy+qz*qz));
	    if (ibin > nbin/2) continue;
	    histbin = nz  + (nbin/2+1)*((ny+nbin/2) + nbin*(nx+nbin/2));
	    sqr = (hist[2*m][histbin]);
	    sqi = (hist[1+2*m][histbin]);	
	    array[ibin][2+m] += sqr*sqr+sqi*sqi;	  
	  }
    }
  }
  
  for (ibin = 0; ibin < size_array_rows; ibin++) 
    array[ibin][1] = counts[ibin];
  
  for (m = 0; m < npairs; m++) 
    for (ibin = 0; ibin < size_array_rows; ibin++) 
      array[ibin][2+m] /= counts[ibin];

}


int ComputeStructureFactorRadial::q_to_bin(double q)
{
  return static_cast<int> ((2*q/sqrt_dq2+1)/2);
}
