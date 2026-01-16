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

#include "fix_population_sensing_stress_divide.h"

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
#include "math_eigen.h"

#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstring>


#include<unistd.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPopulationSensingStressDivide::FixPopulationSensingStressDivide(LAMMPS *lmp, int narg, char **arg) :
  FixPopulationSensing(lmp, narg, arg)
{

  if (force->pair->single_enable != 1)
    error->all(FLERR, "all pair forces must have single_enable to use population/sensing/stressdivide command");    


  minimum_stress_flag = false;
  memory->create(stress,6,"POP_stressdivide:stress");
  memory->create(eigvals,3,"POP_stressdivide:eigvals");
  memory->create(eigvecs,3,3,"POP_stressdivide:eigvecs");
  memory->create(stress_matrix,3,3,"POP_stressdivide:stress_matrix");
}

FixPopulationSensingStressDivide::~FixPopulationSensingStressDivide()
{
  memory->destroy(stress);
  memory->destroy(eigvals);
  memory->destroy(eigvecs);
  memory->destroy(stress_matrix);
}


void FixPopulationSensingStressDivide::init()
{
  FixPopulationSensing::init();
  // need full neighbor list to calculate local stress on whichever
  //   dividing atom needs it
  neighbor->add_request(this,NeighConst::REQ_FULL);
  full_neigh_list = true;

}
void FixPopulationSensingStressDivide::get_division_direction(int i)
{

  // copied below from compute_stress_atom.cpp

  int pairflag = 1;
  int bondflag = 0;
  int angleflag = 0;
  int dihedralflag = 0;
  int improperflag = 0;
  int kspaceflag = 0;
  int fixflag = 0;
  
    // clear local stress array

  int *type = atom->type;
  double **x = atom->x;


  int itype = type[i];

  int *jlist = list->firstneigh[i];
  int jnum = list->numneigh[i];
  int j;

  double xtmp = x[i][0];
  double ytmp = x[i][1];
  double ztmp = x[i][2];
  
  int jtype;
  double delx,dely,delz,rsq;
  double factor_coul,factor_lj;
  double fforce;
  double v[6];

  zero_stress_flag = 1;
  
  for (int index = 0; index < 6; index++) stress[index] = 0.0;

  // add in per-atom contributions from each force


  for (int jj = 0; jj < jnum; jj++) {
    j = jlist[jj];

    factor_coul = force->special_coul[sbmask(j)];
    factor_lj = force->special_lj[sbmask(j)];

    
    j &= NEIGHMASK;
    
    
    jtype = type[j];
    
    delx = xtmp - x[j][0];
    dely = ytmp - x[j][1];
    delz = ztmp - x[j][2];


    rsq = delx*delx + dely*dely + delz*delz;

    if (rsq < force->pair->cutsq[itype][jtype]) {

      if (pairflag && force->pair && force->pair->compute_flag) {

	zero_stress_flag = 0;
	force->pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,
			    fforce);

	if (rsq == 0.0) {
	  fforce = 0.0; // if particles overlap and the potential is very soft, then this is fine...
	}
	if (update->ntimestep == 15449 || update->ntimestep == 15450) {
	  printf("force = %e\n",fforce);
	  printf("factor_lj = %e\n",factor_lj);
	  printf("itype = %d\n",itype);
	  printf("jtype = %d\n",jtype);
	}
	
	
	v[0] = delx*delx*fforce;
	v[1] = dely*dely*fforce;
	v[2] = delz*delz*fforce;
	v[3] = delx*dely*fforce;
	v[4] = delx*delz*fforce;
	v[5] = dely*delz*fforce;

	for (int index = 0; index < 6; index++)
	  stress[index] += 0.5*v[index];
      }
    }

    // Check whether the total stress is nonzero or not.
    // This is necessary for two reasons:
    //  1) in case all forces perfectly balance (unlikely)
    //  2) pair interactions are type dependent, so two atoms
    //     might be within range but if they are of types which
    //     have interactions turned off (which can be the case if
    //     if specified in the input script e.g. using pair_coeff)
    //     the force will be zero even if atoms are close
    bool anystress = false;
    for (int index = 0; index < 6; index++)
      if (std::abs(stress[index])>0.0)
	anystress = true;

    if (!anystress)
      zero_stress_flag = 1;
      

    
  }

  /*
    if (bondflag && force->bond) {

	force->bond->single(btype,rsq,itype,jtype,rsq,force->special_lj[force->pair->sbmask(j)],fforce);
	v[0] = delx*delx*fforce;
	v[1] = dely*dely*fforce;
	v[2] = delz*delz*fforce;
	v[3] = delx*dely*fforce;
	v[4] = delx*delz*fforce;
	v[5] = dely*delz*fforce;

	
	for (int index = 0; index < 6; index++)
	  stress[index] += 0.5*v[index];
      }

  if (angleflag && force->angle) {
    double **vatom = force->angle->vatom;
    for (i = 0; i < nbond; i++)
      for (j = 0; j < 6; j++) stress[i][j] += vatom[i][j];
  }

  if (dihedralflag && force->dihedral) {
    double **vatom = force->dihedral->vatom;
    for (i = 0; i < nbond; i++)
      for (j = 0; j < 6; j++) stress[i][j] += vatom[i][j];
  }

  if (improperflag && force->improper) {
    double **vatom = force->improper->vatom;
    for (i = 0; i < nbond; i++)
      for (j = 0; j < 6; j++) stress[i][j] += vatom[i][j];
  }

  if (kspaceflag && force->kspace && force->kspace->compute_flag) {
    double **vatom = force->kspace->vatom;
    for (i = 0; i < nkspace; i++)
      for (j = 0; j < 6; j++) stress[i][j] += vatom[i][j];
  }

  // add in per-atom contributions from relevant fixes
  // skip if vatom = nullptr
  // possible during setup phase if fix has not initialized its vatom yet
  // e.g. fix ave/spatial defined before fix shake,
  //   and fix ave/spatial uses a per-atom stress from this compute as input

  if (fixflag) {
    for (auto &ifix : modify->get_fix_list())
      if (ifix->virial_peratom_flag && ifix->thermo_virial) {
        double **vatom = ifix->vatom;
        if (vatom)
          for (i = 0; i < nlocal; i++)
            for (j = 0; j < 6; j++) stress[i][j] += vatom[i][j];
      }
  }

      
    }



    
  */

}


void FixPopulationSensingStressDivide::divide(int i, int j)
{

  get_division_direction(i);
  
  double **x = atom->x;
  double **v = atom->v;

  
  double polar,azim;
  double cx,cy,cz,dx,dy,dz;

  cx = x[i][0];
  cy = x[i][1];
  cz = x[i][2];


  if (zero_stress_flag) {
  
    if (domain->dimension == 3) {
      polar = M_PI*rng->uniform();
      azim = 2*M_PI*rng->uniform();
      dx = shift*cos(azim)*sin(polar)/2.;
      dy = shift*sin(azim)*sin(polar)/2.;
      dz = shift*cos(polar)/2;
      
    } else {
      azim = 2*M_PI*rng->uniform();
      dx = shift*cos(azim)/2.;
      dy = shift*sin(azim)/2.;

      dz = 0.0;
    }
  } else {

    if (domain->dimension == 3) {

      stress_matrix[0][0] = stress[0];
      stress_matrix[1][1] = stress[1];
      stress_matrix[2][2] = stress[2];
      stress_matrix[0][1] = stress[3];
      stress_matrix[1][0] = stress[3];
      stress_matrix[0][2] = stress[4];
      stress_matrix[2][0] = stress[4];
      stress_matrix[1][2] = stress[5];
      stress_matrix[2][1] = stress[5];
      
      MathEigen::jacobi3(stress_matrix,eigvals,eigvecs);


    } else {



      // set eigenvector matrix z components to zero
      eigvecs[0][2] = 0;
      eigvecs[1][2] = 0;
      eigvecs[2][2] = 0;
      eigvecs[2][0] = 0;
      eigvecs[2][1] = 0;

      // compute 2D eigenvalues by hand

      double norm;
      double b = - stress[0] - stress[1]; // trace of 2D stress matrix
      double c = stress[0]*stress[1] - stress[3]*stress[3]; // determinant of 2D stress matrix


      printf("stress[0] = %lf, stress[1] = %lf, stress[3] = %lf\n,",
	     stress[0],stress[1],stress[3]);
      printf("b = %e, c = %e\n,",b,c);
      // quadratic formula to solve for eigenvalues
      eigvals[0] = (-b + sqrt(b*b-4*c))/2.;
      eigvals[1] = (-b - sqrt(b*b-4*c))/2.;

      // calculate eigenvectors
      for (int index = 0; index < 2; index++) {
	if (stress[0] - eigvals[index] == 0) {
	  eigvecs[index][0] = 1-index;
	  eigvecs[index][1] = index;
	}
	else {
	  norm
	    = sqrt(stress[3]*stress[3]
		   /((stress[0]-eigvals[index])*(stress[0]-eigvals[index]))
		   +1);
	  
	  eigvecs[index][0] = -stress[3]/(stress[0]-eigvals[index])/norm;
	  eigvecs[index][1] = 1/norm;
	}
      }


    }


    int arg_eig;
    if (minimum_stress_flag)
      arg_eig = get_min_index(eigvals);
    else
      arg_eig = get_max_index(eigvals);
    printf("eigenvector[%d] = %lf, %lf, %lf with eigenvalue %lf on step %lld for atom at position %lf,%lf,%lf\n",
	   arg_eig,eigvecs[arg_eig][0],
	   eigvecs[arg_eig][1],eigvecs[arg_eig][2],
	   eigvals[arg_eig],
	   update->ntimestep,x[i][0],x[i][1],x[i][2]);

    if (arg_eig == 0)
      printf("OTHER eigenvector[%d] = %lf, %lf, %lf with eigenvalue %e on step %lld for atom at position %lf,%lf,%lf\n",
	     1,eigvecs[1][0],
	     eigvecs[1][1],eigvecs[1][2],
	     eigvals[1],
	     update->ntimestep,x[i][0],x[i][1],x[i][2]);

    else if (arg_eig == 1)
      printf("OTHER eigenvector[%d] = %lf, %lf, %lf with eigenvalue %e on step %lld for atom at position %lf,%lf,%lf\n",
	   0,eigvecs[0][0],
	   eigvecs[0][1],eigvecs[0][2],
	   eigvals[0],
	   update->ntimestep,x[i][0],x[i][1],x[i][2]);
    
    
    dx = shift*eigvecs[arg_eig][0]/2.;
    dy = shift*eigvecs[arg_eig][1]/2.;
    dz = shift*eigvecs[arg_eig][2]/2.;
    
  }
  
  atom->type[j] = alivetype;
  
  x[j][0] = cx+dx;
  x[j][1] = cy+dy;
  x[j][2] = cz+dz;
  
  x[i][0] = cx-dx;
  x[i][1] = cy-dy;
  x[i][2] = cz-dz;
  
  v[j][0] = v[i][0]/2.;
  v[j][1] = v[i][1]/2.;
  v[j][2] = v[i][2]/2.;
  
  v[i][0] /= 2;
  v[i][1] /= 2;
  v[i][2] /= 2;
    
  divdeath_array[j][0] = divdeath_array[i][0];
  divdeath_array[j][1] = divdeath_array[i][1];
  return;
}


int FixPopulationSensingStressDivide::get_max_index(double *eigval)
{
  double max = std::abs(eigval[0]);
  int arg_eig = 0;
  for (int index = 1; index < domain->dimension; index++)
    if (std::abs(eigval[index])>=max) {
      max = std::abs(eigval[index]);
      arg_eig = index;
    }
  return arg_eig;
}


int FixPopulationSensingStressDivide::get_min_index(double *eigval)
{
  double min = std::abs(eigval[0]);
  int arg_eig = 0;
  for (int index = 1; index < domain->dimension; index++)
    if (std::abs(eigval[index]) <= min) {
      min = std::abs(eigval[index]);
      arg_eig = index;
    }
  return arg_eig;
}

