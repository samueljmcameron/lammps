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
   Originally modified from USER-CGDNA/fix_nve_dotc_langevin.cpp. 

   Contributing author: Majid Moseyabi (University of Bristol)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "fix_bd_virial_activity.h"
#include "math_extra.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "domain.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBdVirialActivity::FixBdVirialActivity(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  thermo_virial = 1;
  virial_flag = 1;
  time_integrate = 1;
  
  if (narg != 7 && narg != 8)
    error->all(FLERR,"Illegal fix bd/activity command");

  if (domain->dimension != 2)
    error->all(FLERR,"Fix bd/activity requires 2d simulation");

  diff_t = force->numeric(FLERR,arg[3]);
  if (diff_t <= 0.0)
    error->all(FLERR,"Fix bd/activity translational diffusion "
	       "coefficient must be > 0.");

  diff_r = force->numeric(FLERR,arg[4]);
  if (diff_r <= 0.0)
    error->all(FLERR,"Fix bd/activity rotational diffusion "
	       "coefficient must be > 0.");
  
  activity = force->numeric(FLERR,arg[5]);
  if (activity <0) {
    error->all(FLERR,"Illegal fix bd/activity command");
  }
  
  seed = force->inumeric(FLERR,arg[6]);
  if (seed <= 0) error->all(FLERR,"Fix bd/activity seed must be > 0.");

  unwrapcoords_flag = 0;
  if (narg == 8) {
    if (strcmp(arg[7],"unwrapcoords") == 0) {
      unwrapcoords_flag = 1;
    } else {
      error->all(FLERR,"Illegal fix bd/activity command");
    }
  }

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);


}

/* ---------------------------------------------------------------------- */

int FixBdVirialActivity::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

FixBdVirialActivity::~FixBdVirialActivity()
{

  delete random;

}


/* ---------------------------------------------------------------------- */

void FixBdVirialActivity::init()
{

  if (!atom->mu_flag)
    error->all(FLERR,"Fix bd/euler requires atom attributes mu");


  gamma1 = diff_t / force->boltz; 
  gamma2 = sqrt( 24 * diff_t );
  gamma3 = sqrt( 24 * diff_r );
  gamma4 = diff_r / force->boltz; 


}


/* ---------------------------------------------------------------------- */

void FixBdVirialActivity::initial_integrate(int vflag)
{
  double **x = atom->x;
  double **v = atom->v;
  double **mu = atom->mu;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double fx,fy;

  // for virial only
  double vi[6];
  double tmp_x[3];

  imageint *image = atom->image;
  int nlocal = atom->nlocal;

  // energy and virial setup

  if (vflag) v_setup(vflag);
  else evflag = 0;



  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed
  dt = update->dt;
  sqrtdt = sqrt(dt);


  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      // set active forces and add them to the forces
      // force = activity*mu
      fx = activity*mu[i][0];
      fy = activity*mu[i][1];
      f[i][0] += fx;
      f[i][1] += fy;

      // compute virial contribution mu_i \cdot x_i
      if (evflag) {
	if (unwrapcoords_flag) {
	  domain->unmap(x[i],image[i],tmp_x);
	} else {
	  tmp_x[0] = x[i][0];
	  tmp_x[1] = x[i][1];
	  tmp_x[2] = x[i][2];
	}

	vi[0] = fx*tmp_x[0];
	vi[1] = fy*tmp_x[1];
	vi[2] = 0;
	vi[3] = fx*tmp_x[1];
	vi[4] = 0;
	vi[5] = 0;
	v_tally(i, vi);
      }
      
      // update positions and orientations
      
      da = (dt * gamma1 * f[i][0] 
	    +    sqrtdt * gamma2 * (random->uniform()-0.5));
      double dx = da;
      x[i][0] +=  da;
      v[i][0]  =  da/dt;
      da = (dt * gamma1 * f[i][1] 
	    +    sqrtdt * gamma2 * (random->uniform()-0.5));
      x[i][1] +=  da;
      double dy = da;
      v[i][1]  =  da/dt;


      da   =   sqrtdt * gamma3 * (random->uniform()-0.5);
      cosda = cos(da);
      sinda = sin(da);


      da = mu[i][0];
      mu[i][0] =  mu[i][0]*cosda - mu[i][1]*sinda;
      mu[i][1] =  sinda * da     + mu[i][1]*cosda;

    }

}
