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
  
  if (narg != 9 && narg != 10)
    error->all(FLERR,"Illegal fix bd/activity command");

  if (domain->dimension != 2)
    error->all(FLERR,"Fix bd/activity requires 2d simulation");

  t_start = force->numeric(FLERR,arg[3]);
  t_target = t_start;
  t_stop = force->numeric(FLERR,arg[4]);
  diff_t = force->numeric(FLERR,arg[5]);
  if (diff_t <= 0.0)
    error->all(FLERR,"Fix bd/activity translational diffusion "
	       "coefficient must be > 0.");

  diff_r = force->numeric(FLERR,arg[6]);
  if (diff_r <= 0.0)
    error->all(FLERR,"Fix bd/activity rotational diffusion "
	       "coefficient must be > 0.");
  
  activity = force->numeric(FLERR,arg[7]);
  if (activity>0) {
    activity_flag = 1;
  } else if (activity <0) {
    error->all(FLERR,"Illegal fix bd/activity command");
  } else {
    activity_flag = 0;
  }
  seed = force->inumeric(FLERR,arg[8]);
  if (seed <= 0) error->all(FLERR,"Fix bd/activity seed must be > 0.");

  if (narg == 10) {
    if (strcmp(arg[9],"unwrapcoords") == 0) {
      unwrapcoords_flag = 1;
    } else {
      error->all(FLERR,"Illegal fix bd/activity command");
    }
  }

  force_flag = 0;

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);


}

/* ---------------------------------------------------------------------- */

int FixBdVirialActivity::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
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

  //TBD: check if muz is zero for all particles


  // set square root of temperature
  compute_target();

  gamma1 = diff_t / force->boltz; 
  gamma2 = sqrt( 24 * diff_t );
  gamma3 = sqrt( 24 * diff_r );
  gamma4 = diff_r / force->boltz; 

  //  dtv = update->dt;
  //  dtf = 0.5 * update->dt * force->ftm2v;

  //if (strstr(update->integrate_style,"respa"))
  //  step_respa = ((Respa *) update->integrate)->step;

}

/* ----------------------------------------------------------------------
   set current t_target and t_sqrt
------------------------------------------------------------------------- */

void FixBdVirialActivity::compute_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // Only homogeneous temperature supported
  t_target = t_start + delta * (t_stop-t_start);
  tsqrt = sqrt(t_target);

}

/* ---------------------------------------------------------------------- */

void FixBdVirialActivity::initial_integrate(int /* vflag */)
{
  double **x = atom->x;
  double **v = atom->v;
  double **mu = atom->mu;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;


  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA
  dt = update->dt;
  sqrtdt = sqrt(dt);

  // set square root of temperature
  compute_target();


  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      da = (dt * gamma1 * f[i][0] / t_target
	    +    sqrtdt * gamma2 * (random->uniform()-0.5));
      double dx = da;
      x[i][0] +=  da;
      v[i][0]  =  da/dt;
      da = (dt * gamma1 * f[i][1] / t_target
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


/* ----------------------------------------------------------------------
   apply active force, stolen from MISC/fix_efield.cpp 
------------------------------------------------------------------------- */

void FixBdVirialActivity::pre_force(int vflag)
{
  double **f = atom->f;
  double *q = atom->q;
  int *mask = atom->mask;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  double fsum[4];

  // energy and virial setup

  if (vflag) v_setup(vflag);
  else evflag = 0;


  // fsum[0] = "potential energy" for added force
  // fsum[123] = extra force added to atoms

  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;
  force_flag = 0;

  double **x = atom->x;
  double **mu = atom->mu;
  double fx,fy;
  double v[6];

  // constant activity parameter

  double tmp_x[3];

  // charge interactions
  // force = activity*mu, potential energy = F dot x in unwrapped coords

  if (activity_flag) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {

	fx = activity*mu[i][0];
	fy = activity*mu[i][1];
	f[i][0] += fx;
	f[i][1] += fy;

	if (unwrapcoords_flag) {
	  domain->unmap(x[i],image[i],tmp_x);
	} else {
	  tmp_x[0] = x[i][0];
	  tmp_x[1] = x[i][1];
	  tmp_x[2] = x[i][2];
	}
	fsum[0] -= fx*tmp_x[0]+fy*tmp_x[1];
	fsum[1] += fx;
	fsum[2] += fy;
	if (evflag) {
	  v[0] = fx*tmp_x[0];
	  v[1] = fy*tmp_x[1];
	  v[2] = 0;
	  v[3] = fx*tmp_x[1];
	  v[4] = 0;
	  v[5] = 0;
	  v_tally(i, v);
	}
      }
  }
}
