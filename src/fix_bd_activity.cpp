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
#include "fix_bd_activity.h"
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

FixBDActivity::FixBDActivity(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  if (narg != 7) error->all(FLERR,"Illegal fix bd/euler command");

  if (domain->dimension != 2) error->all(FLERR,"Fix bd/euler requires 2d simulation");

  t_start = force->numeric(FLERR,arg[3]);
  t_target = t_start;
  t_stop = force->numeric(FLERR,arg[4]);
  diff = force->numeric(FLERR,arg[5]);
  if (diff <= 0.0) error->all(FLERR,"Fix bd/euler diffusion coefficient must be > 0.0");
  seed = force->inumeric(FLERR,arg[6]);
  if (seed <= 0) error->all(FLERR,"Illegal fix bd/euler command");

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp,seed + comm->me);
}

/* ---------------------------------------------------------------------- */

FixBDActivity::~FixBDActivity()
{

  delete random;

}


/* ---------------------------------------------------------------------- */

void FixBDActivity::init()
{

  if (!atom->mu_flag)
    error->all(FLERR,"Fix bd/euler requires atom attributes mu");

  //TBD: check if muz is zero for all particles


  // set square root of temperature
  compute_target();

  gamma1 = diff / force->boltz; 
  gamma2 = sqrt( 24 * diff );
  gamma3 = sqrt( 24 * 3 * diff );
  gamma4 = 3 * diff / force->boltz; 

  FixNVE::init();
}

/* ----------------------------------------------------------------------
   set current t_target and t_sqrt
------------------------------------------------------------------------- */

void FixBDActivity::compute_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // Only homogeneous temperature supported
  t_target = t_start + delta * (t_stop-t_start);
  tsqrt = sqrt(t_target);

}


/* ---------------------------------------------------------------------- */

void FixBDActivity::initial_integrate(int vflag)
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


