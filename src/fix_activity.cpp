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
   Contributing authors: Carolyn Phillips (U Mich), reservoir energy tally
                         Aidan Thompson (SNL) GJF formulation
                         Charles Sievers & Niels Gronbech-Jensen (UC Davis)
                             updated GJF formulation and included
                             statistically correct 2GJ velocity
------------------------------------------------------------------------- */

#include "fix_activity.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "respa.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;


/* ---------------------------------------------------------------------- */

FixActivity::FixActivity(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), gfactor1(NULL), gfactor2(NULL),
  random(NULL)
{
  if (narg < 8) error->all(FLERR,"Illegal fix activity command");

  global_freq = 1;

  momentum_flag = 0;
  
  fp = force->numeric(FLERR,arg[3])
  Dt = force->numeric(FLERR,arg[4]);
  gamma = force->numeric(FLERR,arg[5]);
  Dr = force->numeric(FLERR,arg[6]);
  seed = force->inumeric(FLERR,arg[7]);

  if (gamma <= 0.0) error->all(FLERR,"Fix activity gamma must be > 0.0");
  if (seed <= 0) error->all(FLERR,"Illegal fix activity command");



  if (narg == 9) {
    if (strcomp(arg[narg-1],"momentum_dynamics") == 0) {
      momentum_flag = 1;
    } else {
      error->all(FLERR,"Illegal fix activity command");
    }
  } else if (narg > 9) {
    error->all(FLERR,"Illegal fix activity command");
  }
  
  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

  // allocate per-type arrays for force prefactors

  gfactor1 = new double[atom->ntypes+1];
  gfactor2 = new double[atom->ntypes+1];


}

/* ---------------------------------------------------------------------- */

FixActivity::~FixActivity()
{
  delete random;
  delete [] gfactor1;
  delete [] gfactor2;

}

/* ---------------------------------------------------------------------- */

int FixActivity::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixActivity::init()
{

  if (!atom->mu_flag)
    error->all(FLERR,"Fix activity requires atom attributes mu");


  // check that all group particles are finite-size

  
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;


  // set force prefactors

  if (!atom->rmass) {
    for (int i = 1; i <= atom->ntypes; i++) {
      gfactor1[i] = -atom->mass[i] / t_period / force->ftm2v;
      gfactor2[i] = sqrt(atom->mass[i]) *
	sqrt(24.0*force->boltz/t_period/update->dt/force->mvv2e) /
	force->ftm2v;
    }
  }

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

}

/* ---------------------------------------------------------------------- */

void FixActivity::setup(int vflag)
{

  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }

}

/* ---------------------------------------------------------------------- */

void FixActivity::post_force(int /*vflag*/)
{
  double *rmass = atom->rmass;

  if (rmass)
    post_force_templated<1>();
  else
    post_force_templated<0>();
  
}

/* ---------------------------------------------------------------------- */

void FixActivity::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ----------------------------------------------------------------------
   modify forces
------------------------------------------------------------------------- */

template < int Tp_RMASS >
void FixActivity::post_force_templated()
{
  double gamma1,gamma2;

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // apply damping and thermostat to atoms in group

  // for Tp_RMASS:
  //   use per-atom masses
  //   else use per-type masses

  double fdrag[3],fran[3],fsum[3],fsumall[3];
  bigint count;
  double fswap;

  double boltz = force->boltz;
  double dt = update->dt;
  double mvv2e = force->mvv2e;
  double ftm2v = force->ftm2v;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (Tp_RMASS) {
        gamma1 = -rmass[i] / t_period / ftm2v;
	gamma2 = sqrt(rmass[i]) * sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
        gamma2 *= tsqrt;
      } else {
        gamma1 = gfactor1[type[i]];
        gamma2 = gfactor2[type[i]] * tsqrt;
      }

      
      fran[0] = gamma2*(random->uniform()-0.5);
      fran[1] = gamma2*(random->uniform()-0.5);
      fran[2] = gamma2*(random->uniform()-0.5);



      fdrag[0] = gamma1*v[i][0];
      fdrag[1] = gamma1*v[i][1];
      fdrag[2] = gamma1*v[i][2];



      f[i][0] += fdrag[0] + fran[0];
      f[i][1] += fdrag[1] + fran[1];
      f[i][2] += fdrag[2] + fran[2];

    }
  }


  // thermostat rotational degrees of freedom

  omega_thermostat();
}



/* ----------------------------------------------------------------------
   thermostat rotational dof via omega
------------------------------------------------------------------------- */

void FixActivity::omega_thermostat()
{
  double gamma1,gamma2;

  double boltz = force->boltz;
  double dt = update->dt;
  double mvv2e = force->mvv2e;
  double ftm2v = force->ftm2v;

  double **torque = atom->torque;
  double **omega = atom->omega;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  // rescale gamma1/gamma2 by 10/3 & sqrt(10/3) for spherical particles
  // does not affect rotational thermosatting
  // gives correct rotational diffusivity behavior

  double tendivthree = 10.0/3.0;
  double tran[3];
  double inertiaone;

  for (int i = 0; i < nlocal; i++) {
    if ((mask[i] & groupbit) && (radius[i] > 0.0)) {
      inertiaone = SINERTIA*radius[i]*radius[i]*rmass[i];
      gamma1 = -tendivthree*inertiaone / t_period / ftm2v;
      gamma2 = sqrt(inertiaone) * sqrt(80.0*boltz/t_period/dt/mvv2e) / ftm2v;
      gamma2 *= tsqrt;
      tran[0] = gamma2*(random->uniform()-0.5);
      tran[1] = gamma2*(random->uniform()-0.5);
      tran[2] = gamma2*(random->uniform()-0.5);
      torque[i][0] += gamma1*omega[i][0] + tran[0];
      torque[i][1] += gamma1*omega[i][1] + tran[1];
      torque[i][2] += gamma1*omega[i][2] + tran[2];
    }
  }
}


void FixActivity::reset_dt()
{
  if (atom->mass) {
    for (int i = 1; i <= atom->ntypes; i++) {
      gfactor2[i] = sqrt(atom->mass[i]) *
	sqrt(24.0*force->boltz/t_period/update->dt/force->mvv2e) /
	force->ftm2v;
      gfactor2[i] *= 1.0;
    }
  }
}


