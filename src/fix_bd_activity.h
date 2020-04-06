/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(bd/activity,FixBDActivity)

#else

#ifndef LMP_FIX_BD_ACTIVITY_H
#define LMP_FIX_BD_ACTIVITY_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixBDActivity : public Fix {
 public:
  FixBDActivity(class LAMMPS *, int, char **);
  virtual ~FixBDActivity();
  void init();
  void initial_integrate(int);
  //void final_integrate();
  void pre_force(int);
  int setmask();
  void setup(int vflag);
 private:
  double dt, sqrtdt;
 protected:
  int seed;
  double t_start,t_stop,t_target,tsqrt;
  double diff;
  double activity;
  double gamma1,gamma2, gamma3, gamma4;
  double cosda, sinda, da, dar;
  class RanMars *random;
  void compute_target();

  int force_flag;
  int activity_flag;
};

}
#endif
#endif

/* ERROR/WARNING messages:

E: Compute bd/euler requires atom style ellipsoid

Self-explanatory.

E: Fix bd/euler requires extended particles

This fix can only be used for particles with a shape setting.

*/
