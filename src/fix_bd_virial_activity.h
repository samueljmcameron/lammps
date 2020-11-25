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

FixStyle(bd/virial/activity,FixBdVirialActivity)

#else

#ifndef LMP_FIX_BD_VIRIAl_ACTIVITY_H
#define LMP_FIX_BD_VIRIAl_ACTIVITY_H

#include "fix_nve.h"

namespace LAMMPS_NS {

class FixBdVirialActivity : public Fix {
 public:
  FixBdVirialActivity(class LAMMPS *, int, char **);
  virtual ~FixBdVirialActivity();
  void init();
  void initial_integrate(int);
  int setmask();

 private:
  double dt, sqrtdt;
  double gamma1,gamma2, gamma3, gamma4;
  double cosda, sinda, da, dar;
  int seed;
  int unwrapcoords_flag;
  double diff_t;
  double diff_r;
  double activity;

 protected:

  class RanMars *random;


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
