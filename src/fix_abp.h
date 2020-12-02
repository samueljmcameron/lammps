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

FixStyle(abp,FixABP)

#else

#ifndef LMP_FIX_ABP_H
#define LMP_FIX_ABP_H

namespace LAMMPS_NS {

class FixABP : public Fix {
 public:
  FixABP(class LAMMPS *, int, char **);
  virtual ~FixABP();
  void init();
  void post_force(int);
  void setup(int);
  int setmask();

 private:
  double selfpropulsionforce;

};

}
#endif
#endif

/* ERROR/WARNING messages:

E: Fix abp requires atom attribute mu

Self-explanatory.

*/
