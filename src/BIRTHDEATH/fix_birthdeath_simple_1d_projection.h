/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */


#ifndef LMP_FIX_BIRTHDEATH_SIMPLE_1D_PROJECTION_H
#define LMP_FIX_BIRTHDEATH_SIMPLE_1D_PROJECTION_H

#include "fix_birthdeath_simple.h"

namespace LAMMPS_NS {

class FixBirthDeathSimple1dProjection : public FixBirthDeathSimple {
 public:
  FixBirthDeathSimple1dProjection(class LAMMPS *, int, char **);
  int setmask() override;
  void post_force(int) override;
  virtual double force(double) = 0;
 protected:
  double omega;
 private:
  virtual void procreate(int,int) override;

};

}    // namespace LAMMPS_NS

#endif
