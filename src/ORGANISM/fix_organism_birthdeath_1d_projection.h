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


#ifndef LMP_FIX_ORGANISM_BIRTHDEATH_1D_PROJECTION_H
#define LMP_FIX_ORGANISM_BIRTHDEATH_1D_PROJECTION_H

#include "fix_organism_birthdeath_logistic.h"

namespace LAMMPS_NS {

class FixOrganismBirthDeath1dProjection : public FixOrganismBirthDeathLogistic {
 public:
  FixOrganismBirthDeath1dProjection(class LAMMPS *, int, char **);
  int setmask() override;
  void post_force(int) override;
  virtual double force(double) = 0;

 private:
  virtual void procreate(int,int) override;

};

}    // namespace LAMMPS_NS

#endif
