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

#ifdef FIX_CLASS
// clang-format off
FixStyle(organism/birthdeath/logistic/ratchet1d,FixOrganismBirthDeathLogisticRatchet1D);
// clang-format on
#else

#ifndef LMP_FIX_ORGANISM_BIRTHDEATH_LOGISTIC_RATCHET1D_H
#define LMP_FIX_ORGANISM_BIRTHDEATH_LOGISTIC_RATCHET1D_H

#include "fix_organism_birthdeath_logistic.h"

namespace LAMMPS_NS {

class FixOrganismBirthDeathLogisticRatchet1D : public FixOrganismBirthDeathLogistic {
 public:
  FixOrganismBirthDeathLogisticRatchet1D(class LAMMPS *, int, char **);
  int setmask() override;
  void post_force(int);

 private:
  virtual void procreate(int,int) override;
  bool first_section(double );
protected:

  double height, first_length, second_length, midvertex;
};

}    // namespace LAMMPS_NS

#endif
#endif
