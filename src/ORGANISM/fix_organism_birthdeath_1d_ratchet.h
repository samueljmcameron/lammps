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
FixStyle(organism/birthdeath/1D/ratchet,FixOrganismBirthDeath1dRatchet);
// clang-format on
#else

#ifndef LMP_FIX_ORGANISM_BIRTHDEATH_1D_RATCHET_H
#define LMP_FIX_ORGANISM_BIRTHDEATH_1D_RATCHET_H

#include "fix_organism_birthdeath_1d_projection.h"

namespace LAMMPS_NS {

class FixOrganismBirthDeath1dRatchet : public FixOrganismBirthDeath1dProjection {
 public:
  FixOrganismBirthDeath1dRatchet(class LAMMPS *, int, char **);
  
  virtual double force(double) override;
  
 private:

  bool first_section(double );

  double height, first_length, second_length, midvertex;
};

}    // namespace LAMMPS_NS

#endif
#endif
