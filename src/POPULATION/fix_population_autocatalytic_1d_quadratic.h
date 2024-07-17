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
FixStyle(population/autocatalytic/1D/quadratic,FixPopulationAutocatalytic1dQuadratic);
// clang-format on
#else

#ifndef LMP_FIX_POPULATION_AUTOCATALYTIC_1D_QUADRATIC_H
#define LMP_FIX_POPULATION_AUTOCATALYTIC_1D_QUADRATIC_H

#include "fix_population_autocatalytic_1d_projection.h"

namespace LAMMPS_NS {

class FixPopulationAutocatalytic1dQuadratic : public FixPopulationAutocatalytic1dProjection {
 public:
  FixPopulationAutocatalytic1dQuadratic(class LAMMPS *, int, char **);

  virtual double force(double) override;
  
 private:

  double springconstant,minloc;
  
};

}    // namespace LAMMPS_NS

#endif
#endif
