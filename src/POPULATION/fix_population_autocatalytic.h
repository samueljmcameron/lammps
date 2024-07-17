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
FixStyle(population/autocatalytic,FixPopulationAutocatalytic);
// clang-format on
#else

#ifndef LMP_FIX_POPULATION_AUTOCATALYTIC_H
#define LMP_FIX_POPULATION_AUTOCATALYTIC_H

#include "fix_population_base.h"

namespace LAMMPS_NS {

class FixPopulationAutocatalytic : public FixPopulationBase {
 public:
  FixPopulationAutocatalytic(class LAMMPS *, int, char **);
  void reset_dt() override;
  

 private:

  virtual void compute_division_and_death_rates() override;

  
protected:
  double divisionrate; // rate of division per cell
  double deathrate; // rate of death per cell divided by number of cells minus 1;

};

}    // namespace LAMMPS_NS

#endif
#endif
