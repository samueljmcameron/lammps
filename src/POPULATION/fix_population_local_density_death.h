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
FixStyle(population/local/density/death,FixPopulationLocalDensityDeath);
// clang-format on
#else

#ifndef LMP_FIX_POPULATION_LOCAL_DENSITY_DEATH_H
#define LMP_FIX_POPULATION_LOCAL_DENSITY_DEATH_H

#include "fix_population_base.h"
#include <random>

namespace LAMMPS_NS {

class FixPopulationLocalDensityDeath : public FixPopulationBase {
 public:
  FixPopulationLocalDensityDeath(class LAMMPS *, int, char **);
  ~FixPopulationLocalDensityDeath() override;
  

 private:

  virtual void compute_division_and_death_rates() override;


  
  double b0; 
  double d0; 


  double xbinsize,ybinsize,zbinsize;
  double binvolume;


  char *id_density;
  class ComputeChunkAtom *density;

  std::vector<double> count_local, count_global;
  
  
protected:
  class NeighList *list;    // neighbor list
  bool full_neigh_list;

};

}    // namespace LAMMPS_NS

#endif
#endif
