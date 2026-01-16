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
FixStyle(population/sensing/stressdivide,FixPopulationSensingStressDivide);
// clang-format on
#else

#ifndef LMP_FIX_POPULATION_SENSING_STRESS_DIVIDE_H
#define LMP_FIX_POPULATION_SENSING_STRESS_DIVIDE_H

#include "fix_population_sensing.h"
#include <random>

namespace LAMMPS_NS {

class FixPopulationSensingStressDivide : public FixPopulationSensing {
 public:
  FixPopulationSensingStressDivide(class LAMMPS *, int, char **);
  ~FixPopulationSensingStressDivide() override;
  void init() override;

 private:

  virtual void divide(int,int) override;
  void get_division_direction(int);

  int get_max_index(double *);
  int get_min_index(double *);
  
  bool minimum_stress_flag;
  int zero_stress_flag;

  double *stress; // components of symmetric local stress tensor of dividing particle 
  double **stress_matrix;     // 3x3 local stress tensor (filled in version of stress)
  double *eigvals; // eigenvalues of local stress tensor
  double **eigvecs; // eigenvectors of local stress tensor;
  
  
protected:
  inline int sbmask(int j) const { return j >> SBBITS & 3; }

};

}    // namespace LAMMPS_NS

#endif
#endif
