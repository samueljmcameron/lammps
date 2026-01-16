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
FixStyle(population/local/neighbor/death,FixPopulationLocalNeighborDeath);
// clang-format on
#else

#ifndef LMP_FIX_POPULATION_LOCAL_NEIGHBOR_DEATH_H
#define LMP_FIX_POPULATION_LOCAL_NEIGHBOR_DEATH_H

#include "fix_population_base.h"
#include <random>

namespace LAMMPS_NS {

class FixPopulationLocalNeighborDeath : public FixPopulationBase {
 public:
  FixPopulationLocalNeighborDeath(class LAMMPS *, int, char **);
  ~FixPopulationLocalNeighborDeath() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  
  int pack_reverse_comm(int, int , double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

 private:

  virtual void compute_division_and_death_rates() override;
  void allocate();
  void coeff(char **);


  
  int **setflag;      // 0/1 = whether each i,j has been set
  double **cg_volume; // coarse-grained volume for local density for i-j types
  double **b0;        // birthrate for i-j types
  double **d0;        // deathrate for i-j types
  double **sigma;     // effective diameter for i-j types
  double **width;     // width of sigmoid function for i-j types

  int cutflag;
  double cutoff_user;
  double mycutneigh;       // user-specified cutoff + neighbor skin

  
  const int ncoeff;
  
  
protected:
  class NeighList *list;    // neighbor list
  bool full_neigh_list;

};

}    // namespace LAMMPS_NS

#endif
#endif
