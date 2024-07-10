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
FixStyle(birthdeath/simple,FixBirthDeathSimple);
// clang-format on
#else

#ifndef LMP_FIX_BIRTHDEATH_SIMPLE_H
#define LMP_FIX_BIRTHDEATH_SIMPLE_H

#include "fix.h"
#include <random>

namespace LAMMPS_NS {

class FixBirthDeathSimple : public Fix {
 public:
  FixBirthDeathSimple(class LAMMPS *, int, char **);
  ~FixBirthDeathSimple() override;
  void reset_dt() override;
  int setmask() override;
  void init() override;
  void post_integrate() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double compute_vector(int) override;

 private:
  bool birth_from_dead(int);
  void create_new_atoms(const std::vector<int> &);
  void delete_dead_atoms();
  void count_vitals();
  
  virtual void procreate(int,int);

  std::mt19937 gen;

  std::uniform_int_distribution<> proc_dist;
  std::uniform_int_distribution<> atom_dist;

  bool exactpoisson;
  
protected:
  class RanMars *rng;
  int alivetype,deadtype;
  int seed;

  double shift;             // radial amount to shift the parent and daughter cells
  int *local_alive_list;
  int *local_dead_list;

  int *global_alive_list;
  int *global_dead_list;

  int atom_swap_nmax;
  int nalive,ndead;
  int localalive,localdead;
  int cleanevery;
  int nspecified_args;

};

}    // namespace LAMMPS_NS

#endif
#endif
