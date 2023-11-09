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
FixStyle(organism/birthdeath/logistic,FixOrganismBirthDeathLogistic);
// clang-format on
#else

#ifndef LMP_FIX_ORGANISM_BIRTHDEATH_LOGISTIC_H
#define LMP_FIX_ORGANISM_BIRTHDEATH_LOGISTIC_H

#include "fix.h"

namespace LAMMPS_NS {

class FixOrganismBirthDeathLogistic : public Fix {
 public:
  FixOrganismBirthDeathLogistic(class LAMMPS *, int, char **);
  ~FixOrganismBirthDeathLogistic() override;
  int setmask() override;
  void init() override;
  void post_integrate() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:
  bool birth_from_dead(int);
  void create_new_atoms(const std::vector<int> &);
  void delete_dead_atoms();
  void count_vitals();
  
  virtual void procreate(int,int);
  
protected:
  class RanMars *rng;
  int alivetype,deadtype;
  int seed;
  
  double birthrate,deathrate;

  double shift;             // radial amount to shift the parent and daughter cells
  int *local_alive_list;
  int *local_dead_list;

  int atom_swap_nmax;
  int nalive,ndead;
  int localalive,localdead;
  int cleanevery;

};

}    // namespace LAMMPS_NS

#endif
#endif
