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


#ifndef LMP_FIX_POPULATION_BASE_H
#define LMP_FIX_POPULATION_BASE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPopulationBase : public Fix {
 public:
  FixPopulationBase(class LAMMPS *, int, char **);
  ~FixPopulationBase() override;
  void reset_dt() override;
  int setmask() override;
  void init() override;
  void post_integrate() override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  double compute_vector(int) override;
  void grow_arrays(int ) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;
  
 private:
  bool recycle_from_dead(int);
  void create_new_atoms(const std::vector<int> &);
  void delete_dead_atoms();
  void count_vitals();

  virtual void compute_division_and_death_rates() = 0;
  virtual void divide(int,int);

  int delete_flag;
  int recycle_flag;
protected:
  enum {NEVER, EVERY, EFFICIENT};
  class RanMars *rng; // used to decide for each atom whether it divides, dies, or does neither
  int alivetype,deadtype; // alivetype == type of atom which we consider alive (similar for deadtype)
  int seed; // rng seed

  int *alive_indices;   // list of indices for alive particles (local to comm->me)
  int *dead_indices;    // list of indices for dead particles (local to comm->me)

  int *nalive_per_proc; // number of alive particles on each processor
  int *ndead_per_proc; // number of dead particles on each processor

  int nalive,ndead;  // total number of alive and dead particles globally
  int localalive,localdead; // totall number of alive and dead particles locally
  int cleanevery; // how often to remove dead particles
  int nspecified_args; // count arguments specified

  double shift; // the total displacement distance between daughter atoms after a division
  double **divdeath_array; // rates of division (first column) and death (second column)

  int total_atoms_from_scratch, total_atoms_recycled, total_atoms_killed, total_atoms_deleted;

};

}    // namespace LAMMPS_NS

#endif
