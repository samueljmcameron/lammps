/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(janus/cut,PairJanusCut);
// clang-format on
#else

#ifndef LMP_PAIR_JANUS_CUT_H
#define LMP_PAIR_JANUS_CUT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairJanusCut : public Pair {
 public:
  PairJanusCut(class LAMMPS *);
  ~PairJanusCut() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  void *extract(const char *, int &) override;

 protected:
  class AtomVecEllipsoid *avec;

  double cut_global;
  double **cut, **cut_sq;

  bool fixed_orientations;
  bool self_interacting;

  double *p1_body_frame, *p2_body_frame, *connector_body_frame;


private:

  double *Fdummy, *Fcumulative, *T_i_dummy, *T_j_dummy, *cross_dummy, *p_cross_r;
  double *rij, *dummy_connector, *dummy_distance;
  double *p1_i, *p1_j, *p2_i, *p2_j, *connector_i, *connector_j;
  
  void Fpp(const double *, const double *, const double *,
	   double *); // Xichen's 2.4

  void Tpp(const double *, const double *, const double *,
	   double *); // Xichen's 2.7
;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
