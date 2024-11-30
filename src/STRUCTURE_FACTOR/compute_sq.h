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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(sq,ComputeSQ);
// clang-format on
#else

#ifndef LMP_COMPUTE_SQ_H
#define LMP_COMPUTE_SQ_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeSQ : public Compute {
 public:
  ComputeSQ(class LAMMPS *, int, char **);
  ~ComputeSQ() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_array() override;

 private:
  int nbin;                // # of sq bins
  int cutflag;             // user cutoff flag
  int npairs;              // # of sq pairs
  double *delta_q;    // bin width and its inverse
  double cutoff_user;      // user-specified cutoff
  double mycutneigh;       // user-specified cutoff + neighbor skin
  int ***sqpair;          // map 2 type pair to sq pair for each histo
  int **nsqpair;          // # of histograms for each type pair
  int *ilo, *ihi, *jlo, *jhi;
  double **hist;       // histogram bins
  double **histall;    // summed histogram bins across all procs

  int *typecount;
  int *icount, *jcount;
  int *duplicates;

  static void callback(int, char *, void *);

  class NeighList *list;    // half neighbor list
  void init_norm();
  bigint natoms_old;
};

}    // namespace LAMMPS_NS

#endif
#endif
