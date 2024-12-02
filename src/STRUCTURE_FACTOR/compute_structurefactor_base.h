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


#ifndef LMP_COMPUTE_STRUCTURE_FACTOR_BASE_H
#define LMP_COMPUTE_STRUCTURE_FACTOR_BASE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeStructureFactorBase : public Compute {
 public:
  ComputeStructureFactorBase(class LAMMPS *, int, char **);
  ~ComputeStructureFactorBase() override;
  void init() override;
  void compute_array() override;
  
 protected:
  virtual void write_qs() = 0;
  int nbin;                // # of sq bins
  int npairs;              // # of sq pairs
  double *delta_q;    // bin width and its inverse
  int ***sqpair;          // map 2 type pair to sq pair for each histo
  int **nsqpair;          // # of histograms for each type pair
  int *ilo, *ihi, *jlo, *jhi;
  double **hist;       // histogram bins
  double **histall;    // summed histogram bins across all procs

  int hist_bins;
  int *typecount;
  int *icount, *jcount;
  int *duplicates;

  static void callback(int, char *, void *);

  void init_norm();
  bigint natoms_old;
};

}    // namespace LAMMPS_NS

#endif
