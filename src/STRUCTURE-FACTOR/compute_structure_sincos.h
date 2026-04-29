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
ComputeStyle(structure/sincos,ComputeStructureSinCos);
// clang-format on
#else


#ifndef LMP_COMPUTE_STRUCTURE_SINCOS_H
#define LMP_COMPUTE_STRUCTURE_SINCOS_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeStructureSinCos : public Compute {
 public:
  ComputeStructureSinCos(class LAMMPS *, int, char **);
  ~ComputeStructureSinCos() override;
  void init() override;
  void compute_array() override;
  
 protected:

  void init_norm();
  int nq_bin;                // # of q bins
  double delta_q;    // bin width and its inverse
  double **hist;       // histogram bins
  double *q_counts;     // number qs in a qbin
  double **histall;    // summed histogram bins across all procs

};

}    // namespace LAMMPS_NS

#endif
#endif
