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
ComputeStyle(structure/sincos/single,ComputeStructureSinCosSingleQ);
// clang-format on
#else


#ifndef LMP_COMPUTE_STRUCTURE_SINCOS_SINGLE_Q_H
#define LMP_COMPUTE_STRUCTURE_SINCOS_SINGLE_Q_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeStructureSinCosSingleQ : public Compute {
 public:
  ComputeStructureSinCosSingleQ(class LAMMPS *, int, char **);
  ~ComputeStructureSinCosSingleQ() override;
  void init() override;
  void compute_vector() override;
  
 protected:

  void init_norm();
  int single_bin;        // which q bin to calculate values at
  int q_counts;
  double delta_q;    // bin width and its inverse
  double binned_q;   // the average value of q within the spherical/circular shell
  double cosval,sinval; // local values of cos(qdotr) and sin(qdotr)
  
};

}    // namespace LAMMPS_NS

#endif
#endif
