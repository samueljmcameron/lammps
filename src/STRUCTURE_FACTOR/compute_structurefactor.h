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
ComputeStyle(structurefactor,ComputeStructureFactor);
// clang-format on
#else

#ifndef LMP_COMPUTE_STRUCTURE_FACTOR_H
#define LMP_COMPUTE_STRUCTURE_FACTOR_H

#include "compute_structurefactor_base.h"

namespace LAMMPS_NS {

class ComputeStructureFactor : public ComputeStructureFactorBase {
 public:
  ComputeStructureFactor(class LAMMPS *, int, char **);
  ~ComputeStructureFactor() override;
  void compute_array() override;
  
 protected:
  void write_qs() override;
  
};

}    // namespace LAMMPS_NS

#endif
#endif
