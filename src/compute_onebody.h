/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef COMPUTE_CLASS

ComputeStyle(onebody,ComputeOneBody)

#else

#ifndef LMP_COMPUTE_ONEBODY_H
#define LMP_COMPUTE_ONEBODY_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeOneBody : public Compute {
 public:
  ComputeOneBody(class LAMMPS *, int, char **);
  ~ComputeOneBody();
  void init();
  void compute_array();

 private:
  int nbin;             // # of bins (distance bins)
  int hist_rows;        // # of different rho_i(r) (i = 0 to atom->ntypes)
  int cutflag;          // user cutoff flag
  double delr,delrinv;  // bin width and its inverse for distance
  double radius;        // radius of container
  double mycutneigh;    // user-specified cutoff + neighbor skin
  double **hist;        // histogram bins
  double **histall;     // summed histogram bins across all procs

  int *typecount;
  int *icount;


  void init_norm();
  bigint natoms_old;

  void set_array(double, double,int);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Compute rdf requires a pair style be defined or cutoff specified

UNDOCUMENTED

E: Compure rdf cutoff exceeds ghost atom range - use comm_modify cutoff command

UNDOCUMENTED

W: Compute rdf cutoff less than neighbor cutoff - forcing a needless neighbor list build

UNDOCUMENTED

U: Compute rdf requires a pair style be defined

Self-explanatory.

*/
