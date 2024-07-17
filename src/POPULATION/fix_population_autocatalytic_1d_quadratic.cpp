/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Sam Cameron
------------------------------------------------------------------------- */

#include "fix_population_autocatalytic_1d_quadratic.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "random_mars.h"
#include "pair.h"
#include "update.h"
#include "modify.h"

#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstring>

#include<unistd.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPopulationAutocatalytic1dQuadratic::FixPopulationAutocatalytic1dQuadratic(LAMMPS *lmp, int narg, char **arg) :
  FixPopulationAutocatalytic1dProjection(lmp, narg, arg)
{

  // required args

  if (strcmp(arg[nspecified_args++], "springconstant") != 0)
    error->all(FLERR, "Need to specify springconstant in fix population/autocatalytic/1D/quadratic command");
  springconstant = utils::numeric(FLERR, arg[nspecified_args++], false, lmp);
  if (strcmp(arg[nspecified_args++], "minloc") != 0)
    error->all(FLERR, "Need to specify springconstant in fix population/autocatalytic/1D/quadratic command");
  minloc = utils::numeric(FLERR, arg[nspecified_args++], false, lmp);
  
}


double FixPopulationAutocatalytic1dQuadratic::force(double x)
{

  return -springconstant*(x-minloc);

}
