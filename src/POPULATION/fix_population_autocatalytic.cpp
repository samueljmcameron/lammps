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

#include "fix_population_autocatalytic.h"

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

FixPopulationAutocatalytic::FixPopulationAutocatalytic(LAMMPS *lmp, int narg, char **arg) :
  FixPopulationBase(lmp, narg, arg)
{

  divisionrate = utils::numeric(FLERR, arg[nspecified_args++], false, lmp);
  
  if (divisionrate < 0)
    error->all(FLERR, "divisionrate must be positive in fix population/autocatalytic command");

  deathrate = utils::numeric(FLERR, arg[nspecified_args++], false, lmp);
  
  if (deathrate < 0)
    error->all(FLERR, "deathrate must be positive in fix population/autocatalytic command");

  shift = utils::numeric(FLERR, arg[nspecified_args++], false, lmp);
  
  if (shift < 0)
    error->all(FLERR, "shift must be positive in fix population/autocatalytic command");

  
}

void FixPopulationAutocatalytic::reset_dt()
{
  
  double dt = update->dt;

  double Nav = divisionrate/deathrate;
  
  if (Nav*divisionrate*dt > 1.0)
    error->all(FLERR,"Birth rate in fix population/autocatalytic command is too large, multiple events will occur for given dt.");    
  else if (Nav*dt*(deathrate*(Nav-1) + divisionrate) > 1.0) 
    error->all(FLERR,"Death rate in fix population/autocatalytic command is too large, multiple events will occur for given dt.");

}

 
/* ---------------------------------------------------------------------- */


void FixPopulationAutocatalytic::compute_division_and_death_rates()
{
  // set divdeath to be zero 
  for (int i = 0; i < atom->nlocal; i++) {
    divdeath_array[i][0] = divisionrate;
    divdeath_array[i][1] = deathrate*(nalive-1);
  }
}

