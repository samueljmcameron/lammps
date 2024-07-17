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

#include "fix_population_autocatalytic_1d_projection.h"

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

FixPopulationAutocatalytic1dProjection::FixPopulationAutocatalytic1dProjection(LAMMPS *lmp, int narg, char **arg) :
  FixPopulationAutocatalytic(lmp, narg, arg)
{
  if (strcmp(arg[nspecified_args++], "omega") != 0)
    error->all(FLERR, "Need to specify omega in fix population/autocatalytic/1D command");
  omega = utils::numeric(FLERR, arg[nspecified_args++], false, lmp);
  
}



/* ---------------------------------------------------------------------- */

int FixPopulationAutocatalytic1dProjection::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}


void FixPopulationAutocatalytic1dProjection::post_force(int /* vflag */)
{
  double **x = atom->x;
  double **f = atom->f;

  for (int i = 0; i < atom->nlocal; i++) {

    if (atom->type[i] == alivetype) {
      f[i][0] += force(x[i][0]);
    }
  }
  
}



void FixPopulationAutocatalytic1dProjection::divide(int i, int j)
{

  double **x = atom->x;
  double **v = atom->v;

  double cx,cy,cz;

  cx = x[i][0];
  cy = x[i][1];
  cz = x[i][2];

  double sigma_x,sigma_hi,sigma_lo;

  
  sigma_x = shift/sqrt(1+force(cx)*force(cx)/(omega*omega));

  sigma_hi = 0.5*sigma_x;
  sigma_lo = 0.5*sigma_x;


  atom->type[j] = alivetype;
	    
  x[j][0] = cx+sigma_hi;
  x[j][1] = cy;
  x[j][2] = cz;
  
  x[i][0] = cx-sigma_lo;
  x[i][1] = cy;
  x[i][2] = cz;
  
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  return;
}
