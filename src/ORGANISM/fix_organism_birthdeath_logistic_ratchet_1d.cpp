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

#include "fix_organism_birthdeath_logistic_ratchet_1d.h"

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

FixOrganismBirthDeathLogisticRatchet1D::FixOrganismBirthDeathLogisticRatchet1D(LAMMPS *lmp, int narg, char **arg) :
  FixOrganismBirthDeathLogistic(lmp, narg, arg)
{

  // required args
  
  double fraction = utils::numeric(FLERR, arg[10], false, lmp);

  if (fraction <= 0 || fraction >= 1) 
    error->all(FLERR, "Invalid ratchet fraction in organism/birthdeath/ratchet1d");
  
  height = utils::numeric(FLERR, arg[11], false, lmp);

  first_length = fraction*domain->prd[0];
  second_length = (1-fraction)*domain->prd[0];
  midvertex = domain->boxlo[0]+first_length;
}


/* ---------------------------------------------------------------------- */

int FixOrganismBirthDeathLogisticRatchet1D::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

void FixOrganismBirthDeathLogisticRatchet1D::post_force(int /* vflag */)
{
  double **x = atom->x;
  double **f = atom->f;

  for (int i = 0; i < atom->nlocal; i++) {

    if (atom->type[i] == alivetype) {
      
      if (first_section(x[i][0])) {
	f[i][0] += height/first_length;
      
      } else {

	f[i][0] += -height/second_length;
      }
      
    }
  }
  
}



bool FixOrganismBirthDeathLogisticRatchet1D::first_section(double x)
{
  
  return ((domain->boxlo[0] <= x) && (x < midvertex)) || x > domain->boxhi[0];
}


void FixOrganismBirthDeathLogisticRatchet1D::procreate(int i, int j)
{

  double **x = atom->x;
  double **v = atom->v;

  
  double polar,azim;
  double cx,cy,cz;

  cx = x[i][0];
  cy = x[i][1];
  cz = x[i][2];

  double sigma_x;
  double slope;

  if (first_section(cx)) {
    slope = -height/first_length;
  } else {
    slope = height/second_length;
  }

  
  sigma_x = shift/sqrt(1+slope*slope);

  

  atom->type[j] = alivetype;
	    
  x[j][0] = cx+sigma_x/2.;
  x[j][1] = cy;
  x[j][2] = cz;
  
  x[i][0] = cx-sigma_x/2.;
  x[i][1] = cy;
  x[i][2] = cz;
  
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  return;
}
