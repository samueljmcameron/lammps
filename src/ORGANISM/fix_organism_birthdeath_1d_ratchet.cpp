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

#include "fix_organism_birthdeath_1d_ratchet.h"

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

FixOrganismBirthDeath1dRatchet::FixOrganismBirthDeath1dRatchet(LAMMPS *lmp, int narg, char **arg) :
  FixOrganismBirthDeath1dProjection(lmp, narg, arg)
{
  // required args
  
  double fraction = utils::numeric(FLERR, arg[10], false, lmp);

  if (fraction <= 0 || fraction >= 1) 
    error->all(FLERR, "Invalid ratchet fraction in organism/birthdeath/1D/ratchet");
  
  height = utils::numeric(FLERR, arg[11], false, lmp);


  first_length = fraction*domain->prd[0];
  second_length = (1-fraction)*domain->prd[0];
  midvertex = domain->boxlo[0]+first_length;
  
}




double FixOrganismBirthDeath1dRatchet::force(double x) {
  
  if (first_section(x)) {
    
    return height/first_length;
    
  } else if (std::abs(x-midvertex) < 1e-14) {
    return 0.0;
    
  } else {
    
    return -height/second_length;
    
  }

}



bool FixOrganismBirthDeath1dRatchet::first_section(double x)
{

  double period = domain->prd[0];

  // get x into periodic position (since PBC aren't applied every step)
  //if (x < domain->boxlo[0]-second_length)
  //  printf("x is far to the left of PBC\n");

  while (x < domain->boxlo[0])
    x += period;

  //if (x >= domain->boxhi[0] + first_length)
  //  printf("x is far to the right of PBC\n");

  while (x >= domain->boxhi[0])
    x -= period;

  return ((domain->boxlo[0] <= x) && (x < midvertex)) || x > domain->boxhi[0];
}
