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

#include "fix_organism_birthdeath_1d_smoothratchet.h"

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

FixOrganismBirthDeath1dSmoothRatchet::FixOrganismBirthDeath1dSmoothRatchet(LAMMPS *lmp, int narg, char **arg) :
  FixOrganismBirthDeath1dProjection(lmp, narg, arg)
{

  // required args
  
  n_pot = utils::inumeric(FLERR, arg[10], false, lmp);
  
  height = utils::numeric(FLERR, arg[11], false, lmp);


  length = domain->prd[0];


  yend = endpoint();
  normfac = gfunc(yend);
  
}


double FixOrganismBirthDeath1dSmoothRatchet::force(double x)
{

  double y = 2*M_PI/length*x + yend -M_PI;
  return -height*M_PI/length*gfunc_deriv(y)/normfac;

}



double FixOrganismBirthDeath1dSmoothRatchet::binom(int n, int k)
{
  if (n == k) return 1.0;
  if (k == 0) return 1.0;
  double output = 1.0;

  for (int i = n; i > n -k; i--)
    output *= i;

  for (int i = 1; i <= k ; i++)
    output /= i;
  
  return output;
}

double FixOrganismBirthDeath1dSmoothRatchet::gfunc(double x)
{

  int n = n_pot;
  double out = 0;
  for (int k = 1; k <= n; k++) 
    out += -binom(2*n,n-k)/(binom(2*n,n)*k)*sin(k*x);

  return out;

}

double FixOrganismBirthDeath1dSmoothRatchet::gfunc_deriv(double x)
{
  int n = n_pot;
  return -pow(2,(2*n-1))/binom(2*n,n)*pow(cos(x/2.),2*n)+0.5;

}

double FixOrganismBirthDeath1dSmoothRatchet::endpoint() {

  int n = n_pot;
  return 2*acos(0.5*(pow(binom(2*n,n),(1./(2*n)))));
  
}
