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

#include "fix_population_local_density_death.h"

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
#include "compute_chunk_atom.h"
#include "pair.h"
#include "update.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstring>

#include<unistd.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPopulationLocalDensityDeath::FixPopulationLocalDensityDeath(LAMMPS *lmp, int narg, char **arg) :
  FixPopulationBase(lmp, narg, arg)
{


  shift = utils::numeric(FLERR, arg[nspecified_args++], false, lmp);
  
  if (shift < 0)
    error->all(FLERR, "shift must be positive in fix population/sensing command");

  xbinsize = utils::numeric(FLERR, arg[nspecified_args++], false, lmp);
  if (xbinsize < 0)
    error->all(FLERR, "xbinsize must be positive integer in fix population/sensing command");
  
  ybinsize = utils::numeric(FLERR, arg[nspecified_args++], false, lmp);
  if (ybinsize < 0)
    error->all(FLERR, "ybinsize must be positive integer in fix population/sensing command");
  
  if (domain->dimension == 3) {
    zbinsize = utils::inumeric(FLERR, arg[nspecified_args++], false, lmp);
    if (zbinsize < 0)
      error->all(FLERR, "zbinsize must be positive integer in fix population/sensing command");
  }
  
  std::string alive_group = std::string(arg[nspecified_args++]);

  b0 = utils::numeric(FLERR, arg[nspecified_args++], false, lmp);
  d0 = utils::numeric(FLERR, arg[nspecified_args++], false, lmp);
  
  id_density = utils::strdup(std::string(id) + "_density");
  if (domain->dimension==2)
    density = dynamic_cast<ComputeChunkAtom *>
      (modify->add_compute(fmt::format("{} {} chunk/atom bin/2d x lower {} y lower {} units box",
				       id_density,alive_group,xbinsize,ybinsize)));
  else
    density = dynamic_cast<ComputeChunkAtom *>
      (modify->add_compute(fmt::format("{} {} chunk/atom bin/3d x lower {} y lower {} z lower {} units box",
				       id_density,alive_group,xbinsize,ybinsize,zbinsize)));


  if (domain->dimension == 2)
    binvolume = xbinsize*ybinsize;
  else
    binvolume = xbinsize*ybinsize*zbinsize;
  
}


FixPopulationLocalDensityDeath::~FixPopulationLocalDensityDeath()
{
  density = nullptr;
}


/* ---------------------------------------------------------------------- */


void FixPopulationLocalDensityDeath::compute_division_and_death_rates()
{



  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double local_density;

  int nchunk;

  // Compute density, store in 
  nchunk = density->setup_chunks();
  density->compute_ichunk();

  int *ichunk = density->ichunk;


  count_local.resize(nchunk);
  count_global.resize(nchunk);

  for (int m = 0; m < nchunk; m++)  {
    count_local[m] = 0.0;
    count_global[m] = 0.0;
  }
  
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && ichunk[i] > 0)
      count_local[ichunk[i]-1]++;


  MPI_Allreduce(count_local.data(),count_global.data(),nchunk,MPI_DOUBLE,MPI_SUM,world);
  
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit & type[i] == alivetype) {
      local_density = count_global[ichunk[i]-1]/binvolume;
      //printf("count of atom %d at ( %lf, %lf ) is %lf, binvolume = %lf\n",
      //	     atom->tag[i],x[i][0],x[i][1],count_global[ichunk[i]-1],binvolume);
      divdeath_array[i][0] = b0;
      divdeath_array[i][1] = d0*local_density;

    }
}
