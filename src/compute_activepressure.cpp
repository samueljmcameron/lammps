/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "compute_pressure.h"
#include <mpi.h>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "modify.h"
#include "fix.h"
#include "force.h"
#include "pair.h"
#include "pair_hybrid.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeActivePressure::ComputeActivePressure(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  vptr(NULL), id_temp(NULL)
{
  if (narg < 3) error->all(FLERR,"Illegal compute activepressure command");
  if (igroup) error->all(FLERR,"Compute activepressure must use group all");

  scalar_flag = 1;
  extscalar = 0;
  timeflag = 1;

  // process optional args

  if (narg == 3) {
    unwrapcoordsflag = 1;
    fixflag = 1;
  } else if (narg ==4) {

    if (strcmp(arg[iarg],"wrapcoords") == 0) {
      unwrapcoordsflag=0;
    } else {
      error->all(FLERR,"Illegal compute activepressure command");
    }

  } else {
    error->all(FLERR,"Illegal compute activepressure command");
  }

  nvirial = 0;
  vptr = NULL;
}

/* ---------------------------------------------------------------------- */

ComputeActivePressure::~ComputeActivePressure()
{
  delete [] vptr;
}

/* ---------------------------------------------------------------------- */

void ComputeActivePressure::init()
{
  boltz = force->boltz;
  nktv2p = force->nktv2p;
  dimension = domain->dimension;

  // detect contributions to virial
  // vptr points to all virial[6] contributions

  delete [] vptr;
  nvirial = 0;
  vptr = NULL;

  if (fixflag)
    for (int i = 0; i < modify->nfix; i++)
      if (modify->fix[i]->thermo_virial) nvirial++;

  if (nvirial) {
    vptr = new double*[nvirial];
    nvirial = 0;
    if (fixflag)
      for (int i = 0; i < modify->nfix; i++)
        if (modify->fix[i]->thermo_virial)
          vptr[nvirial++] = modify->fix[i]->virial;
  }

}

/* ----------------------------------------------------------------------
   compute total pressure, averaged over Pxx, Pyy, Pzz
------------------------------------------------------------------------- */

double ComputeActivePressure::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (update->vflag_global != invoked_scalar)
    error->all(FLERR,"Virial was not tallied on needed timestep");


  if (dimension == 3) {
    inv_volume = 1.0 / (domain->xprd * domain->yprd * domain->zprd);
    virial_compute(3,3);
    scalar = (virial[0] + virial[1] + virial[2]) / 3.0 * inv_volume * nktv2p;
  } else {
    inv_volume = 1.0 / (domain->xprd * domain->yprd);
    virial_compute(2,2);
    scalar = (virial[0] + virial[1]) / 2.0 * inv_volume * nktv2p;
  }

  return scalar;
}


/* ---------------------------------------------------------------------- */

void ComputeActivePressure::virial_compute(int n, int ndiag)
{
  int i,j;
  double v[6],*vcomponent;

  for (i = 0; i < n; i++) v[i] = 0.0;

  // sum contributions to virial from forces and fixes

  for (j = 0; j < nvirial; j++) {
    vcomponent = vptr[j];
    for (i = 0; i < n; i++) v[i] += vcomponent[i];
  }

  // sum virial across procs
  MPI_Allreduce(v,virial,n,MPI_DOUBLE,MPI_SUM,world);

}

/* ---------------------------------------------------------------------- */

void ComputeActivePressure::reset_extra_compute_fix(const char *id_new)
{
  delete [] id_temp;
  int n = strlen(id_new) + 1;
  id_temp = new char[n];
  strcpy(id_temp,id_new);
}
