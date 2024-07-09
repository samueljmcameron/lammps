/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include "pair_birthdeath.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairBirthDeath::PairBirthDeath(LAMMPS *lmp) : Pair(lmp)
{

}

/* ---------------------------------------------------------------------- */

PairBirthDeath::~PairBirthDeath()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(b0);
    memory->destroy(d0);
    memory->destroy(sigma);
    memory->destroy(width);
    memory->destroy(cnum);
  }
}

/* ---------------------------------------------------------------------- */

void PairBirthDeath::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double r_dist,rsq;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *division = atom->division;
  double *death = atom->death;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (i = 0; i < atom->nlocal; i++) 
    division[i] = b0[type[i]][type[i]];
  
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];


    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      r_dist = sqrt(rsq);
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

	division[i] -= b0[itype][jtype]/(1+exp((r_dist-sigma[itype][jtype])/width[itype][jtype]))/cnum[itype][jtype];
	death[i] += d0[itype][jtype]/(1+exp((r_dist-sigma[itype][jtype])/width[itype][jtype]))/cnum[itype][jtype];
	
        if (newton_pair || j < nlocal) {
	  division[j] -= b0[jtype][itype]/(1+exp((r_dist-sigma[jtype][itype])/width[jtype][itype]))/cnum[jtype][itype];
	  death[j] += d0[jtype][itype]/(1+exp((r_dist-sigma[jtype][itype])/width[jtype][itype]))/cnum[jtype][itype];
        }

      }
    }

  }

}


/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBirthDeath::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n, n, "pair:cutsq");

  memory->create(cut, n, n, "pair:cut");
  memory->create(d0, n, n, "pair:d0");
  memory->create(b0, n, n, "pair:b0");
  memory->create(sigma, n, n, "pair:sigma");
  memory->create(width, n, n, "pair:width");
  memory->create(cnum, n, n, "pair:cnum");

}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBirthDeath::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style command");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);
  

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBirthDeath::coeff(int narg, char **arg)
{
  if (narg < 7 || narg > 8) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double b0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double d0_one = utils::numeric(FLERR, arg[3], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[4], false, lmp);
  double width_one = utils::numeric(FLERR, arg[5], false, lmp);
  double cnum_one = utils::numeric(FLERR, arg[6], false, lmp);

  double cut_one = cut_global;

  if (narg == 8) cut_one = utils::numeric(FLERR, arg[7], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      b0[i][j] = b0_one;
      d0[i][j] = d0_one;
      sigma[i][j] = sigma_one;
      width[i][j] = width_one;
      cut[i][j] = cut_one;
      cnum[i][j] = cnum_one;
      setflag[i][j] = 1;
      count++;
    }
  }


  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */
double PairBirthDeath::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    b0[i][j] = mix_distance(b0[i][i],b0[j][j]);
    d0[i][j] = mix_distance(d0[i][i],d0[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    width[i][j] = mix_distance(width[i][i],width[j][j]);
    cnum[i][j] = mix_distance(cnum[i][i],cnum[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  double cutout = cut[i][j];
  return cutout;
}



/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBirthDeath::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&b0[i][j], sizeof(double), 1, fp);
	fwrite(&d0[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
	fwrite(&width[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
	fwrite(&cnum[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBirthDeath::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &b0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &d0[i][j], sizeof(double), 1, fp, nullptr, error);
	  utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
	  utils::sfread(FLERR, &width[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&b0[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&d0[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
	MPI_Bcast(&width[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cnum[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBirthDeath::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBirthDeath::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairBirthDeath::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d %g %g %g %g %g\n",
						  i, b0[i][i],d0[i][i], sigma[i][i],
						  width[i][i],cnum[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairBirthDeath::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g %g %g\n", i, j, b0[i][j], d0[i][j],sigma[i][j],
	      width[i][j],cnum[i][j],cut[i][j]);
}

