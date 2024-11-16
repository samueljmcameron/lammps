// clang-format off
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

#include "pair_janus_cut.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "atom_vec_ellipsoid.h"
#include "math_extra.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairJanusCut::PairJanusCut(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;

}

/* ---------------------------------------------------------------------- */

PairJanusCut::~PairJanusCut()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(cut_sq);
    memory->destroy(p1_body_frame);
    memory->destroy(p2_body_frame);
    memory->destroy(connector_body_frame);


    memory->destroy(Fdummy);
    memory->destroy(Fcumulative);
    memory->destroy(T_i_dummy);
    memory->destroy(T_j_dummy);
    memory->destroy(cross_dummy);
    memory->destroy(p_cross_r);
    memory->destroy(rij);
    memory->destroy(dummy_connector);
    memory->destroy(dummy_distance);

    memory->destroy(connector_i);
    memory->destroy(connector_j);
    memory->destroy(p1_i);
    memory->destroy(p1_j);
    memory->destroy(p2_i);
    memory->destroy(p2_j);    
  }
}

/* ---------------------------------------------------------------------- */

void PairJanusCut::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul;
  double rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double **torque = atom->torque;
  double *quat;
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;

  double Q[3][3];

  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  
  
  
  
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    quat = bonus[ellipsoid[i]].quat;
    MathExtra::quat_to_mat(quat, Q);

    // connector_i is u_i in Xichen's notation
    MathExtra::matvec(Q, connector_body_frame, connector_i);


    if (!fixed_orientations) {
      MathExtra::matvec(Q, p1_body_frame, p1_i);
      MathExtra::matvec(Q, p2_body_frame, p2_i);
    } else {
      for (int dumi = 0; dumi < 3; dumi++) {
	p1_i[dumi] = p1_body_frame[dumi];
	p2_i[dumi] = p2_body_frame[dumi];
      }
    }

    itype = type[i];

    if (self_interacting) {

      for (int dim = 0; dim < 3; dim++)
	T_i_dummy[dim] = 0.0;
      
      
      // first term of equation 2.24
      
      for (int dim = 0; dim < 3; dim++)
	dummy_connector[dim] = -2*connector_i[dim];
      
      // compute the term from equation 2.24
      Tpp(p1_i,p2_i,dummy_connector,T_i_dummy,itype,itype);
      
      
      // second term of equation 2.24
      for (int dim = 0; dim < 3; dim++)
	dummy_connector[dim] = 2*connector_i[dim];
      
      // compute the term from equation 2.24
      Tpp(p2_i,p1_i,dummy_connector,T_i_dummy,itype,itype);
      
      
      // third term of equation 2.24
      Fpp(p1_i,p2_i,dummy_connector,Fdummy,itype,itype);
      
      cross3_addition(dummy_connector,Fdummy,T_i_dummy);
      
      for (int dim = 0; dim < 3; dim++)
	torque[i][dim] += T_i_dummy[dim];
    }
      
    //itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;
      
      rij[0] = xtmp - x[j][0];
      rij[1] = ytmp - x[j][1];
      rij[2] = ztmp - x[j][2];


      quat = bonus[ellipsoid[j]].quat;
      MathExtra::quat_to_mat(quat, Q);
      MathExtra::matvec(Q, connector_body_frame, connector_j);

      if (!fixed_orientations) {
	MathExtra::matvec(Q, p1_body_frame, p1_j);
	MathExtra::matvec(Q, p2_body_frame, p2_j);
      } else {
	for (int dumi = 0; dumi < 3; dumi++) {
	  p1_j[dumi] = p1_body_frame[dumi];
	  p2_j[dumi] = p2_body_frame[dumi];
	}
      }
      
      rsq = MathExtra::dot3(rij,rij);
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        // atom can have both a charge and dipole
        // i,j = charge-charge, dipole-dipole, dipole-charge, or charge-dipole


	for (int dim = 0; dim < 3; dim++) {
	  Fcumulative[dim] = 0.0;
	  T_i_dummy[dim] = 0.0;
	  T_j_dummy[dim] = 0.0;
	}
	
	// laboriously sum over all the four pairs in Xichen's equation 2.5:

	// first pair:
	for (int dim = 0; dim < 3; dim++)
	  dummy_distance[dim] = (connector_j[dim] - connector_i[dim])  + rij[dim];

	// compute the term from equation 2.5
	Fpp(p1_i,p1_j,dummy_distance,Fdummy,itype,jtype);

	// store the force in cumulative force
	for (int dim = 0; dim < 3; dim++)
	  Fcumulative[dim] += Fdummy[dim];

	// compute the term from equation 2.23
	cross3_addition(connector_j,Fdummy,T_j_dummy);

	// compute the term from equation 2.22
	Tpp(p1_j,p1_i,dummy_distance,T_j_dummy,itype,jtype);

	
	// sign change
	for (int dim = 0; dim < 3; dim++)
	  Fdummy[dim] *= -1;
	
	// compute the term from equation 2.20
	cross3_addition(connector_i,Fdummy,T_i_dummy);

	// compute the term from equation 2.19
	Tpp(p1_i,p1_j,dummy_distance,T_i_dummy,itype,jtype);




	

	
	// second pair:
	for (int dim = 0; dim < 3; dim++)
	  dummy_distance[dim] = (connector_j[dim] + connector_i[dim])  + rij[dim];

	// compute the term from equation 2.5
	Fpp(p2_i,p1_j,dummy_distance,Fdummy,itype,jtype);

	
	for (int dim = 0; dim < 3; dim++)
	  Fcumulative[dim] += Fdummy[dim];

	// compute the term from equation 2.20
	cross3_addition(connector_i,Fdummy,T_i_dummy);

	// compute the term from equation 2.19
	Tpp(p2_i,p1_j,dummy_distance,T_i_dummy,itype,jtype);


	// sign change
	for (int dim = 0; dim < 3; dim++)
	  Fdummy[dim] *= -1;
	
	// compute the term from equation 2.23
	cross3_addition(connector_j,Fdummy,T_j_dummy);

	// compute the term from equation 2.22
	Tpp(p1_j,p2_i,dummy_distance,T_j_dummy,itype,jtype);



	
	// third pair:
	for (int dim = 0; dim < 3; dim++)
	  dummy_distance[dim] = -(connector_j[dim] + connector_i[dim])  + rij[dim];

	// compute the term from equation 2.5
	Fpp(p1_i,p2_j,dummy_distance,Fdummy,itype,jtype);

	for (int dim = 0; dim < 3; dim++)
	  Fcumulative[dim] += Fdummy[dim];

	// compute the term from equation 2.23
	cross3_addition(connector_j,Fdummy,T_j_dummy);

	// compute the term from equation 2.22
	Tpp(p2_j,p1_i,dummy_distance,T_j_dummy,itype,jtype);

	// sign change
	for (int dim = 0; dim < 3; dim++)
	  Fdummy[dim] *= -1;

	// compute the term from equation 2.20
	cross3_addition(connector_i,Fdummy,T_i_dummy);

	// compute the term from equation 2.19
	Tpp(p1_i,p2_j,dummy_distance,T_i_dummy,itype,jtype);




	// fourth pair:
	for (int dim = 0; dim < 3; dim++)
	  dummy_distance[dim] = -(connector_j[dim] - connector_i[dim])  + rij[dim];

	// compute the term from equation 2.5
	Fpp(p2_i,p2_j,dummy_distance,Fdummy,itype,jtype);

	for (int dim = 0; dim < 3; dim++)
	  Fcumulative[dim] += Fdummy[dim];

	// compute the term from equation 2.20
	cross3_addition(connector_i,Fdummy,T_i_dummy);

	// compute the term from equation 2.19
	Tpp(p2_i,p2_j,dummy_distance,T_i_dummy,itype,jtype);

	// sign change
	for (int dim = 0; dim < 3; dim++)
	  Fdummy[dim] *= -1;

	// compute the term from equation 2.23
	cross3_addition(connector_j,Fdummy,T_j_dummy);

	// compute the term from equation 2.22
	Tpp(p2_j,p2_i,dummy_distance,T_j_dummy,itype,jtype);

	
	for (int dim = 0; dim < 3; dim++)
	  torque[i][dim] += T_i_dummy[dim];

	for (int dim = 0; dim < 3; dim++)
	  f[i][dim] += Fcumulative[dim];

        if (newton_pair || j < nlocal) {
	  for (int dim = 0; dim < 3; dim++) {
	    f[j][dim] -= Fcumulative[dim];
	    torque[j][dim] += T_j_dummy[dim];
	  }
        }

        if (eflag) {
	  ecoul = 0.0;
	  evdwl = 0.0;
        }

        if (evflag) ev_tally_xyz(i,j,nlocal,newton_pair,
                                 evdwl,ecoul,Fcumulative[0],Fcumulative[1],Fcumulative[2],
				 delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();


  
  
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairJanusCut::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut_lj");
  memory->create(cut_sq,n+1,n+1,"pair:cut_ljsq");
  memory->create(p1_body_frame,3,"pair:p1_body_frame");
  memory->create(p2_body_frame,3,"pair:p2_body_frame");
  memory->create(connector_body_frame,3,"pair:connector_body_frame");


  memory->create(Fdummy,3,"pairtmp:Fdummy");
  memory->create(Fcumulative,3,"pairtmp:Fcumulative");
  memory->create(T_i_dummy,3,"pairtmp:T_i_dummy");
  memory->create(T_j_dummy,3,"pairtmp:T_j_dummy");
  memory->create(cross_dummy,3,"pairtmp:cross_dummy");
  memory->create(p_cross_r,3,"pairtmp:p_cross_r");
  memory->create(rij,3,"pairtmp:r");
  memory->create(dummy_connector,3,"pairtmp:dummy_connector");
  memory->create(dummy_distance,3,"pairtmp:dummy_distance");


  memory->create(connector_i,3,"pairtmp:connector_i");
  memory->create(connector_j,3,"pairtmp:connector_j");

  
  memory->create(p1_i,3,"pairtmp:p1_i");
  memory->create(p1_j,3,"pairtmp:p1_j");
  memory->create(p2_i,3,"pairtmp:p2_i");
  memory->create(p2_j,3,"pairtmp:p2_j");
  

  
  
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairJanusCut::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2)
    error->all(FLERR,"Incorrect args in pair_style command");

  if (strcmp(update->unit_style,"electron") == 0)
    error->all(FLERR,"Cannot (yet) use 'electron' units with dipoles");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut[i][j] = cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairJanusCut::coeff(int narg, char **arg)
{
  if (narg != 16 && narg != 17) {
    printf("%d\n",narg);
    error->all(FLERR,"Incorrect args for pair coefficients");
  }
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double p1_dipole = utils::numeric(FLERR,arg[2],false,lmp);
  if (p1_dipole < 0)
    error->all(FLERR,"p1_dipole magnitude must be > 0");
  double p2_dipole = utils::numeric(FLERR,arg[3],false,lmp);
  if (p2_dipole < 0)
    error->all(FLERR,"p2_dipole magnitude must be > 0");



  double norm;
  p1_body_frame[0] = utils::numeric(FLERR,arg[4],false,lmp);
  p1_body_frame[1] = utils::numeric(FLERR,arg[5],false,lmp);
  p1_body_frame[2] = utils::numeric(FLERR,arg[6],false,lmp);


  norm = 0;
  for (int i = 0; i < 3; i++)
    norm += p1_body_frame[i]*p1_body_frame[i];
  norm = sqrt(norm);
  for (int i = 0; i < 3; i++)
    p1_body_frame[i] = p1_body_frame[i]/norm*p1_dipole;


  p2_body_frame[0] = utils::numeric(FLERR,arg[7],false,lmp);
  p2_body_frame[1] = utils::numeric(FLERR,arg[8],false,lmp);
  p2_body_frame[2] = utils::numeric(FLERR,arg[9],false,lmp);

  norm = 0;
  for (int i = 0; i < 3; i++)
    norm += p2_body_frame[i]*p2_body_frame[i];
  norm = sqrt(norm);
  for (int i = 0; i < 3; i++)
    p2_body_frame[i] = p2_body_frame[i]/norm*p2_dipole;
  

  // connector_distance is d (half the distance between the dipoles) in Xichen's paper
  double connector_distance = utils::numeric(FLERR,arg[10],false,lmp);
  
  connector_body_frame[0] = utils::numeric(FLERR,arg[11],false,lmp);
  connector_body_frame[1] = utils::numeric(FLERR,arg[12],false,lmp);
  connector_body_frame[2] = utils::numeric(FLERR,arg[13],false,lmp);

  norm = 0;
  for (int i = 0; i < 3; i++)
    norm += connector_body_frame[i]*connector_body_frame[i];
  norm = sqrt(norm);
  for (int i = 0; i < 3; i++)
    connector_body_frame[i] = connector_body_frame[i]/norm*connector_distance;


  if (strcmp(arg[14], "false") == 0) 
    fixed_orientations = false;
  else if (strcmp(arg[14], "true") == 0)
    fixed_orientations = true;
  else
    error->all(FLERR, "Illegal pair janus/cut command.");


  if (strcmp(arg[15], "false") == 0) 
    self_interacting = false;
  else if (strcmp(arg[15], "true") == 0)
    self_interacting = true;
  else
    error->all(FLERR, "Illegal pair janus/cut command.");
  
  double cut_one = cut_global;

  
  if (narg == 17) cut_one = utils::numeric(FLERR,arg[16],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairJanusCut::init_style()
{
  avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  if (!avec) error->all(FLERR, "Pair janus/cut requires atom style ellipsoid");

  
  if (!atom->torque_flag)
    error->all(FLERR,"Pair janus/cut requires atom attribute torque");
  if (!atom->ellipsoid_flag)
    error->all(FLERR,"Pair janus/cut requires atom attribute ellipsoid");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairJanusCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }


  cut_sq[i][j] = cut[i][j] * cut[i][j];


  cut_sq[j][i] = cut_sq[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairJanusCut::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairJanusCut::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);

      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairJanusCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairJanusCut::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

void *PairJanusCut::extract(const char *str, int &dim)
{
  dim = 2;
  return nullptr;
}


void PairJanusCut::Fpp(const double *pi, const double *pj, const double *rij,
			     double *F, const int itype, const int jtype)

// NOT A CUMULATIVE SUM OF FORCES
{

  double r = MathExtra::len3(rij);
  double r5 = r*r*r*r*r;
  double r7 = r5*r*r;
  //251024 define local cutoff
  double tmp_cut = cut[itype][jtype];
	

  //140824 adding a smooth cutoff at cut_global, following LAMMPS pair dipole shift.
  double r3 = r*r*r;
  double r4 = r3*r;
  
  //double cut3 = cut_global*cut_global*cut_global;
  double cut3 = tmp_cut*tmp_cut*tmp_cut;
  //double cut4 = cut3*cut_global;      
  double cut4 = cut3*tmp_cut;

  double pref1 = 1-r4/cut4;
  double pref2 = 1-4*r3/cut3+3*r4/cut4;

  double pi_dot_pj = MathExtra::dot3(pi,pj);
  double pi_dot_rij = MathExtra::dot3(pi,rij);
  double pj_dot_rij = MathExtra::dot3(pj,rij);


  for (int dim = 0; dim < 3; dim++) 
    F[dim] = 3./r5*(pref1*(pi_dot_pj-3./(r*r)*pi_dot_rij*pj_dot_rij)*rij[dim] + pref2*(pj_dot_rij*pi[dim] + pi_dot_rij*pj[dim]-2./(r*r)*pi_dot_rij*pj_dot_rij*rij[dim]));	  

  return;

}

void PairJanusCut::Tpp(const double *pi, const double *pj, const double *rij,
		       double *T, const int itype, const int jtype)
// A CUMULATIVE SUM OF TORQUES
{
  double r = MathExtra::len3(rij);
  double r3 = r*r*r;
  double r5 = r3*r*r;
  //251024 define local cutoff
  double tmp_cut = cut[itype][jtype];
  
  //140824 adding a smooth cutoff at cut_global, following LAMMPS pair dipole shift.
  double r4 = r3*r;
  //double cut3 = cut_global*cut_global*cut_global;
  double cut3 = tmp_cut*tmp_cut*tmp_cut;
  //double cut4 = cut3*cut_global;	
  double cut4 = cut3*tmp_cut;
  double pref = 1-4*r3/cut3+3*r4/cut4;

  double pj_dot_rij = MathExtra::dot3(pj,rij);

	

  MathExtra::cross3(pi,pj,cross_dummy);
  MathExtra::cross3(pi,rij,p_cross_r);

  
  for (int dim = 0; dim < 3; dim ++ )
    T[dim] += pref*(-1./r3*cross_dummy[dim] + 3./r5*pj_dot_rij*p_cross_r[dim]);

  return;

}


void PairJanusCut::cross3_addition(const double *v1, const double *v2, double *ans)
{
  ans[0] += v1[1] * v2[2] - v1[2] * v2[1];
  ans[1] += v1[2] * v2[0] - v1[0] * v2[2];
  ans[2] += v1[0] * v2[1] - v1[1] * v2[0];
}
