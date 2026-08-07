/* ----------------------------------------------------------------------
    This is the

    LIGGGHTS(R) - DEM simulation engine
    https://www.cfdem.com

    LIGGGHTS(R) is open-source, distributed under the terms of the GNU
    Public License, version 2 or later.

    This common pair-style plumbing is adapted from the public LAMMPS
    lubrication pair styles.  The unequal-sphere resistance model itself
    is implemented in pair_lubricate_bmpoly.cpp.
------------------------------------------------------------------------- */

#include <stdio.h>
#include <string.h>

#include "pair_lubricate.h"
#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLubricate::PairLubricate(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  comm_forward = 6;
  no_virial_fdotr_compute = 1;

  mu = 0.0;
  flaglog = flagfld = flagHI = flagVF = 0;
  cut_inner_global = cut_global = 0.0;
  cut = cut_inner = NULL;

  vol_P = 0.0;
  R0 = RT0 = RS0 = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Ef[i][j] = 0.0;

  shearing = flagdeform = flagwall = 0;
  wallfix = NULL;
}
/* ---------------------------------------------------------------------- */

PairLubricate::~PairLubricate()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(cut_inner);
  }
}

/* ---------------------------------------------------------------------- */

void PairLubricate::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(cut_inner,n+1,n+1,"pair:cut_inner");

  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;
}

/* ---------------------------------------------------------------------- */

void PairLubricate::settings(int narg, char **arg)
{
  if (narg != 5 && narg != 7)
    error->all(FLERR,"Illegal pair_style command");

  mu = force->numeric(FLERR,arg[0]);
  flaglog = force->inumeric(FLERR,arg[1]);
  flagfld = force->inumeric(FLERR,arg[2]);
  cut_inner_global = force->numeric(FLERR,arg[3]);
  cut_global = force->numeric(FLERR,arg[4]);
  flagHI = flagVF = 1;
  if (narg == 7) {
    flagHI = force->inumeric(FLERR,arg[5]);
    flagVF = force->inumeric(FLERR,arg[6]);
  }

  if (mu <= 0.0)
    error->all(FLERR,"Pair lubricate requires positive viscosity");
  if ((flaglog != 0 && flaglog != 1) ||
      (flagfld != 0 && flagfld != 1) ||
      (flagHI != 0 && flagHI != 1) ||
      (flagVF != 0 && flagVF != 1))
    error->all(FLERR,"Pair lubricate flags must be 0 or 1");
  if (cut_inner_global <= 0.0 || cut_global <= 0.0 ||
      cut_inner_global > cut_global)
    error->all(FLERR,"Pair lubricate requires 0 < cutinner <= cutoff");

  if (flaglog && !flagHI) {
    error->warning(FLERR,"Pair lubricate log terms require 1/r terms; setting flagHI to 1");
    flagHI = 1;
  }

  if (allocated) {
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_inner[i][j] = cut_inner_global;
          cut[i][j] = cut_global;
        }
  }
}

/* ---------------------------------------------------------------------- */

void PairLubricate::coeff(int narg, char **arg)
{
  if (narg != 2 && narg != 4)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double cut_inner_one = cut_inner_global;
  double cut_one = cut_global;
  if (narg == 4) {
    cut_inner_one = force->numeric(FLERR,arg[2]);
    cut_one = force->numeric(FLERR,arg[3]);
  }

  if (cut_inner_one <= 0.0 || cut_one <= 0.0 ||
      cut_inner_one > cut_one)
    error->all(FLERR,"Pair lubricate requires 0 < cutinner <= cutoff");

  int count = 0;
  for (int i = ilo; i <= ihi; i++)
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      cut_inner[i][j] = cut_inner_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }

  if (count == 0)
    error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ---------------------------------------------------------------------- */

void PairLubricate::init_style()
{
  if (!atom->omega_flag)
    error->all(FLERR,"Pair lubricate requires atom attribute omega");
  if (!atom->radius_flag)
    error->all(FLERR,"Pair lubricate requires atom attribute radius");
  if (!atom->torque_flag)
    error->all(FLERR,"Pair lubricate requires atom attribute torque");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair lubricate requires ghost atoms store velocity");

  int bad = 0;
  for (int i = 0; i < atom->nlocal; i++)
    if (atom->radius[i] <= 0.0)
      bad = 1;
  int badall = 0;
  MPI_Allreduce(&bad,&badall,1,MPI_INT,MPI_MAX,world);
  if (badall)
    error->all(FLERR,"Pair lubricate requires positive particle radii");
}

/* ---------------------------------------------------------------------- */

void PairLubricate::init_list(int which, NeighList *ptr)
{
  if (which == 0) list = ptr;
  else error->all(FLERR,"Internal error in PairLubricate::init_list");
}

/* ---------------------------------------------------------------------- */

double PairLubricate::max_radius_for_type(int itype)
{
  double radius_local = 0.0;
  for (int i = 0; i < atom->nlocal; i++)
    if (atom->type[i] == itype && atom->radius[i] > radius_local)
      radius_local = atom->radius[i];

  double radius_global = 0.0;
  MPI_Allreduce(&radius_local,&radius_global,1,MPI_DOUBLE,MPI_MAX,world);
  return radius_global;
}

/* ---------------------------------------------------------------------- */

double PairLubricate::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut_inner[i][j] = mix_distance(cut_inner[i][i],cut_inner[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  if (cut_inner[i][j] <= 0.0 || cut[i][j] <= 0.0 ||
      cut_inner[i][j] > cut[i][j])
    error->all(FLERR,"Pair lubricate has invalid gap cutoffs");

  cut_inner[j][i] = cut_inner[i][j];
  cut[j][i] = cut[i][j];

  double radius_i = max_radius_for_type(i);
  double radius_j = max_radius_for_type(j);

  // cut and cut_inner are surface-to-surface gap distances.  Pair::init()
  // expects init_one() to return the centre-to-centre neighbour cutoff.
  return radius_i + radius_j + cut[i][j];
}

/* ---------------------------------------------------------------------- */

int PairLubricate::pack_comm(int n, int *list, double *buf,
                             int /*pbc_flag*/, int * /*pbc*/)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    buf[m++] = atom->v[j][0];
    buf[m++] = atom->v[j][1];
    buf[m++] = atom->v[j][2];
    buf[m++] = atom->omega[j][0];
    buf[m++] = atom->omega[j][1];
    buf[m++] = atom->omega[j][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairLubricate::unpack_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    atom->v[i][0] = buf[m++];
    atom->v[i][1] = buf[m++];
    atom->v[i][2] = buf[m++];
    atom->omega[i][0] = buf[m++];
    atom->omega[i][1] = buf[m++];
    atom->omega[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void PairLubricate::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cut_inner[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ---------------------------------------------------------------------- */

void PairLubricate::read_restart(FILE *fp, const int major, const int minor)
{
  read_restart_settings(fp,major,minor);
  allocate();

  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) {
      if (comm->me == 0)
        fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (comm->me == 0) {
          fread(&cut_inner[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&cut_inner[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ---------------------------------------------------------------------- */

void PairLubricate::write_restart_settings(FILE *fp)
{
  fwrite(&mu,sizeof(double),1,fp);
  fwrite(&flaglog,sizeof(int),1,fp);
  fwrite(&flagfld,sizeof(int),1,fp);
  fwrite(&cut_inner_global,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&flagHI,sizeof(int),1,fp);
  fwrite(&flagVF,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void PairLubricate::read_restart_settings(FILE *fp, const int /*major*/,
                                           const int /*minor*/)
{
  if (comm->me == 0) {
    fread(&mu,sizeof(double),1,fp);
    fread(&flaglog,sizeof(int),1,fp);
    fread(&flagfld,sizeof(int),1,fp);
    fread(&cut_inner_global,sizeof(double),1,fp);
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&flagHI,sizeof(int),1,fp);
    fread(&flagVF,sizeof(int),1,fp);
  }

  MPI_Bcast(&mu,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&flaglog,1,MPI_INT,0,world);
  MPI_Bcast(&flagfld,1,MPI_INT,0,world);
  MPI_Bcast(&cut_inner_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&flagHI,1,MPI_INT,0,world);
  MPI_Bcast(&flagVF,1,MPI_INT,0,world);
}
