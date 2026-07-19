/* ----------------------------------------------------------------------
    PICI-LIGGGHTS - particle based load balancing support

    Copyright 2026 PICI-LIGGGHTS contributors

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License, version 2 or later.
    The shifted-plane arithmetic follows the GPL LAMMPS balance algorithm.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstring>
#include <vector>
#include "balance.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "lmptype.h"

using namespace LAMMPS_NS;

#define CUT_EPS 1.0e-12

/* ---------------------------------------------------------------------- */

Balance::Balance(LAMMPS *lmp) : Pointers(lmp)
{
  dimstr[0] = '\0';
  ndim = nitermax = 0;
  stopthresh = 1.0;
}

/* ---------------------------------------------------------------------- */

void Balance::shift_setup(const char *dims, int niter, double stop)
{
  ndim = strlen(dims);
  strcpy(dimstr,dims);
  nitermax = niter;
  stopthresh = stop;
}

/* ----------------------------------------------------------------------
   reduced coordinate, independent of orthogonal/triclinic representation
   triclinic atoms are in lamda coordinates while shift() is called
------------------------------------------------------------------------- */

double Balance::reduced_coord(int i, int dim)
{
  if (domain->triclinic) return atom->x[i][dim];
  return (atom->x[i][dim]-domain->boxlo[dim])/domain->prd[dim];
}

/* ---------------------------------------------------------------------- */

int Balance::split_index(double value, int nsplit, double *split)
{
  if (value <= 0.0) return 0;
  if (value >= 1.0) return nsplit-1;

  int lo = 0;
  int hi = nsplit;
  while (lo+1 < hi) {
    int mid = (lo+hi)/2;
    if (value < split[mid]) hi = mid;
    else lo = mid;
  }
  return lo;
}

/* ---------------------------------------------------------------------- */

void Balance::validate_splits(int nsplit, const double *split)
{
  if (!std::isfinite(split[0]) || !std::isfinite(split[nsplit]) ||
      split[0] != 0.0 || split[nsplit] != 1.0)
    error->all(FLERR,"Invalid load balance split endpoints");

  for (int i = 1; i <= nsplit; i++)
    if (!std::isfinite(split[i]) || split[i] <= split[i-1])
      error->all(FLERR,"Load balance generated crossing or zero-width splits");
}

/* ----------------------------------------------------------------------
   cumulative number of atoms below each trial cut
------------------------------------------------------------------------- */

void Balance::cumulative_counts(int dim, const double *cut, int ncut,
                                bigint *counts)
{
  std::vector<bigint> local(ncut,0);
  for (int i = 0; i < atom->nlocal; i++) {
    double value = reduced_coord(i,dim);
    for (int j = 0; j < ncut; j++)
      if (value < cut[j]) local[j]++;
  }
  MPI_Allreduce(&local[0],counts,ncut,MPI_LMP_BIGINT,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   imbalance of the currently installed split arrays
------------------------------------------------------------------------- */

double Balance::imbalance_splits(double &maxload)
{
  const int nprocs = comm->nprocs;
  std::vector<bigint> local(nprocs,0);
  std::vector<bigint> all(nprocs,0);

  for (int i = 0; i < atom->nlocal; i++) {
    int ix = split_index(reduced_coord(i,0),comm->procgrid[0],comm->xsplit);
    int iy = split_index(reduced_coord(i,1),comm->procgrid[1],comm->ysplit);
    int iz = split_index(reduced_coord(i,2),comm->procgrid[2],comm->zsplit);
    int proc = comm->grid2proc[ix][iy][iz];
    local[proc]++;
  }

  MPI_Allreduce(&local[0],&all[0],nprocs,MPI_LMP_BIGINT,MPI_SUM,world);

  bigint ntotal = 0;
  bigint nmax = 0;
  for (int i = 0; i < nprocs; i++) {
    ntotal += all[i];
    if (all[i] > nmax) nmax = all[i];
  }

  maxload = static_cast<double>(nmax);
  if (ntotal == 0) return 0.0;
  return maxload/(static_cast<double>(ntotal)/nprocs);
}

/* ----------------------------------------------------------------------
   imbalance of the present ownership
------------------------------------------------------------------------- */

double Balance::imbalance_factor(double &maxload)
{
  bigint nlocal = atom->nlocal;
  bigint ntotal, nmax;
  MPI_Allreduce(&nlocal,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);
  MPI_Allreduce(&nlocal,&nmax,1,MPI_LMP_BIGINT,MPI_MAX,world);
  maxload = static_cast<double>(nmax);
  if (ntotal == 0) return 0.0;
  return maxload/(static_cast<double>(ntotal)/comm->nprocs);
}

/* ----------------------------------------------------------------------
   density-dependent recursive multisection of selected Cartesian axes
------------------------------------------------------------------------- */

double Balance::shift(int &niter)
{
  bigint nlocal = atom->nlocal;
  bigint ntotal;
  MPI_Allreduce(&nlocal,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);

  niter = 0;
  double maxload;
  if (ntotal == 0) return imbalance_splits(maxload);

  for (int id = 0; id < ndim; id++) {
    int dim = dimstr[id] == 'x' ? 0 : (dimstr[id] == 'y' ? 1 : 2);
    int nsplit = comm->procgrid[dim];
    if (nsplit == 1) continue;

    double *split = dim == 0 ? comm->xsplit :
                    (dim == 1 ? comm->ysplit : comm->zsplit);
    validate_splits(nsplit,split);

    int ncut = nsplit-1;
    std::vector<double> cut(ncut), lower(ncut,0.0), upper(ncut,1.0);
    std::vector<bigint> counts(ncut), nlower(ncut,0), nupper(ncut,ntotal);
    for (int i = 0; i < ncut; i++) cut[i] = split[i+1];

    for (int iter = 0; iter < nitermax; iter++) {
      cumulative_counts(dim,&cut[0],ncut,&counts[0]);

      for (int i = 0; i < ncut; i++) {
        double target = static_cast<double>(ntotal)*(i+1)/nsplit;
        if (static_cast<double>(counts[i]) < target) {
          lower[i] = cut[i];
          nlower[i] = counts[i];
        } else {
          upper[i] = cut[i];
          nupper[i] = counts[i];
        }

        if (nupper[i] > nlower[i]) {
          double fraction = (target-static_cast<double>(nlower[i])) /
            static_cast<double>(nupper[i]-nlower[i]);
          cut[i] = lower[i] + fraction*(upper[i]-lower[i]);
        } else cut[i] = 0.5*(lower[i]+upper[i]);
      }

      // Discrete particle coordinates can yield coincident quantiles.
      // Project them to a strictly ordered set without changing endpoints.

      for (int i = 0; i < ncut; i++) {
        double lo = (i == 0 ? 0.0 : cut[i-1]) + CUT_EPS;
        double hi = 1.0-(ncut-i)*CUT_EPS;
        if (cut[i] < lo) cut[i] = lo;
        if (cut[i] > hi) cut[i] = hi;
      }

      for (int i = 0; i < ncut; i++) split[i+1] = cut[i];
      split[0] = 0.0;
      split[nsplit] = 1.0;
      validate_splits(nsplit,split);
      niter++;

      double imbalance = imbalance_splits(maxload);
      if (imbalance <= stopthresh) return imbalance;
    }
  }

  return imbalance_splits(maxload);
}
