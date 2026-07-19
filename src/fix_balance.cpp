/* ----------------------------------------------------------------------
    PICI-LIGGGHTS - dynamic particle-count load balancing

    Copyright 2026 PICI-LIGGGHTS contributors

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License, version 2 or later.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstring>
#include "fix_balance.h"
#include "atom.h"
#include "balance.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "irregular.h"
#include "lammps.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBalance::FixBalance(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp,narg,arg), balance(NULL), irregular(NULL)
{
  if (narg != 9) error->fix_error(FLERR,this,"illegal command syntax");

  nfreq = force->inumeric(FLERR,arg[3]);
  thresh = force->numeric(FLERR,arg[4]);
  if (strcmp(arg[5],"shift"))
    error->fix_error(FLERR,this,"only the 'shift' style is supported");

  int seen[3] = {0,0,0};
  int ndim = strlen(arg[6]);
  if (ndim < 1 || ndim > 3)
    error->fix_error(FLERR,this,"shift dimensions must contain unique x, y, or z characters");
  for (int i = 0; i < ndim; i++) {
    int dim = arg[6][i] == 'x' ? 0 : (arg[6][i] == 'y' ? 1 :
              (arg[6][i] == 'z' ? 2 : -1));
    if (dim < 0 || seen[dim])
      error->fix_error(FLERR,this,"shift dimensions must contain unique x, y, or z characters");
    if (dim == 2 && domain->dimension == 2)
      error->fix_error(FLERR,this,"cannot balance the z dimension in a 2d simulation");
    seen[dim] = 1;
  }
  strcpy(dimstr,arg[6]);

  nitermax = force->inumeric(FLERR,arg[7]);
  stopthresh = force->numeric(FLERR,arg[8]);
  if (nfreq < 0) error->fix_error(FLERR,this,"Nfreq must be non-negative");
  if (!std::isfinite(thresh) || thresh <= 0.0)
    error->fix_error(FLERR,this,"thresh must be finite and positive");
  if (nitermax <= 0) error->fix_error(FLERR,this,"Niter must be positive");
  if (!std::isfinite(stopthresh) || stopthresh < 1.0)
    error->fix_error(FLERR,this,"stopthresh must be finite and at least 1.0");

  scalar_flag = vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 0;
  extvector = 0;

  force_reneighbor = nfreq > 0;
  pre_exchange_migrate = 1;
  box_change_domain = 1;
  next_reneighbor = -1;

  imbalance = imbalance_before = maxload = 0.0;
  iterations = pending = 0;
  lastcheck = -1;
}

/* ---------------------------------------------------------------------- */

FixBalance::~FixBalance()
{
  delete balance;
  delete irregular;
}

/* ---------------------------------------------------------------------- */

int FixBalance::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBalance::init()
{
  if (domain->is_wedge)
    error->fix_error(FLERR,this,"wedge domains are not supported");
  if (update->whichflag == 2)
    error->fix_error(FLERR,this,"cannot be used during minimization");
  if (!update->integrate_style || strcmp(update->integrate_style,"verlet"))
    error->fix_error(FLERR,this,"requires the non-CUDA Verlet run style");
  if (lmp->cuda)
    error->fix_error(FLERR,this,"does not support CUDA");
  if (modify->n_fixes_style_strict("balance") > 1)
    error->fix_error(FLERR,this,"only one fix balance may be active");
  if (modify->n_fixes_style_strict("couple/pascal") || comm->exchangeEvents)
    error->fix_error(FLERR,this,"is incompatible with PASCAL exchange-event coupling");

  delete balance;
  delete irregular;
  balance = new Balance(lmp);
  balance->shift_setup(dimstr,nitermax,stopthresh);
  irregular = new Irregular(lmp);
  lastcheck = -1;
  schedule_next();
}

/* ----------------------------------------------------------------------
   Initial Verlet setup checks imbalance.  whichflag is zero for write_data
   and write_restart, so those commands cannot cause redistribution.
------------------------------------------------------------------------- */

void FixBalance::setup_pre_exchange()
{
  if (update->whichflag != 1) return;
  check_balance();
  schedule_next();
}

/* ---------------------------------------------------------------------- */

void FixBalance::pre_exchange()
{
  if (update->whichflag != 1) return;
  if (nfreq && update->ntimestep != next_reneighbor) return;
  check_balance();
  schedule_next();
}

/* ---------------------------------------------------------------------- */

void FixBalance::schedule_next()
{
  if (nfreq == 0) {
    next_reneighbor = -1;
    return;
  }
  next_reneighbor = (update->ntimestep/nfreq+1)*nfreq;
}

/* ---------------------------------------------------------------------- */

void FixBalance::check_balance()
{
  if (lastcheck == update->ntimestep) return;
  lastcheck = update->ntimestep;

  imbalance_before = balance->imbalance_factor(maxload);
  imbalance = imbalance_before;
  iterations = 0;
  pending = 0;
  if (imbalance_before > thresh) rebalance();
}

/* ---------------------------------------------------------------------- */

void FixBalance::rebalance()
{
  // Bring atoms and shrink-wrapped bounds current before choosing cuts.

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();

  imbalance = balance->shift(iterations);
  comm->uniform = 0;
  if (domain->triclinic) domain->set_lamda_box();
  domain->set_local_box();

  // Irregular migration preserves all AtomVec and per-atom Fix state and can
  // cross more than one processor plane in a single exchange.

  irregular->migrate_atoms();

  bigint nlocal = atom->nlocal;
  bigint ntotal;
  MPI_Allreduce(&nlocal,&ntotal,1,MPI_LMP_BIGINT,MPI_SUM,world);
  if (ntotal != atom->natoms)
    error->all(FLERR,"Fix balance lost atoms during redistribution");

  if (domain->triclinic) domain->lamda2x(atom->nlocal);
  imbalance = balance->imbalance_factor(maxload);
  pending = 1;
}

/* ----------------------------------------------------------------------
   Called after the normal exchange/border phase; refresh public statistics.
------------------------------------------------------------------------- */

void FixBalance::pre_neighbor()
{
  if (!pending) return;
  imbalance = balance->imbalance_factor(maxload);
  pending = 0;
}

/* ---------------------------------------------------------------------- */

void FixBalance::reset_timestep(bigint)
{
  lastcheck = -1;
  schedule_next();
}

/* ---------------------------------------------------------------------- */

double FixBalance::compute_scalar()
{
  return imbalance;
}

/* ---------------------------------------------------------------------- */

double FixBalance::compute_vector(int i)
{
  if (i == 0) return maxload;
  if (i == 1) return iterations;
  return imbalance_before;
}
