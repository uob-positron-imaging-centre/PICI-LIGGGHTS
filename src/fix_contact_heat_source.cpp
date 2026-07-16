/* ----------------------------------------------------------------------
    This is the

    LIGGGHTS DEM simulation engine, released by
    DCS Computing GmbH, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS is part of CFDEMproject:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

-------------------------------------------------------------------------
    Contributing author for this file:
    OpenAI Codex (PICI-LIGGGHTS contact heat-generation implementation)
------------------------------------------------------------------------- */

#include "fix_contact_heat_source.h"

#include "atom.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_scalar_transport_equation.h"
#include "modify.h"
#include "respa.h"
#include "update.h"
#include "mpi.h"
#include <cmath>
#include <string.h>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixContactHeatSource::FixContactHeatSource(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp,narg,arg),
  fix_contact_heat_source_(NULL),
  nlevels_respa_(0),
  last_reset_(-1)
{
  if (narg != 3)
    error->fix_error(FLERR,this,"does not accept additional arguments");

  peratom_flag = 1;
  size_peratom_cols = 0;
  peratom_freq = 1;

  scalar_flag = 0;
  global_freq = 1;
  extscalar = 1;
}

/* ---------------------------------------------------------------------- */

int FixContactHeatSource::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  mask |= PRE_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixContactHeatSource::post_create()
{
  fix_contact_heat_source_ = static_cast<FixPropertyAtom*>(
      modify->find_fix_property("contactHeatSource","property/atom","scalar",
                                0,0,style,false));

  if (!fix_contact_heat_source_) {
    const char *fixarg[9];
    fixarg[0] = "contactHeatSource";
    fixarg[1] = "all";
    fixarg[2] = "property/atom";
    fixarg[3] = "contactHeatSource";
    fixarg[4] = "scalar";
    fixarg[5] = "no";  // transient: do not write to restart
    fixarg[6] = "no";  // owner values are not needed on ghosts
    fixarg[7] = "yes"; // ghost contributions are reverse communicated
    fixarg[8] = "0.0";
    fix_contact_heat_source_ = modify->add_fix_property_atom(
        9,const_cast<char**>(fixarg),style);
  }

  update_ptrs();
}

/* ---------------------------------------------------------------------- */

void FixContactHeatSource::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa_ = static_cast<Respa*>(update->integrate)->nlevels;

  FixScalarTransportEquation *fix_ste =
      modify->find_fix_scalar_transport_equation("heattransfer");
  if (!fix_ste)
    error->fix_error(FLERR,this,
        "requires fix heat/gran or a heattransfer transportequation/scalar fix");
  fix_ste->validate_capacity_for_contact_heat();

  if (!modify->find_fix_property("Temp","property/atom","scalar",0,0,
                                 style,false))
    error->fix_error(FLERR,this,"requires the per-particle Temp property");

  if (!modify->find_fix_property("heatSource","property/atom","scalar",0,0,
                                 style,false))
    error->fix_error(FLERR,this,"requires the per-particle heatSource property");

  if (!atom->rmass_flag)
    error->fix_error(FLERR,this,"requires per-particle mass");
  for (int i = 0; i < atom->nlocal; i++) {
    if (!std::isfinite(atom->rmass[i]) || atom->rmass[i] <= 0.)
      error->one(FLERR,
          "Contact heat generation requires positive finite particle mass");
  }

  update_ptrs();
}

/* ---------------------------------------------------------------------- */

void FixContactHeatSource::update_ptrs()
{
  if (!fix_contact_heat_source_)
    return;
  vector_atom = fix_contact_heat_source_->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixContactHeatSource::reset_source()
{
  if (last_reset_ == update->ntimestep)
    return;

  update_ptrs();
  const int nall = atom->nlocal + atom->nghost;
  for (int i = 0; i < nall; i++)
    vector_atom[i] = 0.;

  last_reset_ = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixContactHeatSource::setup_pre_force(int vflag)
{
  last_reset_ = -1;
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixContactHeatSource::pre_force(int)
{
  reset_source();
}

/* ---------------------------------------------------------------------- */

void FixContactHeatSource::setup_pre_force_respa(int vflag, int ilevel)
{
  if (ilevel == nlevels_respa_-1) {
    last_reset_ = -1;
    pre_force(vflag);
  }
}

/* ---------------------------------------------------------------------- */

void FixContactHeatSource::pre_force_respa(int vflag, int ilevel, int flag)
{
  if (ilevel == nlevels_respa_-1)
    pre_force(vflag);
  UNUSED(flag);
}

/* ---------------------------------------------------------------------- */

void FixContactHeatSource::enable_tracking()
{
  scalar_flag = 1;
}

/* ---------------------------------------------------------------------- */

void FixContactHeatSource::add_particle_heat(int i, double power)
{
  update_ptrs();
  vector_atom[i] += power;
}

/* ---------------------------------------------------------------------- */

void FixContactHeatSource::reverse_comm()
{
  fix_contact_heat_source_->do_reverse_comm();
}

/* ---------------------------------------------------------------------- */

double *FixContactHeatSource::contact_heat_source()
{
  update_ptrs();
  return vector_atom;
}

/* ---------------------------------------------------------------------- */

double FixContactHeatSource::compute_scalar()
{
  update_ptrs();
  double local_power = 0.;
  for (int i = 0; i < atom->nlocal; i++)
    local_power += vector_atom[i];

  double total_power = 0.;
  MPI_Allreduce(&local_power,&total_power,1,MPI_DOUBLE,MPI_SUM,world);
  return total_power;
}
