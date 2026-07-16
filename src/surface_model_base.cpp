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

#include "surface_model_base.h"

#include "atom.h"
#include "contact_models.h"
#include "error.h"
#include "fix_contact_heat_source.h"
#include "force.h"
#include "modify.h"
#include "settings.h"
#include <cmath>

using namespace LAMMPS_NS;

namespace LIGGGHTS {
namespace ContactModels {

/* ---------------------------------------------------------------------- */

SurfaceModelBase::SurfaceModelBase(LAMMPS *lmp, IContactHistorySetup*,
                                   ContactModelBase *cmb) :
  Pointers(lmp),
  heat_generation_(false),
  heat_generation_tracking_(false),
  is_wall_(cmb ? cmb->is_wall() : false),
  warned_negative_dissipation_(false),
  heat_generation_efficiency_(1.),
  heat_partition_pp_(0.5),
  heat_partition_wall_(0.5),
  fix_contact_heat_source_(NULL)
{}

/* ---------------------------------------------------------------------- */

void SurfaceModelBase::registerContactHeatSettings(Settings &settings)
{
  settings.registerOnOff("heat_generation",heat_generation_,false);
  settings.registerDoubleSetting("heat_generation_efficiency",
                                 heat_generation_efficiency_,1.);
  settings.registerDoubleSetting("heat_partition_pp",heat_partition_pp_,0.5);
  settings.registerDoubleSetting("heat_partition_wall",
                                 heat_partition_wall_,0.5);
  settings.registerOnOff("heat_generation_tracking",
                         heat_generation_tracking_,false);
}

/* ---------------------------------------------------------------------- */

void SurfaceModelBase::postContactHeatSettings()
{
  if (!std::isfinite(heat_generation_efficiency_) ||
      heat_generation_efficiency_ < 0. || heat_generation_efficiency_ > 1.)
    error->all(FLERR,
        "heat_generation_efficiency must be between 0 and 1");
  if (!std::isfinite(heat_partition_pp_) ||
      heat_partition_pp_ < 0. || heat_partition_pp_ > 1.)
    error->all(FLERR,"heat_partition_pp must be between 0 and 1");
  if (!std::isfinite(heat_partition_wall_) ||
      heat_partition_wall_ < 0. || heat_partition_wall_ > 1.)
    error->all(FLERR,"heat_partition_wall must be between 0 and 1");

  if (!heat_generation_ || is_wall_)
    return;

  createContactHeatFix();
}

/* ---------------------------------------------------------------------- */

void SurfaceModelBase::createContactHeatFix()
{
  if (!heat_generation_)
    return;

  fix_contact_heat_source_ = static_cast<FixContactHeatSource*>(
      modify->find_fix_style("contact/heat/source",0));
  if (!fix_contact_heat_source_) {
    const char *fixarg[3];
    fixarg[0] = "totalContactHeatPower";
    fixarg[1] = "all";
    fixarg[2] = "contact/heat/source";
    modify->add_fix(3,const_cast<char**>(fixarg));
    fix_contact_heat_source_ = static_cast<FixContactHeatSource*>(
        modify->find_fix_style("contact/heat/source",0));
  }

  if (!fix_contact_heat_source_)
    error->all(FLERR,"Could not create contact heat source fix");
  if (heat_generation_tracking_)
    fix_contact_heat_source_->enable_tracking();
}

/* ---------------------------------------------------------------------- */

void SurfaceModelBase::depositContactHeat(SurfacesIntersectData &sidata)
{
  if (!heat_generation_ || !sidata.computeflag)
    return;

  double dissipated_power = sidata.P_diss;
  if (!std::isfinite(dissipated_power))
    error->one(FLERR,"Contact dissipated power is not finite");

  if (dissipated_power < 0.) {
    if (dissipated_power < -1.e-12 && !warned_negative_dissipation_) {
      error->warning(FLERR,
          "Negative contact dissipated power detected; contribution ignored");
      warned_negative_dissipation_ = true;
    }
    dissipated_power = 0.;
  }

  const double heat_power = heat_generation_efficiency_ * dissipated_power;
  if (heat_power == 0.)
    return;

  if (sidata.is_wall) {
    fix_contact_heat_source_->add_particle_heat(
        sidata.i,heat_partition_wall_ * heat_power);
    return;
  }

  const bool i_is_first = atom->tag[sidata.i] < atom->tag[sidata.j];
  const double fraction_i = i_is_first ? heat_partition_pp_ :
                                           1. - heat_partition_pp_;
  fix_contact_heat_source_->add_particle_heat(sidata.i,
                                               fraction_i * heat_power);

  if (force->newton_pair || sidata.j < atom->nlocal)
    fix_contact_heat_source_->add_particle_heat(
        sidata.j,(1. - fraction_i) * heat_power);
}

/* ---------------------------------------------------------------------- */

void SurfaceModelBase::finishContactHeatPass()
{
  if (heat_generation_ && !is_wall_ && force->newton_pair)
    fix_contact_heat_source_->reverse_comm();
}

}
}
