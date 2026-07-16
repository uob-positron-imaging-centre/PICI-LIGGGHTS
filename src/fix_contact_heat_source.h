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

#ifdef FIX_CLASS

FixStyle(contact/heat/source,FixContactHeatSource)

#else

#ifndef LMP_FIX_CONTACT_HEAT_SOURCE_H
#define LMP_FIX_CONTACT_HEAT_SOURCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixContactHeatSource : public Fix {
 public:
  FixContactHeatSource(class LAMMPS *, int, char **);

  int setmask();
  void post_create();
  void init();
  void setup_pre_force(int);
  void pre_force(int);
  void setup_pre_force_respa(int,int);
  void pre_force_respa(int,int,int);
  double compute_scalar();

  void enable_tracking();
  void add_particle_heat(int, double);
  void reverse_comm();
  double *contact_heat_source();

 private:
  void update_ptrs();
  void reset_source();

  class FixPropertyAtom *fix_contact_heat_source_;
  int nlevels_respa_;
  bigint last_reset_;
};

}

#endif
#endif
