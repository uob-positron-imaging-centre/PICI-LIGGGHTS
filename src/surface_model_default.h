/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:

    Christoph Kloss (DCS Computing GmbH, Linz)
    Christoph Kloss (JKU Linz)
    Richard Berger (JKU Linz)

    Copyright 2012-     DCS Computing GmbH, Linz
    Copyright 2009-2012 JKU Linz
------------------------------------------------------------------------- */

#ifdef SURFACE_MODEL
SURFACE_MODEL(SURFACE_DEFAULT,default,0)
#else
#ifndef SURFACE_MODEL_DEFAULT_H_
#define SURFACE_MODEL_DEFAULT_H_
#include "contact_models.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix_property_atom.h"
#include "surface_model_base.h"
#include "math_extra.h"
#include "vector_liggghts.h"

namespace LIGGGHTS {
namespace ContactModels
{
  using namespace MathExtra;
  template<>
  class SurfaceModel<SURFACE_DEFAULT> : public SurfaceModelBase
  {
  public:
    SurfaceModel(LAMMPS * lmp, IContactHistorySetup * hsetup, class ContactModelBase * cmb) :
        SurfaceModelBase(lmp, hsetup, cmb),
        elasticpotflag_(false),
        dissipatedflag_(false),
        delta_offset_(-1),
        dissipation_offset_(-1),
        fix_dissipated_(NULL)
    {
      
    }

    inline void registerSettings(Settings& settings)
    {
        settings.registerOnOff("computeElasticPotential", elasticpotflag_, false);
        settings.registerOnOff("computeDissipatedEnergy", dissipatedflag_, false);
    }

    inline void postSettings(IContactHistorySetup * hsetup, ContactModelBase *cmb)
    {
        if (dissipatedflag_)
        {
            if (cmb->is_wall())
            {
                fix_dissipated_ = static_cast<FixPropertyAtom*>(modify->find_fix_property("dissipated_energy_wall", "property/atom", "vector", 0, 0, "dissipated energy"));
                if (!fix_dissipated_)
                    error->one(FLERR, "Could not find dissipated_energy_wall atom property. Ensure that fix calculate/wall_dissipated_energy is before fix wall/gran");
            }
            else
            {
                char * fixarg[15];
                fixarg[0]  = (char*)"dissipated_energy_";
                fixarg[1]  = (char*)"all";
                fixarg[2]  = (char*)"property/atom";
                fixarg[3]  = (char*)"dissipated_energy";
                fixarg[4]  = (char*)"vector";
                fixarg[5]  = (char*)"yes";
                fixarg[6]  = (char*)"yes";
                fixarg[7]  = (char*)"no";
                fixarg[8]  = (char*)"0.0"; // energy
                fixarg[9]  = (char*)"0.0"; // fx
                fixarg[10] = (char*)"0.0"; // fy
                fixarg[11] = (char*)"0.0"; // fz
                fixarg[12] = (char*)"0.0"; // tx
                fixarg[13] = (char*)"0.0"; // ty
                fixarg[14] = (char*)"0.0"; // tz
                fix_dissipated_ = modify->add_fix_property_atom(15, static_cast<char**>(fixarg), "dissipated energy");
            }
        }
        if (cmb->is_wall() && (dissipatedflag_ || elasticpotflag_))
        {
            delta_offset_ = hsetup->add_history_value("delta_0", "1");
            hsetup->add_history_value("delta_1", "1");
            hsetup->add_history_value("delta_2", "1");
            cmb->add_history_offset("delta", delta_offset_);
            if (dissipatedflag_)
            {
                dissipation_offset_ = hsetup->add_history_value("diss_f_0", "1");
                hsetup->add_history_value("diss_f_1", "1");
                hsetup->add_history_value("diss_f_2", "1");
                cmb->add_history_offset("dissipation_force", dissipation_offset_);
            }
        }
    }

    inline void connectToProperties(PropertyRegistry&) {}

    inline bool checkSurfaceIntersect(SurfacesIntersectData & sidata)
    {
        #ifdef SUPERQUADRIC_ACTIVE_FLAG
            sidata.is_non_spherical = false;
        #endif
        
        return true;
    }

    inline void surfacesIntersect(SurfacesIntersectData & sidata, ForceData&, ForceData&)
    {
      #ifdef SUPERQUADRIC_ACTIVE_FLAG
      if(sidata.is_non_spherical)
        error->one(FLERR,"Using default surface model for non-spherical particles!");
      #endif
      const double enx = sidata.en[0];
      const double eny = sidata.en[1];
      const double enz = sidata.en[2];
      const bool sync = sidata.synchronized_verlet && !sidata.is_wall;

      // relative translational velocity
      const double vr1 = sidata.v_i[0] - sidata.v_j[0];
      const double vr2 = sidata.v_i[1] - sidata.v_j[1];
      const double vr3 = sidata.v_i[2] - sidata.v_j[2];

      // relative rotational velocity
      const double deltan = sidata.radsum - sidata.r;
      const double dx = sidata.delta[0];
      const double dy = sidata.delta[1];
      const double dz = sidata.delta[2];
      const double rinv = sidata.rinv;
      double wr1, wr2, wr3;

      const double cri_local = sidata.radi - 0.5 * deltan;
      const double crj_local = sidata.radj - 0.5 * deltan;

      if(sidata.is_wall) {
        // in case of wall contact, r is the contact radius
        const double cr = cri_local;
        wr1 = cr * sidata.omega_i[0] * rinv;
        wr2 = cr * sidata.omega_i[1] * rinv;
        wr3 = cr * sidata.omega_i[2] * rinv;
        sidata.cri = cr;
      } else {
        wr1 = (cri_local * sidata.omega_i[0] + crj_local * sidata.omega_j[0]) * rinv;
        wr2 = (cri_local * sidata.omega_i[1] + crj_local * sidata.omega_j[1]) * rinv;
        wr3 = (cri_local * sidata.omega_i[2] + crj_local * sidata.omega_j[2]) * rinv;
        sidata.cri = cri_local;
        sidata.crj = crj_local;
      }

      if (sync) {
        double nhalf[3] = {sidata.en[0], sidata.en[1], sidata.en[2]};
        const double dt_half = 0.5 * update->dt;
        double dx_mid[3] = {dx - vr1 * dt_half, dy - vr2 * dt_half, dz - vr3 * dt_half};
        const double mid_len = vectorMag3D(dx_mid);
        if (mid_len > 0.0) {
          vectorScalarDiv3D(dx_mid, mid_len);
          vectorCopy3D(dx_mid, nhalf);
        }
        vectorCopy3D(nhalf, sidata.en_half);

        double contact_vec_i[3] = {-sidata.cri * nhalf[0], -sidata.cri * nhalf[1], -sidata.cri * nhalf[2]};
        double rot_i[3];
        vectorCross3D(sidata.omega_i, contact_vec_i, rot_i);

        double rot_j[3] = {0.0, 0.0, 0.0};
        if (!sidata.is_wall) {
          double contact_vec_j[3] = {sidata.crj * nhalf[0], sidata.crj * nhalf[1], sidata.crj * nhalf[2]};
          vectorCross3D(sidata.omega_j, contact_vec_j, rot_j);
        }

        double v_contact[3] = {vr1 + rot_i[0] - rot_j[0], vr2 + rot_i[1] - rot_j[1], vr3 + rot_i[2] - rot_j[2]};
        const double vn_half = vectorDot3D(v_contact, nhalf);
        sidata.vn = vn_half;
        double vt_half[3] = {
          v_contact[0] - vn_half * nhalf[0],
          v_contact[1] - vn_half * nhalf[1],
          v_contact[2] - vn_half * nhalf[2]
        };
        sidata.vtr_half1 = vt_half[0];
        sidata.vtr_half2 = vt_half[1];
        sidata.vtr_half3 = vt_half[2];

        double vtr_full[3] = {vt_half[0], vt_half[1], vt_half[2]};
        const double proj_full = vectorDot3D(vtr_full, sidata.en);
        vtr_full[0] -= proj_full * sidata.en[0];
        vtr_full[1] -= proj_full * sidata.en[1];
        vtr_full[2] -= proj_full * sidata.en[2];

        sidata.vtr1 = vtr_full[0];
        sidata.vtr2 = vtr_full[1];
        sidata.vtr3 = vtr_full[2];
      } else {
        // normal component
        const double vn = vr1 * enx + vr2 * eny + vr3 * enz;
        const double vn1 = vn * enx;
        const double vn2 = vn * eny;
        const double vn3 = vn * enz;

        // tangential component
        const double vt1 = vr1 - vn1;
        const double vt2 = vr2 - vn2;
        const double vt3 = vr3 - vn3;

        const double vtr1 = vt1 - (dz * wr2 - dy * wr3);
        const double vtr2 = vt2 - (dx * wr3 - dz * wr1);
        const double vtr3 = vt3 - (dy * wr1 - dx * wr2);
        sidata.vn = vn;
        sidata.vtr1 = vtr1;
        sidata.vtr2 = vtr2;
        sidata.vtr3 = vtr3;
        sidata.en_half[0] = enx;
        sidata.en_half[1] = eny;
        sidata.en_half[2] = enz;
        sidata.vtr_half1 = vtr1;
        sidata.vtr_half2 = vtr2;
        sidata.vtr_half3 = vtr3;
      }

      sidata.deltan = deltan;
      sidata.wr1 = wr1;
      sidata.wr2 = wr2;
      sidata.wr3 = wr3;
      sidata.P_diss = 0.;
    }

    inline void endSurfacesIntersect(SurfacesIntersectData &sidata,TriMesh *, double * const) {}
    inline void surfacesClose(SurfacesCloseData &scdata, ForceData&, ForceData&) {}
    void beginPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    void endPass(SurfacesIntersectData&, ForceData&, ForceData&){}
    inline void tally_pp(double,int,int,int) {}
    inline void tally_pw(double,int,int,int) {}

  private:
    bool elasticpotflag_;
    bool dissipatedflag_;
    int delta_offset_;
    int dissipation_offset_;
    FixPropertyAtom *fix_dissipated_;
  };
}
}
#endif // SURFACE_MODEL_DEFAULT_H_
#endif
