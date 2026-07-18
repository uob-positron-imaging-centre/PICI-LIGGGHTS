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

    Dominik Werner

    Copyright 2026-     Dominik Werner
------------------------------------------------------------------------- */

#ifndef JKR_SOLVERS_H_
#define JKR_SOLVERS_H_

#include <cmath>
#include <stdio.h>

#include "comm.h"
#include "global_properties.h"
#include "jkr_lookup_table.h"
#include "pointers.h"
#include "property_registry.h"

namespace MODEL_PARAMS
{
  static ScalarProperty* createResolutionJKR(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createScalarProperty(registry, "resolutionJKR", caller);
  }

  static MatrixProperty* createWorkOfAdhesion(PropertyRegistry & registry, const char * caller, bool sanity_checks)
  {
    return createPerTypePairProperty(registry, "workOfAdhesion", caller);
  }
}

namespace LIGGGHTS {
namespace ContactModels {

struct JKRResult
{
  double a_a0;
  double fn_fc;
};

class JKRLookupSolver : protected Pointers
{
public:
  JKRLookupSolver(LAMMPS * lmp) :
    Pointers(lmp),
    resolutionJKR(0.)
  {}

  void connectToProperties(PropertyRegistry & registry)
  {
    registry.registerProperty("resolutionJKR", &MODEL_PARAMS::createResolutionJKR, "model jkr");
    registry.connect("resolutionJKR", resolutionJKR, "model jkr");

    lut.init(resolutionJKR);

    if (comm->me == 0) {
      if (screen)
        fprintf(screen,"JKR table created with:\nResolution = %g \nMaximum(delta_delta_c) = %g\nPoints = %i\n",
                resolutionJKR,lut.get_max_delta_delta_c(),lut.get_nr_of_points());
      if (logfile)
        fprintf(logfile,"JKR table created with:\nResolution = %g \nMaximum(delta_delta_c) = %g\nPoints = %i\n",
                resolutionJKR,lut.get_max_delta_delta_c(),lut.get_nr_of_points());
    }
  }

  JKRResult evaluate(double delta_delta_c)
  {
    JKRResult result;
    result.fn_fc = lut.get_fn_fc(delta_delta_c);
    result.a_a0 = lut.get_a_a0(delta_delta_c);
    return result;
  }

  double calculate_a0(double Wij, double reff, double Yeff)
  {
    return lut.calculate_a0(Wij, reff, Yeff);
  }

  double calculate_delta_c(double a_0, double reff)
  {
    return lut.calculate_delta_c(a_0, reff);
  }

  double calculate_fc(double Wij, double reff)
  {
    return lut.calculate_fc(Wij, reff);
  }

private:
  LUT lut;
  double resolutionJKR;
};

class JKRAnalyticalSolver
{
public:
  JKRAnalyticalSolver(LAMMPS *) {}

  void connectToProperties(PropertyRegistry &) {}

  JKRResult evaluate(double delta_delta_c) const
  {
    JKRResult result;
    if (delta_delta_c < -1.) {
      result.a_a0 = 0.;
      result.fn_fc = 0.;
      return result;
    }

    const double cbrt_six = 1.8171205928321396588912117563272605;
    double z;

    if (delta_delta_c <= 0.) {
      const double discriminant = fmax(0., 1. + delta_delta_c * delta_delta_c * delta_delta_c);
      const double sqrt_discriminant = sqrt(discriminant);
      const double u = cbrt(2. * (1. + sqrt_discriminant) / 9.);
      const double v = cbrt(2. * (1. - sqrt_discriminant) / 9.);
      z = u + v;
    } else {
      const double p = 2. * delta_delta_c / cbrt_six;
      const double q = -4. / 9.;
      const double argument = (-3. * q / (2. * p)) * sqrt(3. / p);
      z = 2. * sqrt(p / 3.) * sinh(asinh(argument) / 3.);
    }

    const double sqrt_z = sqrt(fmax(0., z));
    const double second_radicand = fmax(0., -z + 4. / (3. * sqrt_z));
    const double g = 0.5 * (sqrt_z + sqrt(second_radicand));
    result.a_a0 = g * g;
    result.fn_fc = 4. * (result.a_a0 * result.a_a0 * result.a_a0 -
                         result.a_a0 * g);
    return result;
  }

  double calculate_a0(double Wij, double reff, double Yeff) const
  {
    return pow(9. / 2. * Wij * M_PI * reff * reff / Yeff, 1. / 3.);
  }

  double calculate_delta_c(double a_0, double reff) const
  {
    return a_0 * a_0 / (2. * cbrt(6.) * reff);
  }

  double calculate_fc(double Wij, double reff) const
  {
    return 3. / 2. * M_PI * Wij * reff;
  }
};

}
}

#endif
