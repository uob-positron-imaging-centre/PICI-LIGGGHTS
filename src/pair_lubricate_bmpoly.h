/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lubricate/bmpoly,PairLubricateBmpoly)

#else

#ifndef LMP_PAIR_LUBRICATE_BMPOLY_H
#define LMP_PAIR_LUBRICATE_BMPOLY_H

#include "pair_lubricate.h"

namespace LAMMPS_NS {

class PairLubricateBmpoly : public PairLubricate {
 public:
  PairLubricateBmpoly(class LAMMPS *);
  ~PairLubricateBmpoly() {}
  void compute(int, int);
  void init_style();

 protected:
  void update_fld_coefficients();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Pair lubricate/bmpoly requires newton pair off

Self-explanatory.

E: Pair lubricate/bmpoly requires ghost atoms store velocity

Use the comm_modify vel yes command to enable this.

E: Pair lubricate/bmpoly requires atom style granular or sphere

Self-explanatory.

E: Pair lubricate/bmpoly requires positive particle radii

One of the particles has radius 0.0.

E: Using pair lubricate with inconsistent fix deform remap option

Must use remap v option with fix deform with this pair style.

E: Cannot use multiple fix wall commands with pair lubricate/poly

Self-explanatory.

E: Using pair lubricate/bmpoly with inconsistent fix deform remap option

If fix deform is used, the remap v option is required.

*/
