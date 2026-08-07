/* ----------------------------------------------------------------------
    This is the

    LIGGGHTS(R) - DEM simulation engine
    https://www.cfdem.com

    LIGGGHTS(R) is open-source, distributed under the terms of the GNU
    Public License, version 2 or later.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

#else

#ifndef LMP_PAIR_LUBRICATE_H
#define LMP_PAIR_LUBRICATE_H

#include <stdio.h>

#include "pair.h"

namespace LAMMPS_NS {

class PairLubricate : public Pair {
 public:
  PairLubricate(class LAMMPS *);
  virtual ~PairLubricate();

  virtual void compute(int, int) = 0;
  void settings(int, char **);
  void coeff(int, char **);
  virtual void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);

  void write_restart(FILE *);
  void read_restart(FILE *, const int, const int);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *, const int, const int);

 protected:
  double mu;
  int flaglog;
  int flagfld;
  int flagHI;
  int flagVF;
  double cut_inner_global;
  double cut_global;

  double **cut;
  double **cut_inner;

  double vol_P;
  double R0;
  double RT0;
  double RS0;
  double Ef[3][3];

  int shearing;
  int flagdeform;
  int flagwall;
  class FixWall *wallfix;

  void allocate();
  double max_radius_for_type(int);
};

}

#endif
#endif
