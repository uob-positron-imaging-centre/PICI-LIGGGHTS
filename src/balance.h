/* ----------------------------------------------------------------------
    PICI-LIGGGHTS - particle based load balancing support

    Copyright 2026 PICI-LIGGGHTS contributors

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License, version 2 or later.
------------------------------------------------------------------------- */

#ifndef LMP_BALANCE_H
#define LMP_BALANCE_H

#include "pointers.h"

namespace LAMMPS_NS {

class Balance : protected Pointers {
 public:
  Balance(class LAMMPS *);
  ~Balance() {}

  void shift_setup(const char *, int, double);
  double shift(int &);
  double imbalance_factor(double &);
  double imbalance_splits(double &);

 private:
  char dimstr[4];
  int ndim, nitermax;
  double stopthresh;

  double reduced_coord(int, int);
  int split_index(double, int, double *);
  void cumulative_counts(int, const double *, int, bigint *);
  void validate_splits(int, const double *);
};

}

#endif
