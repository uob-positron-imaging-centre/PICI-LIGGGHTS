/* ----------------------------------------------------------------------
    PICI-LIGGGHTS - dynamic particle-count load balancing

    Copyright 2026 PICI-LIGGGHTS contributors

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License, version 2 or later.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(balance,FixBalance)

#else

#ifndef LMP_FIX_BALANCE_H
#define LMP_FIX_BALANCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBalance : public Fix {
 public:
  FixBalance(class LAMMPS *, int, char **);
  ~FixBalance();

  int setmask();
  void init();
  void setup_pre_exchange();
  void pre_exchange();
  void pre_neighbor();
  void reset_timestep(bigint);
  double compute_scalar();
  double compute_vector(int);

 private:
  int nfreq, nitermax;
  char dimstr[4];
  double thresh, stopthresh;
  double imbalance, imbalance_before, maxload;
  int iterations, pending;
  bigint lastcheck;

  class Balance *balance;
  class Irregular *irregular;

  void check_balance();
  void rebalance();
  void schedule_next();
};

}

#endif
#endif
