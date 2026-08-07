/* ----------------------------------------------------------------------
    This is the

    LIGGGHTS(R) - DEM simulation engine
    https://www.cfdem.com

    LIGGGHTS(R) is open-source, distributed under the terms of the GNU
    General Public License, version 2 or later.

    The unequal-sphere lubrication resistance model follows:

      O. Cheal and C. Ness, "Rheology of dense granular suspensions under
      extensional flow", arXiv:1801.03805v3.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstring>

#include "pair_lubricate_bmpoly.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "fix_deform.h"
#include "fix_wall.h"
#include "force.h"
#include "input.h"
#include "math_const.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "variable.h"

using namespace LAMMPS_NS;
using namespace MathConst;

enum { NO_REMAP, X_REMAP, V_REMAP };
enum { EDGE, CONSTANT, VARIABLE };

/* ---------------------------------------------------------------------- */

PairLubricateBmpoly::PairLubricateBmpoly(LAMMPS *lmp) : PairLubricate(lmp)
{
  no_virial_fdotr_compute = 1;
}

/* ----------------------------------------------------------------------
   Update the background-fluid resistance coefficients.  The volume
   fraction is recalculated when the box or a supported wall moves.
------------------------------------------------------------------------- */

void PairLubricateBmpoly::update_fld_coefficients()
{
  double vol_T;

  if (!flagwall) {
    vol_T = domain->xprd * domain->yprd * domain->zprd;
  } else {
    double wallhi[3],walllo[3];
    walllo[0] = domain->boxlo[0];
    walllo[1] = domain->boxlo[1];
    walllo[2] = domain->boxlo[2];
    wallhi[0] = domain->boxhi[0];
    wallhi[1] = domain->boxhi[1];
    wallhi[2] = domain->boxhi[2];

    for (int m = 0; m < wallfix->nwall; m++) {
      int dim = wallfix->wallwhich[m] / 2;
      int side = wallfix->wallwhich[m] % 2;
      double wallcoord;

      if (wallfix->xstyle[m] == VARIABLE) {
        wallfix->xindex[m] = input->variable->find(wallfix->xstr[m]);
        if (wallfix->xindex[m] < 0)
          error->all(FLERR,"Variable name for lubrication wall does not exist");
        wallcoord = input->variable->compute_equal(wallfix->xindex[m]);
      } else {
        wallcoord = wallfix->coord0[m];
      }

      if (side == 0) walllo[dim] = wallcoord;
      else wallhi[dim] = wallcoord;
    }

    vol_T = (wallhi[0] - walllo[0]) * (wallhi[1] - walllo[1]) *
      (wallhi[2] - walllo[2]);
  }

  if (vol_T <= 0.0)
    error->all(FLERR,"Pair lubricate/bmpoly requires a positive fluid volume");

  double vol_f = 0.0;
  if (flagVF) vol_f = vol_P / vol_T;

  if (flaglog == 0) {
    R0 = 6.0*MY_PI*mu*(1.0 + 2.16*vol_f);
    RT0 = 8.0*MY_PI*mu;
    RS0 = 20.0/3.0*MY_PI*mu *
      (1.0 + 3.33*vol_f + 2.80*vol_f*vol_f);
  } else {
    R0 = 6.0*MY_PI*mu *
      (1.0 + 2.725*vol_f - 6.583*vol_f*vol_f);
    RT0 = 8.0*MY_PI*mu *
      (1.0 + 0.749*vol_f - 2.469*vol_f*vol_f);
    RS0 = 20.0/3.0*MY_PI*mu *
      (1.0 + 3.64*vol_f - 6.95*vol_f*vol_f);
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLubricateBmpoly::init_style()
{
  PairLubricate::init_style();

  if (force->newton_pair)
    error->all(FLERR,"Pair lubricate/bmpoly requires newton pair off");
  if (!atom->sphere_flag)
    error->all(FLERR,"Pair lubricate/bmpoly requires atom style granular or sphere");

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  shearing = flagdeform = flagwall = 0;
  wallfix = NULL;

  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style,"deform") == 0) {
      shearing = flagdeform = 1;
      if (((FixDeform *) modify->fix[i])->remapflag != V_REMAP)
        error->all(FLERR,"Using pair lubricate/bmpoly with inconsistent "
                   "fix deform remap option");
    }

    // The legacy lubrication model can use the volume enclosed by the
    // original fix wall command.  Granular wall fixes have a different
    // interface and are intentionally not cast to FixWall here.
    if (strcmp(modify->fix[i]->style,"wall") == 0) {
      if (flagwall)
        error->all(FLERR,"Cannot use multiple fix wall commands with "
                   "pair lubricate/bmpoly");
      flagwall = 1;
      wallfix = (FixWall *) modify->fix[i];
      if (wallfix->xflag) flagwall = 2;
    }
  }

  double volP = 0.0;
  for (int i = 0; i < atom->nlocal; i++)
    volP += 4.0/3.0*MY_PI*atom->radius[i]*atom->radius[i]*atom->radius[i];
  MPI_Allreduce(&volP,&vol_P,1,MPI_DOUBLE,MPI_SUM,world);

  update_fld_coefficients();

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Ef[i][j] = 0.0;
}

/* ---------------------------------------------------------------------- */

void PairLubricateBmpoly::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  if (flagVF && (flagdeform || flagwall == 2))
    update_fld_coefficients();

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double vxmu2f = force->vxmu2f;

  double wf[3] = {0.0,0.0,0.0};
  double *h_rate = domain->h_rate;
  double *h_ratelo = domain->h_ratelo;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      Ef[i][j] = 0.0;

  if (shearing) {
    Ef[0][0] = h_rate[0]/domain->xprd;
    Ef[1][1] = h_rate[1]/domain->yprd;
    Ef[2][2] = h_rate[2]/domain->zprd;
    Ef[0][1] = Ef[1][0] = 0.5*h_rate[5]/domain->yprd;
    Ef[0][2] = Ef[2][0] = 0.5*h_rate[4]/domain->zprd;
    Ef[1][2] = Ef[2][1] = 0.5*h_rate[3]/domain->zprd;

    wf[0] = -0.5*h_rate[3]/domain->zprd;
    wf[1] =  0.5*h_rate[4]/domain->zprd;
    wf[2] = -0.5*h_rate[5]/domain->yprd;
  }

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    int itype = type[i];
    double radi = radius[i];

    if (flagfld) {
      double vstream[3] = {0.0,0.0,0.0};

      if (shearing) {
        double lamda[3];
        domain->x2lamda(x[i],lamda);
        vstream[0] = h_rate[0]*lamda[0] + h_rate[5]*lamda[1] +
          h_rate[4]*lamda[2] + h_ratelo[0];
        vstream[1] = h_rate[1]*lamda[1] + h_rate[3]*lamda[2] +
          h_ratelo[1];
        vstream[2] = h_rate[2]*lamda[2] + h_ratelo[2];
      }

      f[i][0] -= vxmu2f*R0*radi*(v[i][0] - vstream[0]);
      f[i][1] -= vxmu2f*R0*radi*(v[i][1] - vstream[1]);
      f[i][2] -= vxmu2f*R0*radi*(v[i][2] - vstream[2]);
      torque[i][0] -= vxmu2f*RT0*radi*radi*radi*(omega[i][0] - wf[0]);
      torque[i][1] -= vxmu2f*RT0*radi*radi*radi*(omega[i][1] - wf[1]);
      torque[i][2] -= vxmu2f*RT0*radi*radi*radi*(omega[i][2] - wf[2]);

      if (shearing && vflag_either) {
        double vRS0 = -vxmu2f*RS0*radi*radi*radi;
        v_tally_tensor(i,i,nlocal,newton_pair,
                       vRS0*Ef[0][0],vRS0*Ef[1][1],vRS0*Ef[2][2],
                       vRS0*Ef[0][1],vRS0*Ef[0][2],vRS0*Ef[1][2]);
      }
    }

    if (!flagHI) continue;

    int jnum = numneigh[i];
    int *jlist = firstneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      double delx = x[i][0] - x[j][0];
      double dely = x[i][1] - x[j][1];
      double delz = x[i][2] - x[j][2];
      double rsq = delx*delx + dely*dely + delz*delz;
      if (rsq <= 0.0) continue;

      double r = sqrt(rsq);
      double radj = radius[j];
      double h_sep = r - radi - radj;
      if (h_sep >= cut[itype][type[j]]) continue;
      if (h_sep < cut_inner[itype][type[j]])
        h_sep = cut_inner[itype][type[j]];

      // n points from particle 2 to particle 1, as in the paper.
      double nx = delx/r;
      double ny = dely/r;
      double nz = delz/r;
      double beta = radj/radi;
      double betai = radi/radj;
      double xi = 2.0*h_sep/(radi + radj);
      double loggap = flaglog ? log(1.0/xi) : 0.0;

      double XA11 = 6.0*MY_PI*mu*radi *
        (2.0*beta*beta/(1.0+beta)/(1.0+beta)/(1.0+beta)/xi +
         beta*(1.0+7.0*beta+beta*beta)/
         (5.0*(1.0+beta)*(1.0+beta)*(1.0+beta))*loggap);
      double YA11 = 6.0*MY_PI*mu*radi *
        (4.0*beta*(2.0+beta+2.0*beta*beta)/
         (15.0*(1.0+beta)*(1.0+beta)*(1.0+beta))*loggap);
      double YB11 = -4.0*MY_PI*mu*radi*radi *
        (beta*(4.0+beta)/(5.0*(1.0+beta)*(1.0+beta))*loggap);
      double YB21 = -4.0*MY_PI*mu*radj*radj *
        (betai*(4.0+betai)/(5.0*(1.0+betai)*(1.0+betai))*loggap);
      double YC11 = 8.0*MY_PI*mu*radi*radi*radi *
        (2.0*beta/(5.0*(1.0+beta))*loggap);
      double YC12 = 8.0*MY_PI*mu*radi*radi*radi *
        (beta/(10.0*(1.0+beta))*loggap);

      double vr1 = v[j][0] - v[i][0];
      double vr2 = v[j][1] - v[i][1];
      double vr3 = v[j][2] - v[i][2];
      double vnnr = vr1*nx + vr2*ny + vr3*nz;
      double vn1 = vnnr*nx;
      double vn2 = vnnr*ny;
      double vn3 = vnnr*nz;
      double vt1 = vr1 - vn1;
      double vt2 = vr2 - vn2;
      double vt3 = vr3 - vn3;

      double w1cn1 = omega[i][1]*nz - omega[i][2]*ny;
      double w1cn2 = omega[i][2]*nx - omega[i][0]*nz;
      double w1cn3 = omega[i][0]*ny - omega[i][1]*nx;
      double w2cn1 = omega[j][1]*nz - omega[j][2]*ny;
      double w2cn2 = omega[j][2]*nx - omega[j][0]*nz;
      double w2cn3 = omega[j][0]*ny - omega[j][1]*nx;

      double f1x = vxmu2f*(XA11*vn1 + YA11*vt1 +
                          YB11*w1cn1 + YB21*w2cn1);
      double f1y = vxmu2f*(XA11*vn2 + YA11*vt2 +
                          YB11*w1cn2 + YB21*w2cn2);
      double f1z = vxmu2f*(XA11*vn3 + YA11*vt3 +
                          YB11*w1cn3 + YB21*w2cn3);
      f[i][0] += f1x;
      f[i][1] += f1y;
      f[i][2] += f1z;

      double vrcrossn1 = vr2*nz - vr3*ny;
      double vrcrossn2 = vr3*nx - vr1*nz;
      double vrcrossn3 = vr1*ny - vr2*nx;
      double wnni = omega[i][0]*nx + omega[i][1]*ny + omega[i][2]*nz;
      double wnnj = omega[j][0]*nx + omega[j][1]*ny + omega[j][2]*nz;
      double wti1 = omega[i][0] - wnni*nx;
      double wti2 = omega[i][1] - wnni*ny;
      double wti3 = omega[i][2] - wnni*nz;
      double wtj1 = omega[j][0] - wnnj*nx;
      double wtj2 = omega[j][1] - wnnj*ny;
      double wtj3 = omega[j][2] - wnnj*nz;

      double g1x = vxmu2f*(YB11*vrcrossn1 - YC11*wti1 - YC12*wtj1);
      double g1y = vxmu2f*(YB11*vrcrossn2 - YC11*wti2 - YC12*wtj2);
      double g1z = vxmu2f*(YB11*vrcrossn3 - YC11*wti3 - YC12*wtj3);
      torque[i][0] += g1x;
      torque[i][1] += g1y;
      torque[i][2] += g1z;

      if (vflag_either)
        v_tally_tensor(i,nlocal,nlocal,newton_pair,
                       f1x*delx,f1y*dely,f1z*delz,
                       0.5*(f1x*dely + f1y*delx),
                       0.5*(f1x*delz + f1z*delx),
                       0.5*(f1y*delz + f1z*dely));
    }
  }
}
