/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Tim Linke & Niels Gronbech-Jensen (UC Davis)
------------------------------------------------------------------------- */

#include "fix_nvtmu.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "random_mars.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

FixNVTMU::FixNVTMU(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), random(nullptr), Ne_rand_G(std::nan(""))
{
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix nvtmu", error);

  time_integrate = 1;

  temperature = utils::numeric(FLERR, arg[3], false, lmp);
  t_period = utils::numeric(FLERR, arg[4], false, lmp);
  t_period_Ne = utils::numeric(FLERR, arg[5], false, lmp);
  seed = utils::inumeric(FLERR, arg[6], false, lmp);
  mu_target = utils::numeric(FLERR, arg[7], false, lmp);

  if (t_period <= 0.0) error->all(FLERR, "Fix nvtmu period must be > 0.0");
  if (t_period_Ne <= 0.0) error->all(FLERR, "Fix nvtmu Ne period must be > 0.0");
  if (seed <= 0) error->all(FLERR, "Illegal fix nvtmu command");

  // initialize Marsaglia RNG with processor-unique seed
  random = new RanMars(lmp, seed + comm->me);

  // optional args
  Ne_mass = -1.0;

  int iarg = 8;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "mass_Ne") == 0) {
      Ne_mass = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else error->all(FLERR, "Illegal fix nvtmu command");
  }

  if (Ne_mass < 0.0) Ne_mass = 1.0;
}

/* ---------------------------------------------------------------------- */

FixNVTMU::~FixNVTMU()
{
  delete random;
}

/* ---------------------------------------------------------------------- */

int FixNVTMU::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVTMU::init()
{
  double dt = update->dt;
  fric_exp = exp(-dt / t_period);
  fric_exp_Ne = exp(-dt / t_period_Ne);
}

/* ---------------------------------------------------------------------- */

void FixNVTMU::initial_integrate(int )
{
  double **x = atom->x;
  double **v = atom->v;
  double *Ne = atom->Ne;
  double *Ne_dot = atom->Ne_dot;
  double *mass = atom->mass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double dt = update->dt;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] += 0.5 * dt * v[i][0];
      x[i][1] += 0.5 * dt * v[i][1];
      x[i][2] += 0.5 * dt * v[i][2];
    }
  }

  Ne[0] += 0.5 * dt * Ne_dot[0];
}

void FixNVTMU::final_integrate()
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *Ne = atom->Ne;
  double *Ne_dot = atom->Ne_dot;
  double *dE_dN = atom->dE_dN;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double rand[3];

  double boltz = force->boltz;
  double dt = update->dt;
  double mvv2e = force->mvv2e;
  double ftm2v = force->ftm2v;

  double dtf = 0.5 * dt * ftm2v;
  double dtfm;
  double m, beta;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      m = mass[type[i]];
      beta = sqrt((1 - fric_exp * fric_exp) * boltz * temperature / m / mvv2e);

      dtfm = dtf / m;
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      rand[0] = random->gaussian();
      rand[1] = random->gaussian();
      rand[2] = random->gaussian();

      v[i][0] = fric_exp * v[i][0] + beta * rand[0];
      v[i][1] = fric_exp * v[i][1] + beta * rand[1];
      v[i][2] = fric_exp * v[i][2] + beta * rand[2];

      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      x[i][0] += 0.5 * dt * v[i][0];
      x[i][1] += 0.5 * dt * v[i][1];
      x[i][2] += 0.5 * dt * v[i][2];
    }
  }

  beta = sqrt((1 - fric_exp_Ne * fric_exp_Ne) * boltz * temperature / Ne_mass / mvv2e);
  Ne_rand_G = random->gaussian();

  dtfm = dtf / Ne_mass;
  Ne_dot[0] += dtfm * (-dE_dN[0] + mu_target);
  Ne_dot[0] = fric_exp_Ne * Ne_dot[0] + beta * Ne_rand_G;
  Ne_dot[0] += dtfm * (-dE_dN[0] + mu_target);

  Ne[0] += 0.5 * dt * Ne_dot[0];
}