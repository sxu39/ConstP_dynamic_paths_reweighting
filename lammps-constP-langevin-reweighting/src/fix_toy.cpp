#include "fix_toy.h"
#include "atom.h"
#include "error.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "memory.h"
#include <cmath>
#include <iostream>

using namespace LAMMPS_NS;

FixToy::FixToy(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg) {
  if (narg != 8) error->all(FLERR,"fix toy usage: fix ID group toy omega g Ed H delta");
  omega = utils::numeric(FLERR,arg[3],false,lmp);
  g = utils::numeric(FLERR,arg[4],false,lmp);
  Ed = utils::numeric(FLERR,arg[5],false,lmp);
  H = utils::numeric(FLERR,arg[6],false,lmp);
  delta = utils::numeric(FLERR,arg[7],false,lmp);
  Uext = 0.0;
  // enable reporting of a global scalar energy
  scalar_flag = 1;
  extscalar = 1;
  energy_global_flag = 1;
  thermo_energy = 1;
}

int FixToy::setmask() {
  int mask = 0;
  mask |= FixConst::PRE_FORCE;
  return mask;
}

void FixToy::pre_force(int vflag) {
  Uext = 0.0;
  double **x = atom->x;
  double **f = atom->f;
  double *mass = atom->mass;
  double *Ne = atom->Ne;
  double *dE_dN = atom->dE_dN;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    double m = mass[type[i]];
    double x0 = x[i][0];
    double u = 0.5 * m * omega * omega * x0 * x0 + (*Ne) * (sqrt(2.0) * g * x0 + Ed) + H / (1 + exp((*Ne)/delta)) + H / (1 + exp((-*Ne+1)/delta));
    double f0 = -m * omega * omega * x0 - (*Ne) * sqrt(2.0) * g;
    double dEdN = sqrt(2.0) * g * x0 + Ed - H * exp((*Ne)/delta) / (delta * pow(1 + exp((*Ne)/delta), 2)) + H * exp((-*Ne+1)/delta) / (delta * pow(1 + exp((-*Ne+1)/delta), 2));

    f[i][0] = f0;
    (*dE_dN) = dEdN;
    Uext += u;
  }
}

double FixToy::compute_scalar() {
  return Uext;
}