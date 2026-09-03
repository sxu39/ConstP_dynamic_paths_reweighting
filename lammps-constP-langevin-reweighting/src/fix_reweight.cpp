#include "fix_reweight.h"

#include "fix_nvtmu.h"

#include "atom.h"
#include "comm.h"
#include "compute.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

FixReweight::FixReweight(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix reweight", error);

  int ifix = modify->find_fix(arg[3]);
  if (ifix < 0) error->all(FLERR, "fix reweightmu: could not find fix nvtmu");
  fix_nvtmu = (FixNVTMU *) modify->fix[ifix];
  if (!fix_nvtmu) error->all(FLERR, "fix reweightmu: nvtmu fix pointer is null");

  num_remu = narg - 4;
  memory->create(remu_lst, num_remu, "fix_reweight:remu_lst");
  for (int i = 4; i < narg; i++) {
    remu_lst[i - 4] = utils::numeric(FLERR, arg[i], false, lmp);
  }
  memory->create(reweight_factor_lst, num_remu, "fix_reweight:reweight_factor_lst");
  memory->create(delta_rand_G_lst, num_remu, "fix_reweight:delta_rand_G_lst");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = num_remu;
}

/* ---------------------------------------------------------------------- */

FixReweight::~FixReweight()
{
  memory->destroy(remu_lst);
  memory->destroy(reweight_factor_lst);
  memory->destroy(delta_rand_G_lst);
}

/* ---------------------------------------------------------------------- */

int FixReweight::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixReweight::init()
{
  double mu_set = fix_nvtmu->get_Ne_mu();
  for (int i = 0; i < num_remu; i++) {
    reweight_factor_lst[i] = exp(atom->Ne[0] * (remu_lst[i] - mu_set) / (force->boltz * fix_nvtmu->get_temperature()));
  }

  double ftm2v = force->ftm2v;
  double mvv2e = force->mvv2e;
  double boltz = force->boltz;
  double dt = update->dt;

  double Ne_mass = fix_nvtmu->get_Ne_mass();
  double temperature = fix_nvtmu->get_temperature();
  double fric_exp_Ne = fix_nvtmu->get_fric_exp_Ne();

  double dtfm = 0.5 * dt * ftm2v / Ne_mass;
  double beta = sqrt((1 - fric_exp_Ne * fric_exp_Ne) * boltz * temperature / Ne_mass / mvv2e);
  double delta_mu;

  std::cout << "------------reweight result------------" << std::endl;

  for (int i = 0; i < num_remu; i++) {
    delta_mu = remu_lst[i] - mu_set;
    delta_rand_G_lst[i] = -(1 + fric_exp_Ne) * dtfm * delta_mu / beta;
    std::cout << "delta mu: " << delta_mu << " , delta rand_G: " << delta_rand_G_lst[i] << std::endl;
  }

  std::cout << "----------reweight parameters----------" << std::endl;
  std::cout << "ftm2v: " << ftm2v << std::endl;
  std::cout << "mvv2e: " << mvv2e << std::endl;
  std::cout << "boltz: " << boltz << std::endl;
  std::cout << "dt: " << dt << std::endl;
  std::cout << "Ne: " << atom->Ne[0] << std::endl;
  std::cout << "Ne_mass: " << Ne_mass << std::endl;
  std::cout << "boltz: " << force->boltz << std::endl;
  std::cout << "temperature: " << temperature << std::endl;
  std::cout << "fric_exp_Ne: " << fric_exp_Ne << std::endl;
  std::cout << "---------------------------------------" << std::endl;
}

/* ---------------------------------------------------------------------- */

void FixReweight::end_of_step()
{
  /*
  double ftm2v = force->ftm2v;
  double mvv2e = force->mvv2e;
  double dt = update->dt;
  double dtf = 0.5 * dt * ftm2v;

  double Ne_mass = fix_gjfmu->get_Ne_mass();
  double tsqrt = fix_gjfmu->get_tsqrt();
  double t_period_Ne = fix_gjfmu->get_t_period_Ne();
  double Ne_dot_init = fix_gjfmu->get_Ne_dot_init();
  double Ne_dot_final = atom->Ne_dot[0];
  double Ne_rand_G = fix_gjfmu->get_Ne_rand_G();
  double *dE_dN = atom->dE_dN;
  double gjfc2_Ne = fix_gjfmu->get_gjfc2_Ne();
  double csq_Ne = fix_gjfmu->get_csq_Ne();

  double beta = tsqrt * sqrt(2.0 * dt * Ne_mass * force->boltz / t_period_Ne / mvv2e) / ftm2v;
  double dtfm = dtf / Ne_mass;
  double Ne_dot_init_after, Ne_dot_final_before, Ne_rand_G_new;
  double reweight_factor_check;

  std::cout << "beta: " << beta << std::endl;

  for (int i = 0; i < num_remu; i++) {
    // reweight factor calculation
    Ne_dot_init_after = Ne_dot_init + csq_Ne * dtfm * (-dE_dN[0] + remu_lst[i]);
    Ne_dot_final_before = Ne_dot_final - csq_Ne * dtfm * (-dE_dN[0] + remu_lst[i]);
    Ne_rand_G_new = (Ne_dot_final_before - gjfc2_Ne * Ne_dot_init_after) / ((1 + gjfc2_Ne) * (0.5 * csq_Ne / Ne_mass) * ftm2v * beta);
    reweight_factor_lst[i] = exp(-(Ne_rand_G_new * Ne_rand_G_new - Ne_rand_G * Ne_rand_G) / 2.0);
    std::cout << "delta rand_G: " << Ne_rand_G_new - Ne_rand_G << std::endl;
  }
  */

  double Ne_rand_G = fix_nvtmu->get_Ne_rand_G();
  for (int i = 0; i < num_remu; i++) {
    reweight_factor_lst[i] = exp(-(delta_rand_G_lst[i] * (delta_rand_G_lst[i] + 2 * Ne_rand_G)) / 2.0);
  }
}

/* ---------------------------------------------------------------------- */

double FixReweight::compute_scalar()
{
  return fix_nvtmu->get_Ne_rand_G();
}

/* ---------------------------------------------------------------------- */

double FixReweight::compute_vector(int n)
{
  if (n < 0 || n >= num_remu) {
    error->all(FLERR, "Fix reweight: compute_vector index out of bounds");
  }
  return reweight_factor_lst[n];
}

