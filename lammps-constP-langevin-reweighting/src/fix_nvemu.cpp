// clang-format off
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

#include "fix_nvemu.h"

#include "atom.h"
#include "force.h"
#include "respa.h"
#include "update.h"
#include "error.h"

#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVEMU::FixNVEMU(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), Ne_mass(nullptr)
{
  if (!utils::strmatch(style,"^nvemu/sphere") && narg < 3)
    utils::missing_cmd_args(FLERR, "fix nvemu", error);

  dynamic_group_allow = 1;
  time_integrate = 1;
  ustat_flag = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"mu") == 0) {   //?
      ustat_flag = 1;
      u_target = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else error->all(FLERR,"Unknown fix {} keyword: {}", style, arg[iarg]);
  }

  if (ustat_flag) {  //?
    Ne_mass = new double[1];
    *Ne_mass = 16.0;
  }

}

/* ---------------------------------------------------------------------- */

FixNVEMU::~FixNVEMU()
{
  if (ustat_flag) {  //?
    delete [] Ne_mass;
  }
}

/* ---------------------------------------------------------------------- */

int FixNVEMU::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVEMU::init()
{
  dtv = 1.0 * update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (utils::strmatch(update->integrate_style,"^respa"))
    step_respa = (dynamic_cast<Respa *>(update->integrate))->step;
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVEMU::initial_integrate(int /*vflag*/)
{
  double dtfm;

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *Ne = atom->Ne;  
  double *Ne_dot = atom->Ne_dot;
  double *dE_dN = atom->dE_dN;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

//  int test_id = 152; //?
//  std::cout << "----------nvemu initial_integrate()----------" << std::endl;
//  std::cout << "v before: " << v[test_id][0] << ", " << v[test_id][1] << ", " << v[test_id][2] << std::endl;
//  std::cout << "x before: " << x[test_id][0] << ", " << x[test_id][1] << ", " << x[test_id][2] << std::endl;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }

//  std::cout << "dtfm: " << dtf/mass[type[test_id]] << std::endl;
//  std::cout << "f: " << f[test_id][0] << ", " << f[test_id][1] << ", " << f[test_id][2] << std::endl;

//  std::cout << "v after: " << v[test_id][0] << ", " << v[test_id][1] << ", " << v[test_id][2] << std::endl;
//  std::cout << "x after: " << x[test_id][0] << ", " << x[test_id][1] << ", " << x[test_id][2] << std::endl;

  if (ustat_flag) {
//    std::cout << "Ne before: " << *Ne << std::endl;
//    std::cout << "Ne_dot before: " << *Ne_dot << std::endl;

    dtfm = dtf / *Ne_mass;
    *Ne_dot += dtfm * (-(*dE_dN) + u_target);  //?
    *Ne += dtv * (*Ne_dot); 

//    std::cout << "dtfm: " << dtfm << std::endl;
//    std::cout << "dtv: " << dtv << std::endl;

//    std::cout << "dE_dN: " << *dE_dN << std::endl;
//    std::cout << "u_target: " << u_target << std::endl;

//    std::cout << "Ne after: " << *Ne << std::endl;
//    std::cout << "Ne_dot after: " << *Ne_dot << std::endl;
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEMU::final_integrate()
{
  double dtfm;

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double *Ne = atom->Ne;  
  double *Ne_dot = atom->Ne_dot;
  double *dE_dN = atom->dE_dN;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        //x[i][0] += dtv * v[i][0];
        //x[i][1] += dtv * v[i][1];
        //x[i][2] += dtv * v[i][2];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
  }

  if (ustat_flag) {

//    std::cout << "----------nvemu final_integrate()----------" << std::endl;

//    std::cout << "Ne before: " << *Ne << std::endl;
//    std::cout << "Ne_dot before: " << *Ne_dot << std::endl;

    dtfm = dtf / *Ne_mass;
    *Ne_dot += dtfm * (-(*dE_dN) + u_target);  //?

//    std::cout << "dtfm: " << dtfm << std::endl;
//    std::cout << "dE_dN: " << *dE_dN << std::endl;
//    std::cout << "u_target: " << u_target << std::endl;
    
//    std::cout << "Ne after: " << *Ne << std::endl;
//    std::cout << "Ne_dot after: " << *Ne_dot << std::endl;
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEMU::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVEMU update of v and x
  // all other levels - NVEMU update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVEMU::final_integrate_respa(int ilevel, int /*iloop*/)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVEMU::reset_dt()
{
  dtv = 1.0 * update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
