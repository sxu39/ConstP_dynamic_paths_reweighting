/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(nvtmu,FixNVTMU);
// clang-format on
#else

#ifndef LMP_FIX_NVTMU_H
#define LMP_FIX_NVTMU_H

#include "fix.h"

namespace LAMMPS_NS {

class FixNVTMU : public Fix {
 public:
  FixNVTMU(class LAMMPS *, int, char **);
  ~FixNVTMU() override;
  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;

  // for reweighting
  double get_Ne_mu() { return mu_target; }
  double get_temperature() { return temperature; }
  double get_fric_exp_Ne() { return fric_exp_Ne; }
  double get_Ne_mass() { return Ne_mass; }
  double get_Ne_rand_G() { return Ne_rand_G; }

 protected:
  double temperature, t_period, t_period_Ne, mu_target;
  double fric_exp, fric_exp_Ne;
  double Ne_mass;

  double Ne_rand_G; // for reweighting

  class RanMars *random;
  int seed;
};

}    // namespace LAMMPS_NS

#endif
#endif
