#ifdef FIX_CLASS
// clang-format off
FixStyle(reweight,FixReweight);
// clang-format on
#else

#ifndef LMP_FIX_REWEIGHT_H
#define LMP_FIX_REWEIGHT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixReweight : public Fix {
 public:
  FixReweight(class LAMMPS *, int, char **);
  ~FixReweight() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;

  double compute_scalar() override;
  double compute_vector(int) override;

 protected:
  int num_remu;
  double *remu_lst;
  double *delta_rand_G_lst;
  double *reweight_factor_lst;
  class FixNVTMU *fix_nvtmu;
};

}    // namespace LAMMPS_NS

#endif
#endif
