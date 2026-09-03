#ifdef FIX_CLASS
// clang-format off
FixStyle(toy,FixToy);
// clang-format on
#else

#ifndef LMP_FIX_TOY_H
#define LMP_FIX_TOY_H
#include "fix.h"
namespace LAMMPS_NS {
class FixToy : public Fix {
 public:
  FixToy(class LAMMPS *, int, char **);
  ~FixToy() override = default;
  int setmask() override;
  void pre_force(int) override;
  double compute_scalar() override;

 protected:
  double omega;
  double g;
  double Ed;
  double H;
  double delta;
  double Uext;
};
}
#endif
#endif