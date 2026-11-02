/* -*- c++ -*- ----------------------------------------------------------
   Ewald2: stock LAMMPS Ewald with the reciprocal-space data needed by
   compute stress/atom/resolved exposed as public members
------------------------------------------------------------------------- */

#ifndef LMP_EWALD2_H
#define LMP_EWALD2_H

#include "ewald.h"

namespace LAMMPS_NS {

class Ewald2 : public Ewald {
 public:
  Ewald2(class LAMMPS *lmp) : Ewald(lmp) {}

  using Ewald::kcount;
  using Ewald::kmax;
  using Ewald::kxvecs;
  using Ewald::kyvecs;
  using Ewald::kzvecs;
  using Ewald::ug;
  using Ewald::vg;
  using Ewald::cs;
  using Ewald::sn;
  using KSpace::scale;
};

}    // namespace LAMMPS_NS

#endif
