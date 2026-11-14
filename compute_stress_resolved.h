/* -*- c++ -*- ----------------------------------------------------------
   compute stress/resolved: pair + kspace version of compute stress/tally,
   per-atom and global virial from interactions between group and group2
------------------------------------------------------------------------- */

#ifndef LMP_COMPUTE_STRESS_RESOLVED_H
#define LMP_COMPUTE_STRESS_RESOLVED_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeStressResolved : public Compute {
 public:
  ComputeStressResolved(class LAMMPS *, int, char **);
  ~ComputeStressResolved() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  double compute_scalar() override;
  void compute_peratom() override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;

 private:
  int groupbit2;
  int nmax, kmax_sfac;
  double **stress;
  double *sfac, *sfac_all;
  double virial[6];
  class NeighList *list;
  class Ewald2 *ewald;

  void compute_pair();
  void compute_kspace();

  // sfac offset of an atom: 0 = group only, 2 = group2 only, 4 = both, -1 = neither

  int membership(int m) const
  {
    if (m & groupbit) return (m & groupbit2) ? 4 : 0;
    return (m & groupbit2) ? 2 : -1;
  }
};

}    // namespace LAMMPS_NS

#endif
