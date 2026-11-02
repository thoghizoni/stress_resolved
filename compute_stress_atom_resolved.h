/* -*- c++ -*- ----------------------------------------------------------
   compute stress/atom/resolved: per-atom pair + kspace virial of atom i,
   split by the group (A or B) of its interaction partners
------------------------------------------------------------------------- */

#ifndef LMP_COMPUTE_STRESS_ATOM_RESOLVED_H
#define LMP_COMPUTE_STRESS_ATOM_RESOLVED_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeStressAtomResolved : public Compute {
 public:
  ComputeStressAtomResolved(class LAMMPS *, int, char **);
  ~ComputeStressAtomResolved() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;

 private:
  int groupbitA, groupbitB;
  int nmax, kmax_sfac;
  double **stress;
  double *sfac, *sfac_all;
  class NeighList *list;
  class Ewald2 *ewald;

  void compute_pair();
  void compute_kspace();
};

}    // namespace LAMMPS_NS

#endif
