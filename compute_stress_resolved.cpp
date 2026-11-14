/* ----------------------------------------------------------------------
   compute ID group stress/resolved group2

   same outputs as compute stress/tally, with kspace added:
     per-atom array with 6 columns xx,yy,zz,xy,xz,yz, same sign and units
       as compute stress/atom
     global scalar = trace/3 of the summed virial, as compute stress/tally

   pair i-j with i in group and j in group2 (or vice versa)
     gives 1/2 r_ij (x) F_ij to both i and j
   kspace gives atom i its Ewald virial with the charges it pairs with,
     using S = S_a + S_b + S_c of kspace_style ewald2, where a = group only,
     b = group2 only, c = both
------------------------------------------------------------------------- */

#include "compute_stress_resolved.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "ewald2.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeStressResolved::ComputeStressResolved(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), stress(nullptr), sfac(nullptr), sfac_all(nullptr), list(nullptr),
    ewald(nullptr)
{
  if (narg != 4) error->all(FLERR, "Illegal compute stress/resolved command");

  int igroup2 = group->find(arg[3]);
  if (igroup2 == -1)
    error->all(FLERR, "Could not find compute stress/resolved second group ID {}", arg[3]);
  groupbit2 = group->bitmask[igroup2];

  scalar_flag = 1;
  extscalar = 0;
  peratom_flag = 1;
  size_peratom_cols = 6;
  comm_reverse = 6;

  nmax = 0;
  kmax_sfac = 0;
}

/* ---------------------------------------------------------------------- */

ComputeStressResolved::~ComputeStressResolved()
{
  memory->destroy(stress);
  memory->destroy(sfac);
  memory->destroy(sfac_all);
}

/* ---------------------------------------------------------------------- */

void ComputeStressResolved::init()
{
  if (force->pair == nullptr)
    error->all(FLERR, "Compute stress/resolved requires a pair style");

  // single() of these styles does not return the force of compute():
  // EAM adds the embedding energy to phi before taking its derivative,
  // DPD styles give only the conservative force, or none
  // pair_match() also checks hybrid sub-styles and accelerated variants

  for (const char *pattern : {"^eam", "dpd"}) {
    Pair *bad = force->pair_match(pattern, 0, 1);
    if (bad)
      error->all(FLERR, "Compute stress/resolved does not support pair style {}: its single() "
                 "does not give the force used in the run", force->pair_match_ptr(bad));
  }

  if (force->pair->single_enable == 0)
    error->all(FLERR, "Pair style does not support compute stress/resolved");
  if (force->pair->manybody_flag && (comm->me == 0))
    error->warning(FLERR, "Compute stress/resolved used with many-body pair style {}: single() "
                   "may not include the many-body part of the forces", force->pair_style);

  ewald = nullptr;
  if (force->kspace) {
    ewald = dynamic_cast<Ewald2 *>(force->kspace);
    if (!ewald) error->all(FLERR, "Compute stress/resolved requires kspace_style ewald2");
    if (ewald->slabflag)
      error->all(FLERR, "Compute stress/resolved does not support kspace slab correction");
    if (!ewald->compute_flag)
      error->all(FLERR, "Compute stress/resolved requires kspace_modify compute yes");
  }

  neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputeStressResolved::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

double ComputeStressResolved::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (invoked_peratom != invoked_scalar) compute_peratom();

  // sum accumulated virial across procs

  double vector[6];
  MPI_Allreduce(virial, vector, 6, MPI_DOUBLE, MPI_SUM, world);

  if (domain->dimension == 3)
    scalar = (vector[0] + vector[1] + vector[2]) / 3.0;
  else
    scalar = (vector[0] + vector[1]) / 2.0;

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeStressResolved::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  if (atom->nmax > nmax) {
    memory->destroy(stress);
    nmax = atom->nmax;
    memory->create(stress, nmax, 6, "stress/resolved:stress");
    array_atom = stress;
  }

  const int nlocal = atom->nlocal;
  const int ntotal = force->newton_pair ? nlocal + atom->nghost : nlocal;
  for (int i = 0; i < ntotal; i++)
    for (int j = 0; j < 6; j++) stress[i][j] = 0.0;
  for (int j = 0; j < 6; j++) virial[j] = 0.0;

  compute_pair();
  if (ewald) compute_kspace();

  if (force->newton_pair) comm->reverse_comm(this);

  // convert to stress*volume units = -pressure*volume

  const double nktv2p = -force->nktv2p;
  for (int i = 0; i < nlocal; i++)
    for (int j = 0; j < 6; j++) stress[i][j] *= nktv2p;
}

/* ----------------------------------------------------------------------
   pair virial from Pair::single() on an occasional half list
   same ghost logic as Pair::ev_tally()
------------------------------------------------------------------------- */

void ComputeStressResolved::compute_pair()
{
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, fpair, factor_coul, factor_lj;
  double v[6];
  int *jlist;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  Pair *pair = force->pair;
  double **cutsq = pair->cutsq;

  neighbor->build_one(list);

  int inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      // same pair selection as compute stress/tally

      if (!(((mask[i] & groupbit) && (mask[j] & groupbit2)) ||
            ((mask[i] & groupbit2) && (mask[j] & groupbit))))
        continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];
      if (rsq >= cutsq[itype][jtype]) continue;

      pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);

      fpair *= 0.5;
      v[0] = delx * delx * fpair;
      v[1] = dely * dely * fpair;
      v[2] = delz * delz * fpair;
      v[3] = delx * dely * fpair;
      v[4] = delx * delz * fpair;
      v[5] = dely * delz * fpair;

      if (newton_pair || i < nlocal)
        for (int c = 0; c < 6; c++) {
          stress[i][c] += v[c];
          virial[c] += v[c];
        }
      if (newton_pair || j < nlocal)
        for (int c = 0; c < 6; c++) {
          stress[j][c] += v[c];
          virial[c] += v[c];
        }
    }
  }
}

/* ----------------------------------------------------------------------
   kspace virial: Ewald per-atom virial with S split into S_a, S_b and S_c
   reuses cs/sn computed by Ewald2 for the current positions
------------------------------------------------------------------------- */

void ComputeStressResolved::compute_kspace()
{
  int i, k, c, kx, ky, kz;
  double cypz, sypz, exprl, expim, sre, sim, partial;

  const int kcount = ewald->kcount;
  const int *kxvecs = ewald->kxvecs;
  const int *kyvecs = ewald->kyvecs;
  const int *kzvecs = ewald->kzvecs;
  const double *ug = ewald->ug;
  double **vg = ewald->vg;
  double ***cs = ewald->cs;
  double ***sn = ewald->sn;

  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (kcount > kmax_sfac) {
    memory->destroy(sfac);
    memory->destroy(sfac_all);
    kmax_sfac = kcount;
    memory->create(sfac, 6 * kmax_sfac, "stress/resolved:sfac");
    memory->create(sfac_all, 6 * kmax_sfac, "stress/resolved:sfac_all");
  }

  // partial structure factors: sfac[6k+0,1] = Re,Im S_a  sfac[6k+2,3] = Re,Im S_b
  // sfac[6k+4,5] = Re,Im S_c, with a = group only, b = group2 only, c = both

  for (k = 0; k < 6 * kcount; k++) sfac[k] = 0.0;

  for (k = 0; k < kcount; k++) {
    kx = kxvecs[k];
    ky = kyvecs[k];
    kz = kzvecs[k];
    for (i = 0; i < nlocal; i++) {
      c = membership(mask[i]);
      if (c < 0) continue;
      cypz = cs[ky][1][i] * cs[kz][2][i] - sn[ky][1][i] * sn[kz][2][i];
      sypz = sn[ky][1][i] * cs[kz][2][i] + cs[ky][1][i] * sn[kz][2][i];
      exprl = cs[kx][0][i] * cypz - sn[kx][0][i] * sypz;
      expim = sn[kx][0][i] * cypz + cs[kx][0][i] * sypz;
      sfac[6 * k + c] += q[i] * exprl;
      sfac[6 * k + c + 1] += q[i] * expim;
    }
  }

  MPI_Allreduce(sfac, sfac_all, 6 * kcount, MPI_DOUBLE, MPI_SUM, world);

  // per-atom virial of every local atom with the charges it pairs with:
  // group only -> S_b + S_c, group2 only -> S_a + S_c, both -> S_a + S_b + S_c

  const double qscale = force->qqrd2e * ewald->scale;

  for (i = 0; i < nlocal; i++) {
    c = membership(mask[i]);
    if (c < 0) continue;
    const int usea = (c != 0);
    const int useb = (c != 2);
    double vi[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (k = 0; k < kcount; k++) {
      kx = kxvecs[k];
      ky = kyvecs[k];
      kz = kzvecs[k];
      cypz = cs[ky][1][i] * cs[kz][2][i] - sn[ky][1][i] * sn[kz][2][i];
      sypz = sn[ky][1][i] * cs[kz][2][i] + cs[ky][1][i] * sn[kz][2][i];
      exprl = cs[kx][0][i] * cypz - sn[kx][0][i] * sypz;
      expim = sn[kx][0][i] * cypz + cs[kx][0][i] * sypz;
      sre = sfac_all[6 * k + 4];
      sim = sfac_all[6 * k + 5];
      if (usea) {
        sre += sfac_all[6 * k];
        sim += sfac_all[6 * k + 1];
      }
      if (useb) {
        sre += sfac_all[6 * k + 2];
        sim += sfac_all[6 * k + 3];
      }
      partial = ug[k] * (exprl * sre + expim * sim);
      for (c = 0; c < 6; c++) vi[c] += vg[k][c] * partial;
    }
    for (c = 0; c < 6; c++) {
      stress[i][c] += vi[c] * q[i] * qscale;
      virial[c] += vi[c] * q[i] * qscale;
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputeStressResolved::pack_reverse_comm(int n, int first, double *buf)
{
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++)
    for (int j = 0; j < 6; j++) buf[m++] = stress[i][j];
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeStressResolved::unpack_reverse_comm(int n, int *list, double *buf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    int j = list[i];
    for (int c = 0; c < 6; c++) stress[j][c] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

double ComputeStressResolved::memory_usage()
{
  return (double) nmax * 6 * sizeof(double) + (double) kmax_sfac * 12 * sizeof(double);
}
