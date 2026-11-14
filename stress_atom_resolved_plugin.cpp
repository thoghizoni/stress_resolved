// plugin registration, the single entry point LAMMPS calls on load
#include "lammpsplugin.h"
#include "version.h"

#include "compute_stress_atom_resolved.h"
#include "ewald2.h"

using namespace LAMMPS_NS;

static KSpace *ewald2_creator(LAMMPS *lmp)
{
  return new Ewald2(lmp);
}

static Compute *stress_atom_resolved_creator(LAMMPS *lmp, int narg, char **arg)
{
  return new ComputeStressAtomResolved(lmp, narg, arg);
}

extern "C" void lammpsplugin_init(void *lmp, void *handle, void *regfunc)
{
  lammpsplugin_t plugin;
  lammpsplugin_regfunc register_plugin = (lammpsplugin_regfunc) regfunc;

  plugin.version = LAMMPS_VERSION;
  plugin.style = "kspace";
  plugin.name = "ewald2";
  plugin.info = "Ewald with public reciprocal-space data";
  plugin.author = "Thomaz R. Ghizoni (thoghizoni@gmail.com)";
  plugin.creator.v1 = (lammpsplugin_factory1 *) &ewald2_creator;
  plugin.handle = handle;
  (*register_plugin)(&plugin, lmp);

  plugin.style = "compute";
  plugin.name = "stress/atom/resolved";
  plugin.info = "Pair + kspace stress/tally between two groups";
  plugin.creator.v2 = (lammpsplugin_factory2 *) &stress_atom_resolved_creator;
  (*register_plugin)(&plugin, lmp);
}
