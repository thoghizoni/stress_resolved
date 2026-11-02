# compute stress/atom/resolved

A LAMMPS plugin that splits the per-atom virial stress of every atom into two parts. One part comes from the atom's interaction partners in group A and the other from its partners in group B. It covers both pair and long-range (Ewald) Coulomb interactions.

A typical use is a liquid in contact with a solid. You can separate the stress an atom feels from the liquid from the stress it feels from the wall.

The plugin adds two styles to LAMMPS:

| style | kind | purpose |
|---|---|---|
| `stress/atom/resolved` | compute | the 12-column resolved per-atom stress |
| `ewald2` | kspace | the stock LAMMPS `ewald` with its reciprocal-space data exposed to the compute. It gives exactly the same forces, energies and pressures as `ewald`. |

## What it computes

For each atom *i* in the compute group, a per-atom array with 12 columns:

| columns | content |
|---|---|
| 1–6 | stress of atom *i* from partners *j* in **group A**: xx, yy, zz, xy, xz, yz |
| 7–12 | stress of atom *i* from partners *j* in **group B**: xx, yy, zz, xy, xz, yz |

**Pair part.** Every pair *i–j* within the cutoff gives ½ **r**<sub>ij</sub> ⊗ **F**<sub>ij</sub> to atom *i*. The group of *j* decides which column block it goes into. Atom *j* gets the same amount, in the block chosen by the group of *i*. Forces come from the pair style's `single()` function, so special-bond factors are applied as in the pair style.

**k-space part.** The Ewald per-atom virial of atom *i* is

  W<sub>i</sub> = q<sub>i</sub> Σ<sub>k</sub> u<sub>k</sub> **v**<sub>k</sub> Re[ e<sup>i**k**·**r**<sub>i</sub></sup> S(**k**)* ],  where S(**k**) = Σ<sub>j</sub> q<sub>j</sub> e<sup>i**k**·**r**<sub>j</sub></sup>

The compute splits the structure factor into S = S<sub>A</sub> + S<sub>B</sub>, summing only over charges in A or only over charges in B. W<sub>i</sub> is then evaluated once with S<sub>A</sub> (columns 1–6) and once with S<sub>B</sub> (columns 7–12).

**Units and sign** are the same as `compute stress/atom`: stress × volume, in pressure·volume units (atm·Å³ for `units real`). As with `stress/atom`, the per-atom pressure tensor is the negative of these values divided by a volume.

**Consistency.** If A and B together contain every atom, then for each component, columns *n* + *n*+6 equal `compute stress/atom NULL pair kspace` to machine precision. This holds with newton on and off and on any number of MPI ranks.

### What is not included

- bond, angle, dihedral and improper terms, fixes, and kinetic energy. Only `pair` and `kspace` contributions are included.
- interactions with partners that are in neither A nor B. These are left out of both blocks.

## Requirements

- LAMMPS built as a shared library or executable with the **PLUGIN** and **KSPACE** packages. The plugin must be compiled against the same LAMMPS version and with the same MPI.
- A pair style that implements `single()`. Most pairwise styles do, for example `lj/cut/coul/long`, `born/coul/long`, `buck/coul/long` and `hybrid` / `hybrid/overlay` of such styles. Many-body styles (Tersoff, EAM, ReaxFF, …) do not.
- Long-range Coulomb has to use `kspace_style ewald2`. Other k-space styles (`ewald`, `pppm`, …) are rejected with an error. Without a k-space style the compute gives the pair part only.
- With `ewald2`, slab correction (`kspace_modify slab`) is not supported and `kspace_modify compute` must be `yes` (the default).
- Groups A and B must be **disjoint**. If any atom belongs to both, LAMMPS stops with an error. They do not have to cover all atoms.

Tested with LAMMPS 2 Aug 2023 and OpenMPI on Linux.

## Building

### 1. LAMMPS

The plugin is compiled against the LAMMPS source tree and loaded into a LAMMPS executable built from that same tree. If your LAMMPS doesn't have the PLUGIN and KSPACE packages yet, build it with CMake, for example:

```bash
cd /path/to/lammps
mkdir build && cd build
cmake ../cmake -D PKG_PLUGIN=on -D PKG_KSPACE=on -D BUILD_SHARED_LIBS=on   # plus your usual packages
cmake --build . -j 8
```

`lmp -h` lists the installed packages, so you can check that PLUGIN and KSPACE are there.

### 2. The plugin

You need CMake ≥ 3.16, a C++17 compiler and the MPI compiler wrappers that were used to build LAMMPS. From this folder:

```bash
bash build.sh /path/to/lammps/src
```

This writes `build/stress_atom_resolved.so`.

Optional variables, set in the environment:

| variable | meaning | default |
|---|---|---|
| `LAMMPS_SRC` | the LAMMPS `src` folder, if not given as the argument (it must contain `lammpsplugin.h` and `KSPACE/ewald.h`) | none, required |
| `MPICC`, `MPICXX` | MPI compiler wrappers used to build LAMMPS | `mpicc`, `mpicxx` from `PATH` |
| `BUILD_TYPE` | `Release` (-O3) or `Debug` | `Release` |
| `BUILD_DIR` | output folder | `build` |

Example with a specific MPI installation:

```bash
MPICC=/opt/openmpi/bin/mpicc MPICXX=/opt/openmpi/bin/mpicxx bash build.sh /path/to/lammps/src
```

To use CMake directly:

```bash
cmake -S . -B build -D LAMMPS_SOURCE_DIR=/path/to/lammps/src -D CMAKE_CXX_COMPILER=mpicxx
cmake --build build
```

### 3. Test

Run the example in [`examples/il_silica`](examples/il_silica). It takes about 15 s on 4 cores and checks the plugin against `compute stress/atom`.

### Note on the link flags

The plugin is linked with `-Wl,-z,nodelete`. Without that flag, LAMMPS (2 Aug 2023) crashes at exit: it unloads k-space plugins before deleting the active k-space object, so `ewald2`'s destructor would no longer be in memory. The flag only keeps the plugin code loaded until the process ends and does not affect results.

Both `-z nodelete` and `-rdynamic` are GNU/Linux linker options. On macOS, `CMakeLists.txt` would need to be adapted.

## Usage

```
compute ID group-ID stress/atom/resolved groupA groupB
```

- `ID`, `group-ID`: as for any compute. Atoms outside `group-ID` get zeros, but they still count as partners if they belong to A or B.
- `groupA`, `groupB`: the two disjoint partner groups.

Example: a liquid (types 1–16) against a silica wall (types 17–20). The complete, runnable version is in [`examples/il_silica`](examples/il_silica).

```
plugin load     /path/to/build/stress_atom_resolved.so

units           real
atom_style      full
pair_style      born/coul/long 15.0 15.0
kspace_style    ewald2 1.0e-5          # instead of ewald
read_data       system.lmpsys
...

group           liquid type 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
group           solid  type 17 18 19 20

compute         res all stress/atom/resolved liquid solid

# per-atom output: id, then stress from liquid partners (1-6) and from solid partners (7-12)
dump            d all custom 100 resolved.dump id c_res[*]

# or a global total of each block
compute         tot all reduce sum c_res[*]
thermo_style    custom step temp press c_tot[1] c_tot[7]
```

The array can be used anywhere per-atom arrays are accepted: `dump custom`, `fix ave/atom`, `fix ave/chunk`, `compute reduce`, atom-style variables, and so on.

To check that A + B reproduces LAMMPS's own per-atom stress when A ∪ B is all atoms:

```
compute ref all stress/atom NULL pair kspace
variable d1 atom abs(c_res[1]+c_res[7]-c_ref[1])
```

## Cost

The compute only does work on steps where something asks for it (a dump, thermo, fix ave/…). On those steps it:

- builds an occasional half neighbor list and calls `single()` on every pair within the cutoff
- runs one extra pass over atoms × k-vectors, and one `MPI_Allreduce` of 4 × (number of k-vectors) values

For an ionic liquid in a silica pore (10281 atoms), evaluating it on every step costs about 2.3 times a plain MD step.

## Files

| file | content |
|---|---|
| `compute_stress_atom_resolved.h/.cpp` | the compute |
| `ewald2.h` | `Ewald2`, a subclass of LAMMPS's `Ewald` that exposes k-vectors, `ug`, `vg`, `cs`, `sn` and `scale` |
| `stress_atom_resolved_plugin.cpp` | plugin entry point, registers `ewald2` and `stress/atom/resolved` |
| `CMakeLists.txt` | CMake build of the plugin module |
| `build.sh` | build script with configurable paths |
| `examples/il_silica/` | test system and input that checks the plugin against `compute stress/atom` |
