# Example: ionic liquid in a silica slit pore

This example tests a build of the plugin. The system has 10281 atoms: an ionic liquid (atom types 1–16) between two amorphous silica walls (types 17–20). It uses `pair_style hybrid lj/cut/coul/long + born/coul/long` and `kspace_style ewald2`.

The input defines

```
compute ll liquid stress/atom/resolved liquid
compute ls liquid stress/atom/resolved solid
compute ss solid stress/atom/resolved solid
compute ref all stress/atom NULL pair kspace
```

and checks that, for every atom and component, `ll + ls + ss` = `stress/atom`. This must hold because `liquid` and `solid` are disjoint and together contain all atoms. Every pair and k-space term then falls in exactly one of the three computes, as in `compute stress/tally`.

## Files

| file | content |
|---|---|
| `il_silica.lmp` | LAMMPS input: 200 NVE steps with output every 50 |
| `il_silica.lmpsys` | data file (`atom_style full`, `units real`), including bond, angle and dihedral coefficients |
| `pair_coefficients_long.txt` | pair coefficients, read with `include` |

## Running

Build the plugin first (see the main README). Then, from this folder:

```bash
mpirun -np 4 lmp -in il_silica.lmp
```

Replace `lmp` with the name of your LAMMPS executable (`lmp_mpi`, …). By default the input loads `../../build/stress_atom_resolved.so`. If the plugin is somewhere else, give its path:

```bash
mpirun -np 4 lmp -in il_silica.lmp -var plugin /path/to/stress_atom_resolved.so
```

To test with newton off, add `-var newt off`.

On 4 cores the run takes about 15 s.

## Expected output

```
Step  Temp  PotEng  c_ll  c_ls  c_ss  v_maxref  v_maxdiff  v_reldiff
```

| column | meaning |
|---|---|
| `c_ll`, `c_ls`, `c_ss` | global scalars (as in `stress/tally`) of the liquid–liquid, liquid–solid and solid–solid interactions |
| `v_maxref` | largest \|stress/atom\| component of any atom |
| `v_maxdiff` | largest \|ll + ls + ss − stress/atom\| of any atom and component |
| `v_reldiff` | `v_maxdiff / v_maxref` |

`v_reldiff` should be at rounding level, around 1e-15 (anything below about 1e-10 is fine), at every output step, with newton on and off and on any number of cores. The other values change with the number of cores and the MPI library, because the trajectory depends on floating-point summation order.

Results with LAMMPS 2 Aug 2023 on 4 cores, newton on:

```
Step  v_reldiff
   0  3.66e-15
  50  2.57e-15
 100  4.08e-15
 150  2.87e-15
 200  2.61e-15
```
