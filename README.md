# MTR + snapshot output

This directory is MTR code with **only a
snapshot-output feature added**. The integrator itself (integrator.c, maps.c,
kepler.c, coordinates.c) is identical to original version of MTR coded by David.

## Enabling snapshots

Add the following keys to the input file (without them the code behaves
exactly like `original MTR`, i.e. no snapshots):

```
dt_snap         1.0        # snapshot interval (code time units)
snap_dir        snap       # output directory name (default "snap")
```

- The interval is rounded to global steps: one snapshot every
  `round(dt_snap / timestep)` steps.
- The output directory is created automatically at runtime.
- Snapshots are written under the **current working directory at run time**,
  so run the code from inside the run directory:

```bash
cd <run_dir>            # where input.txt lives
./mtr input.txt
```

`snap000000.dat` is t = 0 (the initial conditions); subsequent files follow
every `dt_snap`.

## Initial conditions (IC_file/)

`IC_file/` contains initial input files for the Duncan et al. (1998)
§6.5 "Accretion of the Moon" problem with N = 50, 100, 200, 400, 800, and
975 particles (`input_N####.txt`), used for the survey. The other
survey cases (rhill multiples, hsub = 4) differ only in the `rhill` / `hsub`
parameter lines — the particle data are identical. See `IC_file/initial.md`
for details on the ICs, units.

## Snapshot format

GPLUM-compatible snapshot format (tab-separated) with three extra orbital-
element columns appended. Positions and velocities are **heliocentric**
(relative to body 0, the central body), and **the central body itself is not
included** in the file. All quantities are in the code's internal units.

### Line 1: header (13 columns)

```
time  n_body  id_next  0 0 0 0 0 0 0 0 0 0
```

- `time` — snapshot time
- `n_body` — number of particles (central body excluded)
- `id_next` — next ID to assign (= n_body)
- remaining 10 columns — placeholders (0) for GPLUM's energy header fields

### Line 2 onward: one line per particle (15 columns)

```
id  mass  r_planet  f  x  y  z  vx  vy  vz  nn  flag  a  e  inc
```

| column | content |
|---|---|
| 1 | id (0-based) |
| 2 | mass |
| 3 | r_planet (unused, always 0) |
| 4 | f (placeholder, always 1) |
| 5–7 | x, y, z (heliocentric position) |
| 8–10 | vx, vy, vz (heliocentric velocity) |
| 11–12 | nn, flag (placeholders, always 0) |
| 13–15 | **a, e, inc** — heliocentric osculating elements computed with mu = G(M_central + m_i); inc in radians; **a < 0 means unbound (ejected)** |

## Changing the output format

All formatting are coded in `write_snapshot()` at the end of `physics.c`;
nothing else needs to be touched.

- **Header line**: the first `fprintf` pair. To drop the GPLUM placeholder
  columns, remove the `for (k = 0; k < 10; k++)` loop.
- **Particle lines**: the large `fprintf` inside the `for (i = 1; i < n; i++)`
  loop. Add or remove columns by editing the format string and its argument
  list together (they must stay in one-to-one correspondence). For example,
  to drop the placeholder columns (r_planet, f, nn, flag) and write only
  `id mass x y z vx vy vz a e inc`, delete the corresponding `%` fields
  and arguments.
- Additional per-particle quantities are available in the loop:
  `dx[]`, `dv[]` (heliocentric), `x[k][i]`, `v[k][i]` (barycentric/DHC-derived
  Cartesian), and anything computable from them.
- The separator is `\t`; change it in the format strings if needed.
