# IC_file — input files for the N-scaling survey

Ready-to-run MTR input files for "The accretion of the Moon" problem
(Duncan, Levison & Lee 1998, §6.5), used for the N = 50–975 survey.

| file | particles (+ Earth as body 0) |
|---|---|
| input_N0050.txt | 50 |
| input_N0100.txt | 100 |
| input_N0200.txt | 200 |
| input_N0400.txt | 400 |
| input_N0800.txt | 800 |
| input_N0975.txt | 975 (full Duncan IC) |

Each file is the first N particles of the Duncan §6.5 initial conditions
(masses 3.2e-7–3.2e-4 M_earth, a in 0.7–1.3 Roche radii, e/inc from a
Rayleigh distribution with mean 0.3 truncated at 0.95 / 50 deg, pericenters
< R_earth removed), with the Earth added as body 0 and the system shifted to
the barycenter.

Units (see the header comments in each file): G equals the value hardcoded
in mtr.h; time/length/velocity numbers equal the problem-unit values
(Roche radius = 1, an orbit at a = 1 has period 2*pi, so dt = 0.2 and
t_end = 6*pi appear literally), and only masses are rescaled by 1/G
(Earth = 3379.38).

Parameters as shipped: flag = 1, hsub = 3, dt = 0.2, t_end = 6*pi,
rmax = 3, rsub = 2, rhill = mutual Hill radius of the two most massive
particles of that subset, dt_snap = 1 (snapshots every problem time unit).
The survey's other cases differ only in the `rhill` line (2x, 4x, 10x, 15x
the shipped value) and the `hsub` line (4) — edit those two lines to
reproduce them; the particle data are identical.