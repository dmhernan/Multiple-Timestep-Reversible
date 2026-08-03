# ============================================================
# MTR N-body Integrator — Input File README and Template
# ============================================================
#
# OVERVIEW
# --------
# This file describes the input format for the MTR (Multiple
# Timestep Reversible) N-body integrator. Copy this file,
# fill in your initial conditions, and run:
#
#   ./mtr input.txt
#
# All quantities are in the following units:
#   Length   : au (astronomical units)
#   Time     : days
#   Mass     : solar masses
#   Velocity : au/day
#
# G = 2.95912208232212840e-4   au^3 day^-2 solar_mass^-1
#
# ============================================================
# PARAMETERS
# ============================================================
#
# n_bodies        [integer]
#   Number of bodies. The first body is the central body
#   (e.g. the star), treated specially in the democratic
#   heliocentric coordinate system.
#
# timestep        [days]
#   Global timestep h_0. Example: 0.01 yr = 3.6525 days.
#   Smaller timesteps give better energy conservation but
#   are slower.
#
# tmax            [days]
#   Total integration time. Example: 100 yr = 36525 days.
#
# output_interval [integer]
#   Write to fcons.txt every this many steps.
#   1 = every step, 10 = every 10th step, etc.
#
# rmax            [dimensionless]
#   Radius threshold for level 1 (the coarsest active level).
#   Pairs with ratio > rmax are not subcycled (stay at level 1).
#   Typical value: 30.0 for free-fall criterion (flag=0),
#   or 3.0 for Hill radius criterion (flag=1).
#
# rsub            [dimensionless]
#   Ratio between consecutive level radius thresholds.
#   rlev[i] = rmax / rsub^i.
#   Typical value: 2.0.
#
# hsub            [integer]
#   Substep multiplier between levels.
#   hlev[i] = h_0 / hsub^i.
#   Typical values: 3 or 4.
#
# levtot          [integer]
#   Maximum number of timestep levels.
#   30 is sufficient for most applications.
#
# flag            [0 or 1]
#   Level assignment criterion:
#   0 = free-fall timescale: ratio = sqrt(r^3 / (G*M)) / h_0
#   1 = Hill radius:         ratio = r / r_Hill
#
# rhill           [au]
#   Hill radius. Only used when flag=1. Set to 0.0 if flag=0.
#
# ============================================================
# INITIAL CONDITIONS
# ============================================================
#
# One line per body, in the format:
#   mass   x   y   z   vx   vy   vz
#
# Units:
#   mass         : solar masses
#   x, y, z      : au
#   vx, vy, vz   : au/day
#
# The first body is the central body (star).
# The centre of mass should be at the origin with zero
# total momentum — the code will subtract the CM if needed,
# but it is best practice to provide CM-corrected ICs.
#
# ============================================================
# TEMPLATE — replace values below with your own
# ============================================================

n_bodies        3
timestep        3.65250000000000000e+00
tmax            3.65250000000000000e+04
output_interval 1
rmax            3.00000000000000000e+01
rsub            2.00000000000000000e+00
hsub            3
levtot          30
flag            0
rhill           0.00000000000000000e+00

# mass                x                   y                   z                   vx                  vy                  vz
1.00000000000000000e+00  0.00000000000000000e+00  0.00000000000000000e+00  0.00000000000000000e+00  0.00000000000000000e+00  0.00000000000000000e+00  0.00000000000000000e+00
1.00000000000000000e-03  1.00000000000000000e+00  0.00000000000000000e+00  0.00000000000000000e+00  0.00000000000000000e+00  1.72069716278813000e-02  0.00000000000000000e+00
1.00000000000000000e-03  2.00000000000000000e+00  0.00000000000000000e+00  0.00000000000000000e+00  0.00000000000000000e+00  1.21668046398827000e-02  0.00000000000000000e+00
