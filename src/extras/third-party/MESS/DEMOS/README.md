We have collected a number of demonstration examples that can serve as
the starting point for your own scripts. The directories serve the
following purposes.

**models/**
  Contains the actual benchmark model data or functios for their generation

**FDM/ Rail/**
  Both demonstrate the use of the *"default"* operators, i.e.,
  the default set of user supplied functions (usfs). While Rail/ contains a
  symmetric system for a FEM semi-discretized heat equation, FDM uses
  a finite difference semi-discretized heat equation with convection
  and is thus non-symmetric.

**DAE2/**
  Demonstrates the use of index-2 Stokes-type DAE systems. By
  default a scalable finite Volume discretization for a Stokes system
  is used. There is also an option to use the FEM model for a Karman
  vortex shedding in a 2d domain. This, however, requires an
  additional download of approx. 270MB.

**TripleChain/**
  Uses the Truhar/Veselic model with three coupled
  mass-spring-damper-chains to demonstrate the use of the structure
  exploiting second order usfs. The directory contains both the case
  for the *"so_1"* usfs and the *"default"* set to demonstrate the
  advantage of the first.

**DAE2_SO/ DAE3_SO/**
  Show the use of the *"dae2_so"* and *"dae3_so"* usfs for
  second order index-2 and index-3 systems.

**RI/**
  Demonstrates the solver for Riccati equations with indefinite quadratic terms.
  
