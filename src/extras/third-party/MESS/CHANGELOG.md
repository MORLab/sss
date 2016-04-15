## version 1.0
Compared to the predecessor lyapack a couple of things have changed.
- The user supplied functions are now managed by an operator manager
- The low rank ADI now has:
  - optimized treatment of E matrices in generalized equations
  - more chioces for shift selection, including completely automatic generation of shifts
  - improved stopping criteria based on low rank factors of the current residual
  - automatic generation of real low rank factors also for complex shifts
- The Newton-Kleinman iteration features:
  - optimized treatment of E matrices in generalized equations
  - improved stopping criteria based on low rank factors of the current residual
  - inexact Newton, line search and Galerkin projection acceleration
- Examples have been extended
- The Riccati iteration for H-infinity Riccati equations was added
- DSPMR has not yet been ported to the new infrastructure
- The SRM routine for balanced truncation is only available for none-DAE systems. Still, DAE versions are included in the corresponding DEMOS.
