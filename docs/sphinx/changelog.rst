.. _changelog:


#########
Changelog
#########

******************
0.3.1 (unreleased)
******************

Internal Changes
================
- Added the ability to update the micromorphic hydra object with a new unknown vector (:pull:`34`). By `Nathan Miller`_.
- Added the calculation of the current stress measures in micromorphic linear elasticity (:pull:`35`). By `Nathan Miller`_.

******************
0.3.0 (01-03-2023)
******************

Release
=======
- Released version (:pull:`33`). By `Nathan Miller`_.

Breaking Changes
================
- Added macros for the setter functions (:pull:`24`). By `Nathan Miller`_.
- Added macros for the getter functions (:pull:`25`). By `Nathan Miller`_.

New Features
============
- Added general setter functions for iteration and previous data (:pull:`23`). By `Nathan Miller`_.
- Added calculation of previous linear elastic stress (:pull:`26`). By `Nathan Miller`_.
- Added an isotropic damage configuration residual (:pull:`32`). By `Nathan Miller`_.

Internal Changes
================
- Copied over micromorphic linear elasticity subroutines to tardigrade hydra (:pull:`17`). By `Nathan Miller`_.
- Added initial micromorphic linear elastic residual (:pull:`18`). By `Nathan Miller`_.
- Added the calculation of the micromorphic linear elastic derived deformation measures (:pull:`19`). By `Nathan Miller`_.
- Added the calculation of the micromorphic linear elastic reference stress measures (:pull:`20`). By `Nathan Miller`_.
- Added the calculation of the Peryzna-based damage and the Jacobians (:pull:`30`). By `Nathan Miller`_.
- Added the calculation of the Peryzna-based damage deformation gradient's Jacobians (:pull:`31`). By `Nathan Miller`_.

******************
0.2.0 (12-11-2023)
******************

Release
=======
- Released version 0.2.0 (:pull:`15`). By `Nathan Miller`_.

Breaking Changes
================
- Changed hydra function calls to be more general (:pull:`1`, :pull:`2`, :pull:`3`, :pull:`4`, :pull:`5`, :pull:`6`). By `Nathan Miller`_.

New Features
============
- Added micromorphic hydra object (:pull:`7`). By `Nathan Miller`_.

Internal Changes
================
- Added decomposition of the micro deformations (:pull:`8`). By `Nathan Miller`_.
- Added the calculation of sub micro configurations (:pull:`9`). By `Nathan Miller`_.
- Added the jacobians of the sub micro configurations w.r.t. the micro configurations (:pull:`10`). By `Nathan Miller`_.
- Added generalization of the computation of the Jacobians of the first configurations (:pull:`11`). By `Nathan Miller`_.
- Added computation of the Jacobian of the first micro-configuration (:pull:`12`). By `Nathan Miller`_.
- Added computation of the gradient of the micro-deformations in their local reference configurations (:pull:`13`). By `Nathan Miller`_.
- Added computation of Jacobian of the gradient of the micro-deformations in their local reference configurations (:pull:`14`). By `Nathan Miller`_.

******************
0.1.2 (12-06-2023)
******************

Breaking Changes
================
- Changed getSubConfiguration to not include the upper bound (:merge:`7`). By `Nathan Miller`_.
- Change project name to tardigrade-hydra (:merge:`17`). by `Nathan Miller`_.

New Features
============
- Added calculation of the gradients of the current and previous F1 configurations (:merge:`11`). By `Nathan Miller`_.
- Added residual class for constructing the residual equations (:merge:`12`). By `Nathan Miller`_.
- Added the initialization of the unknown vector (:merge:`14`). By `Nathan Miller`_.
- Added setting and checking the tolerance of the non-linear solve (:merge:`14`). By `Nathan Miller`_.
- Added setting and checking the tolerance for the line-search of the non-linear solve (:merge:`14`). By `Nathan Miller`_.
- Added the decomposition of the unknown vector and its application to the solution quantities (:merge:`14`). By `Nathan Miller`_.
- Added the solution of the non-linear problem (:merge:`14`). By `Nathan Miller`_.
- Added a linear elastic implementation of a residual for use in testing (:merge:`18`). By `Nathan Miller`_.
- Added the evaluation of hydra to compute the required quantities (:merge:`18`). By `Nathan Miller`_.
- Added a linear viscoelastic implementation of a residual (:merge:`20`). By `Nathan Miller`_.
- Added the residual for a thermal expansion model (:merge:`21`). By `Nathan Miller`_.
- Added the residual for a Peryzna viscoplasticity model (:merge:`24`). By `Nathan Miller`_. 

Internal Changes
================
- Initialized the repository from cpp_stub (:merge:`1`). By `Nathan Miller`_.
- Added getters for the base quantities (:merge:`2`). By `Nathan Miller`_.
- Added additional libraries required for the project to update the environment (:merge:`4`). By `Nathan Miller`_.
- Updated the environment.txt file to reflect the new recipe (:merge:`5`). By `Nathan Miller`_.
- Added the decomposition of the incoming state variable vector into the configurations, state variables
  in the non-linear solve, and additional state variables (:merge:`3`). By `Nathan Miller`_.
- Added function to get a subset of the full deformation gradient (:merge:`6`). By `Nathan Miller`_.
- Added functions to get the part of the sub-configuration preceding and following a given
  configuration (:merge:`7`). By `Nathan Miller`_.
- Generalized the computation of the sub-configuration so that one can use either the current or
  previous configurations (:merge:`8`). By `Nathan Miller`_.
- Added the computation of the previous sub-configurations and the previous preceding and following
  sub-configurations given a configuration (:merge:`8`). By `Nathan Miller`_.
- Added the computation of the gradient of a sub-configuration by all of the configurations (:merge:`9`). By `Nathan Miller`_.
- Added gradients for the preceding and following sub-configurations for the current and previous configurations (:merge:`9`). By `Nathan Miller`_.
- Required >= version 0.5.3 of vector_tools (:merge:`10`). By `Nathan Miller`_.
- Added the construction residual, Jacobian, and other values (:merge:`12`). By `Nathan Miller`_.
- Clean up conda package CI files after ``conda build`` (:issue:`2`, :merge:`15`). By `Sergio Cordova`_.
- Changed the convergence_error type to use standard strings (:merge:`18`). By `Nathan Miller`_.
- Changed the version extraction script (:merge:`19`). By `Nathan Miller`_.
- linearViscoelasticity: Added elastic deformation gradient decomposition to linear viscoelasticity (:merge:`20`). By `Nathan Miller`_.
- linearViscoelasticity: Generalized the decomposition of the current elastic deformation gradient to current and previous (:merge:`20`). By `Nathan Miller`_.
- linearViscoelasticity: Added the decomposition of the additional state variable vector into volumetric and isochoric parts (:merge:`20`). By `Nathan Miller`_.
- linearViscoelasticity: Added the computation of the rate multipliers and the integration alpha parameter (:merge:`20`). By `Nathan Miller`_.
- linearViscoelasticity: Added the construction of the viscoelastic parameter vectors which are able to be parsed by stressTools::linearViscoelasticity (:merge:`20`). By `Nathan Miller`_.
- linearViscoelasticity: Changed the isochoric moduli going into linear viscoelasticity to be 2x the moduli (:merge:`20`). By `Nathan Miller`_.
- linearViscoelasticity: Added the computation of the mean and isochoric viscoelastic PK2 stresses (:merge:`20`). By `Nathan Miller`_.
- linearViscoelasticity: Added the computation of the PK2 stress (:merge:`20`). By `Nathan Miller`_.
- linearViscoelasticity: Added the gradients of the rate multipliers w.r.t. the temperatures (:merge:`20`). By `Nathan Miller`_.
- linearElasticity: Exposed dPK2StressdFe to users through getter-setter functions (:merge:`20`). By `Nathan Miller`_.
- linearElasticity: Changed dPK2dXXX names to dPK2StressdXXX (:merge:`20`). By `Nathan Miller`_.
- linearViscoelasticity: Added the computation of dPK2StressdFe and dPK2StressdT (:merge:`20`). By `Nathan Miller`_.
- linearElasticity: Changed XXXdPK2 names to XXXdPK2Stress (:merge:`20`). By `Nathan Miller`_.
- Updated documentation strings to eliminate all undefined references in the documentation generation (:merge:`21`). By `Nathan Miller`_.
- thermalExpansion: Added the remaining derivatives of the residual (:merge:`21`). By `Nathan Miller`_.
- thermalExpansion: Removed extraneous print statements (:merge:`22`). By `Nathan Miller`_.
