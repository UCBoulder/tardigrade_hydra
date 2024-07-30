.. _changelog:


#########
Changelog
#########

******************
0.5.0 (unreleased)
******************

New Features
============
- Added a data storage set object that will automatically set an object when the destructor is called (:pull:`143`). By `Nathan Miller`_.
- Added a setDataStorage to the macros for data storage definitions (:pull:`144`). By `Nathan Miller`_.

Internal Changes
================
- Working towards removing floatVector in favor of specific sizes (:pull:`142`). By `Nathan Miller`_.
- Removed set_varname and replaced with setDataStorage (:pull:`145`). By `Nathan Miller`_.

******************
0.4.4 (07-12-2024)
******************

Release
=======
- Released version (:pull:`141`). By `Nathan Miller`_.

New Features
============
- Added projection operators to the residualBase (:pull:`133`). By `Nathan Miller`_.
- Allow the user to turn off rank-deficient errors (:pull:`134`). By `Nathan Miller`_.
- Added a Levenberg-Marquardt solve in the case of a convergence failure (:pull:`135`). By `Nathan Miller`_.
- Added the ability to turn on or off applying projections (:pull:`136`). By `Nathan Miller`_.
- Added a projection for the micromorphic Drucker-Prager plasticity (:pull:`136`). By `Nathan Miller`_.

Internal Changes
================
- Set whether to use the projection to default to false (:pull:`137`). By `Nathan Miller`_.
- Automatically turn on the projectors when Levenberg-Marquardt is enabled (:pull:`138`). By `Nathan Miller`_.
- Set hydra to reset the unknown vector to the initial unknown after a failed solve (:pull:`139`). By `Nathan Miller`_.

Bug Fixes
=========
- Reset the iteration number when a re-attempt at a solve is performed (:pull:`138`). By `Nathan Miller`_.
- Fixed the use of the gradient descent flag for the nonlinear solve (:pull:`140`). By `Nathan Miller`_.

******************
0.4.3 (07-12-2024)
******************

Release
=======
- Released version (:pull:`132`). By `Nathan Miller`_.

New Features
============
- Added the gradient of the residual norm (:pull:`125`). By `Nathan Miller`_.
- Added a gradient step alternative to the Armijo-type line search (:pull:`126`). By `Nathan Miller`_.
- Added a test for whether the proposed direction is a descent direction (:pull:`128`). By `Nathan Miller`_.
- Added an automatic switch to gradient descent if the line search algorithm is not in a minimization direction (:pull:`130`). By `Nathan Miller`_.
- Added the adaptive Levenberg-Marquardt regularization parameter (:pull:`131`). By `Nathan Miller`_.
- Added Levenberg-Marquardt steps (:pull:`131`). By `Nathan Miller`_.

Internal Changes
================
- Moved the Armijo-type line search into a separate function (:pull:`126`). By `Nathan Miller`_.
- Added data containers that will be cleared after each nonlinear iteration (:pull:`126`). By `Nathan Miller`_.
- Added setting required data for gradient descent steps (:pull:`127`). by `Nathan Miller`_.
- Allow for the version number to be specified when doing a FetchContent build (:pull:`129`). By `Nathan Miller`_.

******************
0.4.2 (07-11-2024)
******************

Release
=======
- Released version (:pull:`124`). By `Nathan Miller`_.

New Features
============
- Throw a custom convergence error class rather than a nested exception if a failure happens because of the line-search or the Newton loop iterations (:pull:`70`). By `Nathan Miller`_.
- Added a pre-conditioner (jacobian scaling) to try and improve the stability of the Jacobian (:pull:`98`). By `Nathan Miller`_.
- Added a J2 flow isotropic-kinematic hardening viscoplastic model (:pull:`102`). By `Nathan Miller`_.
- Added a mass-change deformation gradient evolution model (:pull:`104`). By `Nathan Miller`_.
- Added the calculation of the total derivative of the unknown vector w.r.t. the additional degrees of freedom (:pull:`104`). By `Nathan Miller`_.
- Added storage for the derivative of the residual w.r.t. the additional dof (:pull:`104`). By `Nathan Miller`_.
- Added the ability to initialize the unknown vector (:pull:`109`). By `Nathan Miller`_.
- Added function that returns the size of the unknown vector (:pull:`109`). By `Nathan Miller`_.
- Generalized the mass-change evolution residual to not be just the mass change rate (:pull:`113`). By `Nathan Miller`_.

Breaking Changes
================
- Changed the micromorphic tools to use the vector Jacobian formulations and changed the micromorphic linear elasticity calculation to use vector Jacobian formulations (:pull:`81`). By `Nathan Miller`_.
- Changed Drucker Prager plasticity to use the vector Jacobian formulations (:pull:`81`). By `Nathan Miller`_.
- Changed hydra and hydraMicromorphic to use vector representations of the configurations and their jacobians (:pull:`82`). By `Nathan Miller`_.
- Added a required input for additionalDOF and previousAdditionalDOF to hydraBase and hydraBaseMicromorphic (:pull:`103`). By `Nathan Miller`_.

Internal Changes
================
- Removed extraneous semicolons (:pull:`69`). By `Nathan Miller`_.
- Changed the Jacobians to use row-major vector fomulation rather than vector of vectors (:pull:`77`). By `Nathan Miller`_.
- Changed the computation of the higher order yield surface to use row-major vector formation rather than vector of vectors (:pull:`78`). By `Nathan Miller`_.
- Updated to use the row-major vector Jacobians for tardigrade_constitutive_tools (:pull:`86`). By `Nathan Miller`_.
- Added definitions for common tensor sizes to the hydra base class (:pull:`87`). By `Nathan Miller`_.
- Changed inverses to fixed size where possible (:pull:`88`). By `Nathan Miller`_.
- Using constexpr instead of const when possible (:pull:`89`). By `Nathan Miller`_.
- Improved the efficiency of hydraBase (:pull:`90`). By `Nathan Miller`_.
- Improved the efficiency of hydraBaseMicromorphic (:pull:`91`). By `Nathan Miller`_.
- Improved the efficiency of tardigradeHydraMicromorphicDruckerPrager (:pull:`92`). By `Nathan Miller`_.
- Moved tardigrade_abaqus_tools.h from the header to the source file for tardigrade_hydra (:pull:`94`). By `Nathan Miller`_.
- Changed fatal error for non-full rank internal Jacobians to convergence errors (:pull:`95`). By `Nathan Miller`_.
- Changed additional fatal error for non-full rank internal Jacobians to convergence errors (:pull:`96`). By `Nathan Miller`_.
- Removed all sayHello tests (:pull:`97`). By `Nathan Miller`_.
- Improved performance of the linear elasticity subroutine (:pull:`99`). By `Nathan Miller`_.
- Using new error_tools check for error function (:pull:`100`). By `Nathan Miller`_.
- Changed Jacobian, dRdF, and dRdD to row-major vectors (:pull:`101`). By `Nathan Miller`_.
- Replaced queries to getUnknownVector purely to get the size of the vector (:pull:`109`). By `Nathan Miller`_.
- Added a better guess for the mass-change residual to improve convergence (:pull:`110`). By `Nathan Miller`_.
- Replaced the trapezoidal evolveF with the exponential map version (:pull:`111`). By `Nathan Miller`_.
- Rolled back exponential integrator for micromorphic (:pull:`114`). By `Nathan Miller`_.
- Added test for a fully directional integration where we know the answer (:pull:`117`). By `Nathan Miller`_.
- Added test for a fully spherical integration where we know the answer (:pull:`118`). By `Nathan Miller`_.
- Added test for when the mass-change rate is zero (:pull:`119`). By `Nathan Miller`_.
- Moved the Newton solve to its own function (:pull:`121`). By `Nathan Miller`_.
- Moved the preconditioned Newton solve to its own function (:pull:`122`). By `Nathan Miller`_.
- Changed the function calls for the Newton solve to a more general LHS and RHS form (:pull:`122`). By `Nathan Miller`_.
- Removed all of the calls to fuzzyEquals for the tests (:pull:`123`). By `Nathan Miller`_.

Bug Fixes
=========
- Corrected bug where the plastic state variable integration parameter was one minus the expected value (:pull:`71`). By `Nathan Miller`_.
- Corrected issue where libxsmm is not being used but was still required to be installed (:pull:`93`). By `Nathan Miller`_.
- Residuals setting initial guesses now force a reset of the current configurations (:pull:`110`). By `Nathan Miller`_.
- Direction vectors of length zero are now handled correctly (:pull:`116`). By `Nathan Miller`_.
- Removed extra whitespace in add_library from CMakeLists file (:pull:`120`). By `Nathan Miller`_.

******************
0.4.1 (01-24-2024)
******************

Release
=======
- Released version (:pull:`68`). By `Nathan Miller`_.

Internal Changes
================
- Removed unused variables (:pull:`67`). By `Nathan Miller`_.

******************
0.4.0 (01-24-2024)
******************

Release
=======
- Released version (:pull:`66`). By `Nathan Miller`_.

New Features
============
- Added setting the stresses and previous stresses for micromorphic linear elasticity (:pull:`54`). By `Nathan Miller`_.
- Added setting dRdT for micromorphic linear elasticity (:pull:`55`). By `Nathan Miller`_.
- Added calculations of the total derivative of the unknown vector (:pull:`57`). By `Nathan Miller`_.
- Added weakened Macaulay brackets (:pull:`62`). By `Nathan Miller`_.
- Added weakened state variable residuals (:pull:`63`). By `Nathan Miller`_.

Internal Changes
================
- Generalized the size of dRdF (:pull:`56`). by `Nathan Miller`_.
- Added the initialization of the unknown vector (:pull:`60`). By `Nathan Miller`_.
- Added dRdT to micromorphic Drucker Prager plasticity (:pull:`61`). By `Nathan Miller`_.
- Simplified the plastic multiplier residuals (:pull:`64`). By `Nathan Miller`_.
- Updated changelog for release (:pull:`65`). By `Nathan Miller`_.

Bug Fixes
=========
- Found problem with lack of generality when computing dRdF (:pull:`58`). By `Nathan Miller`_.
- Found issue with include guards for micromorphic Drucker-Prager plasticity (:pull:`59`). By `Nathan Miller`_.
- Found bug in the state variable residual Jacobians (:pull:`63`). By `Nathan Miller`_.
- Changed the plastic-multiplier residual so that it will attempt to force the plastic multipliers to be positive (:pull:`64`). By `Nathan Miller`_.

******************
0.3.1 (01-19-2024)
******************

Release
=======
- Released version (:pull:`53`). By `Nathan Miller`_.

New Features
============
- Added the micromorphic linear elasticity residual (:pull:`36`). By `Nathan Miller`_.
- Added the micromorphic Drucker Prager plasticity residual (:pull:`52`). By `Nathan Miller`_.

Internal Changes
================
- Added the ability to update the micromorphic hydra object with a new unknown vector (:pull:`34`). By `Nathan Miller`_.
- Added the calculation of the current stress measures in micromorphic linear elasticity (:pull:`35`). By `Nathan Miller`_.
- Initial commit of the micromorphic Drucker-Prager plasticity residual (:pull:`37`). By `Nathan Miller`_.
- Added the calculation of the driving stress for the micromorphic Drucker-Prager plasticity residual (:pull:`38`). By `Nathan Miller`_.
- Added the decomposition of the parameter vector (:pull:`39`). By `Nathan Miller`_.
- Added the extraction of the nonlinear state variables (:pull:`40`). By `Nathan Miller`_.
- Added the calculation of the cohesion (:pull:`41`). By `Nathan Miller`_.
- Added the calculation of the required quantities from the flow potential (:pull:`42`). By `Nathan Miller`_.
- Added the calculation of the jacobians of the strain-like ISV evolution rates (:pull:`43`). By `Nathan Miller`_.
- Added the calculation of the values and Jacobians of the strain-like ISVs (:pull:`44`). By `Nathan Miller`_.
- Moved the calculation of the preceding deformation gradient to its own function (:pull:`46`). By `Nathan Miller`_.
- Added a function to calculate the preceding micro deformation (:pull:`47`). By `Nathan Miller`_.
- Added the plastic velocity gradients for Drucker-Prager plasticity (:pull:`48`). By `Nathan Miller`_.
- Added functions to calculation the updated plastic deformations (:pull:`49`). By `Nathan Miller`_.
- Added updating the plastic deformation measures and their jacobians to the residual object (:pull:`50`). By `Nathan Miller`_.
- Added the residuals and jacobians of the state variables (:pull:`51`). By `Nathan Miller`_.

******************
0.3.0 (01-03-2024)
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
