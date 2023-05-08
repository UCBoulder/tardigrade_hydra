.. _changelog:


#########
Changelog
#########

******************
0.1.0 (unreleased)
******************

Breaking Changes
================
- Changed getSubConfiguration to not include the upper bound (:merge:`7`). By `Nathan Miller`_.

New Features
============
- Added calculation of the gradients of the current and previous F1 configurations (:merge:`11`). By `Nathan Miller`_.

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
