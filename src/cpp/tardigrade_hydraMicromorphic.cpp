/**
  ******************************************************************************
  * \file tardigrade_hydraMicromorphic.cpp
  ******************************************************************************
  * A C++ utility for constructing finite deformation micromorphic constitutive
  * models.
  ******************************************************************************
  */

#include<tardigrade_hydraMicromorphic.h>

namespace tardigradeHydra{

    hydraBaseMicromorphic::hydraBaseMicromorphic( const floatType &time, const floatType &deltaTime,
                                                  const floatType &temperature, const floatType &previousTemperature,
                                                  const floatVector &deformationGradient, const floatVector &previousDeformationGradient,
                                                  const floatVector &microDeformation, const floatVector &previousMicroDeformation,
                                                  const floatVector &gradientMicroDeformation, const floatVector &previousGradientMicroDeformation,
                                                  const floatVector &previousStateVariables, const floatVector &parameters,
                                                  const unsigned int numConfigurations, const unsigned int numNonLinearSolveStateVariables,
                                                  const unsigned int dimension, const floatType tolr, const floatType tola, const unsigned int maxIterations,
                                                  const unsigned int maxLSIterations, const floatType lsAlpha ){ }

}
