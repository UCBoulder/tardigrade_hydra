/**
  ******************************************************************************
  * \file tardigrade_hydraMicromorphic.h
  ******************************************************************************
  * A C++ utility for constructing finite deformation micromorphic constitutive
  * models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_MICROMORPHIC_H
#define TARDIGRADE_HYDRA_MICROMORPHIC_H

#include<tardigrade_hydra.h>

namespace tardigradeHydra{

    //! The base class for hydra framework micromorphic material models
    class hydraBaseMicromorphic : public hydraBase{

        public:
            hydraBaseMicromorphic( const floatType &time, const floatType &deltaTime,
                                   const floatType &temperature, const floatType &previousTemperature,
                                   const floatVector &deformationGradient, const floatVector &previousDeformationGradient,
                                   const floatVector &microDeformation, const floatVector &previousMicroDeformation,
                                   const floatVector &gradientMicroDeformation, const floatVector &previousGradientMicroDeformation,
                                   const floatVector &previousStateVariables, const floatVector &parameters,
                                   const unsigned int numConfigurations, const unsigned int numNonLinearSolveStateVariables,
                                   const unsigned int dimension=3, const unsigned int configuration_unknown_count=27,
                                   const floatType tolr=1e-9, const floatType tola=1e-9, const unsigned int maxIterations=20,
                                   const unsigned int maxLSIterations=5, const floatType lsAlpha=1e-4 );

            //! Get the current micro-deformation tensor
            const floatVector *getMicroDeformation( ){ return &_microDeformation; }

            //! Get the previous micro-deformation tensor
            const floatVector *getPreviousMicroDeformation( ){ return &_previousMicroDeformation; }

            //! Get the current spatial gradient w.r.t. the reference configuration of the micro-deformation tensor
            const floatVector *getGradientMicroDeformation( ){ return &_gradientMicroDeformation; }

            //! Get the previous spatial gradient w.r.t. the reference configuration of the micro-deformation tensor
            const floatVector *getPreviousGradientMicroDeformation( ){ return &_previousGradientMicroDeformation; }

        private:

            floatVector _microDeformation; //!< The current micro-deformation

            floatVector _previousMicroDeformation; //!< The previous micro-deformation

            floatVector _gradientMicroDeformation; //!< The spatial gradient of the micro-deformation w.r.t. the reference coordinates

            floatVector _previousGradientMicroDeformation; //!< The previous spatial gradient of the micro-deformation w.r.t. the reference coordinates

    };

}

#endif
