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
                                   const unsigned int dimension=3, const unsigned int configuration_unknown_count=45,
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

            //! Get a reference to the micro configurations
            const floatMatrix* getMicroConfigurations( ){ return &_microConfigurations.second; }

            //! Get a reference to the inverse micro configurations
            const floatMatrix* getInverseMicroConfigurations( ){ return &_inverseMicroConfigurations.second; }

            //! Get a reference to the previous micro configurations
            const floatMatrix* getPreviousMicroConfigurations( ){ return &_previousMicroConfigurations.second; }

            //! Get a reference to the previous inverse micro configurations
            const floatMatrix* getPreviousInverseMicroConfigurations( ){ return &_previousInverseMicroConfigurations.second; }

            floatVector getSubMicroConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getPrecedingMicroConfiguration( const unsigned int &index );

            floatVector getFollowingMicroConfiguration( const unsigned int &index );

            floatVector getMicroConfiguration( const unsigned int &index );

            floatVector getPreviousSubMicroConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getPreviousPrecedingMicroConfiguration( const unsigned int &index );

            floatVector getPreviousFollowingMicroConfiguration( const unsigned int &index );

            floatVector getPreviousMicroConfiguration( const unsigned int &index );

        protected:
            //Setter functions
            void setMicroConfigurations( const floatMatrix &microConfigurations );

            void setInverseMicroConfigurations( const floatMatrix &inverseMicroConfigurations );

            void setPreviousMicroConfigurations( const floatMatrix &previousMicroConfigurations );

            void setPreviousInverseMicroConfigurations( const floatMatrix &previousInverseMicroConfigurations );

            //Utility functions
            virtual void decomposeStateVariableVector( ) override;

            virtual void decomposeStateVariableVectorMicroConfigurations( );

        private:

            floatVector _microDeformation; //!< The current micro-deformation

            floatVector _previousMicroDeformation; //!< The previous micro-deformation

            floatVector _gradientMicroDeformation; //!< The spatial gradient of the micro-deformation w.r.t. the reference coordinates

            floatVector _previousGradientMicroDeformation; //!< The previous spatial gradient of the micro-deformation w.r.t. the reference coordinates

            dataStorage< floatMatrix > _microConfigurations; //!< The current values of the micro-configurations

            dataStorage< floatMatrix > _inverseMicroConfigurations; //!< The current values of the inverse micro-configurations

            dataStorage< floatMatrix > _previousMicroConfigurations; //!< The previous values of the micro-configurations

            dataStorage< floatMatrix > _previousInverseMicroConfigurations; //!< The previous values of the inverse micro-configurations

    };

}

#endif
