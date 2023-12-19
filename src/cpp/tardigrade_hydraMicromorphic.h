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

            //! Get a reference to the gradient of the micro configurations
            const floatMatrix* getGradientMicroConfigurations( ){ return &_gradientMicroConfigurations.second; }

            //! Get a reference to the previous micro configurations
            const floatMatrix* getPreviousMicroConfigurations( ){ return &_previousMicroConfigurations.second; }

            //! Get a reference to the previous inverse micro configurations
            const floatMatrix* getPreviousInverseMicroConfigurations( ){ return &_previousInverseMicroConfigurations.second; }

            //! Get a reference to the gradient of the previous micro configurations
            const floatMatrix* getPreviousGradientMicroConfigurations( ){ return &_previousGradientMicroConfigurations.second; }

            floatVector getSubMicroConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getPrecedingMicroConfiguration( const unsigned int &index );

            floatVector getFollowingMicroConfiguration( const unsigned int &index );

            floatVector getMicroConfiguration( const unsigned int &index );

            floatVector getPreviousSubMicroConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getPreviousPrecedingMicroConfiguration( const unsigned int &index );

            floatVector getPreviousFollowingMicroConfiguration( const unsigned int &index );

            floatVector getPreviousMicroConfiguration( const unsigned int &index );

            floatMatrix getSubMicroConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatMatrix getPrecedingMicroConfigurationJacobian( const unsigned int &index );

            floatMatrix getFollowingMicroConfigurationJacobian( const unsigned int &index );

            floatMatrix getPreviousSubMicroConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatMatrix getPreviousPrecedingMicroConfigurationJacobian( const unsigned int &index );

            floatMatrix getPreviousFollowingMicroConfigurationJacobian( const unsigned int &index );

            const floatMatrix *getdChi1dChi( );

            const floatMatrix *getdChi1dChin( );

            const floatMatrix *getPreviousdChi1dChi( );

            const floatMatrix *getPreviousdChi1dChin( );

            const floatMatrix *getdGradChi1dFn( );

            const floatMatrix *getdGradChi1dChi( );

            const floatMatrix *getdGradChi1dChin( );

            const floatMatrix *getdGradChi1dGradChi( );

            const floatMatrix *getdGradChi1dGradChin( );

            const floatMatrix *getPreviousdGradChi1dFn( );

            const floatMatrix *getPreviousdGradChi1dChi( );

            const floatMatrix *getPreviousdGradChi1dChin( );

            const floatMatrix *getPreviousdGradChi1dGradChi( );

            const floatMatrix *getPreviousdGradChi1dGradChin( );

        protected:
            //Setter functions
            void setMicroConfigurations( const floatMatrix &microConfigurations );

            void setInverseMicroConfigurations( const floatMatrix &inverseMicroConfigurations );

            void setGradientMicroConfigurations( const floatMatrix &gradientMicroConfigurations );

            void setPreviousMicroConfigurations( const floatMatrix &previousMicroConfigurations );

            void setPreviousInverseMicroConfigurations( const floatMatrix &previousInverseMicroConfigurations );

            void setPreviousGradientMicroConfigurations( const floatMatrix &previousGradientMicroConfigurations );

            void setdGradChi1dFn( const floatMatrix &dGradChi1dFn );

            void setdGradChi1dChi( const floatMatrix &dGradChi1dChi );

            void setdGradChi1dChin( const floatMatrix &dGradChi1dChin );

            void setdGradChi1dGradChi( const floatMatrix &dGradChi1dGradChi );

            void setdGradChi1dGradChin( const floatMatrix &dGradChi1dGradChin );

            void setPreviousdGradChi1dFn( const floatMatrix &previousdGradChi1dFn );

            void setPreviousdGradChi1dChi( const floatMatrix &previousdGradChi1dChi );

            void setPreviousdGradChi1dChin( const floatMatrix &previousdGradChi1dChin );

            void setPreviousdGradChi1dGradChi( const floatMatrix &previousdGradChi1dGradChi );

            void setPreviousdGradChi1dGradChin( const floatMatrix &previousdGradChi1dGradChin );

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

            dataStorage< floatMatrix > _gradientMicroConfigurations; //!< The spatial gradients of the micro-configurations w.r.t. their reference configurations

            dataStorage< floatMatrix > _previousMicroConfigurations; //!< The previous values of the micro-configurations

            dataStorage< floatMatrix > _previousInverseMicroConfigurations; //!< The previous values of the inverse micro-configurations

            dataStorage< floatMatrix > _previousGradientMicroConfigurations; //!< The spatial gradients of the previous micro-configurations w.r.t. their reference configurations

            dataStorage< floatMatrix > _dChi1dChi; //!< The Jacobian of the first micro-configuration w.r.t. the total micro-configuration

            dataStorage< floatMatrix > _dChi1dChin; //!< The Jacobian of the first micro-configuration w.r.t. the remaining micro-configurations

            dataStorage< floatMatrix > _previousdChi1dChi; //!< The Jacobian of the previous first micro-configuration w.r.t. the total micro-configuration

            dataStorage< floatMatrix > _previousdChi1dChin; //!< The Jacobian of the previous first micro-configuration w.r.t. the remaining micro-configurations

            dataStorage< floatMatrix > _dGradChi1dFn; //!< The Jacobian of the first spatial gradient of the micro-configuration w.r.t. the remaining sub-configurations

            dataStorage< floatMatrix > _dGradChi1dChi; //!< The Jacobian of the first spatial gradient of the micro-configuration w.r.t. the total micro-configuration

            dataStorage< floatMatrix > _dGradChi1dChin; //!< The Jacobian of the first spatial gradient of the micro-configuration w.r.t. the remaining sub micro-configurations

            dataStorage< floatMatrix > _dGradChi1dGradChi; //!< The Jacobian of the first spatial gradient of the micro-configuration w.r.t. the total spatial gradient of the micro-configuration

            dataStorage< floatMatrix > _dGradChi1dGradChin; //!< The Jacobian of the first spatial gradient of the micro-configuration w.r.t. the remaining spatial gradients of the sub micro-configurations

            dataStorage< floatMatrix > _previousdGradChi1dFn; //!< The Jacobian of the previous first spatial gradient of the micro-configuration w.r.t. the remaining sub-configurations

            dataStorage< floatMatrix > _previousdGradChi1dChi; //!< The Jacobian of the previous first spatial gradient of the micro-configuration w.r.t. the total micro-configuration

            dataStorage< floatMatrix > _previousdGradChi1dChin; //!< The Jacobian of the previous first spatial gradient of the micro-configuration w.r.t. the remaining sub micro-configurations

            dataStorage< floatMatrix > _previousdGradChi1dGradChi; //!< The Jacobian of the previous first spatial gradient of the micro-configuration w.r.t. the total spatial gradient of the micro-configuration

            dataStorage< floatMatrix > _previousdGradChi1dGradChin; //!< The Jacobian of the previous first spatial gradient of the micro-configuration w.r.t. the remaining spatial gradients of the sub micro-configurations

            void setFirstMicroConfigurationJacobians( );

            void setPreviousFirstMicroConfigurationJacobians( );

            void setFirstGradientMicroConfigurationJacobians( );

            void setPreviousFirstGradientMicroConfigurationJacobians( );

            void setdChi1dChi( const floatMatrix &dChi1dChi );

            void setdChi1dChin( const floatMatrix &dChi1dChin );

            void setPreviousdChi1dChi( const floatMatrix &previousdChi1dChi );

            void setPreviousdChi1dChin( const floatMatrix &previousdChi1dChin );

            void computeGradientMicroConfigurations( const floatVector *data_vector, unsigned int start_index,
                                                     const floatMatrix &configurations, const floatMatrix &microConfigurations,
                                                     const floatVector &gradientMicroConfiguration, floatMatrix &gradientMicroConfigurations );

            void calculateFirstConfigurationGradChi( const floatMatrix &configurations, const floatMatrix &microConfigurations,
                                                     const floatVector &gradientMicroConfiguration, floatMatrix &gradientMicroConfigurations );

            void calculateFirstConfigurationGradChiJacobian( const floatMatrix &configurations, const floatMatrix &microConfigurations,
                                                             const floatVector &gradientMicroConfiguration, const floatMatrix &gradientMicroConfigurations,
                                                             const floatMatrix &dChi1dChi, const floatMatrix &dChi1dChin,
                                                             floatMatrix &dGradChi1dCn,
                                                             floatMatrix &dGradChi1dChi, floatMatrix &dGradChi1dChin,
                                                             floatMatrix &dGradChi1dGradChi, floatMatrix &dGradChi1dGradChin );

    };

    class residualBaseMicromorphic : public residualBase{

        public:

            /*!
             * Base class for micromorphic residuals
             * 
             * \param *_hydra: A pointer to the containing hydra object
             * \param _numEquations: The number of equations the residual defines
             */
            residualBaseMicromorphic( hydraBaseMicromorphic *_hydra, unsigned int _numEquations ) : residualBase( _hydra, _numEquations ), hydra( _hydra ){ }

            hydraBaseMicromorphic *hydra;

    };

}

#endif
