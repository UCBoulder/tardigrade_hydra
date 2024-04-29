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

            hydraBaseMicromorphic( ){ }

            hydraBaseMicromorphic( const floatType &time, const floatType &deltaTime,
                                   const floatType &temperature, const floatType &previousTemperature,
                                   const floatVector &deformationGradient, const floatVector &previousDeformationGradient,
                                   const floatVector &microDeformation, const floatVector &previousMicroDeformation,
                                   const floatVector &gradientMicroDeformation, const floatVector &previousGradientMicroDeformation,
                                   const floatVector &previousStateVariables, const floatVector &parameters,
                                   const unsigned int numConfigurations, const unsigned int numNonLinearSolveStateVariables,
                                   const unsigned int dimension=3, const unsigned int configuration_unknown_count=45,
                                   const floatType tolr=1e-9, const floatType tola=1e-9, const unsigned int maxIterations=20,
                                   const unsigned int maxLSIterations=5, const floatType lsAlpha=1e-4,
                                   const bool use_preconditioner=false, const unsigned int preconditioner_type=0 );

            //! Get the current micro-deformation tensor
            const floatVector *getMicroDeformation( ){ return &_microDeformation; }

            //! Get the previous micro-deformation tensor
            const floatVector *getPreviousMicroDeformation( ){ return &_previousMicroDeformation; }

            //! Get the current spatial gradient w.r.t. the reference configuration of the micro-deformation tensor
            const floatVector *getGradientMicroDeformation( ){ return &_gradientMicroDeformation; }

            //! Get the previous spatial gradient w.r.t. the reference configuration of the micro-deformation tensor
            const floatVector *getPreviousGradientMicroDeformation( ){ return &_previousGradientMicroDeformation; }

            floatVector getSubMicroConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getPrecedingMicroConfiguration( const unsigned int &index );

            floatVector getFollowingMicroConfiguration( const unsigned int &index );

            floatVector getMicroConfiguration( const unsigned int &index );

            floatVector getPreviousSubMicroConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getPreviousPrecedingMicroConfiguration( const unsigned int &index );

            floatVector getPreviousFollowingMicroConfiguration( const unsigned int &index );

            floatVector getPreviousMicroConfiguration( const unsigned int &index );

            floatVector getSubMicroConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getPrecedingMicroConfigurationJacobian( const unsigned int &index );

            floatVector getFollowingMicroConfigurationJacobian( const unsigned int &index );

            floatVector getPreviousSubMicroConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getPreviousPrecedingMicroConfigurationJacobian( const unsigned int &index );

            floatVector getPreviousFollowingMicroConfigurationJacobian( const unsigned int &index );

            const floatVector *getFlatdXdD( ){
                /*!
                 * Get the total derivative of the unknown vector w.r.t. the deformation.
                 *
                 * Pass-through to getFlatdXdF just changing the naming convention
                 */

                return hydraBase::getFlatdXdF( );

            }
        protected:

            //Utility functions
            virtual void initializeUnknownVector( ) override;

            virtual void decomposeUnknownVector( ) override;

            virtual void decomposeUnknownVectorMicroConfigurations( );

            virtual void decomposeStateVariableVector( ) override;

            virtual void decomposeStateVariableVectorMicroConfigurations( );

        private:

            floatVector _microDeformation; //!< The current micro-deformation

            floatVector _previousMicroDeformation; //!< The previous micro-deformation

            floatVector _gradientMicroDeformation; //!< The spatial gradient of the micro-deformation w.r.t. the reference coordinates

            floatVector _previousGradientMicroDeformation; //!< The previous spatial gradient of the micro-deformation w.r.t. the reference coordinates

            void setFirstMicroConfigurationJacobians( );

            void setPreviousFirstMicroConfigurationJacobians( );

            void setFirstGradientMicroConfigurationJacobians( );

            void setPreviousFirstGradientMicroConfigurationJacobians( );

            void computeGradientMicroConfigurations( const floatVector *data_vector, unsigned int start_index,
                                                     const floatVector &configurations, const floatVector &microConfigurations,
                                                     const floatVector &gradientMicroConfiguration, floatVector &gradientMicroConfigurations );

            void calculateFirstConfigurationGradChi( const floatVector &configurations, const floatVector &microConfigurations,
                                                     const floatVector &gradientMicroConfiguration, floatVector &gradientMicroConfigurations );

            void calculateFirstConfigurationGradChiJacobian( const floatVector &configurations, const floatVector &microConfigurations,
                                                             const floatVector &gradientMicroConfiguration, const floatVector &gradientMicroConfigurations,
                                                             const floatVector &dChi1dChi, const floatVector &dChi1dChin,
                                                             floatVector &dGradChi1dCn,
                                                             floatVector &dGradChi1dChi,     floatVector &dGradChi1dChin,
                                                             floatVector &dGradChi1dGradChi, floatVector &dGradChi1dGradChin );

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, microConfigurations,                 floatVector, unexpectedError )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, inverseMicroConfigurations,          floatVector, unexpectedError )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, gradientMicroConfigurations,         floatVector, unexpectedError )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousMicroConfigurations,         floatVector, unexpectedError )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousInverseMicroConfigurations,  floatVector, unexpectedError )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousGradientMicroConfigurations, floatVector, unexpectedError )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dChi1dChi,                           floatVector, setFirstMicroConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dChi1dChin,                          floatVector, setFirstMicroConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdChi1dChi,                   floatVector, setPreviousFirstMicroConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousdChi1dChin,                  floatVector, setPreviousFirstMicroConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGradChi1dFn,                        floatVector, setFirstGradientMicroConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGradChi1dChi,                       floatVector, setFirstGradientMicroConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGradChi1dChin,                      floatVector, setFirstGradientMicroConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGradChi1dGradChi,                   floatVector, setFirstGradientMicroConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dGradChi1dGradChin,                  floatVector, setFirstGradientMicroConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE( private, previousdGradChi1dFn,                 floatVector, setPreviousFirstGradientMicroConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE( private, previousdGradChi1dChi,                floatVector, setPreviousFirstGradientMicroConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE( private, previousdGradChi1dChin,               floatVector, setPreviousFirstGradientMicroConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE( private, previousdGradChi1dGradChi,            floatVector, setPreviousFirstGradientMicroConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE( private, previousdGradChi1dGradChin,           floatVector, setPreviousFirstGradientMicroConfigurationJacobians )

    };

    //! The base class for micromorphic residuals
    class residualBaseMicromorphic : public residualBase{

        public:

            using tardigradeHydra::residualBase::residualBase;

            /*!
             * Base class for micromorphic residuals
             * 
             * \param *_hydra: A pointer to the containing hydra object
             * \param _numEquations: The number of equations the residual defines
             */
            residualBaseMicromorphic( hydraBaseMicromorphic *_hydra, unsigned int _numEquations ) : residualBase( _hydra, _numEquations ), hydra( _hydra ){ }

            hydraBaseMicromorphic *hydra; //!< A pointer to the containing hydra object

            virtual void setdRdD( ){
                /*!
                 * Set the derivative of the residual w.r.t. the deformation.
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The derivative of the residual w.r.t. the deformation is not implemented" ) );

            }

            virtual void setdRdF( ) override {
                /*!
                 * Rename setdRdF to setdRdD because we will use it for all of the deformations
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( setdRdD( ) );

            }

            void setdRdD( const floatVector &dRdD ){
                /*!
                 * Set the derivative of the residual w.r.t. the deformation.
                 *
                 * Pass-through to setdRdF just changing the naming convention
                 * 
                 * \param &dRdD: The derivative of the resdual with respect to the deformation (F, chi, gradChi )
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( residualBase::setdRdF( dRdD ) );

            }

            const floatVector *getdRdD( ){
                /*!
                 * Get the derivative of the residual w.r.t. the deformation.
                 *
                 * Pass-through to getdRdF just changing the naming convention
                 */

                return residualBase::getdRdF( );

            }

    };

}

#endif
