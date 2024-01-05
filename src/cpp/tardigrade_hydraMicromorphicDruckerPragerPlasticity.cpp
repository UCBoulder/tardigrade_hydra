/**
  ******************************************************************************
  * \file tardigrade_hydraMicromorphicDruckerPragerPlasticity.cpp
  ******************************************************************************
  * An implementation of micromorphic linear elasticity using the hydra
  * framework. Used as an example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraMicromorphicDruckerPragerPlasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_micromorphic_tools.h>

namespace tardigradeHydra{

    namespace micromorphicDruckerPragerPlasticity{

        void residual::setMacroDrivingStress( ){
            /*!
             * Set the macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration)
             */

            setDrivingStresses( false );

        }

        void residual::setSymmetricMicroDrivingStress( ){
            /*!
             * Set the symmetric micro driving stress (stress in current configuration of plastic configuration)
             */

            setDrivingStresses( false );

        }

        void residual::setHigherOrderDrivingStress( ){
            /*!
             * Set the higher-order driving stress (stress in current configuration of plastic configuration)
             */

            setDrivingStresses( false );

        }

        void residual::setPreviousMacroDrivingStress( ){
            /*!
             * Set the previous macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration)
             */

            setDrivingStresses( true );

        }

        void residual::setPreviousSymmetricMicroDrivingStress( ){
            /*!
             * Set the previous symmetric micro driving stress (stress in current configuration of plastic configuration)
             */

            setDrivingStresses( true );

        }

        void residual::setPreviousHigherOrderDrivingStress( ){
            /*!
             * Set the previous higher-order driving stress (stress in current configuration of plastic configuration)
             */

            setDrivingStresses( true );

        }

        void residual::setDrivingStresses( const bool isPrevious ){
            /*!
             * Set the driving stresses for the plasticity
             * 
             * We here assume that the driving stresses are in the current configuration of
             * this residual's configuration.
             *
             * We also assume that the stress from hydra is in the reference configuration
             */

            const unsigned int *dim = hydra->getDimension( );

            const floatVector *stress;

            floatVector Fp;

            floatVector chip;

            if ( isPrevious ){

                stress = hydra->getPreviousStress( );

                Fp     = hydra->getPreviousFollowingConfiguration(      ( *getPlasticConfigurationIndex( ) ) - 1 );

                chip   = hydra->getPreviousFollowingMicroConfiguration( ( *getPlasticConfigurationIndex( ) ) - 1 );

            }
            else{

                stress = hydra->getStress( );

                Fp     = hydra->getFollowingConfiguration(      ( *getPlasticConfigurationIndex( ) ) - 1 );

                chip   = hydra->getFollowingMicroConfiguration( ( *getPlasticConfigurationIndex( ) ) - 1 );

            }

            // Extract the stresses from the stress vector
            floatVector PK2Stress(                     stress->begin( ),                           stress->begin( ) + 1 * ( *dim ) * ( *dim ) );

            floatVector referenceSymmetricMicroStress( stress->begin( ) + 1 * ( *dim ) * ( *dim ), stress->begin( ) + 2 * ( *dim ) * ( *dim ) );;

            floatVector referenceHigherOrderStress(    stress->begin( ) + 2 * ( *dim ) * ( *dim ), stress->begin( ) + 2 * ( *dim ) * ( *dim ) + ( *dim ) * ( *dim ) * ( *dim ) );

            // Push the stresses forward to the current configuration of the plastic configuration
            floatVector macroDrivingStress;

            floatVector symmetricMicroDrivingStress;

            floatVector higherOrderDrivingStress;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress, Fp, macroDrivingStress ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceSymmetricMicroStress, Fp, symmetricMicroDrivingStress ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, Fp, chip, higherOrderDrivingStress ) );

            if ( isPrevious ){

                set_previousMacroDrivingStress(          macroDrivingStress );

                set_previousSymmetricMicroDrivingStress( symmetricMicroDrivingStress );

                set_previousHigherOrderDrivingStress(    higherOrderDrivingStress );

            }
            else{

                set_macroDrivingStress(          macroDrivingStress );

                set_symmetricMicroDrivingStress( symmetricMicroDrivingStress );

                set_higherOrderDrivingStress(    higherOrderDrivingStress );

            }

        }

        void residual::setdMacroDrivingStressdMacroStress( ){
            /*!
             * Set the jacobian of the macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration) w.r.t. the macro stress
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdSymmetricMicroDrivingStressdMicroStress( ){
            /*!
             * Set the jacobian of the symmetric micro driving stress (stress in current configuration of plastic configuration) w.r.t. the micro stress
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdHigherOrderDrivingStressdHigherOrderStress( ){
            /*!
             * Set the jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdMacroDrivingStressdF( ){
            /*!
             * Set the jacobian of the macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration) w.r.t. the deformation gradient
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdSymmetricMicroDrivingStressdF( ){
            /*!
             * Set the jacobian of the symmetric micro driving stress (stress in current configuration of plastic configuration) w.r.t. the deformation gradient
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdHigherOrderDrivingStressdF( ){
            /*!
             * Set the jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the deformation gradient
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdHigherOrderDrivingStressdChi( ){
            /*!
             * Set the jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the micro deformation
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdMacroDrivingStressdFn( ){
            /*!
             * Set the jacobian of the macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration) w.r.t. the sub deformation gradients
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdSymmetricMicroDrivingStressdFn( ){
            /*!
             * Set the jacobian of the symmetric micro driving stress (stress in current configuration of plastic configuration) w.r.t. the sub deformation gradients
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdHigherOrderDrivingStressdFn( ){
            /*!
             * Set the jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the sub deformation gradients
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdHigherOrderDrivingStressdChin( ){
            /*!
             * Set the jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the sub micro deformations
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setDrivingStressesJacobians( const bool isPrevious ){
            /*!
             * Set the driving stresses for the plasticity along with the Jacobians
             * 
             * We here assume that the driving stresses are in the current configuration of
             * this residual's configuration.
             *
             * We also assume that the stress from hydra is in the reference configuration
             */

            const unsigned int *dim = hydra->getDimension( );

            const floatVector *stress;

            const floatMatrix dFpdSubFs;

            const floatMatrix *dF1dF;

            const floatMatrix *dF1dFn;

            const floatMatrix dChipdSubChis;

            const floatMatrix *dChi1dChi;

            const floatMatrix *dChi1dChin;

            floatVector Fp;

            floatVector chip;

            if ( isPrevious ){

                stress = hydra->getPreviousStress( );

                Fp     = hydra->getPreviousFollowingConfiguration(      ( *getPlasticConfigurationIndex( ) ) - 1 );

                chip   = hydra->getPreviousFollowingMicroConfiguration( ( *getPlasticConfigurationIndex( ) ) - 1 );

            }
            else{

                stress = hydra->getStress( );

                Fp     = hydra->getFollowingConfiguration(      ( *getPlasticConfigurationIndex( ) ) - 1 );

                chip   = hydra->getFollowingMicroConfiguration( ( *getPlasticConfigurationIndex( ) ) - 1 );

            }

            // Assemble the derivatives of the deformation gradient map
            floatMatrix dFpdF;

            floatMatrix dFpdFn;

            floatMatrix dChipdChi;

            floatMatrix dChipdChin;

            // Extract the stresses from the stress vector
            floatVector PK2Stress(                     stress->begin( ),                           stress->begin( ) + 1 * ( *dim ) * ( *dim ) );

            floatVector referenceSymmetricMicroStress( stress->begin( ) + 1 * ( *dim ) * ( *dim ), stress->begin( ) + 2 * ( *dim ) * ( *dim ) );;

            floatVector referenceHigherOrderStress(    stress->begin( ) + 2 * ( *dim ) * ( *dim ), stress->begin( ) + 2 * ( *dim ) * ( *dim ) + ( *dim ) * ( *dim ) * ( *dim ) );

            // Push the stresses forward to the current configuration of the plastic configuration
            floatVector macroDrivingStress;

            floatVector symmetricMicroDrivingStress;

            floatVector higherOrderDrivingStress;

            floatMatrix dMacrodFp;

            floatMatrix dMacrodPK2;

            floatMatrix dMicrodFp;

            floatMatrix dMicrodSigma;

            floatMatrix dHigherdFp;

            floatMatrix dHigherdChip;

            floatMatrix dHigherdM;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress, Fp, macroDrivingStress,
                                                                                                          dMacrodPK2, dMacrodFp ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceSymmetricMicroStress, Fp, symmetricMicroDrivingStress,
                                                                                                                     dMicrodSigma, dMicrodFp ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, Fp, chip, higherOrderDrivingStress,
                                                                                                                  dHigherdM, dHigherdFp, dHigherdChip ) );

            if ( isPrevious ){

                set_previousMacroDrivingStress(          macroDrivingStress );

                set_previousSymmetricMicroDrivingStress( symmetricMicroDrivingStress );

                set_previousHigherOrderDrivingStress(    higherOrderDrivingStress );

            }
            else{

                set_macroDrivingStress(          macroDrivingStress );

                set_symmetricMicroDrivingStress( symmetricMicroDrivingStress );

                set_higherOrderDrivingStress(    higherOrderDrivingStress );

                set_dMacroDrivingStressdMacroStress( dMacrodPK2 );

                set_dMacroDrivingStressdF( dMacrodF  );

                set_dMacroDrivingStressdFn( dMacrodFn );

                set_dSymmetricMicroDrivingStressdMicroStress( dMicrodSigma );

                set_dSymmetricMicroDrivingStressdF( dMicrodF );

                set_dSymmetricMicroDrivingStressdFn( dMicrodFn );

                set_dHigherOrderDrivingStressdHigherOrderStress( dHigherdM );

                set_dHigherOrderDrivingStressdF( dHigherdF );

                set_dHigherOrderDrivingStressdFn( dHigherdFn );

                set_dHigherOrderDrivingStressdChi( dHigherdChi );

                set_dHigherOrderDrivingStressdChin( dHigherdChin );

            }

        }
        const unsigned int* residual::getPlasticConfigurationIndex( ){
            /*!
             * Get plastic configuration index
             */

            return &_plasticConfigurationIndex;

        }

        const std::vector< unsigned int >* residual::getStateVariableIndices( ){
            /*!
             * Get state variable indices
             */

            return &_stateVariableIndices;

        }

        const floatType* residual::getIntegrationParameter( ){
            /*!
             * Get the integration parameter
             */

            return &_integrationParameter;

        }

    }

}
