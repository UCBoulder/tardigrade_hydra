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

        void residual::setPreviousdMacroDrivingStressdMacroStress( ){
            /*!
             * Set the previous jacobian of the macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration) w.r.t. the macro stress
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdSymmetricMicroDrivingStressdMicroStress( ){
            /*!
             * Set the previous jacobian of the symmetric micro driving stress (stress in current configuration of plastic configuration) w.r.t. the micro stress
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderDrivingStressdHigherOrderStress( ){
            /*!
             * Set the previous jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdMacroDrivingStressdF( ){
            /*!
             * Set the previous jacobian of the macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration) w.r.t. the deformation gradient
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdSymmetricMicroDrivingStressdF( ){
            /*!
             * Set the previous jacobian of the symmetric micro driving stress (stress in current configuration of plastic configuration) w.r.t. the deformation gradient
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderDrivingStressdF( ){
            /*!
             * Set the previous jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the deformation gradient
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderDrivingStressdChi( ){
            /*!
             * Set the previous jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the micro deformation
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdMacroDrivingStressdFn( ){
            /*!
             * Set the previous jacobian of the macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration) w.r.t. the sub deformation gradients
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdSymmetricMicroDrivingStressdFn( ){
            /*!
             * Set the previous jacobian of the symmetric micro driving stress (stress in current configuration of plastic configuration) w.r.t. the sub deformation gradients
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderDrivingStressdFn( ){
            /*!
             * Set the previous jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the sub deformation gradients
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderDrivingStressdChin( ){
            /*!
             * Set the previous jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the sub micro deformations
             */

            setDrivingStressesJacobians( true );

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

            floatMatrix dFpdSubFs;

            const floatMatrix *dF1dF;

            const floatMatrix *dF1dFn;

            floatMatrix dChipdSubChis;

            const floatMatrix *dChi1dChi;

            const floatMatrix *dChi1dChin;

            floatVector Fp;

            floatVector chip;

            if ( isPrevious ){

                stress = hydra->getPreviousStress( );

                dF1dF         = hydra->get_previousdF1dF( );

                dF1dFn        = hydra->get_previousdF1dFn( );

                dFpdSubFs     = hydra->getPreviousFollowingConfigurationJacobian( ( *getPlasticConfigurationIndex( ) ) - 1 );

                dChi1dChi     = hydra->get_previousdChi1dChi( );

                dChi1dChin    = hydra->get_previousdChi1dChin( );

                dChipdSubChis = hydra->getPreviousFollowingMicroConfigurationJacobian( ( *getPlasticConfigurationIndex( ) ) - 1 );

                Fp            = hydra->getPreviousFollowingConfiguration(         ( *getPlasticConfigurationIndex( ) ) - 1 );

                chip          = hydra->getPreviousFollowingMicroConfiguration(    ( *getPlasticConfigurationIndex( ) ) - 1 );

            }
            else{

                stress = hydra->getStress( );

                dF1dF         = hydra->get_dF1dF( );

                dF1dFn        = hydra->get_dF1dFn( );

                dFpdSubFs     = hydra->getFollowingConfigurationJacobian( ( *getPlasticConfigurationIndex( ) ) - 1 );

                dChi1dChi     = hydra->get_dChi1dChi( );

                dChi1dChin    = hydra->get_dChi1dChin( );

                dChipdSubChis = hydra->getFollowingMicroConfigurationJacobian( ( *getPlasticConfigurationIndex( ) ) - 1 );

                Fp            = hydra->getFollowingConfiguration(         ( *getPlasticConfigurationIndex( ) ) - 1 );

                chip          = hydra->getFollowingMicroConfiguration(    ( *getPlasticConfigurationIndex( ) ) - 1 );

            }

            // Assemble the derivatives of the deformation gradient map
            floatMatrix dFpdF(  Fp.size( ), floatVector( hydra->getDeformationGradient( )->size( ), 0 ) );

            floatMatrix dFpdFn( Fp.size( ), floatVector( ( ( *hydra->getNumConfigurations( ) ) - 1 ) * hydra->getDeformationGradient( )->size( ), 0 ) );

            floatMatrix dChipdChi(  chip.size( ), floatVector( hydra->getDeformationGradient( )->size( ), 0 ) );

            floatMatrix dChipdChin( chip.size( ), floatVector( ( ( *hydra->getNumConfigurations( ) ) - 1 ) * hydra->getMicroDeformation( )->size( ), 0 ) );

            for ( unsigned int i = 0; i < ( *dim ) * ( *dim ); i++ ){

                for ( unsigned int j = 0; j < ( *dim ) * ( *dim ); j++ ){

                    for ( unsigned int k = 0; k < ( *dim ) * ( *dim ); k++ ){

                        dFpdF[ i ][ j ] += dFpdSubFs[ i ][ k ] * ( *dF1dF )[ k ][ j ];

                        dChipdChi[ i ][ j ] += dChipdSubChis[ i ][ k ] * ( *dChi1dChi )[ k ][ j ];

                    }

                }

                for ( unsigned int j = 0; j < ( ( *hydra->getNumConfigurations( ) ) - 1 ) * ( *dim ) * ( *dim ); j++ ){

                    dFpdFn[ i ][ j ] += dFpdSubFs[ i ][ j + ( *dim ) * ( *dim ) ];

                    dChipdChin[ i ][ j ] += dChipdSubChis[ i ][ j + ( *dim ) * ( *dim ) ];

                    for ( unsigned int k = 0; k < ( *dim ) * ( *dim ); k++ ){

                        dFpdFn[ i ][ j ] += dFpdSubFs[ i ][ k ] * ( *dF1dFn )[ k ][ j ];

                        dChipdChin[ i ][ j ] += dChipdSubChis[ i ][ k ] * ( *dChi1dChin )[ k ][ j ];

                    }

                }

            }

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

                set_previousdMacroDrivingStressdMacroStress( dMacrodPK2 );

                set_previousdMacroDrivingStressdF( tardigradeVectorTools::dot( dMacrodFp, dFpdF ) );

                set_previousdMacroDrivingStressdFn( tardigradeVectorTools::dot( dMacrodFp, dFpdFn ) );

                set_previousdSymmetricMicroDrivingStressdMicroStress( dMicrodSigma );

                set_previousdSymmetricMicroDrivingStressdF( tardigradeVectorTools::dot( dMicrodFp, dFpdF ) );

                set_previousdSymmetricMicroDrivingStressdFn( tardigradeVectorTools::dot( dMicrodFp, dFpdFn ) );

                set_previousdHigherOrderDrivingStressdHigherOrderStress( dHigherdM );

                set_previousdHigherOrderDrivingStressdF( tardigradeVectorTools::dot( dHigherdFp, dFpdF ) );

                set_previousdHigherOrderDrivingStressdFn( tardigradeVectorTools::dot( dHigherdFp, dFpdFn ) );

                set_previousdHigherOrderDrivingStressdChi( tardigradeVectorTools::dot( dHigherdChip, dChipdChi ) );

                set_previousdHigherOrderDrivingStressdChin( tardigradeVectorTools::dot( dHigherdChip, dChipdChin ) );

            }
            else{

                set_macroDrivingStress(          macroDrivingStress );

                set_symmetricMicroDrivingStress( symmetricMicroDrivingStress );

                set_higherOrderDrivingStress(    higherOrderDrivingStress );

                set_dMacroDrivingStressdMacroStress( dMacrodPK2 );

                set_dMacroDrivingStressdF( tardigradeVectorTools::dot( dMacrodFp, dFpdF ) );

                set_dMacroDrivingStressdFn( tardigradeVectorTools::dot( dMacrodFp, dFpdFn ) );

                set_dSymmetricMicroDrivingStressdMicroStress( dMicrodSigma );

                set_dSymmetricMicroDrivingStressdF( tardigradeVectorTools::dot( dMicrodFp, dFpdF ) );

                set_dSymmetricMicroDrivingStressdFn( tardigradeVectorTools::dot( dMicrodFp, dFpdFn ) );

                set_dHigherOrderDrivingStressdHigherOrderStress( dHigherdM );

                set_dHigherOrderDrivingStressdF( tardigradeVectorTools::dot( dHigherdFp, dFpdF ) );

                set_dHigherOrderDrivingStressdFn( tardigradeVectorTools::dot( dHigherdFp, dFpdFn ) );

                set_dHigherOrderDrivingStressdChi( tardigradeVectorTools::dot( dHigherdChip, dChipdChi ) );

                set_dHigherOrderDrivingStressdChin( tardigradeVectorTools::dot( dHigherdChip, dChipdChin ) );

            }

        }

        void residual::extractMaterialParameters( const parameterVector &parameters ){
            /*!
             * Extract the parameters from the parameter vector
             *
             * :param const std::vector< double > &parameters: The incoming parameter vector
             * :param parameterVector &macroHardeningParameters: The parameters used in the hardening of the macro Strain ISV
             *     (initial cohesion, hardening modulus)
             * :param parameterVector &microHardeningParameters: The parameters used in the hardening of the micro Strain ISV
             *     (initial cohesion, hardening modulus)
             * :param parameterVector &microGradientHardeningParameters: The parameters used in the hardening of the micro Gradient Strain ISV
             *     (initial cohesion, hardening modulus)
             * :param parameterVector &macroFlowParameters: The parameters used in the macro flow direction computation.
             *     (friction angle, beta )
             * :param parameterVector &microFlowParameters: The parameters used in the micro flow direction computation
             *     (friction angle, beta )
             * :param parameterVector &microGradientFlowParameters: The parameters used in the micro Gradient flow direction computation.
             *     (friction angle, beta )
             * :param parameterVector &macroYieldParameters: The parameters used in the macro yielding computation.
             *     (friction angle, beta )
             * :param parameterVector &microYieldParameters: The parameters used in the micro yielding computation
             *     (friction angle, beta )
             * :param parameterVector &microGradientYieldParameters: The parameters used in the micro Gradient yielding computation.
             *     (friction angle, beta )
             */
        
            if ( parameters.size() == 0 ){
        
                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "The parameter vector has a length of zero" ) );
        
            }
        
            unsigned int start = 0;
            unsigned int span;
        
            std::vector< parameterVector > outputs( 9 );
        
            //Extract the material parameters
            for ( unsigned int i = 0; i < outputs.size(); i++ ){
                span = ( unsigned int )std::floor( parameters[ start ]  + 0.5 ); //Extract the span of the parameter set
        
                if ( parameters.size() < start + 1 + span ){
                    std::string outstr = "fparams is not long enough to contain all of the required parameters:\n";
                    outstr +=            "    filling variable " + std::to_string( i ) + "\n";
                    outstr +=            "    size =          "  + std::to_string( parameters.size() ) + "\n";
                    outstr +=            "    required size = "  + std::to_string( start + 1 + span );
        
                    TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( outstr ); )
        
                }
        
                outputs[ i ] = parameterVector( parameters.begin() + start + 1, parameters.begin() + start + 1 + span );
        
                start = start + 1 + span;
            }

            //Set the output values
            set_macroHardeningParameters(         outputs[ 0 ] );
        
            set_microHardeningParameters(         outputs[ 1 ] );
        
            set_microGradientHardeningParameters( outputs[ 2 ] );
        
            set_macroFlowParameters(              outputs[ 3 ] );
        
            set_microFlowParameters(              outputs[ 4 ] );
        
            set_microGradientFlowParameters(      outputs[ 5 ] );
        
            set_macroYieldParameters(             outputs[ 6 ] );
        
            set_microYieldParameters(             outputs[ 7 ] );
        
            set_microGradientYieldParameters(     outputs[ 8 ] );
        
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

        void residual::setPlasticStateVariables( const bool isPrevious ){
            /*!
             * Set the plastic state variables
             * 
             * \param isPrevious: Flag for whether to set the current (false) or previous (true) values of the plastic state variables
             */

            floatVector plasticStateVariables( getStateVariableIndices( )->size( ), 0 );

            const floatVector *nonlinearISVs;

            if ( isPrevious ){

                nonlinearISVs = hydra->get_previousNonLinearSolveStateVariables( );

            }
            else{

                nonlinearISVs = hydra->get_nonLinearSolveStateVariables( );

            }

            for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){

                plasticStateVariables[ ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ] = ( *nonlinearISVs )[ *ind ];

            }

            if ( isPrevious ){

                set_previousPlasticStateVariables( plasticStateVariables );

            }
            else{

                set_plasticStateVariables( plasticStateVariables );

            }

        }

        void residual::setPlasticStateVariables( ){
            /*!
             * Set the plastic state variables
             */

            setPlasticStateVariables( false );

        }

        void residual::setPreviousPlasticStateVariables( ){
            /*!
             * Set the previous plastic state variables
             */

            setPlasticStateVariables( true );

        }

        void residual::setPlasticMultipliers( const bool isPrevious ){
            /*!
             * Set the plastic multipliers
             * 
             * \param isPrevious: Flag for whether to set the current (false) or previous (true) values of the plastic multipliers
             */

            floatVector plasticMultipliers( *getNumPlasticMultipliers( ), 0 );

            const floatVector *nonlinearISVs;

            if ( isPrevious ){

                nonlinearISVs = hydra->get_previousNonLinearSolveStateVariables( );

            }
            else{

                nonlinearISVs = hydra->get_nonLinearSolveStateVariables( );

            }

            for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->begin( ) + *getNumPlasticMultipliers( ); ind++ ){

                plasticMultipliers[ ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ] = ( *nonlinearISVs )[ *ind ];

            }

            if ( isPrevious ){

                set_previousPlasticMultipliers( plasticMultipliers );

            }
            else{

                set_plasticMultipliers( plasticMultipliers );

            }

        }

        void residual::setPlasticMultipliers( ){
            /*!
             * Set the plastic multipliers
             */

            setPlasticMultipliers( false );

        }

        void residual::setPreviousPlasticMultipliers( ){
            /*!
             * Set the previous plastic multipliers
             */

            setPlasticMultipliers( true );

        }

        void residual::setPlasticStrainLikeISVs( const bool isPrevious ){
            /*!
             * Set the plastic strain-like internal state variables
             * 
             * \param isPrevious: Flag for whether to set the current (false) or previous (true) values of the plastic multipliers
             */

            floatVector plasticStrainLikeISVs( getStateVariableIndices( )->size( ) - *getNumPlasticMultipliers( ), 0 );

            const floatVector *nonlinearISVs;

            if ( isPrevious ){

                nonlinearISVs = hydra->get_previousNonLinearSolveStateVariables( );

            }
            else{

                nonlinearISVs = hydra->get_nonLinearSolveStateVariables( );

            }

            for ( auto ind = getStateVariableIndices( )->begin( ) + *getNumPlasticMultipliers( ); ind != getStateVariableIndices( )->end( ); ind++ ){

                plasticStrainLikeISVs[ ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) - *getNumPlasticMultipliers( ) ] = ( *nonlinearISVs )[ *ind ];

            }

            if ( isPrevious ){

                set_previousPlasticStrainLikeISVs( plasticStrainLikeISVs );

            }
            else{

                set_plasticStrainLikeISVs( plasticStrainLikeISVs );

            }

        }

        void residual::setPlasticStrainLikeISVs( ){
            /*!
             * Set the plastic strain-like isvs
             */

            setPlasticStrainLikeISVs( false );

        }

        void residual::setPreviousPlasticStrainLikeISVs( ){
            /*!
             * Set the previous plastic strain-like isvs
             */

            setPlasticStrainLikeISVs( true );

        }

        void residual::setMacroCohesion( ){
            /*!
             * Set the macro cohesion
             */

            setCohesions( false );

        }

        void residual::setMicroCohesion( ){
            /*!
             * Set the micro cohesion
             */

            setCohesions( false );

        }

        void residual::setMicroGradientCohesion( ){
            /*!
             * Set the micro gradient cohesion
             */

            setCohesions( false );

        }

        void residual::setPreviousMacroCohesion( ){
            /*!
             * Set the previous macro cohesion
             */

            setCohesions( true );

        }

        void residual::setPreviousMicroCohesion( ){
            /*!
             * Set the previous macro cohesion
             */

            setCohesions( true );

        }

        void residual::setPreviousMicroGradientCohesion( ){
            /*!
             * Set the micro gradient cohesion
             */

            setCohesions( true );

        }

        void residual::setCohesions( const bool isPrevious ){
            /*!
             * Set the values of the cohesion
             * 
             * \param isPrevious: Flag for whether to compute the current (false) or previous (true) cohesions
             */

            const floatVector *plasticStrainLikeISVs;

            if ( isPrevious ){

                plasticStrainLikeISVs = get_previousPlasticStrainLikeISVs( );

            }
            else{

                plasticStrainLikeISVs = get_plasticStrainLikeISVs( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( get_macroHardeningParameters( )->size( ) != 2 ){
    
                    throw std::runtime_error( "The micro hardening parameters must have a length of 2 rather than " + std::to_string( get_macroHardeningParameters( )->size( ) ) );
    
                }
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( get_microHardeningParameters( )->size( ) != 2 ){
    
                    throw std::runtime_error( "The micro hardening parameters must have a length of 2 rather than " + std::to_string( get_microHardeningParameters( )->size( ) ) );
    
                }
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( get_microGradientHardeningParameters( )->size( ) != 2 ){
    
                    throw std::runtime_error( "The micro hardening parameters must have a length of 2 rather than " + std::to_string( get_microGradientHardeningParameters( )->size( ) ) );
    
                }
            )

            floatType macroCohesion           = ( *get_macroHardeningParameters( ) )[ 0 ] + ( *get_macroHardeningParameters( ) )[ 1 ] * ( *plasticStrainLikeISVs )[ 0 ];

            floatType microCohesion           = ( *get_microHardeningParameters( ) )[ 0 ] + ( *get_microHardeningParameters( ) )[ 1 ] * ( *plasticStrainLikeISVs )[ 1 ];

            floatVector microGradientCohesion = ( *get_microGradientHardeningParameters( ) )[ 0 ] + ( *get_microGradientHardeningParameters( ) )[ 1 ] * floatVector( plasticStrainLikeISVs->begin( ) + 2,
                                                                                                                                                                     plasticStrainLikeISVs->end( ) );

            if ( isPrevious ){

                set_previousMacroCohesion( macroCohesion );

                set_previousMicroCohesion( microCohesion );

                set_previousMicroGradientCohesion( microGradientCohesion );

            }
            else{

                set_macroCohesion( macroCohesion );

                set_microCohesion( microCohesion );

                set_microGradientCohesion( microGradientCohesion );

            }

        }

    }

}
