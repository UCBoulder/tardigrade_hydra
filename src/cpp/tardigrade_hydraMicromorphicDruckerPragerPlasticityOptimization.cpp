/**
  ******************************************************************************
  * \file tardigrade_hydraMicromorphicDruckerPragerPlasticityOptimization.cpp
  ******************************************************************************
  * An implementation of micromorphic Drucker-Prager plasticity based on an
  * optimization approach.
  ******************************************************************************
  */

#include<tardigrade_hydraMicromorphicDruckerPragerPlasticityOptimization.h>

namespace tardigradeHydra{

    namespace micromorphicDruckerPragerPlasticityOptimization{

        void residual::setStateVariableResiduals( ){
            /*!
             * Set the state variable residuals
             *
             * We define these residuals as
             *
             * \f$R = \dot{\bar{\gamma}} f \f$
             *
             * and
             *
             * \f$R = Z^{t+1} - Z^{t+1,\text{trial}}\f$
             *
             * Because we have five plastic multipliers (\f$\gamma\f$), five strain-like state variables (\f$Z\f$)
             * and five yield surfaces (\f$f\f$) we can define ten different equations. The remaining places are for
             * slack variables.
             */

            const floatVector *plasticStrainLikeISVs = get_plasticStrainLikeISVs( );

            const floatVector *updatedPlasticStrainLikeISVs = get_updatedPlasticStrainLikeISVs( );

            const unsigned int numPlasticMultipliers = *getNumPlasticMultipliers( );

            unsigned int numPlasticStrainLikeISVs = plasticStrainLikeISVs->size( );

            const floatVector *plasticMultipliers = get_plasticMultipliers( );

            const floatVector slackVariables( get_plasticStateVariables( )->begin( ) + numPlasticMultipliers + numPlasticStrainLikeISVs,
                                              get_plasticStateVariables( )->end( ) );

            const floatType *macroYield = get_macroYield( );

            const floatType *microYield = get_microYield( );

            const dimVector *microGradientYield = get_microGradientYield( );

            auto residual = get_setDataStorage_stateVariableResiduals( );
            residual.zero( get_plasticStateVariables( )->size( ) );

            // Add the consistency condition
            unsigned int row0 = 0;
            for ( unsigned int i = 0; i < numPlasticMultipliers; i++ ){
                ( *residual.value )[ row0 + i ] = ( *plasticMultipliers )[ i ] * slackVariables[ i ];
            }

            // Set the terms associated with the strain-like ISV evolution
            row0 = numPlasticMultipliers;
            for ( unsigned int i = 0; i < numPlasticStrainLikeISVs; i++ ){

                ( *residual.value )[ row0 + i ] = ( *updatedPlasticStrainLikeISVs )[ i ] - ( *plasticStrainLikeISVs )[ i ];

            }

            // Set the terms associated with the slack variables
            row0 = numPlasticMultipliers + numPlasticStrainLikeISVs;
            ( *residual.value )[ row0 + 0 ] = -( *macroYield ) - slackVariables[ 0 ];

            ( *residual.value )[ row0 + 1 ] = -( *microYield ) - slackVariables[ 1 ];

            for ( auto y = microGradientYield->begin( ); y != microGradientYield->end( ); y++ ){

                unsigned int index = ( unsigned int )( y - microGradientYield->begin( ) );

                ( *residual.value )[ row0 + index + 2 ]
                    = -( *y ) - slackVariables[ index + 2 ];

            }

        }

        void residual::setStateVariableJacobians( ){
            /*!
             * Set the state variable residual jacobians
             *
             * We define these residuals as
             *
             * \f$R = \left\langle f \right\rangle - \dot{\gamma} \left\langle -f \right\rangle\f$
             *
             * and
             *
             * \f$R = Z^{t+1} - Z^{t+1,\text{trial}}\f$
             *
             * Because we have five plastic multipliers (\f$\gamma\f$), five strain-like state variables (\f$Z\f$)
             * and five yield surfaces (\f$f\f$) we can define ten different equations.
             *
             * The residual which includes the yield surfaces is somewhat complex as it contains two separate functions.
             * We may include the ability to weaken the Macaulay bracket to hopefully improve convergence.
             */

            const unsigned int dim = hydra->getDimension( );

            const unsigned int numSecondOrderTensor = hydra->getSOTDimension( );

            const unsigned int numThirdOrderTensor  = hydra->getTOTDimension( );

            unsigned int numConfigurations = *hydra->getNumConfigurations( );

            const unsigned int numConfigurationUnknowns = *hydra->getConfigurationUnknownCount( );

            const floatVector *plasticStrainLikeISVs = get_plasticStrainLikeISVs( );

            const unsigned int numPlasticMultipliers = *getNumPlasticMultipliers( );

            unsigned int numPlasticStrainLikeISVs = plasticStrainLikeISVs->size( );

            const floatVector *plasticMultipliers = get_plasticMultipliers( );

            const floatVector slackVariables( get_plasticStateVariables( )->begin( ) + numPlasticMultipliers + numPlasticStrainLikeISVs,
                                              get_plasticStateVariables( )->end( ) );

            const unsigned int numUnknowns = hydra->getNumUnknowns( );

            const unsigned int numISVs = get_plasticStateVariables( )->size( );

            const secondOrderTensor *dMacroYielddStress                 = get_dMacroYielddStress( );

            const floatVector *dMacroYielddFn                     = get_dMacroYielddFn( );

            const floatVector *dMacroYielddStateVariables         = get_dMacroYielddStateVariables( );

            const secondOrderTensor *dMicroYielddStress                 = get_dMicroYielddStress( );

            const floatVector *dMicroYielddFn                     = get_dMicroYielddFn( );

            const floatVector *dMicroYielddStateVariables         = get_dMicroYielddStateVariables( );

            const floatVector *dMicroGradientYielddStress         = get_dMicroGradientYielddStress( );

            const floatVector *dMicroGradientYielddFn             = get_dMicroGradientYielddFn( );

            const floatVector *dMicroGradientYielddChin           = get_dMicroGradientYielddChin( );

            const floatVector *dMicroGradientYielddStateVariables = get_dMicroGradientYielddStateVariables( );

            const floatVector *dUpdatedPlasticStrainLikeISVsdStateVariables = get_dUpdatedPlasticStrainLikeISVsdStateVariables( );

            auto jacobian = get_setDataStorage_stateVariableJacobians( );
            jacobian.zero( numISVs * numUnknowns );

            unsigned int offset = numConfigurations * numConfigurationUnknowns;
            unsigned int row0   = 0;
            for ( unsigned int i = 0; i < numPlasticMultipliers; i++ ){

                ( *jacobian.value )[ numUnknowns * ( i + row0 ) + i + offset ] = slackVariables[ i ];
                ( *jacobian.value )[ numUnknowns * ( i + row0 ) + i + numPlasticMultipliers + numPlasticStrainLikeISVs + offset ] = ( *plasticMultipliers )[ i ];

            }

            offset = numSecondOrderTensor;
            row0 = numPlasticMultipliers;
            offset = numConfigurations * ( 2 * numSecondOrderTensor + numThirdOrderTensor );
            for ( unsigned int i = 0; i < numPlasticStrainLikeISVs; i++ ){

                ( *jacobian.value )[ numUnknowns * ( i + row0 ) + i + offset + numPlasticMultipliers ] -= 1;

                for ( auto j = getStateVariableIndices( )->begin( ); j != getStateVariableIndices( )->end( ); j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + row0 ) + ( *j ) + offset ] += ( *dUpdatedPlasticStrainLikeISVsdStateVariables )[ numISVs * i + ( unsigned int )( j - getStateVariableIndices( )->begin( ) ) ];

                }

            }

            row0 = numPlasticMultipliers + numPlasticStrainLikeISVs;
            offset = numSecondOrderTensor;
            for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                ( *jacobian.value )[ numUnknowns * ( 0 + row0 ) + j ] = -( *dMacroYielddStress )[ j ];

                ( *jacobian.value )[ numUnknowns * ( 1 + row0 ) + j + offset ] = -( *dMicroYielddStress )[ j ];

            }

            offset = 2 * numSecondOrderTensor;
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < numThirdOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 + row0 ) + j + offset ] = -( *dMicroGradientYielddStress )[ numThirdOrderTensor * i + j ];

                }

            }

            // Sub-Deformation gradient jacobians
            offset = 2 * numSecondOrderTensor + numThirdOrderTensor;
            for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                ( *jacobian.value )[ numUnknowns * ( 0 + row0 ) + j + offset ] = -( *dMacroYielddFn )[ j ];

                ( *jacobian.value )[ numUnknowns * ( 1 + row0 ) + j + offset ] = -( *dMicroYielddFn )[ j ];

            }

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 + row0 ) + j + offset ] = -( *dMicroGradientYielddFn )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

            }

            // Sub-Micro deformation jacobians
            offset = 2 * numSecondOrderTensor + numThirdOrderTensor + ( numConfigurations - 1 ) * numSecondOrderTensor;
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 + row0 ) + j + offset ] = -( *dMicroGradientYielddChin )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

            }

            // State Variable Jacobians
            offset = numConfigurations * ( 2 * numSecondOrderTensor + numThirdOrderTensor );
            for ( unsigned int j = 0; j < ( numPlasticMultipliers + numPlasticStrainLikeISVs ); j++ ){

                ( *jacobian.value )[ numUnknowns * ( 0 + row0 ) + j + offset ] += -( *dMacroYielddStateVariables )[ j ];

                ( *jacobian.value )[ numUnknowns * ( 1 + row0 ) + j + offset ] += -( *dMicroYielddStateVariables )[ j ];

            }

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < ( numPlasticMultipliers + numPlasticStrainLikeISVs ); j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 + row0 ) + j + offset ] += -( *dMicroGradientYielddStateVariables )[ numISVs * i + j ];

                }

            }

            offset = numConfigurations * ( 2 * numSecondOrderTensor + numThirdOrderTensor ) + numPlasticMultipliers + numPlasticStrainLikeISVs;
            for ( unsigned int i = 0; i < 5; i++ ){

                ( *jacobian.value )[ numUnknowns * ( i + row0 ) + i + offset ] -= 1;

            }

        }

        void residual::setdStateVariableResidualsdD( ){
            /*!
             * Set the state variable residuals derivatives w.r.t. the deformation measures
             *
             * We define these residuals as
             *
             * \f$R = \dot{\gamma} f \f$
             *
             * and
             *
             * \f$R = Z^{t+1} - Z^{t+1,\text{trial}}\f$
             *
             * Because we have five plastic multipliers (\f$\gamma\f$), five strain-like state variables (\f$Z\f$)
             * and five yield surfaces (\f$f\f$) we can define ten different equations.
             *
             * The residual which includes the yield surfaces is somewhat complex as it contains two separate functions.
             * We may include the ability to weaken the Macaulay bracket to hopefully improve convergence.
             */

            const unsigned int numConfigurationUnknowns = *hydra->getConfigurationUnknownCount( );

            const unsigned int dim = hydra->getDimension( );

            const unsigned int numSecondOrderTensor = hydra->getSOTDimension( );

            const unsigned int numPlasticMultipliers = *getNumPlasticMultipliers( );

            unsigned int numPlasticStrainLikeISVs = get_plasticStrainLikeISVs( )->size( );

            const secondOrderTensor *dMacroYielddF                      = get_dMacroYielddF( );

            const secondOrderTensor *dMicroYielddF                      = get_dMicroYielddF( );

            const thirdOrderTensor  *dMicroGradientYielddF              = get_dMicroGradientYielddF( );

            const thirdOrderTensor  *dMicroGradientYielddChi            = get_dMicroGradientYielddChi( );

            auto dRdD = get_setDataStorage_dStateVariableResidualsdD( );
            dRdD.zero( get_plasticStateVariables( )->size( ) * numConfigurationUnknowns );

            unsigned int offset = 0;
            unsigned int row0   = numPlasticMultipliers + numPlasticStrainLikeISVs;

            for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                ( *dRdD.value )[ numConfigurationUnknowns * ( 0 + row0 ) + j + offset ] = -( *dMacroYielddF )[ j ];

                ( *dRdD.value )[ numConfigurationUnknowns * ( 1 + row0 ) + j + offset ] = -( *dMicroYielddF )[ j ];

            }

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                    ( *dRdD.value )[ numConfigurationUnknowns * ( i + 2 + row0 ) + j + offset ] = -( *dMicroGradientYielddF )[ numSecondOrderTensor * i + j ];

                }

            }

            // Micro deformation jacobians
            offset = numSecondOrderTensor;
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                    ( *dRdD.value )[ numConfigurationUnknowns * ( i + 2 + row0 ) + j + offset ] = -( *dMicroGradientYielddChi )[ numSecondOrderTensor * i + j ];

                }

            }

        }

        void residual::setConstraints( ){
            /*!
             * Set the values of the constraints which are defined as
             * 
             * c = -f - s
             * 
             * Where \f$f\f$ is the yield surface and \f$s\f$ is the slack variable.
             */

            const unsigned int numPlasticMultipliers = *getNumPlasticMultipliers( );

            unsigned int numPlasticStrainLikeISVs = get_plasticStrainLikeISVs( )->size( );

            const floatVector *plasticMultipliers = get_plasticMultipliers( );

            const floatVector slackVariables( get_plasticStateVariables( )->begin( ) + numPlasticMultipliers + numPlasticStrainLikeISVs,
                                              get_plasticStateVariables( )->end( ) );

            const unsigned int numSlackVariables = slackVariables.size( );

            auto constraints = get_setDataStorage_constraints( );
            constraints.zero( numPlasticMultipliers + numSlackVariables );

            // Set the positivity constraints on the plastic multipliers
            for ( unsigned int i = 0; i < numPlasticMultipliers; i++ ){
                ( *constraints.value )[ i ] = ( *plasticMultipliers )[ i ];
            }

            // Set the positivity constraints on the slack variables
            for ( unsigned int i = 0; i < numSlackVariables; i++ ){
                ( *constraints.value )[ i + numPlasticMultipliers ] = slackVariables[ i ];
            }

        }

        void residual::setConstraintJacobians( ){
            /*!
             * Set the values of the constraint Jacobians which are defined as
             * 
             * c = -f - s
             * 
             * Where \f$f\f$ is the yield surface and \f$s\f$ is the slack variable.
             */

            const unsigned int numUnknowns = hydra->getNumUnknowns( );

            const unsigned int numConfigurations = *hydra->getNumConfigurations( );

            const unsigned int numConfigurationUnknowns = *hydra->getConfigurationUnknownCount( );

            const unsigned int numPlasticMultipliers = *getNumPlasticMultipliers( );

            unsigned int numPlasticStrainLikeISVs = get_plasticStrainLikeISVs( )->size( );

            const floatVector slackVariables( get_plasticStateVariables( )->begin( ) + numPlasticMultipliers + numPlasticStrainLikeISVs,
                                              get_plasticStateVariables( )->end( ) );

            const unsigned int numSlackVariables = slackVariables.size( );

            auto jacobian = get_setDataStorage_constraintJacobians( );
            jacobian.zero( ( numPlasticMultipliers + numSlackVariables ) * numUnknowns );

            unsigned int offset = numConfigurations * numConfigurationUnknowns;
            for ( unsigned int i = 0; i < numPlasticMultipliers; i++ ){
                ( *jacobian.value )[ numUnknowns * i + i + offset ] = 1;
            }

            unsigned int row0 = numPlasticMultipliers;
            offset = numConfigurations * numConfigurationUnknowns + numPlasticMultipliers + numPlasticStrainLikeISVs;
            for ( unsigned int i = 0; i < numSlackVariables; i++ ){
                ( *jacobian.value )[ numUnknowns * ( i + row0 ) + i + offset ] = 1;
            }

        }

        void residual::suggestInitialIterateValues( std::vector< unsigned int >   &indices,
                                                    std::vector< floatType > &values ){

            /*!
             * Function which is called which allows the residual to suggest initial values for given
             * configurations. This is called when the unknown vector is being initialized. If more than
             * one residual attempts to set the initial vector the last residual will override all of the others.
             *
             * After the initial iterate has been suggested, the iteration data is cleared so that the residual
             * starts the iteration in a clean state.
             * 
             * \param &indices: The indices of the unknown vector to set
             * \param &values:  The values to be set in the unknown vector
             */

            constexpr unsigned int dim = 3;

            const unsigned int numPlasticMultipliers   = *getNumPlasticMultipliers( );

            const unsigned int numPlasticStrainLikeISVs = get_plasticStrainLikeISVs( )->size( );

            const unsigned int offset = numPlasticMultipliers + numPlasticStrainLikeISVs;

            indices = std::vector< unsigned int >( getStateVariableIndices( )->begin( ) + offset,
                                                   getStateVariableIndices( )->end( ) );

            values  = std::vector< floatType >( 5, 0 );

            values[ 0 ] = std::fmax( 0, -( *get_macroYield( ) ) );

            if ( ( *get_macroYield( ) ) > 0 ){

                // Perturb the plastic multiplier
                floatType delta = std::fmax( ( *get_plasticMultipliers( ) )[ 0 ], ( *hydra->getRelativeTolerance( ) ) * ( *get_macroYield( ) ) + ( *hydra->getAbsoluteTolerance( ) ) );
                indices.push_back( ( *getStateVariableIndices( ) )[ 0 ] );
                values.push_back( delta );

            }

            values[ 1 ] = std::fmax( 0, -( *get_microYield( ) ) );

            if ( ( *get_microYield( ) ) > 0 ){

                // Perturb the plastic multiplier
                floatType delta = std::fmax( ( *get_plasticMultipliers( ) )[ 1 ], ( *hydra->getRelativeTolerance( ) ) * ( *get_microYield( ) ) + ( *hydra->getAbsoluteTolerance( ) ) );
                indices.push_back( ( *getStateVariableIndices( ) )[ 1 ] );
                values.push_back( delta );

            }

            for ( unsigned int i = 0; i < dim; i++ ){

                values[ i + 2 ] = std::fmax( 0, -( *get_microGradientYield( ) )[ i ] );

                if ( ( *get_microGradientYield( ) )[ i ] > 0 ){
    
                    // Perturb the plastic multiplier
                    floatType delta = std::fmax( ( *get_plasticMultipliers( ) )[ i + 2 ], ( *hydra->getRelativeTolerance( ) ) * ( *get_microGradientYield( ) )[ i ] + ( *hydra->getAbsoluteTolerance( ) ) );
                    indices.push_back( ( *getStateVariableIndices( ) )[ i + 2 ] );
                    values.push_back( delta );
    
                }

            }

            indices += ( *hydra->getNumConfigurations( ) ) * ( *hydra->getConfigurationUnknownCount( ) );

        }

    }

}
