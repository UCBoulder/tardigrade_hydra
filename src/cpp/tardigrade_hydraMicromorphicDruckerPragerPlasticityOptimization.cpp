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

            const floatVector *plasticMultipliers = get_plasticMultipliers( );

            const floatVector *plasticStrainLikeISVs = get_plasticStrainLikeISVs( );

            const floatVector *updatedPlasticStrainLikeISVs = get_updatedPlasticStrainLikeISVs( );

            const unsigned int numPlasticMultipliers = *getNumPlasticMultipliers( );

            unsigned int numPlasticStrainLikeISVs = plasticStrainLikeISVs->size( );

            const floatType *macroYield = get_macroYield( );

            const floatType *microYield = get_microYield( );

            const dimVector *microGradientYield = get_microGradientYield( );

            auto residual = get_setDataStorage_stateVariableResiduals( );
            residual.zero( get_plasticStateVariables( )->size( ) );

            // Set the terms associated with the yield surface
            ( *residual.value )[ 0 ] = ( *plasticMultipliers )[ 0 ] * ( *macroYield );

            ( *residual.value )[ 1 ] = ( *plasticMultipliers )[ 1 ] * ( *microYield );

            for ( auto y = microGradientYield->begin( ); y != microGradientYield->end( ); y++ ){

                unsigned int index = ( unsigned int )( y - microGradientYield->begin( ) );

                ( *residual.value )[ index + 2 ]
                    = ( *plasticMultipliers )[ index + 2 ] * ( *y );

            }

            // Set the terms associated with the strain-like ISV evolution
            for ( unsigned int i = 0; i < numPlasticStrainLikeISVs; i++ ){

                ( *residual.value )[ numPlasticMultipliers + i ] = ( *updatedPlasticStrainLikeISVs )[ i ] - ( *plasticStrainLikeISVs )[ i ];

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

            const floatVector *plasticMultipliers = get_plasticMultipliers( );

            const floatVector *plasticStrainLikeISVs = get_plasticStrainLikeISVs( );

            const unsigned int numPlasticMultipliers = *getNumPlasticMultipliers( );

            unsigned int numPlasticStrainLikeISVs = plasticStrainLikeISVs->size( );

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

            const floatType *macroYield = get_macroYield( );

            const floatType *microYield = get_microYield( );

            const dimVector *microGradientYield = get_microGradientYield( );

            auto jacobian = get_setDataStorage_stateVariableJacobians( );
            jacobian.zero( numISVs * numUnknowns );

            unsigned int offset = numSecondOrderTensor;
            for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                ( *jacobian.value )[ numUnknowns * 0 + j ] = ( *plasticMultipliers )[ 0 ] * ( *dMacroYielddStress )[ j ];

                ( *jacobian.value )[ numUnknowns * 1 + j + offset ] = ( *plasticMultipliers )[ 1 ] * ( *dMicroYielddStress )[ j ];

            }

            offset = 2 * numSecondOrderTensor;
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < numThirdOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 ) + j + offset ] = ( *plasticMultipliers )[ i + 2 ] * ( *dMicroGradientYielddStress )[ numThirdOrderTensor * i + j ];

                }

            }

            // Sub-Deformation gradient jacobians
            offset = 2 * numSecondOrderTensor + numThirdOrderTensor;
            for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                ( *jacobian.value )[ numUnknowns * 0 + j + offset ] = ( *plasticMultipliers )[ 0 ] * ( *dMacroYielddFn )[ j ];

                ( *jacobian.value )[ numUnknowns * 1 + j + offset ] = ( *plasticMultipliers )[ 1 ] * ( *dMicroYielddFn )[ j ];

            }

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 ) + j + offset ] = ( *plasticMultipliers )[ i + 2 ] * ( *dMicroGradientYielddFn )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

            }

            // Sub-Micro deformation jacobians
            offset = 2 * numSecondOrderTensor + numThirdOrderTensor + ( numConfigurations - 1 ) * numSecondOrderTensor;
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < ( numConfigurations - 1 ) * numSecondOrderTensor; j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 ) + j + offset ] = ( *plasticMultipliers )[ i + 2 ] * ( *dMicroGradientYielddChin )[ ( numConfigurations - 1 ) * numSecondOrderTensor * i + j ];

                }

            }


            // State Variable Jacobians
            offset = numConfigurations * ( 2 * numSecondOrderTensor + numThirdOrderTensor );

            ( *jacobian.value )[ numUnknowns * 0 + offset + 0 ] += *macroYield;

            ( *jacobian.value )[ numUnknowns * 1 + offset + 1 ] += *microYield;

            for ( unsigned int i = 0; i < dim; i++ ){

                ( *jacobian.value )[ numUnknowns * ( i + 2 ) + offset + i + 2 ] += ( *microGradientYield )[ i ];

            }

            for ( unsigned int j = 0; j < ( numPlasticMultipliers + numPlasticStrainLikeISVs ); j++ ){

                ( *jacobian.value )[ numUnknowns * 0 + j + offset ] += ( *plasticMultipliers )[ 0 ] * ( *dMacroYielddStateVariables )[ j ];

                ( *jacobian.value )[ numUnknowns * 1 + j + offset ] += ( *plasticMultipliers )[ 1 ] * ( *dMicroYielddStateVariables )[ j ];

            }

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < ( numPlasticMultipliers + numPlasticStrainLikeISVs ); j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + 2 ) + j + offset ] += ( *plasticMultipliers )[ i + 2 ] * ( *dMicroGradientYielddStateVariables )[ numISVs * i + j ];

                }

            }

            unsigned int row0 = numPlasticMultipliers;
            offset = numConfigurations * ( 2 * numSecondOrderTensor + numThirdOrderTensor );
            for ( unsigned int i = 0; i < numPlasticStrainLikeISVs; i++ ){

                ( *jacobian.value )[ numUnknowns * ( i + row0 ) + i + offset + numPlasticMultipliers ] -= 1;

                for ( auto j = getStateVariableIndices( )->begin( ); j != getStateVariableIndices( )->end( ); j++ ){

                    ( *jacobian.value )[ numUnknowns * ( i + row0 ) + ( *j ) + offset ] += ( *dUpdatedPlasticStrainLikeISVsdStateVariables )[ numISVs * i + ( unsigned int )( j - getStateVariableIndices( )->begin( ) ) ];

                }

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

            const floatVector *plasticMultipliers = get_plasticMultipliers( );

            const unsigned int numPlasticMultipliers = *getNumPlasticMultipliers( );

            const secondOrderTensor *dMacroYielddF                      = get_dMacroYielddF( );

            const secondOrderTensor *dMicroYielddF                      = get_dMicroYielddF( );

            const thirdOrderTensor  *dMicroGradientYielddF              = get_dMicroGradientYielddF( );

            const thirdOrderTensor  *dMicroGradientYielddChi            = get_dMicroGradientYielddChi( );

            const floatType *macroYield = get_macroYield( );

            const floatType *microYield = get_microYield( );

            const dimVector *microGradientYield = get_microGradientYield( );

            auto dRdD = get_setDataStorage_dStateVariableResidualsdD( );
            dRdD.zero( get_plasticStateVariables( )->size( ) * numConfigurationUnknowns );

            unsigned int offset = 0;
            for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                ( *dRdD.value )[ numConfigurationUnknowns * 0 + j + offset ] = ( *plasticMultipliers )[ 0 ] * ( *dMacroYielddF )[ j ];

                ( *dRdD.value )[ numConfigurationUnknowns * 1 + j + offset ] = ( *plasticMultipliers )[ 1 ] * ( *dMicroYielddF )[ j ];

            }

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                    ( *dRdD.value )[ numConfigurationUnknowns * ( i + 2 ) + j + offset ] = ( *plasticMultipliers )[ i + 2 ] * ( *dMicroGradientYielddF )[ numSecondOrderTensor * i + j ];

                }

            }

            // Micro deformation jacobians
            offset = numSecondOrderTensor;
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < numSecondOrderTensor; j++ ){

                    ( *dRdD.value )[ numConfigurationUnknowns * ( i + 2 ) + j + offset ] = ( *plasticMultipliers )[ i + 2 ] * ( *dMicroGradientYielddChi )[ numSecondOrderTensor * i + j ];

                }

            }

        }

    }

}
