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

    }

}
