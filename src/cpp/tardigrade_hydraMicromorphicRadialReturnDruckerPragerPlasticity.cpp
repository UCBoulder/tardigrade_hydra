/**
  ******************************************************************************
  * \file tardigrade_hydraMicromorphicDruckerPragerPlasticityRadialReturn.cpp
  ******************************************************************************
  * An implementation of micromorphic Drucker-Prager plasticity based on an
  * radial return mapping approach.
  ******************************************************************************
  */

#include<tardigrade_hydraMicromorphicDruckerPragerPlasticityRadialReturn.h>

namespace tardigradeHydra{

    namespace micromorphicDruckerPragerPlasticityRadialReturn{

        void residual::setStateVariableResiduals( ){
            /*!
             * Set the state variable residuals
             *
             * We define these residuals as
             *
             * \f$R = -( c^{\text{trial}} - c ) + \dot{\bar{\gamma}} \frac{\partial f}{\partial Z} \f$
             *
             * and
             *
             * \f$R = f\f$ if the constraint is in the active set and \f$ R = \dot{\bar{\gamma}} \f$ if it isn't
             */

            auto trialMacroCohesion = get_trialMacroCohesion( );

            auto trialMicroCohesion = get_trialMicroCohesion( );

            auto trialMicroGradientCohesion = get_trialMicroGradientCohesion( );

            auto macroCohesion = get_macroCohesion( );

            auto microCohesion = get_microCohesion( );

            auto microGradientCohesion = get_microGradientCohesion( );

            auto numPlasticMultipliers = *getNumPlasticMultipliers( );

            auto plasticMultipliers = get_plasticMultipliers( );

            auto macroYield = get_macroYield( );

            auto microYield = get_microYield( );

            auto microGradientYield = get_microGradientYield( );

            auto dMacroYielddStateVariables = get_dMacroYielddStateVariables( );

            auto dMicroYielddStateVariables = get_dMicroYielddStateVariables( );

            auto dMicroGradientYielddStateVariables = get_dMicroGradientYielddStateVariables( );

            auto activeConstraintSet = get_activeConstraintSet( );

            auto residual = get_setDataStorage_stateVariableResiduals( );
            residual.zero( get_plasticStateVariables( )->size( ) );

            // Add the trial to computed cohesion comparison
            ( *residual.value )[ 0 ] = -( ( *trialMacroCohesion ) - ( *macroCohesion ) );

            ( *residual.value )[ 1 ] = -( ( *trialMicroCohesion ) - ( *microCohesion ) );

            for (
                unsigned int i = 0;
                i < ( const unsigned int )( std::end( *microGradientCohesion ) - std::begin( *microGradientCohesion ) );
                ++i
            ){

                ( *residual.value )[ 2 + i ] = -( ( std::begin( *trialMicroGradientCohesion ) + i ) - ( std::begin( microGradientCohesion ) + i ) );

            }

            // Add the constraint equations
            for ( auto v = std::begin( *activeConstraintSet ); v != std::end( *activeConstraintSet ); ++v ){

                ( *residual.value )[ 0 ] += ( *( std::begin( *plasticMultipliers ) + 0 ) ) * ( *( std::begin( *dMacroYielddStateVariables ) + ( *v ) ) );

                ( *residual.value )[ 1 ] += ( *( std::begin( *plasticMultipliers ) + 1 ) ) * ( *( std::begin( *dMicroYielddStateVariables ) + ( *v ) ) );

                for (
                    unsigned int i = 0;
                    i < ( const unsigned int )( std::end( *microGradientCohesion ) - std::begin( *microGradientCohesion ) );
                    ++i
                ){

                    ( *residual.value )[ 2 + i ] += ( *( std::begin( *plasticMultipliers ) + 2 + i ) ) * ( *( std::begin( *dMicroGradientYielddStateVariables ) + numPlasticStrainLikeISVS * i + ( *v ) ) );

                }

            }

            // Add the yield conditions

        }

    }

}
