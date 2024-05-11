/**
  ******************************************************************************
  * \file tardigrade_hydraPeryznaJ2Viscoplasticity.cpp
  ******************************************************************************
  * An implementation of peryznaJ2Viscoplasticity using the hydra framework.
  * Used as an example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraPeryznaJ2Viscoplasticity.h>

#include<tardigrade_stress_tools.h>

namespace tardigradeHydra{

    namespace peryznaJ2Viscoplasticity{

        void residual::setYieldFunction( const bool isPrevious ){
            /*!
             * Set the value of the yield function
             * 
             * \param isPrevious: Flag for whether this is the previous timestep
             */

            const floatVector* drivingStress;

            const floatVector* stateVariables;

            const floatVector* yieldParameters;

            floatType yieldFunction;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_previousDrivingStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_previousStateVariables( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_drivingStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_stateVariables( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( yieldParameters = get_yieldParameters( ) );

            floatType vonMises;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::calculateVonMisesStress( *drivingStress, vonMises ) );

            yieldFunction = std::fabs( vonMises - ( *yieldParameters )[ 2 ] * ( *stateVariables )[ 1 ] ) - ( *yieldParameters )[ 1 ] * ( *stateVariables )[ 0 ] - ( *yieldParameters )[ 0 ];

            if ( isPrevious ){

                set_previousYieldFunction( yieldFunction );

            }
            else{

                set_yieldFunction( yieldFunction );

            }

        }

        void residual::setYieldFunctionDerivatives( const bool isPrevious ){
            /*!
             * Set the value of the yield function derivatives
             * 
             * \param isPrevious: Flag for whether this is the previous timestep
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector* drivingStress;

            const floatVector* dDrivingStressdCauchyStress;

            const floatVector* dDrivingStressdF;

            const floatVector* dDrivingStressdSubFs;

            const floatVector* stateVariables;

            const floatVector* yieldParameters;

            floatType yieldFunction;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = get_dPreviousDrivingStressdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = get_dPreviousDrivingStressdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = get_dPreviousDrivingStressdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_previousDrivingStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_previousStateVariables( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = get_dDrivingStressdCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = get_dDrivingStressdF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = get_dDrivingStressdSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_drivingStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_stateVariables( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( yieldParameters = get_yieldParameters( ) );

            floatType vonMises;

            floatVector dVonMisesdDrivingStress( sot_dim, 0 );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::calculateVonMisesStress( *drivingStress, vonMises, dVonMisesdDrivingStress ) );

            floatType term1 = vonMises - ( *yieldParameters )[ 2 ] * ( *stateVariables )[ 1 ];

            yieldFunction = std::fabs( term1 ) - ( *yieldParameters )[ 1 ] * ( *stateVariables )[ 0 ] - ( *yieldParameters )[ 0 ];

            floatVector dYieldFunctiondDrivingStress = sgn( term1 ) * dVonMisesdDrivingStress;

            floatVector dYieldFunctiondStateVariables = { -( *yieldParameters )[ 1 ], -sgn( term1 ) * ( *yieldParameters )[ 2 ] };

            floatVector dYieldFunctiondCauchyStress = tardigradeVectorTools::matrixMultiply( dYieldFunctiondDrivingStress, *dDrivingStressdCauchyStress, 1, sot_dim, sot_dim, sot_dim );

            floatVector dYieldFunctiondF = tardigradeVectorTools::matrixMultiply( dYieldFunctiondDrivingStress, *dDrivingStressdF, 1, sot_dim, sot_dim, sot_dim );

            floatVector dYieldFunctiondSubFs = tardigradeVectorTools::matrixMultiply( dYieldFunctiondDrivingStress, *dDrivingStressdSubFs, 1, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            if ( isPrevious ){

                set_previousYieldFunction( yieldFunction );

                set_dPreviousYieldFunctiondPreviousCauchyStress( dYieldFunctiondCauchyStress );

                set_dPreviousYieldFunctiondPreviousF( dYieldFunctiondF );

                set_dPreviousYieldFunctiondPreviousSubFs( dYieldFunctiondSubFs );

                set_dPreviousYieldFunctiondPreviousStateVariables( dYieldFunctiondStateVariables );

            }
            else{

                set_yieldFunction( yieldFunction );

                set_dYieldFunctiondCauchyStress( dYieldFunctiondCauchyStress );

                set_dYieldFunctiondF( dYieldFunctiondF );

                set_dYieldFunctiondSubFs( dYieldFunctiondSubFs );

                set_dYieldFunctiondStateVariables( dYieldFunctiondStateVariables );

            }

        }

        void residual::decomposeParameters( const floatVector &parameters ){
            /*!
             * Decompose the incoming parameter vector
             * 
             * \param &parameters: The incoming parameter vector. We assume a
             *     Peryzna viscoplastic driving stress with a Drucker-Prager
             *     yield and flow potential surface. The parameters are
             *     [ n, q0, C1, C2, Tref, Y, Ei, Ek, hi0, hi1, hk0, hk1 ] where n is the
             *     Peryzna exponent, q0 is the initial drag stress,
             *     C1, C2, and Tref are the temperature effect
             *     parameters, Y is the yield stress for the yield equation,
             *     Ei is the isotroping hardening modulus, Ek is the kinematic hardening
             *     modulus, hi0 and hi1 are the initial and linear hardening moduli for
             *     the isotropic hardening, and hk0 and hk1 are the initial and linear hardening
             *     modulus for the kinematic hardening.
             */

            constexpr unsigned int expectedSize = 12;

            TARDIGRADE_ERROR_TOOLS_CHECK( parameters.size( ) == expectedSize, "The parameters vector is not the correct length.\n  parameters: " + std::to_string( parameters.size( ) ) + "\n  required:   " + std::to_string( expectedSize ) + "\n" );

            set_peryznaParameters( { parameters[ 0 ] } );

            set_dragStressParameters( { parameters[ 1 ] } );

            set_thermalParameters( { parameters[ 2 ], parameters[ 3 ], parameters[ 4 ] } );

            set_yieldParameters( { parameters[ 5 ], parameters[ 6 ], parameters[ 7 ] } );

            set_hardeningParameters( { parameters[ 8 ], parameters[ 9 ], parameters[ 10 ], parameters[ 11 ] } );

        }

    }

}
