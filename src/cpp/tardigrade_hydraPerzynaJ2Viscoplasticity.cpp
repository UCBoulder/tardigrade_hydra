/**
  ******************************************************************************
  * \file tardigrade_hydraPerzynaJ2Viscoplasticity.cpp
  ******************************************************************************
  * An implementation of perzynaJ2Viscoplasticity using the hydra framework.
  * Used as an example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraPerzynaJ2Viscoplasticity.h>

#include<tardigrade_stress_tools.h>

namespace tardigradeHydra{

    namespace perzynaJ2Viscoplasticity{

        void residual::setSignTerm( ){
            /*!
             * Set the sign term for kinematic hardening
             */

            setYieldFunction( false );

        }

        void residual::setPreviousSignTerm( ){
            /*!
             * Set the previous value of the sign term for kinematic hardening
             */

            setYieldFunction( true );

        }

        void residual::setYieldFunction( const bool isPrevious ){
            /*!
             * Set the value of the yield function
             * 
             * \param isPrevious: Flag for whether this is the previous timestep
             */

            const floatVector* drivingStress;

            const floatVector* stateVariables;

            const floatVector* yieldParameters;

            setDataStorageBase< floatType > yieldFunction;

            setDataStorageBase< floatType > signTerm;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_previousDrivingStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_previousStateVariables( ) );

                yieldFunction = get_setDataStorage_previousYieldFunction( );

                signTerm      = get_setDataStorage_previousSignTerm( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_drivingStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_stateVariables( ) );

                yieldFunction = get_setDataStorage_yieldFunction( );

                signTerm      = get_setDataStorage_signTerm( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( yieldParameters = get_yieldParameters( ) );

            floatType vonMises;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::calculateVonMisesStress( *drivingStress, vonMises ) );

            floatType term1 = vonMises - ( *yieldParameters )[ 2 ] * ( *stateVariables )[ 1 ];

            *yieldFunction.value = std::fabs( term1 ) - ( *yieldParameters )[ 1 ] * ( *stateVariables )[ 0 ] - ( *yieldParameters )[ 0 ];

            *signTerm.value      = sgn( term1 );

        }

        void residual::setYieldFunctionDerivatives( const bool isPrevious ){
            /*!
             * Set the value of the yield function derivatives
             * 
             * \param isPrevious: Flag for whether this is the previous timestep
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector* drivingStress;

            const floatVector* dDrivingStressdCauchyStress;

            const floatVector* dDrivingStressdF;

            const floatVector* dDrivingStressdSubFs;

            const floatVector* stateVariables;

            const floatVector* yieldParameters;

            setDataStorageBase< floatType > yieldFunction;

            setDataStorageBase< floatType > signTerm;

            setDataStorageBase< secondOrderTensor > dYieldFunctiondCauchyStress;

            setDataStorageBase< secondOrderTensor > dYieldFunctiondF;

            setDataStorageBase< floatVector > dYieldFunctiondSubFs;

            setDataStorageBase< floatVector > dYieldFunctiondStateVariables;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = get_dPreviousDrivingStressdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = get_dPreviousDrivingStressdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = get_dPreviousDrivingStressdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_previousDrivingStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_previousStateVariables( ) );

                yieldFunction = get_setDataStorage_previousYieldFunction( );

                signTerm      = get_setDataStorage_previousSignTerm( );

                dYieldFunctiondCauchyStress   = get_setDataStorage_dPreviousYieldFunctiondPreviousCauchyStress( );

                dYieldFunctiondF              = get_setDataStorage_dPreviousYieldFunctiondPreviousF( );

                dYieldFunctiondSubFs          = get_setDataStorage_dPreviousYieldFunctiondPreviousSubFs( );

                dYieldFunctiondStateVariables = get_setDataStorage_dPreviousYieldFunctiondPreviousStateVariables( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = get_dDrivingStressdCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = get_dDrivingStressdF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = get_dDrivingStressdSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_drivingStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_stateVariables( ) );

                yieldFunction = get_setDataStorage_yieldFunction( );

                signTerm      = get_setDataStorage_signTerm( );

                dYieldFunctiondCauchyStress   = get_setDataStorage_dYieldFunctiondCauchyStress( );

                dYieldFunctiondF              = get_setDataStorage_dYieldFunctiondF( );

                dYieldFunctiondSubFs          = get_setDataStorage_dYieldFunctiondSubFs( );

                dYieldFunctiondStateVariables = get_setDataStorage_dYieldFunctiondStateVariables( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( yieldParameters = get_yieldParameters( ) );

            floatType vonMises;

            floatVector dVonMisesdDrivingStress( sot_dim, 0 );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::calculateVonMisesStress( *drivingStress, vonMises, dVonMisesdDrivingStress ) );

            floatType term1 = vonMises - ( *yieldParameters )[ 2 ] * ( *stateVariables )[ 1 ];

            *yieldFunction.value = std::fabs( term1 ) - ( *yieldParameters )[ 1 ] * ( *stateVariables )[ 0 ] - ( *yieldParameters )[ 0 ];

            *signTerm.value = sgn( term1 );

            floatVector dYieldFunctiondDrivingStress = ( *signTerm.value ) * dVonMisesdDrivingStress;

            *dYieldFunctiondStateVariables.value = { -( *yieldParameters )[ 1 ], -( *signTerm.value ) * ( *yieldParameters )[ 2 ] };

            auto map_dYieldFunctiondDrivingStress = getFixedSizeMatrixMap< floatType,       1, sot_dim >( dYieldFunctiondDrivingStress.data( ) );
            auto map_dDrivingStressdCauchyStress  = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dDrivingStressdCauchyStress->data( ) );
            auto map_dDrivingStressdF             = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dDrivingStressdF->data( ) );

            auto map_dDrivingStressdSubFs         = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dDrivingStressdSubFs->data( ), ( num_configs - 1 ) * sot_dim );

            auto map_dYieldFunctiondCauchyStress = dYieldFunctiondCauchyStress.zeroMap< floatType, 1, sot_dim >( );

            auto map_dYieldFunctiondF            = dYieldFunctiondF.zeroMap< floatType, 1, sot_dim >( );

            auto map_dYieldFunctiondSubFs        = dYieldFunctiondSubFs.zeroMap< floatType, 1 >( ( num_configs - 1 ) * sot_dim );

            map_dYieldFunctiondCauchyStress = ( map_dYieldFunctiondDrivingStress * map_dDrivingStressdCauchyStress ).eval( );

            map_dYieldFunctiondF            = ( map_dYieldFunctiondDrivingStress * map_dDrivingStressdF ).eval( );

            map_dYieldFunctiondSubFs        = ( map_dYieldFunctiondDrivingStress * map_dDrivingStressdSubFs ).eval( );

        }

        void residual::setDragStress( const bool isPrevious ){
            /*!
             * Set the viscous drag stress
             */

            setDataStorageBase< floatType > dragStress;

            if ( isPrevious ){
                dragStress = get_setDataStorage_previousDragStress( );
            }
            else{
                dragStress = get_setDataStorage_dragStress( );
            }

            *dragStress.value = ( *get_dragStressParameters( ) )[ 0 ];

        }

        void residual::setDragStressDerivatives( const bool isPrevious ){
            /*!
             * Set the derivatives of the drag stress
             */

            setDataStorageBase< floatType > dragStress;
            setDataStorageBase< floatVector > dDragStressdStateVariables;

            if ( isPrevious ){
                dragStress = get_setDataStorage_previousDragStress( );
                dDragStressdStateVariables = get_setDataStorage_dPreviousDragStressdPreviousStateVariables( );
            }
            else{
                dragStress = get_setDataStorage_dragStress( );
                dDragStressdStateVariables = get_setDataStorage_dDragStressdStateVariables( );
            }

            *dragStress.value = ( *get_dragStressParameters( ) )[ 0 ];
            *dDragStressdStateVariables.value = floatVector( get_stateVariables( )->size( ), 0 );

        }
        void residual::setHardeningFunction( const bool isPrevious ){
            /*!
             * Set the value of the hardening function
             * 
             * \param &isPrevious: Flag for whether to compute the values for the
             *     previous timestep
             */

            const floatVector *stateVariables;

            const floatVector *hardeningParameters;

            const floatType *sign_term;

            setDataStorageBase< floatVector > hardeningFunction;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_previousStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( sign_term = get_previousSignTerm( ) );

                hardeningFunction = get_setDataStorage_previousHardeningFunction( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_stateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( sign_term = get_signTerm( ) );

                hardeningFunction = get_setDataStorage_hardeningFunction( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( hardeningParameters = get_hardeningParameters( ) );

            *hardeningFunction.value = { ( ( *hardeningParameters )[ 0 ] + ( *hardeningParameters )[ 1 ] * ( *stateVariables )[ 0 ] ),
                                         ( ( *hardeningParameters )[ 2 ] + ( *hardeningParameters )[ 3 ] * ( *stateVariables )[ 0 ] ) * ( *sign_term ) };

        }

        void residual::setHardeningFunctionDerivatives( const bool isPrevious ){
            /*!
             * Set the value of the derivatives of the hardening function
             * 
             * \param &isPrevious: Flag for whether to compute the values for the
             *     previous timestep
             */

            const floatVector *stateVariables;

            const floatVector *hardeningParameters;

            const floatType *sign_term;

            setDataStorageBase< floatVector > hardeningFunction;

            setDataStorageBase< floatVector > dHardeningFunctiondStateVariables;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_previousStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( sign_term = get_previousSignTerm( ) );

                hardeningFunction = get_setDataStorage_previousHardeningFunction( );

                dHardeningFunctiondStateVariables = get_setDataStorage_dPreviousHardeningFunctiondPreviousStateVariables( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_stateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( sign_term = get_signTerm( ) );

                hardeningFunction = get_setDataStorage_hardeningFunction( );

                dHardeningFunctiondStateVariables = get_setDataStorage_dHardeningFunctiondStateVariables( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( hardeningParameters = get_hardeningParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( hardeningParameters = get_hardeningParameters( ) );

            *hardeningFunction.value = { ( ( *hardeningParameters )[ 0 ] + ( *hardeningParameters )[ 1 ] * ( *stateVariables )[ 0 ] ),
                                         ( ( *hardeningParameters )[ 2 ] + ( *hardeningParameters )[ 3 ] * ( *stateVariables )[ 0 ] ) * ( *sign_term ) };

            *dHardeningFunctiondStateVariables.value = { ( *hardeningParameters )[ 1 ], 0., ( *hardeningParameters )[ 3 ] * ( *sign_term ), 0 };

        }

        void residual::decomposeParameters( const floatVector &parameters ){
            /*!
             * Decompose the incoming parameter vector
             * 
             * \param &parameters: The incoming parameter vector. We assume a
             *     Perzyna viscoplastic driving stress with a Drucker-Prager
             *     yield and flow potential surface. The parameters are
             *     [ n, q0, C1, C2, Tref, Y, Ei, Ek, hi0, hi1, hk0, hk1 ] where n is the
             *     Perzyna exponent, q0 is the initial drag stress,
             *     C1, C2, and Tref are the temperature effect
             *     parameters, Y is the yield stress for the yield equation,
             *     Ei is the isotroping hardening modulus, Ek is the kinematic hardening
             *     modulus, hi0 and hi1 are the initial and linear hardening moduli for
             *     the isotropic hardening, and hk0 and hk1 are the initial and linear hardening
             *     modulus for the kinematic hardening.
             */

            constexpr unsigned int expectedSize = 12;

            TARDIGRADE_ERROR_TOOLS_CHECK( parameters.size( ) == expectedSize, "The parameters vector is not the correct length.\n  parameters: " + std::to_string( parameters.size( ) ) + "\n  required:   " + std::to_string( expectedSize ) + "\n" );

            set_perzynaParameters( { parameters[ 0 ] } );

            set_dragStressParameters( { parameters[ 1 ] } );

            set_thermalParameters( { parameters[ 2 ], parameters[ 3 ], parameters[ 4 ] } );

            set_yieldParameters( { parameters[ 5 ], parameters[ 6 ], parameters[ 7 ] } );

            set_flowParameters( { 0., 0. } );

            set_hardeningParameters( { parameters[ 8 ], parameters[ 9 ], parameters[ 10 ], parameters[ 11 ] } );

        }

        void residual::addParameterizationInfo( std::string &parameterization_info ){
            /*!
             * Add the information about the parameterization to the incoming string
             * 
             * \param &parameterization_info: The incoming parameterization string
             */

            parameterization_info += "class: tardigradeHydra::perzynaJ2Viscoplasticity::residual\n\n";
            parameterization_info += "name,                          description,       units, current value\n";
            parameterization_info += "   n,         the Perzyna exponential term,        none, " + std::to_string(    ( *get_perzynaParameters( ) )[ 0 ] ) + "\n";
            parameterization_info += "  q0,                      the drag stress,      stress, " + std::to_string( ( *get_dragStressParameters( ) )[ 0 ] ) + "\n";
            parameterization_info += "  C1,                 the WLF C1 parameter,        none, " + std::to_string(    ( *get_thermalParameters( ) )[ 0 ] ) + "\n";
            parameterization_info += "  C2,                 the WLF C2 parameter, temperature, " + std::to_string(    ( *get_thermalParameters( ) )[ 1 ] ) + "\n";
            parameterization_info += "Tref,        the WLF reference temperature, temperature, " + std::to_string(    ( *get_thermalParameters( ) )[ 2 ] ) + "\n";
            parameterization_info += "   Y,                 initial yield stress,      stress, " + std::to_string(      ( *get_yieldParameters( ) )[ 0 ] ) + "\n";
            parameterization_info += "  Ei,          isotropic hardening modulus,      stress, " + std::to_string(      ( *get_yieldParameters( ) )[ 1 ] ) + "\n";
            parameterization_info += "  Ek,          kinematic hardening modulus,      stress, " + std::to_string(      ( *get_yieldParameters( ) )[ 2 ] ) + "\n";
            parameterization_info += " hi0, isotropic isv initial evolution rate,        none, " + std::to_string(  ( *get_hardeningParameters( ) )[ 0 ] ) + "\n";
            parameterization_info += " hi1,  isotropic isv linear evolution rate,        none, " + std::to_string(  ( *get_hardeningParameters( ) )[ 1 ] ) + "\n";
            parameterization_info += " hk0, kinematic isv initial evolution rate,        none, " + std::to_string(  ( *get_hardeningParameters( ) )[ 2 ] ) + "\n";
            parameterization_info += " hk1,  kinematic isv linear evolution rate,        none, " + std::to_string(  ( *get_hardeningParameters( ) )[ 3 ] ) + "\n";

        }

    }

}
