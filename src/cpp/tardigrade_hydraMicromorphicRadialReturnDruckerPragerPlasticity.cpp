/**
  ******************************************************************************
  * \file tardigrade_hydraMicromorphicRadialReturnDruckerPragerPlasticity.cpp
  ******************************************************************************
  * An implementation of micromorphic Drucker-Prager plasticity based on an
  * radial return mapping approach.
  ******************************************************************************
  */

#include<tardigrade_hydraMicromorphicRadialReturnDruckerPragerPlasticity.h>

#include<tardigrade_micromorphic_tools.h>

namespace tardigradeHydra{

    namespace micromorphicRadialReturnDruckerPragerPlasticity{

        void residual::projectSuggestedX( std::vector< floatType > &trialX,    
                                          const std::vector< floatType > &Xp ){
            /*!
             * Project the suggested X vector if required
             * 
             * \param &trialX: The proposed unknown vector
             * \param &Xp: The previous unknown vector
             */
        }

        void residual::setActiveConstraints( ){
            /*!
             * Set which constraints are active
             */

            auto activeConstraints = get_setDataStorage_activeConstraints( );

            activeConstraints.value->resize( 5 );

	    // Check the yield stress constraint
            ( *activeConstraints.value )[ 0 ] = ( ( *get_macroYield( ) ) > 0. );

            ( *activeConstraints.value )[ 1 ] = ( ( *get_microYield( ) ) > 0. );

            for ( auto v = std::begin( *get_microGradientYield( ) ); v != std::end( *get_microGradientYield( ) ); ++v ){

                ( *activeConstraints.value )[ ( v - std::begin( *get_microGradientYield( ) ) ) + 2 ] = ( ( *v ) > 0. );

            }

        }

        void residual::updateActiveConstraints( ){
            /*!
             * Update the active constraints based on the plastic multipliers
             */
        
            get_activeConstraints( ); // Initializes the active constraints if they haven't been yet
            auto activeConstraints = get_setDataStorage_activeConstraints( );
        
            TARDIGRADE_ERROR_TOOLS_CHECK( activeConstraints.value->size( ) == get_plasticMultipliers( )->size( ), "The active constraints must be the same size as the plastic multipliers\n  activeConstraints size : " + std::to_string( activeConstraints.value->size( ) ) + "\n  plasticMultipliers size: " + std::to_string( get_plasticMultipliers( )->size( ) ) + "\n" )
        
            for ( auto v = activeConstraints.begin( ); v != activeConstraints.end( ); ++v ){
        
                *v = ( *v ) && ( ( *get_plasticMultipliers( ) )[ v - activeConstraints.begin( ) ] > *getPlasticMultiplierTolerance( ) );
        
            }
        
        }

        void residual::setStateVariableResiduals( ){
            /*!
             * Set the state variable residuals
             *
             * We define these residuals as
             *
             * \f$R = Z^{\text{update}} - Z \f$
             *
             * and
             *
             * \f$R = f\f$ if the constraint is in the active set and \f$ R = \dot{\bar{\gamma}} \f$ if it isn't
             */

            auto updatedPlasticStrainLikeISVs = get_updatedPlasticStrainLikeISVs( );

            auto plasticStrainLikeISVs = get_plasticStrainLikeISVs( );

            auto plasticMultipliers = get_plasticMultipliers( );

            auto residual = get_setDataStorage_stateVariableResiduals( );

            auto num_ISVS = get_plasticStateVariables( )->size( );

            auto num_plastic_state_variables = ( const unsigned int )( std::end( *plasticStrainLikeISVs ) - std::begin( *plasticStrainLikeISVs ) );

            residual.zero( num_ISVS );

            // Set the state variable residuals
            for ( unsigned int i = 0; i < num_plastic_state_variables; ++i ){

                ( *residual.value )[ i ] = ( *updatedPlasticStrainLikeISVs )[ i ] - ( *plasticStrainLikeISVs )[ i ];

            }

            // Assemble the yield surface values
            std::array< floatType, 5 > yieldSurfaceValues = {
                *get_macroYield( ), *get_microYield( ),
                ( *get_microGradientYield( ) )[ 0 ], ( *get_microGradientYield( ) )[ 1 ], ( *get_microGradientYield( ) )[ 2 ]
            };

            // Set the constraints
            for ( auto v = std::begin( *get_activeConstraints( ) ); v != std::end( *get_activeConstraints( ) ); ++v ){

                if ( *v ){
                    ( *residual.value )[ num_plastic_state_variables + ( unsigned int )( v - std::begin( *get_activeConstraints( ) ) ) ] += yieldSurfaceValues[ ( unsigned int )( v - std::begin( *get_activeConstraints( ) ) )  ];
                }
                else{
                    ( *residual.value )[ num_plastic_state_variables + ( unsigned int )( v - std::begin( *get_activeConstraints( ) ) ) ] += ( *plasticMultipliers )[ ( unsigned int )( v - std::begin( *get_activeConstraints( ) ) ) ];
                }

            }

        }

        void residual::setStateVariableJacobians( ){
            /*!
             * Set the Jacobians of the state variable residuals
             *
             * We define these residuals as
             *
             * \f$R = Z^{\text{update}} - Z \f$
             *
             * and
             *
             * \f$R = f\f$ if the constraint is in the active set and \f$ R = \dot{\bar{\gamma}} \f$ if it isn't
             */

            auto dim = hydra->getDimension( );

            auto sot_dim = dim * dim;

            auto tot_dim = dim * dim * dim;

            auto plasticStrainLikeISVs = get_plasticStrainLikeISVs( );

            auto activeConstraints = get_activeConstraints( );

            auto plasticMultipliers = get_plasticMultipliers( );

            auto jacobian = get_setDataStorage_stateVariableJacobians( );

            auto num_ISVS = get_plasticStateVariables( )->size( );

            auto num_unknowns = hydra->getNumUnknowns( );

            auto num_plastic_multipliers     = ( const unsigned int )( std::end( *plasticMultipliers )    - std::begin( *plasticMultipliers ) );

            auto num_plastic_state_variables = ( const unsigned int )( std::end( *plasticStrainLikeISVs ) - std::begin( *plasticStrainLikeISVs ) );

            auto num_configurations = *hydra->getNumConfigurations( );

            jacobian.zero( num_ISVS * num_unknowns );

            // Set the state variable jacobians
            unsigned int offset = num_configurations * ( 2 * sot_dim + tot_dim );

            for ( unsigned int i = 0; i < num_plastic_state_variables; ++i ){

                ( *jacobian.value )[ num_unknowns * i + i + offset + num_plastic_multipliers ] -= 1;

                for ( unsigned int j = 0; j < num_ISVS; ++j ){

                    ( *jacobian.value )[ num_unknowns * i + j + offset ] += ( *get_dUpdatedPlasticStrainLikeISVsdStateVariables( ) )[ num_ISVS * i + j ];

                }

            }

            // Add the active constraint contributions

            // Macro yielding
            unsigned int row = 0;
            if ( ( *activeConstraints )[ 0 ] ){

                row = num_plastic_state_variables;

                // Stress Jacobians
                unsigned int col;
                for ( auto v = std::begin( *get_dMacroYielddStress( ) ); v != std::end( *get_dMacroYielddStress( ) ); ++v ){
                    
                    col = ( unsigned int )( v - std::begin( *get_dMacroYielddStress( ) ) );
                    ( *jacobian.value )[ num_unknowns * row + col ] += *v;

                }

                // Deformation Jacobians
                for ( auto v = std::begin ( *get_dMacroYielddFn( ) ); v != std::end( *get_dMacroYielddFn( ) ); ++v ){

                    col = ( unsigned int )( v - std::begin( *get_dMacroYielddFn( ) ) ) + 2 * sot_dim + tot_dim;
                    ( *jacobian.value )[ num_unknowns * row + col ] += *v;

                }

                // State variable Jacobians
                for ( auto v = std::begin ( *get_dMacroYielddStateVariables( ) ); v != std::end( *get_dMacroYielddStateVariables( ) ); ++v ){

                    col = ( unsigned int )( v - std::begin( *get_dMacroYielddStateVariables( ) ) ) + num_configurations * ( 2 * sot_dim + tot_dim );
                    ( *jacobian.value )[ num_unknowns * row + col ] += *v;

                }

            }

            // Micro yielding
            if ( ( *activeConstraints )[ 1 ] ){

                row = num_plastic_state_variables + 1;

                // Stress Jacobians
                unsigned int col;
                for ( auto v = std::begin( *get_dMicroYielddStress( ) ); v != std::end( *get_dMicroYielddStress( ) ); ++v ){
                    
                    col = sot_dim + ( unsigned int )( v - std::begin( *get_dMicroYielddStress( ) ) );
                    ( *jacobian.value )[ num_unknowns * row + col ] += *v;

                }

                // Deformation Jacobians
                for ( auto v = std::begin ( *get_dMicroYielddFn( ) ); v != std::end( *get_dMicroYielddFn( ) ); ++v ){

                    col = ( unsigned int )( v - std::begin( *get_dMicroYielddFn( ) ) ) + 2 * sot_dim + tot_dim;
                    ( *jacobian.value )[ num_unknowns * row + col ] += *v;

                }

                // State variable Jacobians
                for ( auto v = std::begin ( *get_dMicroYielddStateVariables( ) ); v != std::end( *get_dMicroYielddStateVariables( ) ); ++v ){

                    col = ( unsigned int )( v - std::begin( *get_dMicroYielddStateVariables( ) ) ) + num_configurations * ( 2 * sot_dim + tot_dim );
                    ( *jacobian.value )[ num_unknowns * row + col ] += *v;

                }

            }

            // Micro gradient yielding
            for ( unsigned int i = 0; i < 3; ++i ){

                if ( ( *activeConstraints )[ 2 + i ] ){

                    row = num_plastic_state_variables + 2 + i;

                    // Stress Jacobians
                    unsigned int col;
                    for (
                        auto v = std::begin( *get_dMicroGradientYielddStress( ) ) + tot_dim * i;
                        v != std::begin( *get_dMicroGradientYielddStress( ) ) + tot_dim * ( i + 1 ); ++v
                    ){
                        
                        col = 2 * sot_dim + ( unsigned int )( v - std::begin( *get_dMicroGradientYielddStress( ) ) - tot_dim  * i );
                        ( *jacobian.value )[ num_unknowns * row + col ] += *v;

                    }

                    // Deformation Jacobians
                    for (
                        auto v = std::begin ( *get_dMicroGradientYielddFn( ) ) + sot_dim * ( num_configurations - 1 ) * i;
                        v != std::begin( *get_dMicroGradientYielddFn( ) ) + sot_dim * ( num_configurations - 1 ) * ( i + 1 ); ++v
                    ){

                        col = ( unsigned int )( v - std::begin( *get_dMicroGradientYielddFn( ) ) - sot_dim * ( num_configurations - 1 ) * i ) + 2 * sot_dim + tot_dim;
                        ( *jacobian.value )[ num_unknowns * row + col ] += *v;

                    }

                    for (
                        auto v = std::begin ( *get_dMicroGradientYielddChin( ) ) + sot_dim * ( num_configurations - 1 ) * i;
                        v != std::begin( *get_dMicroGradientYielddChin( ) ) + sot_dim * ( num_configurations - 1 ) * ( i + 1 ); ++v ){

                        col = ( unsigned int )( v - std::begin( *get_dMicroGradientYielddChin( ) ) - sot_dim * ( num_configurations - 1 ) * i ) + 2 * sot_dim + tot_dim + sot_dim * ( num_configurations - 1 );
                        ( *jacobian.value )[ num_unknowns * row + col ] += *v;

                    }

                    // State variable Jacobians
                    for (
                        auto v = std::begin ( *get_dMicroGradientYielddStateVariables( ) ) + num_ISVS * i;
                        v != std::begin( *get_dMicroGradientYielddStateVariables( ) ) + num_ISVS * ( i + 1 ); ++v ){

                        col = ( unsigned int )( v - std::begin( *get_dMicroGradientYielddStateVariables( ) ) - num_ISVS * i ) + num_configurations * ( 2 * sot_dim + tot_dim );
                        ( *jacobian.value )[ num_unknowns * row + col ] += *v;

                    }

                }

            }

            // Set the inactive constraints
            for ( auto v = std::begin( *activeConstraints ); v != std::end( *activeConstraints ); ++v ){

                row = num_plastic_state_variables + ( unsigned int )( v - std::begin( *activeConstraints ) );

                if ( !( *v ) ){

                    offset = num_configurations * ( 2 * sot_dim + tot_dim );
                    ( *jacobian.value )[ num_unknowns * row + ( unsigned int )( v - std::begin( *activeConstraints ) ) + offset ] += 1;

                }

            }

        }

        void residual::setdStateVariableResidualsdD( ){
            /*!
             * Set the Jacobians of the state variable residuals w.r.t. the deformation measure
             *
             * We define these residuals as
             *
             * \f$R = Z^{\text{update}} - Z \f$
             *
             * and
             *
             * \f$R = f\f$ if the constraint is in the active set and \f$ R = \dot{\bar{\gamma}} \f$ if it isn't
             */

            auto dim = hydra->getDimension( );

            auto sot_dim = dim * dim;

            auto tot_dim = dim * dim * dim;

            auto plasticStrainLikeISVs = get_plasticStrainLikeISVs( );

            auto activeConstraints = get_activeConstraints( );

            auto jacobian = get_setDataStorage_dStateVariableResidualsdD( );

            auto num_ISVS = get_plasticStateVariables( )->size( );

            auto num_plastic_state_variables = ( const unsigned int )( std::end( *plasticStrainLikeISVs ) - std::begin( *plasticStrainLikeISVs ) );

            auto num_cols = ( 2 * sot_dim + tot_dim );

            jacobian.zero( num_ISVS * num_cols );

            // Add the active constraint contributions

            // Macro yielding
            unsigned int row = 0;
            unsigned int col = 0;
            if ( ( *activeConstraints )[ 0 ] ){

                row = num_plastic_state_variables;

                // Deformation Jacobians
                for ( auto v = std::begin ( *get_dMacroYielddF( ) ); v != std::end( *get_dMacroYielddF( ) ); ++v ){

                    col = ( unsigned int )( v - std::begin( *get_dMacroYielddF( ) ) );
                    ( *jacobian.value )[ num_cols * row + col ] += *v;

                }

            }

            // Micro yielding
            if ( ( *activeConstraints )[ 1 ] ){

                row = num_plastic_state_variables + 1;

                // Deformation Jacobians
                for ( auto v = std::begin ( *get_dMicroYielddF( ) ); v != std::end( *get_dMicroYielddF( ) ); ++v ){

                    col = ( unsigned int )( v - std::begin( *get_dMicroYielddF( ) ) );
                    ( *jacobian.value )[ num_cols * row + col ] += *v;

                }

            }

            // Micro gradient yielding
            for ( unsigned int i = 0; i < 3; ++i ){

                if ( ( *activeConstraints )[ 2 + i ] ){

                    row = num_plastic_state_variables + 2 + i;

                    // Deformation Jacobians
                    for (
                        auto v = std::begin ( *get_dMicroGradientYielddF( ) ) + sot_dim * i;
                        v != std::begin( *get_dMicroGradientYielddF( ) ) + sot_dim * ( i + 1 ); ++v
                    ){

                        col = ( unsigned int )( v - std::begin( *get_dMicroGradientYielddF( ) ) - sot_dim * i );
                        ( *jacobian.value )[ num_cols * row + col ] += *v;

                    }

                    for (
                        auto v = std::begin ( *get_dMicroGradientYielddChi( ) ) + sot_dim * i;
                        v != std::begin( *get_dMicroGradientYielddChi( ) ) + sot_dim * ( i + 1 ); ++v ){

                        col = ( unsigned int )( v - std::begin( *get_dMicroGradientYielddChi( ) ) - sot_dim * i ) + sot_dim;
                        ( *jacobian.value )[ num_cols * row + col ] += *v;

                    }

                }

            }

        }

        void residual::correctResiduals( ){
            /*!
             * Correct the residuals to be consistent with the new active set
             */

            setStateVariableResiduals( );

            setResidual( );

            floatVector * globalResidual = hydra->getMutableResidual( );

            TARDIGRADE_ERROR_TOOLS_CHECK( globalResidual, "The global residual isn't mutable" );

            unsigned int offset = hydra->getCurrentResidualOffset( );

            std::copy(
                std::begin( *getResidual( ) ),
                std::end( *getResidual( ) ),
                std::begin( *globalResidual ) + offset
            );

        }

        void residual::successfulNLStep( ){
            /*!
             * Function that runs after a successful nonlinear step
             */

            updateActiveConstraints( );
            correctResiduals( );            

        }

        void residual::preNLSolve( ){
            /*!
             * Function that runs prior to a nonlinear solve
             */

            setActiveConstraints( );

        }

        void residual::postNLSolve( ){
            /*!
             * Function that runs after a nonlinear solve
             */

            // Assemble the yield surface values
            std::array< floatType, 5 > yieldSurfaceValues = {
                *get_macroYield( ), *get_microYield( ),
                ( *get_microGradientYield( ) )[ 0 ], ( *get_microGradientYield( ) )[ 1 ], ( *get_microGradientYield( ) )[ 2 ]
            };

            for ( auto v = std::begin( yieldSurfaceValues ); v != std::end( yieldSurfaceValues ); ++v ){

                if ( ( *v ) > *getYieldTolerance( ) ){

                    throw convergence_error( "A yield surface is larger than the yield surface tolerance" );

                }

            }

        }

    }

}
