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

        void residual::setDeltaIntegratedPlasticMultipliers( ){
            /*!
             * Set the net change in the integrated plastic multipliers
             * using trapezoidal integration
             * 
             * 
             */

            auto plasticMultipliers         = get_plasticMultipliers( );

            auto previousPlasticMultipliers = get_previousPlasticMultipliers( );

            auto deltaIntegratedPlasticMultipliers = get_setDataStorage_deltaIntegratedPlasticMultipliers( );

            auto num_plastic_multipliers = ( unsigned int )( std::end( *plasticMultipliers ) - std::begin( *plasticMultipliers ) );

            deltaIntegratedPlasticMultipliers.zero( num_plastic_multipliers );

            for ( unsigned int i = 0; i < num_plastic_multipliers; ++i ){

                ( *deltaIntegratedPlasticMultipliers.value )[ i ] = ( *hydra->getDeltaTime( ) ) * ( ( 1 - ( *getIntegrationParameter( ) ) ) * ( *previousPlasticMultipliers )[ i ] + ( *getIntegrationParameter( ) ) * ( *plasticMultipliers )[ i ] );

            }

        }

        void residual::setdDeltaIntegratedPlasticMultipliersdPlasticMultipliers( ){
            /*!
             * Set the derivative of the net change in the integrated plastic multipliers
             * using trapezoidal integration w.r.t. the plastic multipliers
             * 
             * This just returns the diagonal of the Jacobian since it's the only non-zero part
             */

            auto plasticMultipliers         = get_plasticMultipliers( );

            auto dDeltaIntegratedPlasticMultipliersdPlasticMultipliers = get_setDataStorage_dDeltaIntegratedPlasticMultipliersdPlasticMultipliers( );

            auto num_plastic_multipliers = ( unsigned int )( std::end( *plasticMultipliers ) - std::begin( *plasticMultipliers ) );

            dDeltaIntegratedPlasticMultipliersdPlasticMultipliers.zero( num_plastic_multipliers );

            std::fill(
                dDeltaIntegratedPlasticMultipliersdPlasticMultipliers.begin( ),
                dDeltaIntegratedPlasticMultipliersdPlasticMultipliers.end( ),
                ( *hydra->getDeltaTime( ) ) * ( *getIntegrationParameter( ) )
            );

        }

        void residual::setActiveConstraints( ){
            /*!
             * Set which constraints are active
             */

            auto activeConstraints   = get_setDataStorage_activeConstraints( );

            activeConstraints.value->resize( 5 );

            ( *activeConstraints.value )[ 0 ] = ( *get_macroYield( ) ) > 0;

            ( *activeConstraints.value )[ 1 ] = ( *get_microYield( ) ) > 0;

            for ( auto v = std::begin( *get_microGradientYield( ) ); v != std::end( *get_microGradientYield( ) ); ++v ){

                ( *activeConstraints.value )[ ( v - std::begin( *get_microGradientYield( ) ) ) + 2 ] = ( *v ) > 0;

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

            auto activeConstraints = get_activeConstraints( );

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
            for ( auto v = std::begin( *activeConstraints ); v != std::end( *activeConstraints ); ++v ){

                if ( *v ){
                    ( *residual.value )[ num_plastic_state_variables + ( unsigned int )( v - std::begin( *activeConstraints ) ) ] += yieldSurfaceValues[ ( unsigned int )( v - std::begin( *activeConstraints ) )  ];
                }
                else{
                    ( *residual.value )[ num_plastic_state_variables + ( unsigned int )( v - std::begin( *activeConstraints ) ) ] += ( *plasticMultipliers )[ ( unsigned int )( v - std::begin( *activeConstraints ) ) ];
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

        void residual::setDrivingStresses( const bool isPrevious ){
            /*!
             * Set the driving stresses
             * 
             * \param isPrevious: Flag for if the previous values are to be set (currently throws not implemented error)
             */

            if ( isPrevious ){

                tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::setDrivingStresses( isPrevious );

            }
            else{

                tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::setDrivingStresses(
                    isPrevious, *get_setDataStorage_followingDeformationGradient( ).value, *get_setDataStorage_followingMicroDeformation( ).value
                );
            }

        }

        void residual::setDrivingStressesJacobians( const bool isPrevious ){
            /*!
             * Set the driving stresses and jacobians
             * 
             * \param isPrevious: Flag for if the previous values are to be set (currently throws not implemented error)
             */

            if ( isPrevious ){

                tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::setDrivingStressesJacobians( isPrevious );

            }
            else{

                tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::setDrivingStressesJacobians(
                    isPrevious,
                    *get_setDataStorage_followingDeformationGradient( ).value, *get_setDataStorage_followingMicroDeformation( ).value,
                    *get_setDataStorage_dFfollowdF( ).value,                   *get_setDataStorage_dFfollowdFn( ).value,
                    *get_setDataStorage_dChifollowdChi( ).value,               *get_setDataStorage_dChifollowdChin( ).value
                );

            }

        }

        void residual::setFollowingDeformationGradient( ){
            /*!
             * Set the following deformation gradient
             */

            setDrivingStresses( false );

        }

        void residual::setdFfollowdF( ){
            /*!
             * Set the derivative of the following deformation gradient w.r.t. the deformation gradient
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdFfollowdFn( ){
            /*!
             * Set the derivative of the following deformation gradient w.r.t. the sub deformation gradients
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setFollowingMicroDeformation( ){
            /*!
             * Set the following micro deformation
             */

            setDrivingStresses( false );

        }

        void residual::setdChifollowdChi( ){
            /*!
             * Set the derivative of the following micro deformation w.r.t. the micro deformation
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdChifollowdChin( ){
            /*!
             * Set the derivative of the following micro deformation w.r.t. the sub micro deformations
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setMacroYieldStressGradient( ){
            /*!
             * Set the partial derivative of the macro yield function w.r.t. the stress
             */

            auto dim = hydra->getDimension( );

            auto sot_dim = dim * dim;

            auto macroCohesion = get_macroCohesion( );

            auto macroYieldParameters = get_macroYieldParameters( );

            auto macroDrivingStress = get_macroDrivingStress( );

            auto precedingDeformationGradient = get_precedingDeformationGradient( );

            auto macroYieldStressGradient = get_setDataStorage_macroYieldStressGradient( );

            floatType tempYield;

            floatVector dMacroYielddDrivingStress;

            floatType dMacroYielddCohesion;

            floatVector dMacroYielddPrecedingF;

            floatVector dDrivingStressdStress;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation(
                    *macroDrivingStress, *macroCohesion, *precedingDeformationGradient,
                    ( *macroYieldParameters )[ 0 ], ( *macroYieldParameters )[ 1 ],
                    tempYield, dMacroYielddDrivingStress, dMacroYielddCohesion, dMacroYielddPrecedingF
                )

            );

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeMicromorphicTools::dCauchyStressdPK2Stress( *get_followingDeformationGradient( ), dDrivingStressdStress );
            )

            macroYieldStressGradient.zero( sot_dim );
            for ( unsigned int IJ = 0; IJ < sot_dim; ++IJ ){
                for ( unsigned int KL = 0; KL < sot_dim; ++KL ){
                    ( *macroYieldStressGradient.value )[ KL ] += dMacroYielddDrivingStress[ IJ ] * dDrivingStressdStress[ sot_dim * IJ + KL ];
                }
            }

        }

        void residual::setdMacroYieldStressGradientdStress( ){
            /*!
             * Set the jacobian of the macro yield stress gradient w.r.t. the stress
             */

            setMacroYieldStressGradientJacobians( );

        }

        void residual::setdMacroYieldStressGradientdFn( ){
            /*!
             * Set the jacobian of the macro yield stress gradient w.r.t. the sub deformation gradients
             */

            setMacroYieldStressGradientJacobians( );

        }

        void residual::setdMacroYieldStressGradientdF( ){
            /*!
             * Set the jacobian of the macro yield stress gradient w.r.t. the deformation gradient
             */

            setdMacroYieldStressGradientdD( );

        }

        void residual::setMacroYieldStressGradientJacobians( ){
            /*!
             * Set the jacobians of the gradient of the macro yield function w.r.t. the stress
             */

            auto dim = hydra->getDimension( );

            auto sot_dim = dim * dim;

            auto macroCohesion = get_macroCohesion( );

            auto macroYieldParameters = get_macroYieldParameters( );

            auto macroDrivingStress = get_macroDrivingStress( );

            auto precedingDeformationGradient = get_precedingDeformationGradient( );

            auto macroYieldStressGradient = get_setDataStorage_macroYieldStressGradient( );

            auto dMacroYieldStressGradientdStress = get_setDataStorage_dMacroYieldStressGradientdStress( );

            auto dMacroYieldStressGradientdFn = get_setDataStorage_dMacroYieldStressGradientdFn( );

            floatType tempYield;

            floatVector dMacroYielddDrivingStress;

            floatType dMacroYielddCohesion;

            floatVector dMacroYielddPrecedingF;

            floatVector dDrivingStressdStress;

            floatVector d2MacroYielddDrivingStress2;

            floatVector d2MacroYielddDrivingStressdPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation(
                    *macroDrivingStress, *macroCohesion, *precedingDeformationGradient,
                    ( *macroYieldParameters )[ 0 ], ( *macroYieldParameters )[ 1 ],
                    tempYield, dMacroYielddDrivingStress, dMacroYielddCohesion, dMacroYielddPrecedingF,
                    d2MacroYielddDrivingStress2, d2MacroYielddDrivingStressdPrecedingF
                )

            );

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeMicromorphicTools::dCauchyStressdPK2Stress( *get_followingDeformationGradient( ), dDrivingStressdStress );
            )

            floatVector dMacroYieldStressGradientdDrivingStress( sot_dim * sot_dim, 0 );

            macroYieldStressGradient.zero( sot_dim );
            dMacroYieldStressGradientdStress.zero( sot_dim * sot_dim );
            for ( unsigned int IJ = 0; IJ < sot_dim; ++IJ ){

                for ( unsigned int KL = 0; KL < sot_dim; ++KL ){

                    ( *macroYieldStressGradient.value )[ KL ] += dMacroYielddDrivingStress[ IJ ] * dDrivingStressdStress[ sot_dim * IJ + KL ];

                    for ( unsigned int MN = 0; MN < sot_dim; ++MN ){

                        dMacroYieldStressGradientdDrivingStress[ sot_dim * KL + MN ] += d2MacroYielddDrivingStress2[ sot_dim * IJ + MN ] * dDrivingStressdStress[ sot_dim * IJ + KL ];

                    }

                }

            }

            for ( unsigned int IJ = 0; IJ < sot_dim; ++IJ ){
                for ( unsigned int KL = 0; KL < sot_dim; ++KL ){
                    for ( unsigned int MN = 0; MN < sot_dim; ++MN ){
                        ( *dMacroYieldStressGradientdStress.value )[ sot_dim * IJ + MN ] += dMacroYieldStressGradientdDrivingStress[ sot_dim * IJ + KL ] * dDrivingStressdStress[ sot_dim * KL + MN ];
                    }
                }
            }

        }

        void residual::setdMacroYieldStressGradientdD( ){
            /*!
             * Set the Jacobians of the macro yield stress gradient w.r.t. the fundamental deformation measures
             */
        }

        void residual::setMicroYieldStressGradient( ){
            /*!
             * Set the partial derivative of the micro yield function w.r.t. the stress
             */

            auto dim = hydra->getDimension( );

            auto sot_dim = dim * dim;

            auto microCohesion = get_microCohesion( );

            auto microYieldParameters = get_microYieldParameters( );

            auto microDrivingStress = get_symmetricMicroDrivingStress( );

            auto precedingDeformationGradient = get_precedingDeformationGradient( );

            auto microYieldStressGradient = get_setDataStorage_microYieldStressGradient( );

            floatType tempYield;

            floatVector dMicroYielddDrivingStress;

            floatType dMicroYielddCohesion;

            floatVector dMicroYielddPrecedingF;

            floatVector dDrivingStressdStress;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation(
                    *microDrivingStress, *microCohesion, *precedingDeformationGradient,
                    ( *microYieldParameters )[ 0 ], ( *microYieldParameters )[ 1 ],
                    tempYield, dMicroYielddDrivingStress, dMicroYielddCohesion, dMicroYielddPrecedingF
                )
            );

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeMicromorphicTools::dSymmetricMicroStressdReferenceSymmetricMicroStress( *get_followingDeformationGradient( ), dDrivingStressdStress );
            )

            microYieldStressGradient.zero( sot_dim );
            for ( unsigned int IJ = 0; IJ < sot_dim; ++IJ ){
                for ( unsigned int KL = 0; KL < sot_dim; ++KL ){
                    ( *microYieldStressGradient.value )[ KL ] += dMicroYielddDrivingStress[ IJ ] * dDrivingStressdStress[ sot_dim * IJ + KL ];
                }
            }

        }

        void residual::setMicroGradientYieldStressGradient( ){
            /*!
             * Set the partial derivative of the micro gradient yield function w.r.t. the stress
             */

            auto dim = hydra->getDimension( );

            auto tot_dim = dim * dim * dim;

            auto microGradientCohesion = get_microGradientCohesion( );

            auto microGradientYieldParameters = get_microGradientYieldParameters( );

            auto microGradientDrivingStress = get_higherOrderDrivingStress( );

            auto precedingDeformationGradient = get_precedingDeformationGradient( );

            auto microGradientYieldStressGradient = get_setDataStorage_microGradientYieldStressGradient( );

            dimVector tempVectorYield;

            floatVector dMicroGradientYielddDrivingStress;

            floatVector dMicroGradientYielddCohesion;

            floatVector dMicroGradientYielddPrecedingF;

            floatVector dDrivingStressdStress;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation(
                    *microGradientDrivingStress, *microGradientCohesion, *precedingDeformationGradient,
                    ( *microGradientYieldParameters )[ 0 ], ( *microGradientYieldParameters )[ 1 ],
                    tempVectorYield, dMicroGradientYielddDrivingStress, dMicroGradientYielddCohesion, dMicroGradientYielddPrecedingF
                )
            );

            TARDIGRADE_ERROR_TOOLS_CATCH(
                tardigradeMicromorphicTools::dHigherOrderStressdReferenceHigherOrderStress( *get_followingDeformationGradient( ), *get_followingMicroDeformation( ), dDrivingStressdStress );
            )

            microGradientYieldStressGradient.zero( dim * tot_dim );
            for ( unsigned int I = 0; I < dim; ++I ){
                for ( unsigned int JKL = 0; JKL < tot_dim; ++JKL ){
                    for ( unsigned int MNO = 0; MNO < tot_dim; ++MNO ){
                        ( *microGradientYieldStressGradient.value )[ tot_dim * I + MNO ] += dMicroGradientYielddDrivingStress[ tot_dim * I + JKL ] * dDrivingStressdStress[ tot_dim * JKL + MNO ];
                    }
                }
            }

        }

    }

}
