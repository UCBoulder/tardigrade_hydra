/**
  ******************************************************************************
  * \file tardigrade_hydraPeryznaViscodamage.cpp
  ******************************************************************************
  * An implementation of peryznaViscodamage using the hydra framework. Used
  * as an example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraPeryznaViscodamage.h>
#include<tardigrade_constitutive_tools.h>

namespace tardigradeHydra{

    namespace peryznaViscodamage{


        void residual::setDamage( ){
            /*!
             * Get the value of the damage from the evolved state variables
             */

            set_damage( ( *get_plasticStateVariables( ) )[ damageISVIndex ] );

        }

        void residual::setDamageJacobians( const bool withPrevious){
            /*!
             * Set the damage along with the jacobians of the damage.
             * 
             * \param withPrevious: Flag for whether to include the derivatives w.r.t. the previous values.
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            if ( withPrevious ){

                set_dDamagedPreviousCauchyStress( tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdPreviousCauchyStress( ), num_isvs, sot_dim, damageISVIndex ) );

                set_dDamagedPreviousF( tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdPreviousF( ), num_isvs, sot_dim, damageISVIndex ) );

                set_dDamagedPreviousSubFs( tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdPreviousSubFs( ), num_isvs, ( num_configs - 1 ) * sot_dim, damageISVIndex ) );

                set_dDamagedPreviousT( ( *get_dPlasticStateVariablesdPreviousT( ) )[ damageISVIndex ] );

                set_dDamagedPreviousStateVariables( tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdPreviousStateVariables( ), num_isvs, num_isvs, damageISVIndex ) );

            }

            set_damage( ( *get_plasticStateVariables( ) )[ damageISVIndex ] );

            set_dDamagedCauchyStress( tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdCauchyStress( ), num_isvs, sot_dim, damageISVIndex ) );

            set_dDamagedF( tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdF( ), num_isvs, sot_dim, damageISVIndex ) );

            set_dDamagedSubFs( tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdSubFs( ), num_isvs, ( num_configs - 1 ) * sot_dim, damageISVIndex ) );

            set_dDamagedT( ( *get_dPlasticStateVariablesdT( ) )[ damageISVIndex ] );

            set_dDamagedStateVariables( tardigradeVectorTools::getRow( *get_dPlasticStateVariablesdStateVariables( ), num_isvs, num_isvs, damageISVIndex ) );

        }

        void residual::setDamageJacobians( ){
            /*!
             * Set the damage along with the jacobians of the damage.
             * This routine will not set the jacobians w.r.t. the previous values.
             */

            setDamageJacobians( false );

        }

        void residual::setAllDamageJacobians( ){
            /*!
             * Set the damage along with the jacobians of the damage.
             * This routine will set the jacobians w.r.t. the previous values.
             */

            setDamageJacobians( true );

        }

        void residual::setDamageDeformationGradient( ){
            /*!
             * Set the deformation gradient of the damage
             * 
             * We define the damage through the strain as
             * 
             * \f$ \bf{E}^d = D \bf{E}^r \f$
             * 
             * where \f$\bf{E}^d\f$ is the Green-Lagrange damage strain,
             * \f$D\f$ is the damage and \f$\bf{E}^r\f$ is the recoverable
             * strain and is defined as
             * 
             * \f$ \bf{E}^r = \bf{E}^e + \bf{E}^d \f$ where \f$\bf{E}^e\f$
             * is the elastic strain. We can solve for the damage deformation
             * gradient using the definition of the Green-Lagrange strain and
             * knowledge of the elastic configuration.
             */

            const unsigned int dim = hydra->getDimension( );
            const unsigned int sot_dim = hydra->getSOTDimension( );
            const unsigned int elastic_config_index = *getElasticConfigurationIndex( );

            // Get the elastic deformation gradient
            floatVector Fe = floatVector( hydra->get_configurations( )->begin( ) + sot_dim * elastic_config_index,
                                          hydra->get_configurations( )->begin( ) + sot_dim * ( elastic_config_index + 1 ) );
    
            // Compute the elastic Green-Lagrange strain
            floatVector Ee;
    
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::computeGreenLagrangeStrain( Fe, Ee ) );
    
            // Compute the damage strain
            floatVector Ed = ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) * Ee;
    
            // Compute the square root to solve for the damage deformation gradient
            floatVector eye( sot_dim );
            tardigradeVectorTools::eye( eye );
    
            floatVector Fd;

            TARDIGRADE_ERROR_TOOLS_CATCH( Fd = tardigradeVectorTools::matrixSqrt( 2.0 * Ed + eye, dim ) );
    
            set_damageDeformationGradient( Fd );

        }

        void residual::setDamageDeformationGradientJacobians( ){
            /*!
             * Set the jacobians for the damage deformation gradient
             */

            setDamageDeformationGradientJacobians( false );

        }

        void residual::setAllDamageDeformationGradientJacobians( ){
            /*!
             * Set all of the jacobians for the damage deformation gradient
             */

            setDamageDeformationGradientJacobians( true );

        }

        void residual::setDamageDeformationGradientJacobians( const bool withPrevious ){
            /*!
             * Set the deformation gradient of the damage
             * 
             * We define the damage through the strain as
             * 
             * \f$ \bf{E}^d = D \bf{E}^r \f$
             * 
             * where \f$\bf{E}^d\f$ is the Green-Lagrange damage strain,
             * \f$D\f$ is the damage and \f$\bf{E}^r\f$ is the recoverable
             * strain and is defined as
             * 
             * \f$ \bf{E}^r = \bf{E}^e + \bf{E}^d \f$ where \f$\bf{E}^e\f$
             * is the elastic strain. We can solve for the damage deformation
             * gradient using the definition of the Green-Lagrange strain and
             * knowledge of the elastic configuration.
             * 
             * \param withPrevious: Flag for whether to set the Jacobians w.r.t. the previous unknowns
             */

            const unsigned int dim = hydra->getDimension( );

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int elastic_config_index = *getElasticConfigurationIndex( );
    
            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            // Get the elastic deformation gradient
            floatVector Fe = floatVector( hydra->get_configurations( )->begin( ) + sot_dim * elastic_config_index,
                                          hydra->get_configurations( )->begin( ) + sot_dim * ( elastic_config_index + 1 ) );
    
            // Compute the elastic Green-Lagrange strain
            floatVector Ee;

            floatVector dEedFe;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::computeGreenLagrangeStrain( Fe, Ee, dEedFe ) );

            floatVector dFedF( sot_dim * sot_dim, 0 );

            floatVector dFedSubFs( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            if ( elastic_config_index == 0 ){

                dFedF     = *hydra->get_dF1dF( );

                dFedSubFs = *hydra->get_dF1dFn( );

            }
            else{

                for ( unsigned int i = 0; i < sot_dim; i++ ){

                    dFedSubFs[ ( num_configs - 1 ) * sot_dim * i + i + elastic_config_index - 1 ] = 1;

                }

            }
 
            // Compute the damage strain
            floatVector Ed = ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) * Ee;

            floatVector dEddD = 1 / ( 1 - ( *get_damage( ) ) ) * ( 1 + ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) ) * Ee;

            floatVector dEddF = ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) * tardigradeVectorTools::matrixMultiply( dEedFe, dFedF, sot_dim, sot_dim, sot_dim, sot_dim );
 
            floatVector dEddSubFs = ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) * tardigradeVectorTools::matrixMultiply( dEedFe, dFedSubFs, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );
 
            // Compute the square root to solve for the damage deformation gradient
            floatVector eye( sot_dim );
            tardigradeVectorTools::eye( eye );
    
            floatVector Fd;

            floatMatrix _dAdFe; //A = 2.0 * Ed + eye
                                //
            floatVector dAdFe; //A = 2.0 * Ed + eye

            TARDIGRADE_ERROR_TOOLS_CATCH( Fd = tardigradeVectorTools::matrixSqrt( 2.0 * Ed + eye, dim, _dAdFe ) );

            dAdFe = tardigradeVectorTools::appendVectors( _dAdFe );

            floatVector dFddEd;

            TARDIGRADE_ERROR_TOOLS_CATCH( dFddEd = 2 * tardigradeVectorTools::inverse( dAdFe, sot_dim, sot_dim ) );

            floatVector dFddD = tardigradeVectorTools::matrixMultiply( dFddEd, dEddD, sot_dim, sot_dim, sot_dim, 1 );

            floatVector dFddF = tardigradeVectorTools::matrixMultiply( dFddEd, dEddF, sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dFddSubFs = tardigradeVectorTools::matrixMultiply( dFddEd, dEddSubFs, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            if ( withPrevious ){

                set_dDamageDeformationGradientdPreviousCauchyStress( tardigradeVectorTools::matrixMultiply( dFddD, *get_dDamagedPreviousCauchyStress( ), sot_dim, 1, 1, sot_dim ) );
    
                set_dDamageDeformationGradientdPreviousF( tardigradeVectorTools::matrixMultiply( dFddD, *get_dDamagedPreviousF( ), sot_dim, 1, 1, sot_dim ) );
    
                set_dDamageDeformationGradientdPreviousSubFs( tardigradeVectorTools::matrixMultiply( dFddD, *get_dDamagedPreviousSubFs( ), sot_dim, 1, 1, ( num_configs - 1 ) * sot_dim ) );
    
                set_dDamageDeformationGradientdPreviousT( dFddD * ( *get_dDamagedPreviousT( ) ) );
    
                set_dDamageDeformationGradientdPreviousStateVariables( tardigradeVectorTools::matrixMultiply( dFddD, *get_dDamagedPreviousStateVariables( ), sot_dim, 1, 1, num_isvs ) );

            }

            set_damageDeformationGradient( Fd );

            set_dDamageDeformationGradientdCauchyStress( tardigradeVectorTools::matrixMultiply( dFddD, *get_dDamagedCauchyStress( ), sot_dim, 1, 1, sot_dim ) );

            set_dDamageDeformationGradientdF( tardigradeVectorTools::matrixMultiply( dFddD, *get_dDamagedF( ), sot_dim, 1, 1, sot_dim ) + dFddF );

            set_dDamageDeformationGradientdSubFs( tardigradeVectorTools::matrixMultiply( dFddD, *get_dDamagedSubFs( ), sot_dim, 1, 1, ( num_configs - 1 ) * sot_dim ) + dFddSubFs );

            set_dDamageDeformationGradientdT( dFddD * ( *get_dDamagedT( ) ) );

            set_dDamageDeformationGradientdStateVariables( tardigradeVectorTools::matrixMultiply( dFddD, *get_dDamagedStateVariables( ), sot_dim, 1, 1, num_isvs ) );

        }

        void residual::setStateVariableEvolutionRates( const bool isPrevious ){
            /*!
             * Set the state variable evolution rates including the damage
             * 
             * \param &isPrevious: Whether to compute the previous values or not
             */

            floatVector evolutionRates;

            tardigradeHydra::peryznaViscoplasticity::residual::setStateVariableEvolutionRates( isPrevious );

            if ( isPrevious ){

                evolutionRates = { ( *get_previousStateVariableEvolutionRates( ) )[ 0 ], *get_previousPlasticMultiplier( ) };

                set_previousStateVariableEvolutionRates( evolutionRates );

            }
            else{

                evolutionRates = { ( *get_stateVariableEvolutionRates( ) )[ 0 ], *get_plasticMultiplier( ) };

                set_stateVariableEvolutionRates( evolutionRates );

            }

        }

        void residual::setStateVariableEvolutionRateDerivatives( const bool isPrevious ){
            /*!
             * Set the jacobians of the state variable evolution rates including the damage
             * 
             * \param &isPrevious: Whether to compute the previous values or not
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            floatVector tempJac;

            tardigradeHydra::peryznaViscoplasticity::residual::setStateVariableEvolutionRateDerivatives( isPrevious );

            const floatVector *stateVariableEvolutionRates;

            const floatType   *plasticMultiplier;

            const floatVector *dPlasticMultiplierdCauchyStress;

            setDataStorageBase< floatVector > evolutionRates;

            setDataStorageBase< floatVector > dStateVariableEvolutionRatesdCauchyStress;

            if ( isPrevious ){

                stateVariableEvolutionRates = get_previousStateVariableEvolutionRates( );

                plasticMultiplier           = get_previousPlasticMultiplier( );

                evolutionRates              = get_setDataStorage_previousStateVariableEvolutionRates( );

                dPlasticMultiplierdCauchyStress = get_dPreviousPlasticMultiplierdPreviousCauchyStress( );

                dStateVariableEvolutionRatesdCauchyStress = get_setDataStorage_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( );

            }
            else{

                stateVariableEvolutionRates = get_stateVariableEvolutionRates( );

                plasticMultiplier           = get_plasticMultiplier( );

                evolutionRates              = get_setDataStorage_stateVariableEvolutionRates( );

                dPlasticMultiplierdCauchyStress = get_dPlasticMultiplierdCauchyStress( );

                dStateVariableEvolutionRatesdCauchyStress = get_setDataStorage_dStateVariableEvolutionRatesdCauchyStress( );

            }

            *evolutionRates.value = { ( *stateVariableEvolutionRates )[ 0 ], *plasticMultiplier };

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                ( *dStateVariableEvolutionRatesdCauchyStress.value )[ sot_dim + i ] = ( *dPlasticMultiplierdCauchyStress )[ i ];

            }

            if ( isPrevious ){

                // Set the derivatives w.r.t. the previous deformation gradient
                tempJac = *get_dPreviousStateVariableEvolutionRatesdPreviousF( );
                for ( unsigned int i = 0; i < sot_dim; i++ ){
                    tempJac[ sot_dim + i ] = ( *get_dPreviousPlasticMultiplierdPreviousF( ) )[ i ];
                }

                set_dPreviousStateVariableEvolutionRatesdPreviousF( tempJac );

                // Set the derivatives w.r.t. the previous sub-deformation gradients
                tempJac = *get_dPreviousStateVariableEvolutionRatesdPreviousSubFs( );
                for ( unsigned int i = 0; i < ( num_configs - 1 ) * sot_dim; i++ ){
                    tempJac[ ( num_configs - 1 ) * sot_dim + i ] = ( *get_dPreviousPlasticMultiplierdPreviousSubFs( ) )[ i ];
                }

                set_dPreviousStateVariableEvolutionRatesdPreviousSubFs( tempJac );

                // Set the derivatives w.r.t. the previous temperature
                floatVector tempJacVec = tardigradeVectorTools::appendVectors( { *get_dPreviousStateVariableEvolutionRatesdPreviousT( ),
                                                                               { *get_dPreviousPlasticMultiplierdPreviousT( ) } } );

                set_dPreviousStateVariableEvolutionRatesdPreviousT( tempJacVec );

                // Set the derivatives w.r.t. the previous state variables
                tempJac = { ( *get_dPreviousStateVariableEvolutionRatesdPreviousStateVariables( ) )[ 0 ], 0.,
                            ( *get_dPreviousPlasticMultiplierdPreviousStateVariables( ) )[ 0 ], 0. };

                set_dPreviousStateVariableEvolutionRatesdPreviousStateVariables( tempJac );

            }
            else{

                // Set the derivatives w.r.t. the deformation gradient
                tempJac = *get_dStateVariableEvolutionRatesdF( );
                for ( unsigned int i = 0; i < sot_dim; i++ ){
                    tempJac[ sot_dim + i ] = ( *get_dPlasticMultiplierdF( ) )[ i ];
                }

                set_dStateVariableEvolutionRatesdF( tempJac );

                // Set the derivatives w.r.t. the sub-deformation gradients
                tempJac = *get_dStateVariableEvolutionRatesdSubFs( );
                for ( unsigned int i = 0; i < ( num_configs - 1 ) * sot_dim; i++ ){
                    tempJac[ ( num_configs - 1 ) * sot_dim + i ] = ( *get_dPlasticMultiplierdSubFs( ) )[ i ];
                }

                set_dStateVariableEvolutionRatesdSubFs( tempJac );

                // Set the derivatives w.r.t. the temperature
                floatVector tempJacVec = tardigradeVectorTools::appendVectors( { *get_dStateVariableEvolutionRatesdT( ),
                                                                               { *get_dPlasticMultiplierdT( ) } } );

                set_dStateVariableEvolutionRatesdT( tempJacVec );

                // Set the derivatives w.r.t. the state variables
                tempJac = { ( *get_dStateVariableEvolutionRatesdStateVariables( ) )[ 0 ], 0.,
                            ( *get_dPlasticMultiplierdStateVariables( ) )[ 0 ], 0. };

                set_dStateVariableEvolutionRatesdStateVariables( tempJac );

            }

        }

        void residual::setResidual( ){
            /*!
             * Set the residual vector
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            auto residual = get_setDataStorage_residual( );
            residual.zero( *getNumEquations( ) );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                ( *residual.value )[ i ] = hydra->getConfiguration( *getDamageConfigurationIndex( ) )[ i ] - ( *get_damageDeformationGradient( ) )[ i ];

            }

            for ( unsigned int i = 0; i < num_isvs; i++ ){

                ( *residual.value )[ sot_dim + i ] = ( *get_stateVariables( ) )[ i ] - ( *get_plasticStateVariables( ) )[ i ];

            }

        }

        void residual::setJacobian( ){
            /*!
             * Set the Jacobian matrix
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int damage_configuration_index = *getDamageConfigurationIndex( );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            const unsigned int num_unknowns = hydra->getNumUnknowns( );

            auto jacobian = get_setDataStorage_jacobian( );
            jacobian.zero( *getNumEquations( ) * num_unknowns );

            // Jacobians of the damage deformation gradient
            for ( unsigned int i = 0; i < sot_dim; i++ ){
                unsigned int row = i;

                // Jacobians w.r.t. the Cauchy stress
                for ( unsigned int j = 0; j < sot_dim; j++ ){
                    unsigned int col = j;

                    ( *jacobian.value )[ num_unknowns * row + col ] -= ( *get_dDamageDeformationGradientdCauchyStress( ) )[ sot_dim * i + j ];

                }

                // Jacobians w.r.t. the sub configurations
                ( *jacobian.value )[ num_unknowns * row + sot_dim * damage_configuration_index + i ] += 1;

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){
                    unsigned int col = sot_dim + j;

                    ( *jacobian.value )[ num_unknowns * row + col ] -= ( *get_dDamageDeformationGradientdSubFs( ) )[ ( num_configs - 1 ) * sot_dim * i + j ];

                }

                // Jacobians w.r.t. the state variables
                for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){
                    unsigned int col = sot_dim + ( num_configs - 1 ) * sot_dim + *ind;

                    ( *jacobian.value )[ num_unknowns * row + col ] -= ( *get_dDamageDeformationGradientdStateVariables( ) )[ num_isvs * i + ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ];

                }

            }

            // Jacobians of the damage hardening state variables
            for ( unsigned int i = 0; i < num_isvs; i++ ){
                unsigned int row = sot_dim + i;

                // Jacobians w.r.t. the Cauchy stress
                for ( unsigned int j = 0; j < sot_dim; j++ ){
                    unsigned int col = j;

                    ( *jacobian.value )[ num_unknowns * row + col ] -= ( *get_dPlasticStateVariablesdCauchyStress( ) )[ sot_dim * i + j ];

                }

                // Jacobians w.r.t. the sub configurations
                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){
                    unsigned int col = sot_dim + j;

                    ( *jacobian.value )[ num_unknowns * row + col ] -= ( *get_dPlasticStateVariablesdSubFs( ) )[ ( num_configs - 1 ) * sot_dim * i + j ];

                }

                // Jacobians w.r.t. the state variables
                ( *jacobian.value )[ num_unknowns * row + sot_dim + ( num_configs - 1 ) * sot_dim + ( *getStateVariableIndices( ) )[ i ] ] += 1;
                for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){
                    unsigned int col = sot_dim + ( num_configs - 1 ) * sot_dim + *ind;

                    ( *jacobian.value )[ num_unknowns * row + col ] -= ( *get_dPlasticStateVariablesdStateVariables( ) )[ num_isvs * i + ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ];

                }

            }

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_isvs = get_dPlasticStateVariablesdT( )->size( );

            auto dRdT = get_setDataStorage_dRdT( );

            dRdT.zero( *getNumEquations( ) );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                ( *dRdT.value )[ i ] -= ( *get_dDamageDeformationGradientdT( ) )[ i ];

            }

            for ( unsigned int i = 0; i < num_isvs; i++ ){

                ( *dRdT.value )[ sot_dim + i ] = -( *get_dPlasticStateVariablesdT( ) )[ i ];

            }

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            auto dRdF = get_setDataStorage_dRdF( );
            dRdF.zero( ( *getNumEquations( ) ) * sot_dim );

            *dRdF.value = tardigradeVectorTools::appendVectors( { *get_dDamageDeformationGradientdF( ),
                                                                  *get_dPlasticStateVariablesdF( ) } );

            std::transform( dRdF.value->cbegin( ), dRdF.value->cend( ), dRdF.value->begin( ), std::negate<floatType>( ) );

        }

        void residual::decomposeParameters( const floatVector &parameters ){
            /*!
             * Decompose the incoming parameter vector
             * 
             * \param &parameters: The incoming parameter vector
             * 
             * The form has the same interpretation as the base viscoplastic case
             * except we have a different method of hardening for the damage state
             * variable which we account for here.
             */

            tardigradeHydra::peryznaViscoplasticity::residual::decomposeParameters( parameters );

            //Setting the contribution of damage to the calculation of the drag stress to zero
            set_dragStressParameters( { parameters[ 1 ], parameters[  2 ], 0. } );

            //Setting the contribution of damage to the hardening of the damage state variable to zero
            set_hardeningParameters(  { parameters[ 9 ], parameters[ 10 ], 0. } );

        }

    }

}
