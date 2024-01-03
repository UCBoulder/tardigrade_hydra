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

            if ( withPrevious ){

                set_dDamagedPreviousCauchyStress( ( *get_dPlasticStateVariablesdPreviousCauchyStress( ) )[ damageISVIndex ] );

                set_dDamagedPreviousF( ( *get_dPlasticStateVariablesdPreviousF( ) )[ damageISVIndex ] );

                set_dDamagedPreviousSubFs( ( *get_dPlasticStateVariablesdPreviousSubFs( ) )[ damageISVIndex ] );

                set_dDamagedPreviousT( ( *get_dPlasticStateVariablesdPreviousT( ) )[ damageISVIndex ] );

                set_dDamagedPreviousStateVariables( ( *get_dPlasticStateVariablesdPreviousStateVariables( ) )[ damageISVIndex ] );

            }

            set_damage( ( *get_plasticStateVariables( ) )[ damageISVIndex ] );

            set_dDamagedCauchyStress( ( *get_dPlasticStateVariablesdCauchyStress( ) )[ damageISVIndex ] );

            set_dDamagedF( ( *get_dPlasticStateVariablesdF( ) )[ damageISVIndex ] );

            set_dDamagedSubFs( ( *get_dPlasticStateVariablesdSubFs( ) )[ damageISVIndex ] );

            set_dDamagedT( ( *get_dPlasticStateVariablesdT( ) )[ damageISVIndex ] );

            set_dDamagedStateVariables( ( *get_dPlasticStateVariablesdStateVariables( ) )[ damageISVIndex ] );

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

            const unsigned int *dim = hydra->getDimension( );
    
            // Get the elastic deformation gradient
            floatVector Fe = ( *hydra->get_configurations( ) )[ *getElasticConfigurationIndex( ) ];
    
            // Compute the elastic Green-Lagrange strain
            floatVector Ee;
    
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeGreenLagrangeStrain( Fe, Ee ) );
    
            // Compute the damage strain
            floatVector Ed = ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) * Ee;
    
            // Compute the square root to solve for the damage deformation gradient
            floatVector eye( ( *dim ) * ( *dim ) );
            tardigradeVectorTools::eye( eye );
    
            floatVector Fd;

            TARDIGRADE_ERROR_TOOLS_CATCH( Fd = tardigradeVectorTools::matrixSqrt( 2.0 * Ed + eye, *dim ) );
    
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

            const unsigned int *dim = hydra->getDimension( );
    
            // Get the elastic deformation gradient
            floatVector Fe = ( *hydra->get_configurations( ) )[ *getElasticConfigurationIndex( ) ];
    
            // Compute the elastic Green-Lagrange strain
            floatVector Ee;

            floatMatrix dEedFe;
    
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeGreenLagrangeStrain( Fe, Ee, dEedFe ) );

            floatMatrix dFedF( Fe.size( ), floatVector( hydra->getDeformationGradient( )->size( ), 0 ) );;

            floatMatrix dFedSubFs( Fe.size( ), floatVector( ( *hydra->getNumConfigurations( ) ) * hydra->getDeformationGradient( )->size( ), 0 ) ) ;

            if ( ( *getElasticConfigurationIndex( ) ) == 0 ){

                dFedF     = *hydra->get_dF1dF( );

                dFedSubFs = *hydra->get_dF1dFn( );

            }
            else{

                for ( unsigned int i = 0; i < hydra->getDeformationGradient( )->size( ); i++ ){

                    dFedSubFs[ i ][ i + ( *getElasticConfigurationIndex( ) ) - 1 ] = 1;

                }

            }
 
            // Compute the damage strain
            floatVector Ed = ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) * Ee;

            floatVector dEddD = 1 / ( 1 - ( *get_damage( ) ) ) * ( 1 + ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) ) * Ee;

            floatMatrix dEddF = ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) * tardigradeVectorTools::dot( dEedFe, dFedF );
 
            floatMatrix dEddSubFs = ( *get_damage( ) ) / ( 1 - ( *get_damage( ) ) ) * tardigradeVectorTools::dot( dEedFe, dFedSubFs );
 
            // Compute the square root to solve for the damage deformation gradient
            floatVector eye( ( *dim ) * ( *dim ) );
            tardigradeVectorTools::eye( eye );
    
            floatVector Fd;

            floatMatrix dAdFe; //A = 2.0 * Ed + eye

            TARDIGRADE_ERROR_TOOLS_CATCH( Fd = tardigradeVectorTools::matrixSqrt( 2.0 * Ed + eye, *dim, dAdFe ) );

            floatMatrix dFddEd;

            TARDIGRADE_ERROR_TOOLS_CATCH( dFddEd = tardigradeVectorTools::inflate( 2 * tardigradeVectorTools::inverse( tardigradeVectorTools::appendVectors( dAdFe ),
                                                                                                                       ( *dim ) * ( *dim ), ( *dim ) * ( *dim ) ),
                                                                                   ( *dim ) * ( *dim ), ( * dim ) * ( *dim ) ) );

            floatVector dFddD = tardigradeVectorTools::dot( dFddEd, dEddD );

            floatMatrix dFddF = tardigradeVectorTools::dot( dFddEd, dEddF );

            floatMatrix dFddSubFs = tardigradeVectorTools::dot( dFddEd, dEddSubFs );

            if ( withPrevious ){

                set_dDamageDeformationGradientdPreviousCauchyStress( tardigradeVectorTools::dyadic( dFddD, *get_dDamagedPreviousCauchyStress( ) ) );
    
                set_dDamageDeformationGradientdPreviousF( tardigradeVectorTools::dyadic( dFddD, *get_dDamagedPreviousF( ) ) );
    
                set_dDamageDeformationGradientdPreviousSubFs( tardigradeVectorTools::dyadic( dFddD, *get_dDamagedPreviousSubFs( ) ) );
    
                set_dDamageDeformationGradientdPreviousT( dFddD * ( *get_dDamagedPreviousT( ) ) );
    
                set_dDamageDeformationGradientdPreviousStateVariables( tardigradeVectorTools::dyadic( dFddD, *get_dDamagedPreviousStateVariables( ) ) );

            }

            set_damageDeformationGradient( Fd );

            set_dDamageDeformationGradientdCauchyStress( tardigradeVectorTools::dyadic( dFddD, *get_dDamagedCauchyStress( ) ) );

            set_dDamageDeformationGradientdF( tardigradeVectorTools::dyadic( dFddD, *get_dDamagedF( ) ) + dFddF );

            set_dDamageDeformationGradientdSubFs( tardigradeVectorTools::dyadic( dFddD, *get_dDamagedSubFs( ) ) + dFddSubFs );

            set_dDamageDeformationGradientdT( dFddD * ( *get_dDamagedT( ) ) );

            set_dDamageDeformationGradientdStateVariables( tardigradeVectorTools::dyadic( dFddD, *get_dDamagedStateVariables( ) ) );

        }

        void residual::setStateVariableEvolutionRates( const bool isPrevious ){
            /*!
             * Set the state variable evolution rates including the damage
             * 
             * \param &isPrevious:Whether to compute the previous values or not
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
             * \param &isPrevious:Whether to compute the previous values or not
             */


            floatVector evolutionRates;

            floatMatrix tempJac( 2 );

            tardigradeHydra::peryznaViscoplasticity::residual::setStateVariableEvolutionRateDerivatives( isPrevious );

            if ( isPrevious ){

                evolutionRates = { ( *get_previousStateVariableEvolutionRates( ) )[ 0 ], *get_previousPlasticMultiplier( ) };

                set_previousStateVariableEvolutionRates( evolutionRates );

                // Set the derivatives w.r.t. the previous Cauchy stress
                tempJac[ 0 ] = ( *get_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( ) )[ 0 ];
                tempJac[ 1 ] = *get_dPreviousPlasticMultiplierdPreviousCauchyStress( );

                set_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( tempJac );

                // Set the derivatives w.r.t. the previous deformation gradient
                tempJac[ 0 ] = ( *get_dPreviousStateVariableEvolutionRatesdPreviousF( ) )[ 0 ];
                tempJac[ 1 ] = *get_dPreviousPlasticMultiplierdPreviousF( );

                set_dPreviousStateVariableEvolutionRatesdPreviousF( tempJac );

                // Set the derivatives w.r.t. the previous sub-deformation gradients
                tempJac[ 0 ] = ( *get_dPreviousStateVariableEvolutionRatesdPreviousSubFs( ) )[ 0 ];
                tempJac[ 1 ] = *get_dPreviousPlasticMultiplierdPreviousSubFs( );

                set_dPreviousStateVariableEvolutionRatesdPreviousSubFs( tempJac );

                // Set the derivatives w.r.t. the previous temperature
                floatVector tempJacVec = tardigradeVectorTools::appendVectors( { *get_dPreviousStateVariableEvolutionRatesdPreviousT( ),
                                                                               { *get_dPreviousPlasticMultiplierdPreviousT( ) } } );

                set_dPreviousStateVariableEvolutionRatesdPreviousT( tempJacVec );

                // Set the derivatives w.r.t. the previous state variables
                tempJac[ 0 ] = ( *get_dPreviousStateVariableEvolutionRatesdPreviousStateVariables( ) )[ 0 ];
                tempJac[ 1 ] = *get_dPreviousPlasticMultiplierdPreviousStateVariables( );

                set_dPreviousStateVariableEvolutionRatesdPreviousStateVariables( tempJac );

            }
            else{

                evolutionRates = { ( *get_stateVariableEvolutionRates( ) )[ 0 ], *get_plasticMultiplier( ) };

                set_stateVariableEvolutionRates( evolutionRates );

                // Set the derivatives w.r.t. the previous Cauchy stress
                tempJac[ 0 ] = ( *get_dStateVariableEvolutionRatesdCauchyStress( ) )[ 0 ];
                tempJac[ 1 ] = *get_dPlasticMultiplierdCauchyStress( );

                set_dStateVariableEvolutionRatesdCauchyStress( tempJac );

                // Set the derivatives w.r.t. the previous deformation gradient
                tempJac[ 0 ] = ( *get_dStateVariableEvolutionRatesdF( ) )[ 0 ];
                tempJac[ 1 ] = *get_dPlasticMultiplierdF( );

                set_dStateVariableEvolutionRatesdF( tempJac );

                // Set the derivatives w.r.t. the previous sub-deformation gradients
                tempJac[ 0 ] = ( *get_dStateVariableEvolutionRatesdSubFs( ) )[ 0 ];
                tempJac[ 1 ] = *get_dPlasticMultiplierdSubFs( );

                set_dStateVariableEvolutionRatesdSubFs( tempJac );

                // Set the derivatives w.r.t. the previous temperature
                floatVector tempJacVec = tardigradeVectorTools::appendVectors( { *get_dStateVariableEvolutionRatesdT( ),
                                                                               { *get_dPlasticMultiplierdT( ) } } );

                set_dStateVariableEvolutionRatesdT( tempJacVec );

                // Set the derivatives w.r.t. the previous state variables
                tempJac[ 0 ] = ( *get_dStateVariableEvolutionRatesdStateVariables( ) )[ 0 ];
                tempJac[ 1 ] = *get_dPlasticMultiplierdStateVariables( );

                set_dStateVariableEvolutionRatesdStateVariables( tempJac );

            }

        }

        void residual::setResidual( ){
            /*!
             * Set the residual vector
             */

            std::cout << "initializing \n";
            floatVector residual( *getNumEquations( ), 0 );

            for ( unsigned int i = 0; i < get_damageDeformationGradient( )->size( ); i++ ){

                residual[ i ] = hydra->getConfiguration( *getDamageConfigurationIndex( ) )[ i ] - ( *get_damageDeformationGradient( ) )[ i ];

            }

            std::cout << "residual 1:\n"; tardigradeVectorTools::print( residual );

            std::cout << "stateVariables:\n"; tardigradeVectorTools::print( *get_stateVariables( ) );
            std::cout << "plasticStateVariables:\n"; tardigradeVectorTools::print( *get_plasticStateVariables( ) );

            residual[ get_damageDeformationGradient( )->size( ) + 0 ] = ( *get_stateVariables( ) )[ 0 ] - ( *get_plasticStateVariables( ) ) [ 0 ]; //Evolution of the damage hardening state variable

            residual[ get_damageDeformationGradient( )->size( ) + 1 ] = ( *get_stateVariables( ) )[ 1 ] - ( *get_damage( ) ); //Evolution of the damage

            setResidual( residual );

        }

        void residual::setJacobian( ){
            /*!
             * Set the Jacobian matrix
             */

            floatMatrix jacobian( *getNumEquations( ), floatVector( hydra->getUnknownVector( )->size( ), 0 ) );

           // Jacobians of the damage deformation gradient
            for ( unsigned int i = 0; i < get_damageDeformationGradient( )->size( ); i++ ){
                unsigned int row = i;

                // Jacobians w.r.t. the previous Cauchy stress
                for ( unsigned int j = 0; j < getStress( )->size( ); j++ ){
                    unsigned int col = j;

                    jacobian[ row ][ col ] -= ( *get_dDamageDeformationGradientdCauchyStress( ) )[ i ][ j ];

                }

                // Jacobians w.r.t. the sub configurations
                jacobian[ row ][ hydra->getDeformationGradient( )->size( ) * ( *getDamageConfigurationIndex( ) ) + i ] += 1;

                for ( unsigned int j = 0; j < ( *get_dDamageDeformationGradientdSubFs( ) )[ i ].size( ); j++ ){
                    unsigned int col = ( *get_dDamageDeformationGradientdCauchyStress( ) )[ i ].size( ) + j;

                    jacobian[ row ][ col ] -= ( *get_dDamageDeformationGradientdSubFs( ) )[ i ][ j ];

                }

                // Jacobians w.r.t. the state variables
                for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){
                    unsigned int col = ( *get_dDamageDeformationGradientdCauchyStress( ) )[ i ].size( ) + ( *get_dDamageDeformationGradientdSubFs( ) )[ i ].size( ) + *ind;

                    jacobian[ row ][ col ] -= ( *get_dDamageDeformationGradientdStateVariables( ) )[ i ][ ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ];

                }

            }

            // Jacobians of the damage hardening state variable
            for ( unsigned int i = 0; i < 1; i++ ){
                unsigned int row = get_damageDeformationGradient( )->size( ) + i;

                // Jacobians w.r.t. the previous Cauchy stress
                for ( unsigned int j = 0; j < ( *get_dPlasticStateVariablesdCauchyStress( ) )[ i ].size( ); j++ ){
                    unsigned int col = j;

                    jacobian[ row ][ col ] -= ( *get_dPlasticStateVariablesdCauchyStress( ) )[ i ][ j ];

                }

                // Jacobians w.r.t. the sub configurations
                for ( unsigned int j = 0; j < ( *get_dPlasticStateVariablesdSubFs( ) )[ i ].size( ); j++ ){
                    unsigned int col = ( *get_dDamageDeformationGradientdCauchyStress( ) )[ i ].size( ) + j;

                    jacobian[ row ][ col ] -= ( *get_dPlasticStateVariablesdSubFs( ) )[ i ][ j ];

                }

                // Jacobians w.r.t. the state variables
                jacobian[ row ][ ( *get_dPlasticStateVariablesdCauchyStress( ) )[ i ].size( ) + ( *get_dPlasticStateVariablesdSubFs( ) )[ i ].size( ) + ( *getStateVariableIndices( ) )[ i ] ] += 1;
                for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){
                    unsigned int col = ( *get_dPlasticStateVariablesdCauchyStress( ) )[ i ].size( ) + ( *get_dPlasticStateVariablesdSubFs( ) )[ i ].size( ) + *ind;

                    jacobian[ row ][ col ] += ( *get_dPlasticStateVariablesdStateVariables( ) )[ i ][ ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ];

                }

            }

            //Jacobians of the damage
            for ( unsigned int i = 1; i < 2; i++ ){
                unsigned int row = get_damageDeformationGradient( )->size( ) + i;

                // Jacobians w.r.t. the previous Cauchy stress
                for ( unsigned int j = 0; j < ( *get_dDamagedCauchyStress( ) ).size( ); j++ ){
                    unsigned int col = j;

                    jacobian[ row ][ col ] -= ( *get_dDamagedCauchyStress( ) )[ j ];

                }

                // Jacobians w.r.t. the sub configurations
                for ( unsigned int j = 0; j < ( *get_dDamagedSubFs( ) ).size( ); j++ ){
                    unsigned int col = ( *get_dDamagedCauchyStress( ) ).size( ) + j;

                    jacobian[ row ][ col ] -= ( *get_dDamagedSubFs( ) )[ j ];

                }

                // Jacobians w.r.t. the state variables
                jacobian[ row ][ ( *get_dDamagedCauchyStress( ) ).size( ) + ( *get_dDamagedSubFs( ) ).size( ) + ( *getStateVariableIndices( ) )[ i ] ] += 1;
                for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){
                    unsigned int col = ( *get_dDamagedCauchyStress( ) ).size( ) + ( *get_dDamagedSubFs( ) ).size( ) + *ind;

                    jacobian[ row ][ col ] -= ( *get_dDamagedStateVariables( ) )[ ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ];

                }

            }

            setJacobian( jacobian );

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            floatVector dRdT( *getNumEquations( ), 0 );

            for ( unsigned int i = 0; i < get_damageDeformationGradient( )->size( ); i++ ){

                dRdT[ i ] -= ( *get_dDamageDeformationGradientdT( ) )[ i ];

            }

            dRdT[ get_damageDeformationGradient( )->size( ) ] = -( *get_dPlasticStateVariablesdT( ) )[ 0 ];

            dRdT[ get_damageDeformationGradient( )->size( ) + 1 ] = -( *get_dDamagedT( ) );

            setdRdT( dRdT );

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient
             */

            floatMatrix dRdF( *getNumEquations( ), floatVector( hydra->getDeformationGradient( )->size( ), 0 ) );

            for ( unsigned int i = 0; i < get_damageDeformationGradient( )->size( ); i++ ){

                for ( unsigned int j = 0; j < ( *get_dDamageDeformationGradientdF( ) )[ i ].size( ); j++ ){

                    dRdF[ i ][ j ] -= ( *get_dDamageDeformationGradientdF( ) )[ i ][ j ];

                }

            }

            dRdF[ get_damageDeformationGradient( )->size( ) ] = -( *get_dPlasticStateVariablesdF( ) )[ 0 ];

            dRdF[ get_damageDeformationGradient( )->size( ) + 1 ] = -( *get_dDamagedF( ) );

            setdRdF( dRdF );

        }

    }

}
