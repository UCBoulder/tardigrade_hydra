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
             * Set the value of the damage
             */

            const floatType   *damageMultiplier;

            const floatType   *previousDamageMultiplier;

            TARDIGRADE_ERROR_TOOLS_CATCH( damageMultiplier = get_plasticMultiplier( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousDamageMultiplier = get_previousPlasticMultiplier( ) );

            floatVector previousDamage( get_previousStateVariables( )->begin( ) + 1,
                                        get_previousStateVariables( )->begin( ) + 2 ); //The first state variable is the hardening variable and the second is the damage

            floatVector deltaDamage;

            floatVector updatedDamage;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolution( *hydra->getDeltaTime( ),       previousDamage,
                                                                                                       { *previousDamageMultiplier }, { *damageMultiplier },
                                                                                                       deltaDamage, updatedDamage, 1 - ( *getIntegrationParameter( ) ) ) );

            set_damage( updatedDamage[ 0 ] );

        }

        void residual::setDamageDerivatives( const bool withPrevious){
            /*!
             * Set the damage along with the derivatives of the damage.
             * 
             * \param withPrevious: Flag for whether to include the derivatives w.r.t. the previous values.
             */

            const floatType   *damageMultiplier;

            const floatType   *previousDamageMultiplier;

            TARDIGRADE_ERROR_TOOLS_CATCH( damageMultiplier = get_plasticMultiplier( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousDamageMultiplier = get_previousPlasticMultiplier( ) );

            floatVector previousDamage( get_previousStateVariables( )->begin( ) + 1,
                                        get_previousStateVariables( )->begin( ) + 2 ); //The first state variable is the hardening variable and the second is the damage

            floatVector deltaDamage;

            floatVector updatedDamage;

            floatMatrix dDdDdot;

            if ( withPrevious ){

                floatMatrix dDdDdotp;

                TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolution( *hydra->getDeltaTime( ),       previousDamage,
                                                                                                           { *previousDamageMultiplier }, { *damageMultiplier },
                                                                                                           deltaDamage, updatedDamage, dDdDdot, dDdDdotp, 1 - ( *getIntegrationParameter( ) ) ) );

//                set_dDamagedPreviousCauchyStress( dDdDdotp[ 0 ][ 0 ] * ( *get_previousdPlasticMultiplierdCauchyStress( ) ) );
//    
//                set_dDamagedPreviousF( dDdDdotp[ 0 ][ 0 ] * ( *get_previousdPlasticMultiplierdF( ) ) );
//    
//                set_dDamagedPreviousSubFs( dDdDdotp[ 0 ][ 0 ] * ( *get_previousdPlasticMultiplierdSubFs( ) ) );
//    
//                set_dDamagedPreviousT( dDdDdotp[ 0 ][ 0 ] * ( *get_previousdPlasticMultiplierdT( ) ) );
//    
//                set_dDamagedPreviousStateVariables( dDdDdotp[ 0 ][ 0 ] * ( *get_previousdPlasticMultiplierdStateVariables( ) ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolution( *hydra->getDeltaTime( ),       previousDamage,
                                                                                                           { *previousDamageMultiplier }, { *damageMultiplier },
                                                                                                           deltaDamage, updatedDamage, dDdDdot, 1 - ( *getIntegrationParameter( ) ) ) );

            }

            set_damage( updatedDamage[ 0 ] );

            set_dDamagedCauchyStress( dDdDdot[ 0 ][ 0 ] * ( *get_dPlasticMultiplierdCauchyStress( ) ) );

            set_dDamagedF( dDdDdot[ 0 ][ 0 ] * ( *get_dPlasticMultiplierdF( ) ) );

            set_dDamagedSubFs( dDdDdot[ 0 ][ 0 ] * ( *get_dPlasticMultiplierdSubFs( ) ) );

            set_dDamagedT( dDdDdot[ 0 ][ 0 ] * ( *get_dPlasticMultiplierdT( ) ) );

            set_dDamagedStateVariables( dDdDdot[ 0 ][ 0 ] * ( *get_dPlasticMultiplierdStateVariables( ) ) );

        }

        void residual::setDamageDerivatives( ){
            /*!
             * Set the damage along with the derivatives of the damage.
             * This routine will not set the derivatives w.r.t. the previous values.
             */

            setDamageDerivatives( false );

        }

        void residual::setAllDamageDerivatives( ){
            /*!
             * Set the damage along with the derivatives of the damage.
             * This routine will set the derivatives w.r.t. the previous values.
             */

            setDamageDerivatives( true );

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

    }

}
