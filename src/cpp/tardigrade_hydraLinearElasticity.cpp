/**
  ******************************************************************************
  * \file tardigrade_hydraLinearElasticity.cpp
  ******************************************************************************
  * An implementation of linear elasticity using the hydra framework. Used as an
  * example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraLinearElasticity.h>
#include<tardigrade_constitutive_tools.h>

namespace tardigradeHydra{

    namespace linearElasticity{

        void residual::decomposeParameterVector( const floatVector &parameters ){
            /*!
             * Decompose the parameter vector
             * 
             * \param &parameters: The paramter vector. Assumed to be a vector of length 2 which defines lambda and mu.
             */
    
            if ( parameters.size( ) != 2 ){
    
                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "Parameter vector is expected to have a length of 2 but has a length of " + std::to_string( parameters.size( ) ) ) );
    
            }
    
            setLambda( parameters[ 0 ] );
    
            setMu( parameters[ 1 ] );
    
        }
   
        void residual::setEe( const bool isPrevious ){
            /*!
             * Set the value of the elastic Green-Lagrange strain
             * 
             * \param isPrevious: Flag for whether to compute the strain (false) or the previous strain (true)
             */

            floatVector Fe;
    
            floatMatrix dEedFe;

            if ( isPrevious ){

                Fe = ( *hydra->get_previousConfigurations( ) )[ 0 ];

            }
            else{

                Fe = ( *hydra->get_configurations( ) )[ 0 ];

            }

            floatVector Ee;
   
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeGreenLagrangeStrain( Fe, Ee, dEedFe ) );

            if ( isPrevious ){
    
                set_previousEe( Ee );
    
                set_previousdEedFe( dEedFe );
            }
            else{

                set_Ee( Ee );
    
                set_dEedFe( dEedFe );

            }

        }
 
        void residual::setEe( ){
            /*!
             * Set the current value of the elastic Green-Lagrange strain
             */

            TARDIGRADE_ERROR_TOOLS_CATCH( setEe( false ) )
    
        }
    
        void residual::setdEedFe( ){
            /*!
             * Set the gradient of the elastic Green-Lagrange strain w.r.t. the elastic deformation gradient
             * 
             * Default assumption is that this happens when the Green-Lagrange strain is computed; 
             */
    
            TARDIGRADE_ERROR_TOOLS_CATCH( setEe( false ) );
    
        }
    
        void residual::setPreviousEe( ){
            /*!
             * Set the previous value of the elastic Green-Lagrange strain
             */

            TARDIGRADE_ERROR_TOOLS_CATCH( setEe( true ) )
    
        }
    
        void residual::setPreviousdEedFe( ){
            /*!
             * Set the previous gradient of the elastic Green-Lagrange strain w.r.t. the elastic deformation gradient
             * 
             * Default assumption is that this happens when the Green-Lagrange strain is computed; 
             */
    
            TARDIGRADE_ERROR_TOOLS_CATCH( setEe( true ) );
    
        }
    
        void residual::setPK2Stress( ){
            /*!
             * Compute the Second Piola-Kirchhoff stress
             */
    
            floatVector eye( get_Ee( )->size( ), 0 );
            tardigradeVectorTools::eye( eye );
    
            floatVector PK2Stress = ( *getLambda( ) ) * tardigradeVectorTools::trace( *get_Ee( ) ) * eye + 2 * ( *getMu( ) ) * ( *get_Ee( ) );
    
            set_PK2Stress( PK2Stress );
    
        }
    
        void residual::setdPK2StressdEe( ){
            /*!
             * Compute the gradient of the PK2 stress w.r.t. the elastic Green-Lagrange strain
             */
    
            floatVector eye( get_Ee( )->size( ), 0 );
            tardigradeVectorTools::eye( eye );
    
            floatMatrix EYE = tardigradeVectorTools::eye< floatType >( get_Ee( )->size( ) );
    
            floatMatrix dPK2StressdEe = ( *getLambda( ) ) * tardigradeVectorTools::dyadic( eye, eye ) + 2 * ( *getMu( ) ) * EYE;
    
            set_dPK2StressdEe( dPK2StressdEe );
    
        }
    
        void residual::setdPK2StressdFe( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the elastic
             * deformation gradient
             */

             floatMatrix dPK2StressdFe = tardigradeVectorTools::dot( *get_dPK2StressdEe( ), *get_dEedFe( ) );

             set_dPK2StressdFe( dPK2StressdFe );

        }

        void residual::setStress( ){
            /*!
             * Set the Cauchy stress
             */
    
            floatVector Fe = ( *hydra->get_configurations( ) )[ 0 ];
    
            floatMatrix dFedF = ( *hydra->get_dF1dF( ) );
    
            floatMatrix dFedFn = ( *hydra->get_dF1dFn( ) );
    
            // Compute the gradient of the PK2 stress w.r.t. the elastic deformation gradient
            floatMatrix dPK2StressdFe = *get_dPK2StressdFe( );

            // Compute the Second Piola-Kirchhoff stress and it's gradients
            floatMatrix dPK2StressdF = tardigradeVectorTools::dot( *get_dPK2StressdFe( ), dFedF );
    
            floatMatrix dPK2StressdFn = tardigradeVectorTools::dot( *get_dPK2StressdFe( ), dFedFn );
    
            // Map the PK2 stress to the current configuration
            floatVector cauchyStress;
            floatMatrix dCauchyStressdPK2Stress;
            floatMatrix dCauchyStressdFe;
    
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::pushForwardPK2Stress( *get_PK2Stress( ), Fe, cauchyStress, dCauchyStressdPK2Stress, dCauchyStressdFe ) );
    
            setStress( cauchyStress );
    
            floatMatrix dCauchyStressdF  = tardigradeVectorTools::dot( dCauchyStressdPK2Stress, dPK2StressdF )
                                         + tardigradeVectorTools::dot( dCauchyStressdFe, dFedF );
    
            floatMatrix dCauchyStressdFn = tardigradeVectorTools::dot( dCauchyStressdPK2Stress, dPK2StressdFn )
                                         + tardigradeVectorTools::dot( dCauchyStressdFe, dFedFn );

            set_dCauchyStressdPK2Stress( dCauchyStressdPK2Stress );   

            set_dCauchyStressdF( dCauchyStressdF );
    
            set_dCauchyStressdFn( dCauchyStressdFn );
    
        }
    
        void residual::setdCauchyStressdPK2Stress( ){
            /*!
             * Set the derivative of the computed Cauchy stress w.r.t. the second Piola-Kirchoff stress (this is a partial derivative generally)
             */
    
            setStress( );
    
        }

        void residual::setdCauchyStressdF( ){
            /*!
             * Set the derivative of the computed Cauchy stress w.r.t. F (this is a partial derivative generally)
             */
    
            setStress( );
    
        }

        void residual::setdCauchyStressdFn( ){
            /*!
             * Set the derivative of the computed Cauchy stress w.r.t. the configurations solved for in the non-linear solve
             * (this is a partial derivative generally)
             */
    
            setStress( );
    
        }
    
        void residual::setResidual( ){
            /*!
             * Set the residual value
             */
    
            const floatVector *cauchyStress = getStress( );

            TARDIGRADE_ERROR_TOOLS_CATCH( setResidual( *cauchyStress - *hydra->getStress( ) ) );
    
        }
    
        void residual::setJacobian( ){
            /*!
             * Set the Jacobian value
             */
    
            const unsigned int *dim = hydra->getDimension( );
    
            // Form the Jacobian
            floatMatrix jacobian = floatMatrix( getStress( )->size( ), floatVector( hydra->getUnknownVector( )->size( ), 0 ) );
    
            for ( unsigned int i = 0; i < ( *dim ); i++ ){
    
                for ( unsigned int j = 0; j < ( *dim ); j++ ){
    
                    jacobian[ ( *dim ) * i + j ][ ( *dim ) * i + j ] = -1;
    
                    for ( unsigned int I = 0; I < ( ( *hydra->getNumConfigurations( ) ) - 1 ) * ( *dim ) * ( *dim ); I++ ){
    
                        jacobian[ ( *dim ) * i + j ][ getStress( )->size( ) + I ] = ( *get_dCauchyStressdFn( ) )[ ( *dim ) * i + j ][ I ];
    
                    }
    
                }
    
            }
    
            setJacobian( jacobian );
    
        }
    
        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */
    
            setdRdT( floatVector( *getNumEquations( ), 0 ) );
    
        }
    
        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient
             */
    
            setdRdF( *get_dCauchyStressdF( ) );
    
        }

    }

}
