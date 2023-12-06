/**
  ******************************************************************************
  * \file tardigrade-hydraLinearElasticity.cpp
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
    
        void residual::setEe( ){
            /*!
             * Set the current value of the elastic Green-Lagrange strain
             */
    
            floatVector Fe = ( *hydra->getConfigurations( ) )[ 0 ];
    
            floatVector Ee;
   
            floatMatrix dEedFe;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeGreenLagrangeStrain( Fe, Ee, dEedFe ) );
    
            setEe( Ee );
    
            setdEedFe( dEedFe );
    
        }
    
        void residual::setEe( const floatVector &Ee ){
            /*!
             * Set the elastic Green-Lagrange strain
             * 
             * \param &Ee: The elastic Green-Lagrange strain
             */
    
            _Ee.second = Ee;
    
            _Ee.first = true;
    
            addIterationData( &_Ee );
    
        } 
    
        const floatVector* residual::getEe( ){
            /*!
             * Get the elastic Green-Lagrange strain
             */
    
            if ( !_Ee.first ){
    
                TARDIGRADE_ERROR_TOOLS_CATCH( setEe( ) );
    
            }
    
            return &_Ee.second;
    
        }
    
        void residual::setdEedFe( ){
            /*!
             * Set the gradient of the elastic Green-Lagrange strain w.r.t. the elastic deformation gradient
             * 
             * Default assumption is that this happens when the Green-Lagrange strain is computed; 
             */
    
            TARDIGRADE_ERROR_TOOLS_CATCH( getEe( ) );
    
        }
    
        void residual::setdEedFe( const floatMatrix &dEedFe ){
            /*!
             * Set the gradient of the elastic Green-Lagrange strain w.r.t. the elastic deformation gradient
             * 
             * \param &dEedFe: The derivative of the elastic Green-Lagrange strain w.r.t. the elastic deformation gradient
             */
    
            _dEedFe.second = dEedFe;
    
            _dEedFe.first = true;
    
            addIterationData( &_dEedFe );
    
        }
    
        const floatMatrix* residual::getdEedFe( ){
            /*!
             * Get the elastic Green-Lagrange strain
             */
    
            if ( !_dEedFe.first ){
    
                TARDIGRADE_ERROR_TOOLS_CATCH( setdEedFe( ) );
    
            }
    
            return &_dEedFe.second;
    
        }
    
        void residual::setPK2Stress( ){
            /*!
             * Compute the Second Piola-Kirchhoff stress
             */
    
            floatVector eye( getEe( )->size( ), 0 );
            tardigradeVectorTools::eye( eye );
    
            floatVector PK2Stress = ( *getLambda( ) ) * tardigradeVectorTools::trace( *getEe( ) ) * eye + 2 * ( *getMu( ) ) * ( *getEe( ) );
    
            setPK2Stress( PK2Stress );
    
        }
    
        void residual::setPK2Stress( const floatVector &PK2Stress ){
            /*!
             * Compute the Second Piola-Kirchhoff stress
             * 
             * \param &PK2Stress: The Second Piola-Kirchhoff stress
             */
    
            _PK2Stress.second = PK2Stress;
    
            _PK2Stress.first = true;
    
            addIterationData( &_PK2Stress );
    
        }
    
        const floatVector* residual::getPK2Stress( ){
            /*!
             * Get the Second Piola-Kirchhoff stress
             */
    
            if ( !_PK2Stress.first ){
    
                TARDIGRADE_ERROR_TOOLS_CATCH( setPK2Stress( ) );
    
            }
    
            return &_PK2Stress.second;
    
        }
    
        void residual::setdPK2StressdEe( ){
            /*!
             * Compute the gradient of the PK2 stress w.r.t. the elastic Green-Lagrange strain
             */
    
            floatVector eye( getEe( )->size( ), 0 );
            tardigradeVectorTools::eye( eye );
    
            floatMatrix EYE = tardigradeVectorTools::eye< floatType >( getEe( )->size( ) );
    
            floatMatrix dPK2StressdEe = ( *getLambda( ) ) * tardigradeVectorTools::dyadic( eye, eye ) + 2 * ( *getMu( ) ) * EYE;
    
            setdPK2StressdEe( dPK2StressdEe );
    
        }
    
        void residual::setdPK2StressdEe( const floatMatrix &dPK2StressdEe ){
            /*!
             * Set the gradient of the PK2 stress w.r.t. the elastic Green-Lagrange strain
             *
             * \param &dPK2StressdEe: The gradient of the Second Piola-Kirchhoff stress w.r.t. the elastic Green-Lagrange strain
             */
    
            _dPK2StressdEe.second = dPK2StressdEe;
    
            _dPK2StressdEe.first = true;
    
            addIterationData( &_dPK2StressdEe );
    
        }
    
        const floatMatrix* residual::getdPK2StressdEe( ){
            /*!
             * Get the gradient of the Second Piola-Kirchhoff stress w.r.t. the elastic Green-Lagrange strain
             */
    
            if ( !_dPK2StressdEe.first ){
    
                TARDIGRADE_ERROR_TOOLS_CATCH( setdPK2StressdEe( ) );
    
            }
    
            return &_dPK2StressdEe.second;
    
        }

        void residual::setdPK2StressdFe( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the elastic
             * deformation gradient
             */

             floatMatrix dPK2StressdFe = tardigradeVectorTools::dot( *getdPK2StressdEe( ), *getdEedFe( ) );

             setdPK2StressdFe( dPK2StressdFe );

        }

        void residual::setdPK2StressdFe( const floatMatrix &dPK2StressdFe ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the elastic
             * deformation gradient
             * 
             * \param &dPK2StressdFe: The derivative of the second Piola-Kirchhoff stress w.r.t.
             *     the elastic deformation gradient
             */

            _dPK2StressdFe.second = dPK2StressdFe;

            _dPK2StressdFe.first = true;

            addIterationData( &_dPK2StressdFe );

        }

        const floatMatrix* residual::getdPK2StressdFe( ){
            /*!
             * Get the derivative of the second Piola-Kirchhoff stress w.r.t. the elastic
             * deformation gradient
             */

            if ( !_dPK2StressdFe.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPK2StressdFe( ) )

            }

            return &_dPK2StressdFe.second;

        }

        void residual::setCauchyStress( ){
            /*!
             * Set the Cauchy stress
             */
    
            floatVector Fe = ( *hydra->getConfigurations( ) )[ 0 ];
    
            floatMatrix dFedF = ( *hydra->getdF1dF( ) );
    
            floatMatrix dFedFn = ( *hydra->getdF1dFn( ) );
    
            // Compute the gradient of the PK2 stress w.r.t. the elastic deformation gradient
            floatMatrix dPK2StressdFe = *getdPK2StressdFe( );

            // Compute the Second Piola-Kirchhoff stress and it's gradients
            floatMatrix dPK2StressdF = tardigradeVectorTools::dot( *getdPK2StressdFe( ), dFedF );
    
            floatMatrix dPK2StressdFn = tardigradeVectorTools::dot( *getdPK2StressdFe( ), dFedFn );
    
            // Map the PK2 stress to the current configuration
            floatVector cauchyStress;
            floatMatrix dCauchyStressdPK2Stress;
            floatMatrix dCauchyStressdFe;
    
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::pushForwardPK2Stress( *getPK2Stress( ), Fe, cauchyStress, dCauchyStressdPK2Stress, dCauchyStressdFe ) );
    
            setCauchyStress( cauchyStress );
    
            floatMatrix dCauchyStressdF  = tardigradeVectorTools::dot( dCauchyStressdPK2Stress, dPK2StressdF )
                                         + tardigradeVectorTools::dot( dCauchyStressdFe, dFedF );
    
            floatMatrix dCauchyStressdFn = tardigradeVectorTools::dot( dCauchyStressdPK2Stress, dPK2StressdFn )
                                         + tardigradeVectorTools::dot( dCauchyStressdFe, dFedFn );

            setdCauchyStressdPK2Stress( dCauchyStressdPK2Stress );   

            setdCauchyStressdF( dCauchyStressdF );
    
            setdCauchyStressdFn( dCauchyStressdFn );
    
        }
    
        void residual::setdCauchyStressdPK2Stress( ){
            /*!
             * Set the derivative of the computed Cauchy stress w.r.t. the second Piola-Kirchoff stress (this is a partial derivative generally)
             */
    
            setCauchyStress( );
    
        }

        void residual::setdCauchyStressdF( ){
            /*!
             * Set the derivative of the computed Cauchy stress w.r.t. F (this is a partial derivative generally)
             */
    
            setCauchyStress( );
    
        }

        void residual::setdCauchyStressdPK2Stress( const floatMatrix &dCauchyStressdPK2Stress ){
            /*!
             * Set the partial derivative of the Cauchy stress w.r.t. the second Piola-Kirchhoff stress
             * 
             * \param &dCauchyStressdPK2Stress: The partial derivative of the Cauchy stress w.r.t. the second Piola-Kirchhoff stress
             */
    
            _dCauchyStressdPK2Stress.second = dCauchyStressdPK2Stress;
    
            _dCauchyStressdPK2Stress.first = true;
    
            addIterationData( &_dCauchyStressdPK2Stress );
    
        }

        void residual::setdCauchyStressdF( const floatMatrix &dCauchyStressdF ){
            /*!
             * Set the partial derivative of the Cauchy stress w.r.t. the deformation gradient
             * 
             * \param &dCauchyStressdF: The partial derivative of the Cauchy stress w.r.t. the deformation gradient
             */
    
            _dCauchyStressdF.second = dCauchyStressdF;
    
            _dCauchyStressdF.first = true;
    
            addIterationData( &_dCauchyStressdF );
    
        }

        const floatMatrix* residual::getdCauchyStressdPK2Stress( ){
            /*!
             * Get the derivative of the Cauchy stress w.r.t. the second Piola-Kirchhoff stress
             */

            if ( !_dCauchyStressdPK2Stress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdCauchyStressdPK2Stress( ) );

            }

            return &_dCauchyStressdPK2Stress.second;

        }

        const floatMatrix* residual::getdCauchyStressdF( ){
            /*!
             * Get the derivative of the Cauchy stress w.r.t. the deformation gradient
             */

            if ( !_dCauchyStressdF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdCauchyStressdF( ) );

            }

            return &_dCauchyStressdF.second;

        }

        void residual::setdCauchyStressdFn( ){
            /*!
             * Set the derivative of the computed Cauchy stress w.r.t. the configurations solved for in the non-linear solve
             * (this is a partial derivative generally)
             */
    
            setCauchyStress( );
    
        }
    
        void residual::setdCauchyStressdFn( const floatMatrix &dCauchyStressdFn ){
            /*!
             * Set the partial derivative of the Cauchy stress w.r.t. the configurations solve for in the non-linear solve
             * 
             * \param &dCauchyStressdFn: The partial derivative of the Cauchy stress w.r.t. the configurations solve for
             *     in the non-linear solve
             */
    
            _dCauchyStressdFn.second = dCauchyStressdFn;
    
            _dCauchyStressdFn.first = true;
    
            addIterationData( &_dCauchyStressdFn );
    
        }

        const floatMatrix* residual::getdCauchyStressdFn( ){
            /*!
             * Get the derivative of the Cauchy stress w.r.t. the configurations solve for in the non-linear solve
             */

            if ( !_dCauchyStressdFn.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdCauchyStressdFn( ) );

            }

            return &_dCauchyStressdFn.second;

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
    
                        jacobian[ ( *dim ) * i + j ][ getStress( )->size( ) + I ] = ( *getdCauchyStressdFn( ) )[ ( *dim ) * i + j ][ I ];
    
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
    
            setdRdF( *getdCauchyStressdF( ) );
    
        }

    }

}
