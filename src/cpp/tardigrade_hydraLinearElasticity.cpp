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

        void residual::setFe( const bool isPrevious ){
            /*!
             * Set the value of the elastic deformation gradient
             * 
             * \param isPrevious: Flag for whether to set the current (false) or previous (true) elastic deformation gradient
             */

            const unsigned int dim = 3;
            const unsigned int sot_dim = dim * dim;

            if ( isPrevious ){

                set_previousFe( floatVector( hydra->get_previousConfigurations( )->begin( ),
                                             hydra->get_previousConfigurations( )->begin( ) + sot_dim ) );

            }
            else{

                set_Fe( floatVector( hydra->get_configurations( )->begin( ),
                                     hydra->get_configurations( )->begin( ) + sot_dim ) );

            }

        }

        void residual::setFe( ){
            /*!
             * Set the value of the elastic deformation gradient
             */

            setFe( false );

        }

        void residual::setPreviousFe( ){
            /*!
             * Set the value of the previous elastic deformation gradient
             */

            setFe( true );

        }

        void residual::setFeDerivatives( const bool isPrevious ){
            /*!
             * Set the value of the derivatives of the elastic strain
             */

            if ( isPrevious ){

                set_previousdFedF( *hydra->get_previousdF1dF( ) );

                set_previousdFedFn( *hydra->get_previousdF1dFn( ) );

            }
            else{

                set_dFedF(  *hydra->get_dF1dF( ) );

                set_dFedFn( *hydra->get_dF1dFn( ) );

            }

        }

        void residual::setdFedF( ){
            /*!
             * Set the value of the derivative of the elastic deformation gradient w.r.t. the deformation gradient
             */

            setFeDerivatives( false );

        }

        void residual::setdFedFn( ){
            /*!
             * Set the value of the derivative of the elastic deformation gradient w.r.t. the sub-deformation gradients
             */

            setFeDerivatives( false );

        }

        void residual::setPreviousdFedF( ){
            /*!
             * Set the value of the previous derivative of the elastic deformation gradient w.r.t. the deformation gradient
             */

            setFeDerivatives( true );

        }

        void residual::setPreviousdFedFn( ){
            /*!
             * Set the value of the previous derivative of the elastic deformation gradient w.r.t. the sub-deformation gradients
             */

            setFeDerivatives( true );

        }

        void residual::setEe( const bool isPrevious ){
            /*!
             * Set the value of the elastic Green-Lagrange strain
             * 
             * \param isPrevious: Flag for whether to compute the strain (false) or the previous strain (true)
             */

            const floatVector *Fe;
    
            if ( isPrevious ){

                Fe = get_previousFe( );

            }
            else{

                Fe = get_Fe( );

            }

            floatVector Ee;
   
            floatVector dEedFe;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeGreenLagrangeStrain( *Fe, Ee, dEedFe ) );

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
    
        void residual::setPK2Stress( const bool isPrevious ){
            /*!
             * Compute the Second Piola-Kirchhoff stress
             * 
             * \param isPrevious: Flag for whether to compute the current (false) or previous (true) PK2 stress
             */
    
            floatVector eye( get_Ee( )->size( ), 0 );
            tardigradeVectorTools::eye( eye );

            const floatVector *Ee;

            if ( isPrevious ){

                Ee = get_previousEe( );

            }
            else{

                Ee = get_Ee( );

            }
            floatVector PK2Stress = ( *getLambda( ) ) * tardigradeVectorTools::trace( *Ee ) * eye + 2 * ( *getMu( ) ) * ( *Ee );

            if ( isPrevious ){

                set_previousPK2Stress( PK2Stress );

            }
            else{
 
                set_PK2Stress( PK2Stress );

            }
    
        }
    
        void residual::setPK2Stress( ){
            /*!
             * Compute the Second Piola-Kirchhoff stress
             */

            setPK2Stress( false );
    
        }

        void residual::setPreviousPK2Stress( ){
            /*!
             * Compute the previous Second Piola-Kirchhoff stress
             */

            setPK2Stress( true );
    
        }
    
        void residual::setdPK2StressdEe( const bool isPrevious ){
            /*!
             * Compute the gradient of the PK2 stress w.r.t. the elastic Green-Lagrange strain
             * 
             * \param isPrevious: Flag for whether to compute the current (false) or previous (true) value
             */
   
            const unsigned int dim = *hydra->getDimension( );
            const unsigned int sot_dim = dim * dim;
            const unsigned int fot_dim = sot_dim * sot_dim;
 
            floatVector eye( sot_dim, 0 );
            tardigradeVectorTools::eye( eye );
    
            floatVector EYE( fot_dim, 0 );
            tardigradeVectorTools::eye( EYE );
    
            floatVector dPK2StressdEe = ( *getLambda( ) ) * tardigradeVectorTools::matrixMultiply( eye, eye, sot_dim, 1, 1, sot_dim ) + 2 * ( *getMu( ) ) * EYE;
   
            set_dPK2StressdEe( dPK2StressdEe );

            set_previousdPK2StressdEe( dPK2StressdEe );
    
        }
    
        void residual::setdPK2StressdEe( ){
            /*!
             * Compute the gradient of the PK2 stress w.r.t. the elastic Green-Lagrange strain
             */

            setdPK2StressdEe( false );
    
        }
    
        void residual::setPreviousdPK2StressdEe( ){
            /*!
             * Compute the previous gradient of the PK2 stress w.r.t. the elastic Green-Lagrange strain
             */

            setdPK2StressdEe( true );
    
        }
    
        void residual::setdPK2StressdFe( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the elastic
             * deformation gradient
             */

             unsigned int sot_dim = ( *hydra->getDimension( ) ) * ( *hydra->getDimension( ) );

             floatVector dPK2StressdFe = tardigradeVectorTools::matrixMultiply( *get_dPK2StressdEe( ), *get_dEedFe( ), sot_dim, sot_dim, sot_dim, sot_dim );

             set_dPK2StressdFe( dPK2StressdFe );

        }

        void residual::setdPK2StressdPreviousFe( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the previous elastic
             * deformation gradient
             */

             floatVector dPK2StressdPreviousFe( get_PK2Stress( )->size( ) * get_previousEe( )->size( ), 0 );

             set_dPK2StressdPreviousFe( dPK2StressdPreviousFe );

        }

        void residual::setPreviousdPK2StressdFe( ){
            /*!
             * Set the derivative of the previous second Piola-Kirchhoff stress w.r.t. the elastic
             * deformation gradient
             */

             unsigned int sot_dim = ( *hydra->getDimension( ) ) * ( *hydra->getDimension( ) );

             floatVector previousdPK2StressdFe = tardigradeVectorTools::matrixMultiply( *get_previousdPK2StressdEe( ), *get_previousdEedFe( ), sot_dim, sot_dim, sot_dim, sot_dim );

             set_previousdPK2StressdFe( previousdPK2StressdFe );

        }

        void residual::setCauchyStress( const bool isPrevious ){
            /*!
             * Set the Cauchy stress
             * 
             * \param isPrevious: Whether to compute the current (false) or previous (true) Cauchy stress
             */

            unsigned int dim = *hydra->getDimension( );

            unsigned int sot_dim = dim * dim;

            unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector *Fe;

            const floatVector *dFedF;

            const floatVector *dFedFn;

            const floatVector *PK2Stress;

            const floatVector *dPK2StressdFe;

            if ( isPrevious ){

                Fe            = get_previousFe( );

                dFedF         = get_previousdFedF( );

                dFedFn        = get_previousdFedFn( );

                PK2Stress     = get_previousPK2Stress( );

                dPK2StressdFe = get_previousdPK2StressdFe( );

            }
            else{

                Fe            = get_Fe( );

                dFedF         = get_dFedF( );

                dFedFn        = get_dFedFn( );

                PK2Stress     = get_PK2Stress( );

                dPK2StressdFe = get_dPK2StressdFe( );

            }

            // Compute the Second Piola-Kirchhoff stress and it's gradients
            floatVector dPK2StressdF = tardigradeVectorTools::matrixMultiply( *dPK2StressdFe, *dFedF, sot_dim, sot_dim, sot_dim, sot_dim );
    
            floatVector dPK2StressdFn = tardigradeVectorTools::matrixMultiply( *dPK2StressdFe, *dFedFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );
    
            // Map the PK2 stress to the current configuration
            floatVector cauchyStress;
            floatVector dCauchyStressdPK2Stress;
            floatVector dCauchyStressdFe;
 
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::pushForwardPK2Stress( *PK2Stress, *Fe, cauchyStress, dCauchyStressdPK2Stress, dCauchyStressdFe ) );
    
            floatVector dCauchyStressdF  = tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, dPK2StressdF, sot_dim, sot_dim, sot_dim, sot_dim )
                                         + tardigradeVectorTools::matrixMultiply( dCauchyStressdFe, *dFedF, sot_dim, sot_dim, sot_dim, sot_dim );
    
            floatVector dCauchyStressdFn = tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, dPK2StressdFn, sot_dim, sot_dim, sot_dim, sot_dim * ( num_configs - 1 ) )
                                         + tardigradeVectorTools::matrixMultiply( dCauchyStressdFe, *dFedFn, sot_dim, sot_dim, sot_dim, sot_dim * ( num_configs - 1 ) );

            if ( isPrevious ){

                set_previousCauchyStress( cauchyStress );
    
                set_previousdCauchyStressdPK2Stress( dCauchyStressdPK2Stress );   

                set_previousdCauchyStressdF( dCauchyStressdF );
    
                set_previousdCauchyStressdFn( dCauchyStressdFn );

            }
            else{

                set_cauchyStress( cauchyStress );
    
                set_dCauchyStressdPK2Stress( dCauchyStressdPK2Stress );   

                set_dCauchyStressdF( dCauchyStressdF );
    
                set_dCauchyStressdFn( dCauchyStressdFn );

                floatVector dPK2StressdPreviousF  = tardigradeVectorTools::matrixMultiply( *get_dPK2StressdPreviousFe( ), *get_previousdFedF( ), sot_dim, sot_dim, sot_dim, sot_dim );

                floatVector dPK2StressdPreviousFn = tardigradeVectorTools::matrixMultiply( *get_dPK2StressdPreviousFe( ), *get_previousdFedFn( ), sot_dim, sot_dim, sot_dim, sot_dim * ( num_configs - 1 ) );

                set_dCauchyStressdPreviousF( tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, dPK2StressdPreviousF, sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dCauchyStressdPreviousFn( tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, dPK2StressdPreviousFn, sot_dim, sot_dim, sot_dim, sot_dim * ( num_configs - 1 ) ) );

            }

        }

        void residual::setCauchyStress( ){
            /*!
             * Set the value of the Cauchy stress
             */

            setCauchyStress( false );

        }

        void residual::setPreviousCauchyStress( ){
            /*!
             * Set the previous value of the Cauchy stress
             */

            setCauchyStress( true );

        }

        void residual::setStress( ){
            /*!
             * Set the stress
             *
             * Currently uses the Cauchy stress
             */

            setCauchyStress( false );

            setStress( *get_cauchyStress( ) );

        }
    
        void residual::setPreviousStress( ){
            /*!
             * Set the previous stress
             *
             * Currently uses the Cauchy stress
             */

            setCauchyStress( true );

            setPreviousStress( *get_previousCauchyStress( ) );

        }
    
        void residual::setdCauchyStressdPK2Stress( ){
            /*!
             * Set the derivative of the computed Cauchy stress w.r.t. the second Piola-Kirchoff stress (this is a partial derivative generally)
             */
    
            setCauchyStress( false );
    
        }

        void residual::setdCauchyStressdF( ){
            /*!
             * Set the derivative of the computed Cauchy stress w.r.t. F (this is a partial derivative generally)
             */
    
            setCauchyStress( false );
    
        }

        void residual::setdCauchyStressdFn( ){
            /*!
             * Set the derivative of the computed Cauchy stress w.r.t. the configurations solved for in the non-linear solve
             * (this is a partial derivative generally)
             */
    
            setCauchyStress( false );
    
        }
    
        void residual::setdCauchyStressdPreviousF( ){
            /*!
             * Set the derivative of the computed Cauchy stress w.r.t. the previous value of F (this is a partial derivative generally)
             */
    
            setCauchyStress( false );
    
        }

        void residual::setdCauchyStressdPreviousFn( ){
            /*!
             * Set the derivative of the computed Cauchy stress w.r.t. the previous configurations solved for in the non-linear solve
             * (this is a partial derivative generally)
             */
    
            setCauchyStress( false );
    
        }
    
        void residual::setPreviousdCauchyStressdPK2Stress( ){
            /*!
             * Set the previous derivative of the computed Cauchy stress w.r.t. the second Piola-Kirchoff stress (this is a partial derivative generally)
             */
    
            setCauchyStress( true );
    
        }

        void residual::setPreviousdCauchyStressdF( ){
            /*!
             * Set the previous derivative of the computed Cauchy stress w.r.t. F (this is a partial derivative generally)
             */
    
            setCauchyStress( true );
    
        }

        void residual::setPreviousdCauchyStressdFn( ){
            /*!
             * Set the previous derivative of the computed Cauchy stress w.r.t. the configurations solved for in the non-linear solve
             * (this is a partial derivative generally)
             */
    
            setCauchyStress( true );
    
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

            const unsigned int dim = *hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int num_unknown_config_vars = ( num_configs - 1 ) * sot_dim;

            // Form the Jacobian
            floatMatrix jacobian = floatMatrix( sot_dim, floatVector( hydra->getUnknownVector( )->size( ), 0 ) );

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < dim; j++ ){

                    jacobian[ dim * i + j ][ dim * i + j ] = -1;

                    for ( unsigned int I = 0; I < num_unknown_config_vars; I++ ){

                        jacobian[ dim * i + j ][ getStress( )->size( ) + I ] = ( *get_dCauchyStressdFn( ) )[ dim * num_unknown_config_vars * i + num_unknown_config_vars * j + I ];

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

            const unsigned int dim = *hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            setdRdF( tardigradeVectorTools::inflate( *get_dCauchyStressdF( ), sot_dim, sot_dim ) );

        }

    }

}
