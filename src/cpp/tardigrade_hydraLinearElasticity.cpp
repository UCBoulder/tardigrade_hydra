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

            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;

            if ( isPrevious ){

                auto previousFe = get_setDataStorage_previousFe( );

                *previousFe.value = secondOrderTensor( hydra->get_previousConfigurations( )->begin( ),
                                                       hydra->get_previousConfigurations( )->begin( ) + sot_dim );

            }
            else{

                auto Fe = get_setDataStorage_Fe( );
                
                *Fe.value = secondOrderTensor( hydra->get_configurations( )->begin( ),
                                               hydra->get_configurations( )->begin( ) + sot_dim );

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

                auto previousdFedF = get_setDataStorage_previousdFedF( );

                auto previousdFedFn = get_setDataStorage_previousdFedFn( );

                *previousdFedF.value  = *hydra->get_previousdF1dF( );

                *previousdFedFn.value = *hydra->get_previousdF1dFn( );

            }
            else{

                auto dFedF = get_setDataStorage_dFedF( );

                auto dFedFn = get_setDataStorage_dFedFn( );

                *dFedF.value  = *hydra->get_dF1dF( );

                *dFedFn.value = *hydra->get_dF1dFn( );

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

            const secondOrderTensor *Fe;
   
            setDataStorageBase< secondOrderTensor > Ee;

            setDataStorageBase< fourthOrderTensor > dEedFe;

            if ( isPrevious ){

                Fe = get_previousFe( );

                Ee = get_setDataStorage_previousEe( );

                dEedFe = get_setDataStorage_previousdEedFe( );

            }
            else{

                Fe = get_Fe( );

                Ee = get_setDataStorage_Ee( );

                dEedFe = get_setDataStorage_dEedFe( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeGreenLagrangeStrain( *Fe, *Ee.value, *dEedFe.value ) );

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
   
            constexpr unsigned int dim = 3;

            setDataStorageBase< secondOrderTensor > PK2Stress;

            const secondOrderTensor *Ee;

            if ( isPrevious ){

                Ee = get_previousEe( );

                PK2Stress = get_setDataStorage_previousPK2Stress( );

            }
            else{

                Ee = get_Ee( );

                PK2Stress = get_setDataStorage_PK2Stress( );

            }

            *PK2Stress.value = 2 * ( *getMu( ) ) * ( *Ee );

            floatType trace_Ee = tardigradeVectorTools::trace( *Ee );            

            for ( unsigned int i = 0; i < dim; i++ ){ ( *PK2Stress.value )[ dim * i + i ] += ( *getLambda( ) ) * trace_Ee; }

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
  
            const unsigned int dim     = hydra->getDimension( ); 
            const unsigned int fot_dim = hydra->getFOTDimension( );

            auto dPK2StressdEe = get_setDataStorage_dPK2StressdEe( );

            auto previousdPK2StressdEe = get_setDataStorage_previousdPK2StressdEe( );

            dPK2StressdEe.zero( fot_dim );
            for ( unsigned int i = 0; i < dim; i++ ){
                for ( unsigned int j = 0; j < dim; j++ ){
                    ( *dPK2StressdEe.value )[ dim * dim * ( dim * i + j ) + ( dim * i + j ) ] += 2 * ( *getMu( ) );
                    ( *dPK2StressdEe.value )[ dim * dim * dim * i + dim * dim * i + dim * j + j ] += ( *getLambda( ) );
                }
            }

            *previousdPK2StressdEe.value = *dPK2StressdEe.value;

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

             const unsigned int sot_dim = hydra->getSOTDimension( );

             const unsigned int fot_dim = hydra->getFOTDimension( );

             auto dPK2StressdFe = get_setDataStorage_dPK2StressdFe( );

             dPK2StressdFe.zero( fot_dim );

             for ( unsigned int I = 0; I < sot_dim; I++ ){

                 for ( unsigned int J = 0; J < sot_dim; J++ ){

                     for ( unsigned int K = 0; K < sot_dim; K++ ){

                         ( *dPK2StressdFe.value )[ sot_dim * I + K ] += ( *get_dPK2StressdEe( ) )[ sot_dim * I + J ] * ( *get_dEedFe( ) )[ sot_dim * J + K ];

                     }

                 }

             }

        }

        void residual::setdPK2StressdPreviousFe( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress w.r.t. the previous elastic
             * deformation gradient
             */


            auto dPK2StressdPreviousFe = get_setDataStorage_dPK2StressdPreviousFe( );

            dPK2StressdPreviousFe.zero( hydra->getFOTDimension( ) );

        }

        void residual::setPreviousdPK2StressdFe( ){
            /*!
             * Set the derivative of the previous second Piola-Kirchhoff stress w.r.t. the elastic
             * deformation gradient
             */

             auto previousdPK2StressdFe = get_setDataStorage_previousdPK2StressdFe( );

             previousdPK2StressdFe.zero( hydra->getFOTDimension( ) );

             const unsigned int sot_dim = hydra->getSOTDimension( );

             for ( unsigned int I = 0; I < sot_dim; I++ ){

                 for ( unsigned int J = 0; J < sot_dim; J++ ){

                     for ( unsigned int K = 0; K < sot_dim; K++ ){

                         ( *previousdPK2StressdFe.value )[ sot_dim * I + K ] += ( *get_previousdPK2StressdEe( ) )[ sot_dim * I + J ] * ( *get_previousdEedFe( ) )[ sot_dim * J + K ];

                     }

                 }

             }

        }

        void residual::setCauchyStress( const bool isPrevious ){
            /*!
             * Set the Cauchy stress
             * 
             * \param isPrevious: Whether to compute the current (false) or previous (true) Cauchy stress
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            constexpr unsigned int tot_dim = sot_dim * dim;

            constexpr unsigned int fot_dim = tot_dim * dim;

            unsigned int num_configs = *hydra->getNumConfigurations( );

            const secondOrderTensor *Fe;

            const fourthOrderTensor *dFedF;

            const floatVector *dFedFn;

            const secondOrderTensor *PK2Stress;

            const fourthOrderTensor *dPK2StressdFe;

            setDataStorageBase< secondOrderTensor > cauchyStress;

            setDataStorageBase< fourthOrderTensor > dCauchyStressdPK2Stress;

            setDataStorageBase< fourthOrderTensor > dCauchyStressdF;

            setDataStorageBase< floatVector > dCauchyStressdFn;

            if ( isPrevious ){

                Fe            = get_previousFe( );

                dFedF         = get_previousdFedF( );

                dFedFn        = get_previousdFedFn( );

                PK2Stress     = get_previousPK2Stress( );

                dPK2StressdFe = get_previousdPK2StressdFe( );

                cauchyStress            = get_setDataStorage_previousCauchyStress( );

                dCauchyStressdPK2Stress = get_setDataStorage_previousdCauchyStressdPK2Stress( );

                dCauchyStressdF         = get_setDataStorage_previousdCauchyStressdF( );

                dCauchyStressdFn        = get_setDataStorage_previousdCauchyStressdFn( );

            }
            else{

                Fe            = get_Fe( );

                dFedF         = get_dFedF( );

                dFedFn        = get_dFedFn( );

                PK2Stress     = get_PK2Stress( );

                dPK2StressdFe = get_dPK2StressdFe( );

                cauchyStress             = get_setDataStorage_cauchyStress( );

                dCauchyStressdPK2Stress  = get_setDataStorage_dCauchyStressdPK2Stress( );

                dCauchyStressdF          = get_setDataStorage_dCauchyStressdF( );

                dCauchyStressdFn         = get_setDataStorage_dCauchyStressdFn( );

            }

            // Compute the Second Piola-Kirchhoff stress and it's gradients
            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dPK2StressdFe( dPK2StressdFe->data( ), sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dFedF( dFedF->data( ), sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, -1, Eigen::RowMajor > > map_dFedFn( dFedFn->data( ), sot_dim, sot_dim * ( num_configs - 1 ) );

            fourthOrderTensor dPK2StressdF( fot_dim, 0 );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dPK2StressdF( dPK2StressdF.data( ), sot_dim, sot_dim );

            floatVector dPK2StressdFn( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            Eigen::Map< Eigen::Matrix< floatType, sot_dim, -1, Eigen::RowMajor > > map_dPK2StressdFn( dPK2StressdFn.data( ), sot_dim, sot_dim * ( num_configs - 1 ) );

            map_dPK2StressdF = ( map_dPK2StressdFe * map_dFedF ).eval( );
  
            map_dPK2StressdFn = ( map_dPK2StressdFe * map_dFedFn ).eval( );  

            // Map the PK2 stress to the current configuration
            fourthOrderTensor dCauchyStressdFe;
 
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::pushForwardPK2Stress( *PK2Stress, *Fe, *cauchyStress.value, *dCauchyStressdPK2Stress.value, dCauchyStressdFe ) );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dCauchyStressdPK2Stress( dCauchyStressdPK2Stress.value->data( ), sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dCauchyStressdFe( dCauchyStressdFe.data( ), sot_dim, sot_dim );

            dCauchyStressdF.zero( fot_dim );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dCauchyStressdF( dCauchyStressdF.value->data( ), sot_dim, sot_dim );
            map_dCauchyStressdF = ( map_dCauchyStressdPK2Stress * map_dPK2StressdF + map_dCauchyStressdFe * map_dFedF ).eval( );

            dCauchyStressdFn.zero( sot_dim * ( num_configs - 1 ) * sot_dim );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim, -1, Eigen::RowMajor > > map_dCauchyStressdFn( dCauchyStressdFn.value->data( ), sot_dim, sot_dim * ( num_configs - 1 ) );
            map_dCauchyStressdFn = ( map_dCauchyStressdPK2Stress * map_dPK2StressdFn + map_dCauchyStressdFe * map_dFedFn ).eval( );

            if ( !isPrevious ){

                auto dCauchyStressdPreviousF  = get_setDataStorage_dCauchyStressdPreviousF( );

                auto dCauchyStressdPreviousFn = get_setDataStorage_dCauchyStressdPreviousFn( );

                Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dPK2StressdPreviousFe( get_dPK2StressdPreviousFe( )->data( ), sot_dim, sot_dim );

                Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_previousdFedF( get_previousdFedF( )->data( ), sot_dim, sot_dim );

                Eigen::Map< const Eigen::Matrix< floatType, sot_dim, -1, Eigen::RowMajor > > map_previousdFedFn( get_previousdFedFn( )->data( ), sot_dim, sot_dim * ( num_configs - 1 ) );

                fourthOrderTensor dPK2StressdPreviousF( fot_dim, 0 );
                Eigen::Map< Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dPK2StressdPreviousF( dPK2StressdPreviousF.data( ), sot_dim, sot_dim );

                floatVector dPK2StressdPreviousFn( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );
                Eigen::Map< Eigen::Matrix< floatType, sot_dim, -1, Eigen::RowMajor > > map_dPK2StressdPreviousFn( dPK2StressdPreviousFn.data( ), sot_dim, sot_dim * ( num_configs - 1 ) );

                map_dPK2StressdPreviousF  = ( map_dPK2StressdPreviousFe * map_previousdFedF  ).eval( );

                map_dPK2StressdPreviousFn = ( map_dPK2StressdPreviousFe * map_previousdFedFn ).eval( );

                dCauchyStressdPreviousF.zero( fot_dim );
                Eigen::Map< Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dCauchyStressdPreviousF( dCauchyStressdPreviousF.value->data( ), sot_dim, sot_dim );

                dCauchyStressdPreviousFn.zero( sot_dim * ( num_configs - 1 ) * sot_dim );
                Eigen::Map< Eigen::Matrix< floatType, sot_dim, -1, Eigen::RowMajor > > map_dCauchyStressdPreviousFn( dCauchyStressdPreviousFn.value->data( ), sot_dim, sot_dim * ( num_configs - 1 ) );

                map_dCauchyStressdPreviousF  = ( map_dCauchyStressdPK2Stress * map_dPK2StressdPreviousF  ).eval( );

                map_dCauchyStressdPreviousFn = ( map_dCauchyStressdPK2Stress * map_dPK2StressdPreviousFn ).eval( );

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

            auto stress = get_setDataStorage_stress( );

            *stress.value = *get_cauchyStress( );

        }
    
        void residual::setPreviousStress( ){
            /*!
             * Set the previous stress
             *
             * Currently uses the Cauchy stress
             */

            setCauchyStress( true );

            auto previousStress = get_setDataStorage_previousStress( );

            *previousStress.value = *get_previousCauchyStress( );

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

            auto residual = get_setDataStorage_residual( );

            const secondOrderTensor *stress = getStress( );

            TARDIGRADE_ERROR_TOOLS_CATCH( *residual.value = *stress - *hydra->getStress( ) );
    
        }
    
        void residual::setJacobian( ){
            /*!
             * Set the Jacobian value
             */

            const unsigned int dim = hydra->getDimension( );

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int num_unknown_config_vars = ( num_configs - 1 ) * sot_dim;

            const unsigned int num_unknowns = hydra->getNumUnknowns( );

            // Form the Jacobian
            auto jacobian = get_setDataStorage_jacobian( );

            jacobian.zero( sot_dim * num_unknowns );

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < dim; j++ ){

                    ( *jacobian.value )[ num_unknowns * dim * i + num_unknowns * j + dim * i + j ] = -1;

                    for ( unsigned int I = 0; I < num_unknown_config_vars; I++ ){

                        ( *jacobian.value )[ num_unknowns * dim * i + num_unknowns * j + getStress( )->size( ) + I ] = ( *get_dCauchyStressdFn( ) )[ dim * num_unknown_config_vars * i + num_unknown_config_vars * j + I ];

                    }

                }

            }

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            auto dRdT = get_setDataStorage_dRdT( );

            dRdT.zero( *getNumEquations( ) );

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient
             */

            auto dRdF = get_setDataStorage_dRdF( );

            *dRdF.value = *get_dCauchyStressdF( );

        }

    }

}
