/**
  ******************************************************************************
  * \file tardigrade_hydraLinearTestMaterial.cpp
  ******************************************************************************
  * An implementation of a linear test material using the hydra framework. Used
  * as an example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraLinearTestMaterial.h>
#include<tardigrade_constitutive_tools.h>

namespace tardigradeHydra{

    namespace linearTestMaterial{

        void residual::decomposeParameterVector( const floatVector &parameters ){
            /*!
             * Decompose the parameter vector
             * 
             * \param &parameters: The paramter vector. Assumed to be a vector of length 2 which defines lambda and mu.
             */
    
            if ( parameters.size( ) != 2 ){
    
                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "Parameter vector is expected to have a length of 2 but has a length of " + std::to_string( parameters.size( ) ) ) );
    
            }
    
        }

        void residual::setStress( ){
            /*!
             * Set the stress
             *
             * Currently uses the Cauchy stress
             */

//            setCauchyStress( false );
//
//            auto stress = get_setDataStorage_stress( );
//
//            *stress.value = *get_cauchyStress( );

        }
    
        void residual::setPreviousStress( ){
            /*!
             * Set the previous stress
             *
             * Currently uses the Cauchy stress
             */

//            setCauchyStress( true );
//
//            auto previousStress = get_setDataStorage_previousStress( );
//
//            *previousStress.value = *get_previousCauchyStress( );

        }
    
        void residual::setResidual( ){
            /*!
             * Set the residual value
             */

            auto residual = get_setDataStorage_residual( );

//            const secondOrderTensor *stress = getStress( );
//
 //           TARDIGRADE_ERROR_TOOLS_CATCH( *residual.value = *stress - *hydra->getStress( ) );
    
        }
    
        void residual::setJacobian( ){
            /*!
             * Set the Jacobian value
             */

//            const unsigned int dim = hydra->getDimension( );
//
//            const unsigned int sot_dim = hydra->getSOTDimension( );
//
//            const unsigned int num_configs = *hydra->getNumConfigurations( );
//
//            const unsigned int num_unknown_config_vars = ( num_configs - 1 ) * sot_dim;
//
//            const unsigned int num_unknowns = hydra->getNumUnknowns( );
//
//            // Form the Jacobian
//            auto jacobian = get_setDataStorage_jacobian( );
//
//            jacobian.zero( sot_dim * num_unknowns );
//
//            for ( unsigned int i = 0; i < dim; i++ ){
//
//                for ( unsigned int j = 0; j < dim; j++ ){
//
//                    ( *jacobian.value )[ num_unknowns * dim * i + num_unknowns * j + dim * i + j ] = -1;
//
//                    for ( unsigned int I = 0; I < num_unknown_config_vars; I++ ){
//
//                        ( *jacobian.value )[ num_unknowns * dim * i + num_unknowns * j + getStress( )->size( ) + I ] = ( *get_dCauchyStressdFn( ) )[ dim * num_unknown_config_vars * i + num_unknown_config_vars * j + I ];
//
//                    }
//
//                }
//
//            }

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

//            auto dRdF = get_setDataStorage_dRdF( );
//
//            *dRdF.value = *get_dCauchyStressdF( );

        }

    }

}
