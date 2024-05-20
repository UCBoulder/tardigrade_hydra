/**
  ******************************************************************************
  * \file tardigrade_hydraMassChange.h
  ******************************************************************************
  * An implementation of the mass-change residual. Used as the basis for more
  * complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraMassChange.h>
#include<tardigrade_constitutive_tools.h>

namespace tardigradeHydra{

    namespace massChange{

        void residual::setResidual( ){
            /*!
             * Set the value of the residual
             * 
             * Defined as the residual's computed thermal deformation gradient minus the value stored in hydra's configurations.
             */

            const unsigned int massChangeConfigurationIndex = *getMassChangeConfigurationIndex( );

            setResidual( *get_massChangeDeformationGradient( ) - floatVector( hydra->get_configurations( )->begin( ) +   massChangeConfigurationIndex * 9,
                                                                              hydra->get_configurations( )->begin( ) + ( massChangeConfigurationIndex + 1 ) * 9 ) );

        }

        void residual::setJacobian( ){
            /*!
             * Set the values of the jacobian
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int num_unknowns = hydra->getUnknownVector( )->size( );

            floatVector jacobian( *getNumEquations( ) * num_unknowns, 0 );

            for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                jacobian[ num_unknowns * i + sot_dim * ( *getMassChangeConfigurationIndex( ) ) + i ] = -1;

            }

            setJacobian( jacobian );

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            setdRdT( floatVector( sot_dim, 0 ) );

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient
             */

            TARDIGRADE_ERROR_TOOLS_CHECK( false, "not implemented" );
            setdRdF( floatVector( *getNumEquations( ) * hydra->getDeformationGradient( )->size( ), 0 ) );

        }

        void residual::decomposeParameters( const floatVector &parameters ){
            /*!
             * Decompose the parameter vector
             *
             * \param &parameters: The parameter vector. Assumed to be
             *     of the form ( d ) where d is the mixing parameter between
             *     a spherical and directional response (0 <= d <= 1).
             */

            TARDIGRADE_ERROR_TOOLS_CHECK( parameters.size( ) == 1, "The parameter vector must have a length of 1 rather than " + std::to_string( parameters.size( ) ) );

            set_massDirectionMixingParameter( parameters[ 0 ] );

        }

    }

}
