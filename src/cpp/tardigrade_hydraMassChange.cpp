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

        void residual::decomposeAdditionalDOF( ){
            /*!
             * Decompose the additional DOF vectors
             */

            TARDIGRADE_ERROR_TOOLS_CHECK( hydra->getAdditionalDOF( )->size( ) >= 5, "The additional DOF vector is of size " + std::to_string( hydra->getAdditionalDOF( )->size( ) ) + " which is less than the required size of 5" );

            TARDIGRADE_ERROR_TOOLS_CHECK( hydra->getPreviousAdditionalDOF( )->size( ) >= 5, "The previous additional DOF vector is of size " + std::to_string( hydra->getPreviousAdditionalDOF( )->size( ) ) + " which is less than the required size of 5" );

            set_density( ( *hydra->getAdditionalDOF( ) )[ 0 ] );

            set_massChangeRate( ( *hydra->getAdditionalDOF( ) )[ 1 ] );

            set_massChangeRateGradient( floatVector( *hydra->getAdditionalDOF( )->begin( ) + 2,
                                                     *hydra->getAdditionalDOF( )->begin( ) + 5 ) );

            set_previousDensity( ( *hydra->getPreviousAdditionalDOF( ) )[ 0 ] );

            set_previousMassChangeRate( ( *hydra->getPreviousAdditionalDOF( ) )[ 1 ] );

            set_previousMassChangeRateGradient( floatVector( *hydra->getPreviousAdditionalDOF( )->begin( ) + 2,
                                                             *hydra->getPreviousAdditionalDOF( )->begin( ) + 5 ) );

        }

        void residual::setMassChangeVelocityGradientTrace( const bool isPrevious ){
            /*!
             * Set the mass-change velocity gradient trace \f$ \left( \ell_{\bar{I}\bar{I}}^{A} \right) \f$
             * 
             * \param isPrevious: Flag for whether to compute the current or previous value
             */

            const floatType *density;

            const floatType *mass_change_rate;

            if ( isPrevious ){

                density = get_previousDensity( );

                mass_change_rate = get_previousMassChangeRate( );

            }
            else{

                density = get_density( );

                mass_change_rate = get_massChangeRate( );

            }

            floatType massChangeVelocityGradientTrace = ( *mass_change_rate ) / ( *density );

            if ( isPrevious ){

                set_previousMassChangeVelocityGradientTrace( massChangeVelocityGradientTrace );

            }
            else{

                set_massChangeVelocityGradientTrace( massChangeVelocityGradientTrace );

            }

        }

        void residual::setMassChangeVelocityGradientTraceDerivatives( const bool isPrevious ){
            /*!
             * Set the derivatives mass-change velocity gradient trace \f$ \left( \ell_{\bar{I}\bar{I}}^{A} \right) \f$
             * 
             * \param isPrevious: Flag for whether to compute the current or previous value
             */

            TARDIGRADE_ERROR_TOOLS_CHECK( false, "not implemented" );

        }

        void residual::setMassChangeVelocityGradientTrace( ){
            /*!
             * Set the mass-change velocity gradient trace \f$ \left( \ell_{\bar{I}\bar{I}}^{A} \right) \f$
             */

            setMassChangeVelocityGradientTrace( false );

        }

        void residual::setdMassChangeVelocityGradientTracedDensity( ){
            /*!
             * Set the derivative of the mass-change velocity gradient trace \f$ \left( \ell_{\bar{I}\bar{I}}^{A} \right) \f$ w.r.t. the density
             */

            setMassChangeVelocityGradientTraceDerivatives( false );

        }

        void residual::setdMassChangeVelocityGradientTracedMassChangeRate( ){
            /*!
             * Set the derivative of the mass-change velocity gradient trace \f$ \left( \ell_{\bar{I}\bar{I}}^{A} \right) \f$ w.r.t. the mass change rate
             */

            setMassChangeVelocityGradientTraceDerivatives( false );

        }

        void residual::setPreviousMassChangeVelocityGradientTrace( ){
            /*!
             * Set the previous mass-change velocity gradient trace \f$ \left( \ell_{\bar{I}\bar{I}}^{A} \right) \f$
             */

            setMassChangeVelocityGradientTrace( true );

        }

        void residual::setdPreviousMassChangeVelocityGradientTracedPreviousDensity( ){
            /*!
             * Set the derivative of the previous mass-change velocity gradient trace \f$ \left( \ell_{\bar{I}\bar{I}}^{A} \right) \f$ w.r.t. previous density
             */

            setMassChangeVelocityGradientTraceDerivatives( true );

        }

        void residual::setdPreviousMassChangeVelocityGradientTracedPreviousMassChangeRate( ){
            /*!
             * Set the derivative of the previous mass-change velocity gradient trace \f$ \left( \ell_{\bar{I}\bar{I}}^{A} \right) \f$ w.r.t. the previous mass change rate
             */

            setMassChangeVelocityGradientTraceDerivatives( true );

        }

        void residual::setMassChangeVelocityGradientDirection( ){

            TARDIGRADE_ERROR_TOOLS_CHECK( false, "not implemented" );

        }

        void residual::setMassChangeVelocityGradient( ){

            TARDIGRADE_ERROR_TOOLS_CHECK( false, "not implemented" );

        }

        void residual::setMassChangeDeformationGradient( ){

            TARDIGRADE_ERROR_TOOLS_CHECK( false, "not implemented" );

        }

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
