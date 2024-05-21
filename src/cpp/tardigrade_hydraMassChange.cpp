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

            set_massChangeRateGradient( floatVector( hydra->getAdditionalDOF( )->begin( ) + 2,
                                                     hydra->getAdditionalDOF( )->begin( ) + 5 ) );

            set_previousDensity( ( *hydra->getPreviousAdditionalDOF( ) )[ 0 ] );

            set_previousMassChangeRate( ( *hydra->getPreviousAdditionalDOF( ) )[ 1 ] );

            set_previousMassChangeRateGradient( floatVector( hydra->getPreviousAdditionalDOF( )->begin( ) + 2,
                                                             hydra->getPreviousAdditionalDOF( )->begin( ) + 5 ) );

        }

        void residual::setMassChangeVelocityGradientTrace( const bool &isPrevious ){
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

        void residual::setMassChangeVelocityGradientTraceDerivatives( const bool &isPrevious ){
            /*!
             * Set the derivatives mass-change velocity gradient trace \f$ \left( \ell_{\bar{I}\bar{I}}^{A} \right) \f$
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

            floatType dMassChangeVelocityGradientTracedDensity = - ( *mass_change_rate ) / ( ( *density ) * ( *density ) );

            floatType dMassChangeVelocityGradientTracedMassChangeRate = 1. / ( *density );

            if ( isPrevious ){

                set_previousMassChangeVelocityGradientTrace( massChangeVelocityGradientTrace );

                set_dPreviousMassChangeVelocityGradientTracedPreviousDensity( dMassChangeVelocityGradientTracedDensity );

                set_dPreviousMassChangeVelocityGradientTracedPreviousMassChangeRate( dMassChangeVelocityGradientTracedMassChangeRate );

            }
            else{

                set_massChangeVelocityGradientTrace( massChangeVelocityGradientTrace );

                set_dMassChangeVelocityGradientTracedDensity( dMassChangeVelocityGradientTracedDensity );

                set_dMassChangeVelocityGradientTracedMassChangeRate( dMassChangeVelocityGradientTracedMassChangeRate );

            }

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

        void residual::setDirectionVector( const bool &isPrevious ){
            /*!
             * Set the direction vector
             * 
             * \param &isPrevious: Flag for whether this is the previous or current direction
             */

            const unsigned int dim = hydra->getDimension( );

            const floatVector *massChangeRateGradient;

            if ( isPrevious ){

                massChangeRateGradient = get_previousMassChangeRateGradient( );

            }
            else{

                massChangeRateGradient = get_massChangeRateGradient( );

            }

            floatType normMassChangeRateGradient = tardigradeVectorTools::l2norm( *massChangeRateGradient );

            floatVector directionVector( dim, 0 );

            if ( std::isfinite( normMassChangeRateGradient ) ){

                directionVector = ( *massChangeRateGradient ) / normMassChangeRateGradient;

            }

            if ( isPrevious ){

                set_previousDirectionVector( directionVector );

            }
            else{

                set_directionVector( directionVector );

            }

        }

        void residual::setDirectionVector( ){
            /*!
             * Set the direction vector
             */

            setDirectionVector( false );

        }

        void residual::setPreviousDirectionVector( ){
            /*!
             * Set the previous direction vector
             */

            setDirectionVector( true );

        }

        void residual::setDirectionVectorDerivatives( const bool &isPrevious ){
            /*!
             * Set the direction vector derivatives
             * 
             * \param &isPrevious: Flag for whether this is the previous or current direction
             */

            const unsigned int dim = hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const floatVector *massChangeRateGradient;

            if ( isPrevious ){

                massChangeRateGradient = get_previousMassChangeRateGradient( );

            }
            else{

                massChangeRateGradient = get_massChangeRateGradient( );

            }

            floatType normMassChangeRateGradient = tardigradeVectorTools::l2norm( *massChangeRateGradient );

            floatVector directionVector( dim, 0 );

            floatVector dDirectionVectordMassChangeRateGradient( sot_dim, 0 );

            if ( std::isfinite( normMassChangeRateGradient ) ){

                directionVector = ( *massChangeRateGradient ) / normMassChangeRateGradient;

                for ( unsigned int i = 0; i < dim; i++ ){

                    dDirectionVectordMassChangeRateGradient[ dim * i + i ] += 1 / normMassChangeRateGradient;

                }

                for ( unsigned int i = 0; i < dim; i++ ){

                    for ( unsigned int j = 0; j < dim; j++ ){

                        dDirectionVectordMassChangeRateGradient[ dim * i + j ] -= directionVector[ i ] * directionVector[ j ] / normMassChangeRateGradient;

                    }

                }

            }

            if ( isPrevious ){

                set_previousDirectionVector( directionVector );

                set_dPreviousDirectionVectordPreviousMassChangeRateGradient( dDirectionVectordMassChangeRateGradient );

            }
            else{

                set_directionVector( directionVector );

                set_dDirectionVectordMassChangeRateGradient( dDirectionVectordMassChangeRateGradient );

            }


        }

        void residual::setdDirectionVectordMassChangeRateGradient( ){
            /*!
             * Set the direction vector derivatives
             */

            setDirectionVectorDerivatives( false );

        }

        void residual::setdPreviousDirectionVectordPreviousMassChangeRateGradient( ){
            /*!
             * Set the previous direction vector derivatives
             */

            setDirectionVectorDerivatives( true );

        }

        void residual::setMassChangeVelocityGradient( const bool &isPrevious ){
            /*!
             * Set the mass-change velocity gradient
             * 
             * \param &isPrevious: Flag for whether this is the previous or current direction
             */

            const unsigned int dim = hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const floatType *massDirectionMixingParameter = get_massDirectionMixingParameter( );

            const floatType *velocityGradientTrace = get_massChangeVelocityGradientTrace( );

            const floatVector *directionVector;

            if ( isPrevious ){

                velocityGradientTrace = get_previousMassChangeVelocityGradientTrace( );

                directionVector = get_previousDirectionVector( );

            }
            else{

                velocityGradientTrace = get_massChangeVelocityGradientTrace( );

                directionVector = get_directionVector( );

            }

            floatVector velocityGradient( sot_dim, 0 );

            floatType a = ( *velocityGradientTrace ) / ( 3. - 2. * ( *massDirectionMixingParameter ) );

            if ( tardigradeVectorTools::l2norm( *directionVector ) > 0.5 ){

                for ( unsigned int i = 0; i < dim; i++ ){

                    velocityGradient[ dim * i + i ] += a * ( 1 - *massDirectionMixingParameter );

                }

                for ( unsigned int i = 0; i < dim; i++ ){

                    for ( unsigned int j = 0; j < dim; j++ ){

                        velocityGradient[ dim * i + j ] += a * ( *massDirectionMixingParameter ) * ( *directionVector )[ i ] * ( *directionVector )[ j ];

                    }

                }

            }
            else{

                for ( unsigned int i = 0; i < dim; i++ ){

                    velocityGradient[ dim * i + i ] += *velocityGradientTrace;

                }

            }

            if ( isPrevious ){

                set_previousMassChangeVelocityGradient( velocityGradient );

            }
            else{

                set_massChangeVelocityGradient( velocityGradient );

            }

        }

        void residual::setMassChangeVelocityGradientDerivatives( const bool &isPrevious ){
            /*!
             * Set the mass-change velocity gradient derivatives
             * 
             * \param &isPrevious: Flag for whether this is the previous or current direction
             */

            const unsigned int dim = hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int tot_dim = sot_dim * dim;

            const floatType *massDirectionMixingParameter = get_massDirectionMixingParameter( );

            const floatType *velocityGradientTrace;

            const floatType *dVelocityGradientTracedDensity;

            const floatType *dVelocityGradientTracedMassChangeRate;

            const floatVector *directionVector;

            if ( isPrevious ){

                dVelocityGradientTracedDensity = get_dPreviousMassChangeVelocityGradientTracedPreviousDensity( );

                dVelocityGradientTracedMassChangeRate = get_dPreviousMassChangeVelocityGradientTracedPreviousMassChangeRate( );

                velocityGradientTrace = get_previousMassChangeVelocityGradientTrace( );

                directionVector = get_previousDirectionVector( );

            }
            else{

                dVelocityGradientTracedDensity = get_dMassChangeVelocityGradientTracedDensity( );

                dVelocityGradientTracedMassChangeRate = get_dMassChangeVelocityGradientTracedMassChangeRate( );

                velocityGradientTrace = get_massChangeVelocityGradientTrace( );

                directionVector = get_directionVector( );

            }

            floatVector velocityGradient( sot_dim, 0 );

            floatVector dVelocityGradientdDensity( sot_dim, 0 );

            floatVector dVelocityGradientdMassChangeRate( sot_dim, 0 );

            floatVector dVelocityGradientdMassChangeRateGradient( tot_dim, 0 );

            floatType a = ( *velocityGradientTrace ) / ( 3. - 2. * ( *massDirectionMixingParameter ) );

            if ( tardigradeVectorTools::l2norm( *directionVector ) > 0.5 ){

                for ( unsigned int i = 0; i < dim; i++ ){

                    velocityGradient[ dim * i + i ] += a * ( 1 - *massDirectionMixingParameter );

                }

                for ( unsigned int i = 0; i < dim; i++ ){

                    for ( unsigned int j = 0; j < dim; j++ ){

                        velocityGradient[ dim * i + j ] += a * ( *massDirectionMixingParameter ) * ( *directionVector )[ i ] * ( *directionVector )[ j ];

                    }

                }

            }
            else{

                for ( unsigned int i = 0; i < dim; i++ ){

                    velocityGradient[ dim * i + i ] += *velocityGradientTrace;

                    dVelocityGradientdDensity[ dim * i + i ] += ( *dVelocityGradientTracedDensity );

                    dVelocityGradientdMassChangeRate[ dim * i + i ] += ( *dVelocityGradientTracedMassChangeRate );

                }

            }

            if ( isPrevious ){

                set_previousMassChangeVelocityGradient( velocityGradient );

                set_dPreviousMassChangeVelocityGradientdPreviousDensity( dVelocityGradientdDensity );

                set_dPreviousMassChangeVelocityGradientdPreviousMassChangeRate( dVelocityGradientdMassChangeRate );

                set_dPreviousMassChangeVelocityGradientdPreviousMassChangeRateGradient( dVelocityGradientdMassChangeRateGradient );

            }
            else{

                set_massChangeVelocityGradient( velocityGradient );

                set_dMassChangeVelocityGradientdDensity( dVelocityGradientdDensity );

                set_dMassChangeVelocityGradientdMassChangeRate( dVelocityGradientdMassChangeRate );

                set_dMassChangeVelocityGradientdMassChangeRateGradient( dVelocityGradientdMassChangeRateGradient );

            }

        }

        void residual::setMassChangeVelocityGradient( ){
            /*!
             * Set the value of the mass-change velocity gradient
             */

            setMassChangeVelocityGradient( false );

        }

        void residual::setdMassChangeVelocityGradientdDensity( ){
            /*!
             * Set the derivative of the mass-change velocity gradient w.r.t. the density
             */

            setMassChangeVelocityGradientDerivatives( false );

        }

        void residual::setdMassChangeVelocityGradientdMassChangeRate( ){
            /*!
             * Set the derivative of the mass-change velocity gradient w.r.t. the mass change rate
             */

            setMassChangeVelocityGradientDerivatives( false );

        }

        void residual::setdMassChangeVelocityGradientdMassChangeRateGradient( ){
            /*!
             * Set the derivative of the mass-change velocity gradient w.r.t. the mass change rate gradient
             */

            setMassChangeVelocityGradientDerivatives( false );

        }

        void residual::setPreviousMassChangeVelocityGradient( ){
            /*!
             * Set the value of the previous mass-change velocity gradient
             */

            setMassChangeVelocityGradient( true );

        }

        void residual::setdPreviousMassChangeVelocityGradientdPreviousDensity( ){
            /*!
             * Set the derivative of the previous mass-change velocity gradient w.r.t. the previous density
             */

            setMassChangeVelocityGradientDerivatives( true );

        }

        void residual::setdPreviousMassChangeVelocityGradientdPreviousMassChangeRate( ){
            /*!
             * Set the derivative of the previous mass-change velocity gradient w.r.t. the previous mass change rate
             */

            setMassChangeVelocityGradientDerivatives( true );

        }

        void residual::setdPreviousMassChangeVelocityGradientdPreviousMassChangeRateGradient( ){
            /*!
             * Set the derivative of the previous mass-change velocity gradient w.r.t. the previous mass change rate gradient
             */

            setMassChangeVelocityGradientDerivatives( true );

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
