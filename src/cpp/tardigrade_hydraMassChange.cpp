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

            set_directionVector( floatVector( hydra->getAdditionalDOF( )->begin( ) + 2,
                                                     hydra->getAdditionalDOF( )->begin( ) + 5 ) );

            set_previousDensity( ( *hydra->getPreviousAdditionalDOF( ) )[ 0 ] );

            set_previousMassChangeRate( ( *hydra->getPreviousAdditionalDOF( ) )[ 1 ] );

            set_previousDirectionVector( floatVector( hydra->getPreviousAdditionalDOF( )->begin( ) + 2,
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

        void residual::setUnitDirectionVector( const bool &isPrevious ){
            /*!
             * Set the direction vector
             * 
             * \param &isPrevious: Flag for whether this is the previous or current direction
             */

            const unsigned int dim = hydra->getDimension( );

            const floatVector *directionVector;

            if ( isPrevious ){

                directionVector = get_previousDirectionVector( );

            }
            else{

                directionVector = get_directionVector( );

            }

            floatType normDirectionVector = tardigradeVectorTools::l2norm( *directionVector );

            floatVector unitDirectionVector( dim, 0 );

            if ( std::isfinite( 1. / normDirectionVector ) ){

                unitDirectionVector = ( *directionVector ) / normDirectionVector;

            }

            if ( isPrevious ){

                set_previousUnitDirectionVector( unitDirectionVector );

            }
            else{

                set_unitDirectionVector( unitDirectionVector );

            }

        }

        void residual::setUnitDirectionVector( ){
            /*!
             * Set the direction vector
             */

            setUnitDirectionVector( false );

        }

        void residual::setPreviousUnitDirectionVector( ){
            /*!
             * Set the previous direction vector
             */

            setUnitDirectionVector( true );

        }

        void residual::setUnitDirectionVectorDerivatives( const bool &isPrevious ){
            /*!
             * Set the direction vector derivatives
             * 
             * \param &isPrevious: Flag for whether this is the previous or current direction
             */

            const unsigned int dim = hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const floatVector *directionVector;

            if ( isPrevious ){

                directionVector = get_previousDirectionVector( );

            }
            else{

                directionVector = get_directionVector( );

            }

            floatType normDirectionVector = tardigradeVectorTools::l2norm( *directionVector );

            floatVector unitDirectionVector( dim, 0 );

            floatVector dUnitDirectionVectordDirectionVector( sot_dim, 0 );

            if ( std::isfinite( 1. / normDirectionVector ) ){

                unitDirectionVector = ( *directionVector ) / normDirectionVector;

                for ( unsigned int i = 0; i < dim; i++ ){

                    dUnitDirectionVectordDirectionVector[ dim * i + i ] += 1 / normDirectionVector;

                }

                for ( unsigned int i = 0; i < dim; i++ ){

                    for ( unsigned int j = 0; j < dim; j++ ){

                        dUnitDirectionVectordDirectionVector[ dim * i + j ] -= unitDirectionVector[ i ] * unitDirectionVector[ j ] / normDirectionVector;

                    }

                }

            }

            if ( isPrevious ){

                set_previousUnitDirectionVector( unitDirectionVector );

                set_dPreviousUnitDirectionVectordPreviousDirectionVector( dUnitDirectionVectordDirectionVector );

            }
            else{

                set_unitDirectionVector( unitDirectionVector );

                set_dUnitDirectionVectordDirectionVector( dUnitDirectionVectordDirectionVector );

            }


        }

        void residual::setdUnitDirectionVectordDirectionVector( ){
            /*!
             * Set the direction vector derivatives
             */

            setUnitDirectionVectorDerivatives( false );

        }

        void residual::setdPreviousUnitDirectionVectordPreviousDirectionVector( ){
            /*!
             * Set the previous direction vector derivatives
             */

            setUnitDirectionVectorDerivatives( true );

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

            const floatVector *unitDirectionVector;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradientTrace = get_previousMassChangeVelocityGradientTrace( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( unitDirectionVector = get_previousUnitDirectionVector( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradientTrace = get_massChangeVelocityGradientTrace( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( unitDirectionVector = get_unitDirectionVector( ) );

            }

            floatVector velocityGradient( sot_dim, 0 );

            floatType a = ( *velocityGradientTrace ) / ( 3. - 2. * ( *massDirectionMixingParameter ) );

            if ( tardigradeVectorTools::l2norm( *unitDirectionVector ) > 0.5 ){

                for ( unsigned int i = 0; i < dim; i++ ){

                    velocityGradient[ dim * i + i ] += a * ( 1 - *massDirectionMixingParameter );

                }

                for ( unsigned int i = 0; i < dim; i++ ){

                    for ( unsigned int j = 0; j < dim; j++ ){

                        velocityGradient[ dim * i + j ] += a * ( *massDirectionMixingParameter ) * ( *unitDirectionVector )[ i ] * ( *unitDirectionVector )[ j ];

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

            const floatVector *unitDirectionVector;

            const floatVector *dUnitDirectionVectordDirectionVector;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dVelocityGradientTracedDensity = get_dPreviousMassChangeVelocityGradientTracedPreviousDensity( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dVelocityGradientTracedMassChangeRate = get_dPreviousMassChangeVelocityGradientTracedPreviousMassChangeRate( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradientTrace = get_previousMassChangeVelocityGradientTrace( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dUnitDirectionVectordDirectionVector = get_dPreviousUnitDirectionVectordPreviousDirectionVector( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( unitDirectionVector = get_previousUnitDirectionVector( ) )

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dVelocityGradientTracedDensity = get_dMassChangeVelocityGradientTracedDensity( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dVelocityGradientTracedMassChangeRate = get_dMassChangeVelocityGradientTracedMassChangeRate( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradientTrace = get_massChangeVelocityGradientTrace( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dUnitDirectionVectordDirectionVector = get_dUnitDirectionVectordDirectionVector( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( unitDirectionVector = get_unitDirectionVector( ) )

            }

            floatVector velocityGradient( sot_dim, 0 );

            floatVector dVelocityGradientdDensity( sot_dim, 0 );

            floatVector dVelocityGradientdMassChangeRate( sot_dim, 0 );

            floatVector dVelocityGradientdDirectionVector( tot_dim, 0 );

            floatType a = ( *velocityGradientTrace ) / ( 3. - 2. * ( *massDirectionMixingParameter ) );

            floatType dadDensity = ( *dVelocityGradientTracedDensity ) / ( 3. - 2. * ( *massDirectionMixingParameter ) );

            floatType dadMassChangeRate = ( *dVelocityGradientTracedMassChangeRate ) / ( 3. - 2. * ( *massDirectionMixingParameter ) );

            if ( tardigradeVectorTools::l2norm( *unitDirectionVector ) > 0.5 ){

                for ( unsigned int i = 0; i < dim; i++ ){

                    velocityGradient[ dim * i + i ] += a * ( 1 - *massDirectionMixingParameter );

                    dVelocityGradientdDensity[ dim * i + i ] += dadDensity * ( 1 - *massDirectionMixingParameter );

                    dVelocityGradientdMassChangeRate[ dim * i + i ] += dadMassChangeRate * ( 1 - *massDirectionMixingParameter );

                }

                for ( unsigned int i = 0; i < dim; i++ ){

                    for ( unsigned int j = 0; j < dim; j++ ){

                        velocityGradient[ dim * i + j ] += a * ( *massDirectionMixingParameter ) * ( *unitDirectionVector )[ i ] * ( *unitDirectionVector )[ j ];

                        dVelocityGradientdDensity[ dim * i + j ] += dadDensity * ( *massDirectionMixingParameter ) * ( *unitDirectionVector )[ i ] * ( *unitDirectionVector )[ j ];

                        dVelocityGradientdMassChangeRate[ dim * i + j ] += dadMassChangeRate * ( *massDirectionMixingParameter ) * ( *unitDirectionVector )[ i ] * ( *unitDirectionVector )[ j ];

                        for ( unsigned int k = 0; k < dim; k++ ){

                            dVelocityGradientdDirectionVector[ dim * dim * i + dim * j + k ]
                                += a * ( *massDirectionMixingParameter ) * ( ( *dUnitDirectionVectordDirectionVector )[ dim * i + k ] * ( *unitDirectionVector )[ j ]
                                                                           + ( *unitDirectionVector )[ i ] * ( *dUnitDirectionVectordDirectionVector )[ dim * j + k ] );

                        }

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

                set_dPreviousMassChangeVelocityGradientdPreviousDirectionVector( dVelocityGradientdDirectionVector );

            }
            else{

                set_massChangeVelocityGradient( velocityGradient );

                set_dMassChangeVelocityGradientdDensity( dVelocityGradientdDensity );

                set_dMassChangeVelocityGradientdMassChangeRate( dVelocityGradientdMassChangeRate );

                set_dMassChangeVelocityGradientdDirectionVector( dVelocityGradientdDirectionVector );

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

        void residual::setdMassChangeVelocityGradientdDirectionVector( ){
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

        void residual::setdPreviousMassChangeVelocityGradientdPreviousDirectionVector( ){
            /*!
             * Set the derivative of the previous mass-change velocity gradient w.r.t. the previous mass change rate gradient
             */

            setMassChangeVelocityGradientDerivatives( true );

        }

        void residual::setPrecedingDeformationGradient( const bool &isPrevious ){
            /*!
             * Set the preceding deformation gradient
             *
             * \param &isPrevious: Flag for whether to set the current (false) or previous (true) value
             */

            if ( isPrevious ){

                set_previousPrecedingDeformationGradient( hydra->getPreviousPrecedingConfiguration( *getMassChangeConfigurationIndex( ) ) );

            }
            else{

                set_precedingDeformationGradient( hydra->getPrecedingConfiguration( *getMassChangeConfigurationIndex( ) ) );

            }

        }

        void residual::setPrecedingDeformationGradientDerivatives( const bool &isPrevious ){
            /*!
             * Set the derivatives of the preceding deformation gradient
             *
             * \param &isPrevious Flag for whether to set the current (false) or previous (true) values
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector *dF1dF;

            const floatVector *dF1dFn;

            floatVector dpFdFs;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dF = hydra->get_previousdF1dF( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dFn = hydra->get_previousdF1dFn( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dpFdFs = hydra->getPreviousPrecedingConfigurationJacobian( *getMassChangeConfigurationIndex( ) ) )

                set_previousPrecedingDeformationGradient( hydra->getPreviousPrecedingConfiguration( *getMassChangeConfigurationIndex( ) ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dF = hydra->get_dF1dF( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dFn = hydra->get_dF1dFn( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dpFdFs = hydra->getPrecedingConfigurationJacobian( *getMassChangeConfigurationIndex( ) ) )

                set_precedingDeformationGradient( hydra->getPrecedingConfiguration( *getMassChangeConfigurationIndex( ) ) );

            }

            floatVector dpFdF( sot_dim * sot_dim, 0 );

            floatVector dpFdFn( sot_dim * sot_dim * ( num_configs - 1 ), 0 );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int j = 0; j < sot_dim; j++ ){

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        dpFdF[ sot_dim * i + k ] += dpFdFs[ num_configs * sot_dim * i + j ] * ( *dF1dF )[ sot_dim * j + k ];

                    }

                }

            }

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){

                    dpFdFn[ ( num_configs - 1 ) * sot_dim * i + j ] += dpFdFs[ num_configs * sot_dim * i + j + sot_dim ];

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        dpFdFn[ ( num_configs - 1 ) * sot_dim * i + j ] += dpFdFs[ num_configs * sot_dim * i + k ] * ( *dF1dFn )[ ( num_configs - 1 ) * sot_dim * k + j ];

                    }

                }

            }

            if ( isPrevious ){

                set_dPreviousPrecedingDeformationGradientdPreviousDeformationGradient( dpFdF );

                set_dPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients( dpFdFn );

            }
            else{

                set_dPrecedingDeformationGradientdDeformationGradient( dpFdF );

                set_dPrecedingDeformationGradientdSubDeformationGradients( dpFdFn );

            }

        }

        void residual::setPrecedingDeformationGradient( ){
            /*!
             * Set the value of the preceding deformation gradient
             */

            setPrecedingDeformationGradient( false );

        }

        void residual::setPreviousPrecedingDeformationGradient( ){
            /*!
             * Set the value of the previous preceding deformation gradient
             */

            setPrecedingDeformationGradient( true );

        }

        void residual::setdPrecedingDeformationGradientdDeformationGradient( ){
            /*!
             * Set the derivative of the preceding deformation gradient w.r.t. the total deformation gradient
             */

            setPrecedingDeformationGradientDerivatives( false );

        }

        void residual::setdPrecedingDeformationGradientdSubDeformationGradients( ){
            /*!
             * Set the derivative of the preceding deformation gradient w.r.t. the sub-deformation gradients
             */

            setPrecedingDeformationGradientDerivatives( false );

        }

        void residual::setdPreviousPrecedingDeformationGradientdPreviousDeformationGradient( ){
            /*!
             * Set the derivative of the previous preceding deformation gradient w.r.t. the previous total deformation gradient
             */

            setPrecedingDeformationGradientDerivatives( true );

        }

        void residual::setdPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients( ){
            /*!
             * Set the derivative of the previous preceding deformation gradient w.r.t. the previous sub-deformation gradients
             */

            setPrecedingDeformationGradientDerivatives( true );

        }

        void residual::setMassChangeIntermediateVelocityGradient( const bool &isPrevious ){
            /*!
             * Set the velocity gradient in the intermediate configuration
             *
             * \param &isPrevious: Flag for whether this is being computed for the current or previous timestep
             */

            const floatVector *velocityGradient;

            const floatVector *precedingDeformationGradient;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradient = get_previousMassChangeVelocityGradient( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingDeformationGradient = get_previousPrecedingDeformationGradient( ) )

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradient = get_massChangeVelocityGradient( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingDeformationGradient = get_precedingDeformationGradient( ) )

            }

            floatVector intermediateVelocityGradient;

            tardigradeConstitutiveTools::pullBackVelocityGradient( *velocityGradient, *precedingDeformationGradient, intermediateVelocityGradient );

            if ( isPrevious ){

                set_previousMassChangeIntermediateVelocityGradient( intermediateVelocityGradient );

            }
            else{

                set_massChangeIntermediateVelocityGradient( intermediateVelocityGradient );

            }

        }

        void residual::setMassChangeIntermediateVelocityGradientDerivatives( const bool &isPrevious ){
            /*!
             * Set the derivatives of the velocity gradient in the intermediate configuration
             *
             * \param &isPrevious: Flag for whether this is being computed for the current or previous timestep
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            constexpr unsigned int tot_dim = sot_dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector *velocityGradient;

            const floatVector *precedingDeformationGradient;

            const floatVector *dLdRho;

            const floatVector *dLdC;

            const floatVector *dLdGradC;

            const floatVector *dPFdF;

            const floatVector *dPFdFn;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dLdRho = get_dPreviousMassChangeVelocityGradientdPreviousDensity( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dLdC = get_dPreviousMassChangeVelocityGradientdPreviousMassChangeRate( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dLdGradC = get_dPreviousMassChangeVelocityGradientdPreviousDirectionVector( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dPFdF = get_dPreviousPrecedingDeformationGradientdPreviousDeformationGradient( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dPFdFn = get_dPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradient = get_previousMassChangeVelocityGradient( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingDeformationGradient = get_previousPrecedingDeformationGradient( ) )

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dLdRho = get_dMassChangeVelocityGradientdDensity( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dLdC = get_dMassChangeVelocityGradientdMassChangeRate( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dLdGradC = get_dMassChangeVelocityGradientdDirectionVector( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dPFdF = get_dPrecedingDeformationGradientdDeformationGradient( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dPFdFn = get_dPrecedingDeformationGradientdSubDeformationGradients( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradient = get_massChangeVelocityGradient( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingDeformationGradient = get_precedingDeformationGradient( ) )

            }

            floatVector intermediateVelocityGradient;

            floatVector dILdL, dILdPF;

            tardigradeConstitutiveTools::pullBackVelocityGradient( *velocityGradient, *precedingDeformationGradient, intermediateVelocityGradient, dILdL, dILdPF );

            floatVector dILdRho( sot_dim, 0 );

            floatVector dILdC( sot_dim, 0 );

            floatVector dILdGradC( tot_dim, 0 );

            floatVector dILdF( sot_dim * sot_dim, 0 );

            floatVector dILdFn( ( num_configs - 1 ) * sot_dim * sot_dim, 0 );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int j = 0; j < sot_dim; j++ ){

                    dILdRho[ i ] += dILdL[ sot_dim * i + j ] * ( *dLdRho )[ j ];

                    dILdC[ i ] += dILdL[ sot_dim * i + j ] * ( *dLdC )[ j ];

                    for ( unsigned int k = 0; k < dim; k++ ){

                        dILdGradC[ dim * i + k ] += dILdL[ sot_dim * i + j ] * ( *dLdGradC )[ dim * j + k ];

                    }

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        dILdF[ sot_dim * i + k ] += dILdPF[ sot_dim * i + j ] * ( *dPFdF )[ sot_dim * j + k ];

                    }

                    for ( unsigned int k = 0; k < ( num_configs - 1 ) * sot_dim; k++ ){

                        dILdFn[ ( num_configs - 1 ) * sot_dim * i + k ] += dILdPF[ sot_dim * i + j ] * ( *dPFdFn )[ ( num_configs - 1 ) * sot_dim * j + k ];

                    }

                }

            }

            if ( isPrevious ){

                set_previousMassChangeIntermediateVelocityGradient( intermediateVelocityGradient );

                set_dPreviousMassChangeIntermediateVelocityGradientdPreviousDensity( dILdRho );

                set_dPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeRate( dILdC );

                set_dPreviousMassChangeIntermediateVelocityGradientdPreviousDirectionVector( dILdGradC );

                set_dPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient( dILdF );

                set_dPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients( dILdFn );

            }
            else{

                set_massChangeIntermediateVelocityGradient( intermediateVelocityGradient );

                set_dMassChangeIntermediateVelocityGradientdDensity( dILdRho );

                set_dMassChangeIntermediateVelocityGradientdMassChangeRate( dILdC );

                set_dMassChangeIntermediateVelocityGradientdDirectionVector( dILdGradC );

                set_dMassChangeIntermediateVelocityGradientdDeformationGradient( dILdF );

                set_dMassChangeIntermediateVelocityGradientdSubDeformationGradients( dILdFn );

            }

        }

        void residual::setMassChangeIntermediateVelocityGradient( ){
            /*!
             * Set the current intermediate velocity gradient
             */

            setMassChangeIntermediateVelocityGradient( false );

        }

        void residual::setPreviousMassChangeIntermediateVelocityGradient( ){
            /*!
             * Set the previous intermediate velocity gradient
             */

            setMassChangeIntermediateVelocityGradient( true );

        }

        void residual::setdMassChangeIntermediateVelocityGradientdDensity( ){
            /*!
             * Set the derivative of the current intermediate velocity gradient w.r.t. the density
             */

            setMassChangeIntermediateVelocityGradientDerivatives( false );

        }

        void residual::setdMassChangeIntermediateVelocityGradientdMassChangeRate( ){
            /*!
             * Set the derivative of the current intermediate velocity gradient w.r.t. the mass change rate
             */

            setMassChangeIntermediateVelocityGradientDerivatives( false );

        }

        void residual::setdMassChangeIntermediateVelocityGradientdDirectionVector( ){
            /*!
             * Set the derivative of the current intermediate velocity gradient w.r.t. the mass change rate gradient
             */

            setMassChangeIntermediateVelocityGradientDerivatives( false );

        }

        void residual::setdMassChangeIntermediateVelocityGradientdDeformationGradient( ){
            /*!
             * Set the derivative of the current intermediate velocity gradient w.r.t. the deformation gradient
             */

            setMassChangeIntermediateVelocityGradientDerivatives( false );

        }

        void residual::setdMassChangeIntermediateVelocityGradientdSubDeformationGradients( ){
            /*!
             * Set the derivative of the current intermediate velocity gradient w.r.t. the sub-deformation gradients
             */

            setMassChangeIntermediateVelocityGradientDerivatives( false );

        }

        void residual::setdPreviousMassChangeIntermediateVelocityGradientdPreviousDensity( ){
            /*!
             * Set the derivative of the previous intermediate velocity gradient w.r.t. the previous density
             */

            setMassChangeIntermediateVelocityGradientDerivatives( true );

        }

        void residual::setdPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeRate( ){
            /*!
             * Set the derivative of the previous intermediate velocity gradient w.r.t. the previous mass change rate
             */

            setMassChangeIntermediateVelocityGradientDerivatives( true );

        }

        void residual::setdPreviousMassChangeIntermediateVelocityGradientdPreviousDirectionVector( ){
            /*!
             * Set the derivative of the previous intermediate velocity gradient w.r.t. the previous mass change rate gradient
             */

            setMassChangeIntermediateVelocityGradientDerivatives( true );

        }

        void residual::setdPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient( ){
            /*!
             * Set the derivative of the previous intermediate velocity gradient w.r.t. the previous deformation gradient
             */

            setMassChangeIntermediateVelocityGradientDerivatives( true );

        }

        void residual::setdPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients( ){
            /*!
             * Set the derivative of the previous intermediate velocity gradient w.r.t. the previous sub-deformation gradients
             */

            setMassChangeIntermediateVelocityGradientDerivatives( true );

        }

        void residual::setMassChangeDeformationGradient( ){
            /*!
             * Set the mass-change deformation gradient
             */

            const floatVector *intermediateVelocityGradient = get_massChangeIntermediateVelocityGradient( );

            const floatVector *previousIntermediateVelocityGradient = get_previousMassChangeIntermediateVelocityGradient( );

            const floatVector previousMassChangeDeformationGradient = hydra->getPreviousConfiguration( *getMassChangeConfigurationIndex( ) );

            floatVector massChangeDeformationGradient;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::evolveFExponentialMap( *hydra->getDeltaTime( ), previousMassChangeDeformationGradient,
                                                                                              *previousIntermediateVelocityGradient, *intermediateVelocityGradient,
                                                                                               massChangeDeformationGradient,
                                                                                              *getIntegrationParameter( ) ) )

            set_massChangeDeformationGradient( massChangeDeformationGradient );

        }

        void residual::setMassChangeDeformationGradientDerivatives( const bool &computePrevious ){
            /*!
             * Compute the derivatives of the mass-change deformation gradient
             *
             * \param &computePrevious: Compute the gradients w.r.t. previous values
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            constexpr unsigned int tot_dim = sot_dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector *intermediateVelocityGradient = get_massChangeIntermediateVelocityGradient( );

            const floatVector *dLdRho   = get_dMassChangeIntermediateVelocityGradientdDensity( );

            const floatVector *dLdC     = get_dMassChangeIntermediateVelocityGradientdMassChangeRate( );

            const floatVector *dLdGradC = get_dMassChangeIntermediateVelocityGradientdDirectionVector( );

            const floatVector *dLdF     = get_dMassChangeIntermediateVelocityGradientdDeformationGradient( );

            const floatVector *dLdFn    = get_dMassChangeIntermediateVelocityGradientdSubDeformationGradients( );

            const floatVector *previousIntermediateVelocityGradient = get_previousMassChangeIntermediateVelocityGradient( );

            const floatVector previousMassChangeDeformationGradient = hydra->getPreviousConfiguration( *getMassChangeConfigurationIndex( ) );

            floatVector massChangeDeformationGradient;

            floatVector dFmdL;

            if ( computePrevious ){

                const floatVector *dLpdRho   = get_dPreviousMassChangeIntermediateVelocityGradientdPreviousDensity( );

                const floatVector *dLpdC     = get_dPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeRate( );

                const floatVector *dLpdGradC = get_dPreviousMassChangeIntermediateVelocityGradientdPreviousDirectionVector( );

                const floatVector *dLpdF     = get_dPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient( );

                const floatVector *dLpdFn    = get_dPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients( );

                floatVector dFmdFp;

                floatVector dFmdLp;

                TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::evolveFExponentialMap( *hydra->getDeltaTime( ), previousMassChangeDeformationGradient,
                                                                                                  *previousIntermediateVelocityGradient, *intermediateVelocityGradient,
                                                                                                  massChangeDeformationGradient,
                                                                                                  dFmdL, dFmdFp, dFmdLp,
                                                                                                  *getIntegrationParameter( ) ) )

                floatVector dFmdPreviousRho( sot_dim, 0 );

                floatVector dFmdPreviousC( sot_dim, 0 );

                floatVector dFmdPreviousGradC( tot_dim, 0 );

                floatVector dFmdPreviousF( sot_dim * sot_dim, 0 );

                floatVector dFmdPreviousFn( ( num_configs - 1 ) * sot_dim * sot_dim, 0 );

                for ( unsigned int i = 0; i < sot_dim; i++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        dFmdPreviousRho[ i ] += dFmdLp[ sot_dim * i + j ] * ( *dLpdRho )[ j ];

                        dFmdPreviousC[ i ] += dFmdLp[ sot_dim * i + j ] * ( *dLpdC )[ j ];

                        dFmdPreviousFn[ ( num_configs - 1 ) * sot_dim * i + j + ( ( *getMassChangeConfigurationIndex( ) ) - 1 ) * sot_dim ]
                            += dFmdFp[ sot_dim * i + j ];

                        for ( unsigned int k = 0; k < dim; k++ ){

                            dFmdPreviousGradC[ dim * i + k ] += dFmdLp[ sot_dim * i + j ] * ( *dLpdGradC )[ dim * j + k ];

                        }

                        for ( unsigned int k = 0; k < sot_dim; k++ ){

                            dFmdPreviousF[ sot_dim * i + k ] += dFmdLp[ sot_dim * i + j ] * ( *dLpdF )[ sot_dim * j + k ];

                        }

                        for ( unsigned int k = 0; k < ( num_configs - 1 ) * sot_dim; k++ ){

                            dFmdPreviousFn[ ( num_configs - 1 ) * sot_dim * i + k ] += dFmdLp[ sot_dim * i + j ] * ( *dLpdFn )[ ( num_configs - 1 ) * sot_dim * j + k ];

                        }

                    }

                }

                set_dMassChangeDeformationGradientdPreviousDensity( dFmdPreviousRho );

                set_dMassChangeDeformationGradientdPreviousMassChangeRate( dFmdPreviousC );

                set_dMassChangeDeformationGradientdPreviousDirectionVector( dFmdPreviousGradC );

                set_dMassChangeDeformationGradientdPreviousDeformationGradient( dFmdPreviousF );

                set_dMassChangeDeformationGradientdPreviousSubDeformationGradients( dFmdPreviousFn );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::evolveFExponentialMap( *hydra->getDeltaTime( ), previousMassChangeDeformationGradient,
                                                                                                  *previousIntermediateVelocityGradient, *intermediateVelocityGradient,
                                                                                                  massChangeDeformationGradient,
                                                                                                  dFmdL,
                                                                                                  *getIntegrationParameter( ) ) )

            }

            floatVector dFmdRho( sot_dim, 0 );

            floatVector dFmdC( sot_dim, 0 );

            floatVector dFmdGradC( tot_dim, 0 );

            floatVector dFmdF( sot_dim * sot_dim, 0 );

            floatVector dFmdFn( ( num_configs - 1 ) * sot_dim * sot_dim, 0 );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int j = 0; j < sot_dim; j++ ){

                    dFmdRho[ i ] += dFmdL[ sot_dim * i + j ] * ( *dLdRho )[ j ];

                    dFmdC[ i ] += dFmdL[ sot_dim * i + j ] * ( *dLdC )[ j ];

                    for ( unsigned int k = 0; k < dim; k++ ){

                        dFmdGradC[ dim * i + k ] += dFmdL[ sot_dim * i + j ] * ( *dLdGradC )[ dim * j + k ];

                    }

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        dFmdF[ sot_dim * i + k ] += dFmdL[ sot_dim * i + j ] * ( *dLdF )[ sot_dim * j + k ];

                    }

                    for ( unsigned int k = 0; k < ( num_configs - 1 ) * sot_dim; k++ ){

                        dFmdFn[ ( num_configs - 1 ) * sot_dim * i + k ] += dFmdL[ sot_dim * i + j ] * ( *dLdFn )[ ( num_configs - 1 ) * sot_dim * j + k ];

                    }

                }

            }

            set_massChangeDeformationGradient( massChangeDeformationGradient );

            set_dMassChangeDeformationGradientdDensity( dFmdRho );

            set_dMassChangeDeformationGradientdMassChangeRate( dFmdC );

            set_dMassChangeDeformationGradientdDirectionVector( dFmdGradC );

            set_dMassChangeDeformationGradientdDeformationGradient( dFmdF );

            set_dMassChangeDeformationGradientdSubDeformationGradients( dFmdFn );

        }

        void residual::setdMassChangeDeformationGradientdDensity( ){
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the density
             */

            setMassChangeDeformationGradientDerivatives( false );

        }

        void residual::setdMassChangeDeformationGradientdMassChangeRate( ){
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the mass change rate
             */

            setMassChangeDeformationGradientDerivatives( false );

        }

        void residual::setdMassChangeDeformationGradientdDirectionVector( ){
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the mass change rate gradient
             */

            setMassChangeDeformationGradientDerivatives( false );

        }

        void residual::setdMassChangeDeformationGradientdDeformationGradient( ){
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the deformation gradient
             */

            setMassChangeDeformationGradientDerivatives( false );

        }

        void residual::setdMassChangeDeformationGradientdSubDeformationGradients( ){
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the sub deformation gradients
             */

            setMassChangeDeformationGradientDerivatives( false );

        }

        void residual::setdMassChangeDeformationGradientdPreviousDensity( ){
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the previous density
             */

            setMassChangeDeformationGradientDerivatives( true );

        }

        void residual::setdMassChangeDeformationGradientdPreviousMassChangeRate( ){
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the previous mass change rate
             */

            setMassChangeDeformationGradientDerivatives( true );

        }

        void residual::setdMassChangeDeformationGradientdPreviousDirectionVector( ){
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the previous mass change rate gradient
             */

            setMassChangeDeformationGradientDerivatives( true );

        }

        void residual::setdMassChangeDeformationGradientdPreviousDeformationGradient( ){
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the previous deformation gradient
             */

            setMassChangeDeformationGradientDerivatives( true );

        }

        void residual::setdMassChangeDeformationGradientdPreviousSubDeformationGradients( ){
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the previous sub deformation gradients
             */

            setMassChangeDeformationGradientDerivatives( true );

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

            const unsigned int num_unknowns = hydra->getNumUnknowns( );

            const unsigned int num_equations = *getNumEquations( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            floatVector jacobian( num_equations * num_unknowns, 0 );

            const floatVector *dFmdFn = get_dMassChangeDeformationGradientdSubDeformationGradients( );

            for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                jacobian[ num_unknowns * i + sot_dim * ( *getMassChangeConfigurationIndex( ) ) + i ] += -1;

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){

                    jacobian[ num_unknowns * i + j + sot_dim ] += ( *dFmdFn )[ ( num_configs - 1 ) * sot_dim * i + j ];

                }

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

            setdRdF( *get_dMassChangeDeformationGradientdDeformationGradient( ) );

        }

        void residual::setdRdAdditionalDOF( ){
            /*!
             * Set the additional derivatives
             */

            const unsigned int dim = hydra->getDimension( );

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int num_equations = *getNumEquations( );

            const unsigned int num_additional_dof = hydra->getAdditionalDOF( )->size( );

            const floatVector *dMassChangeDeformationdDensity = get_dMassChangeDeformationGradientdDensity( );

            const floatVector *dMassChangeDeformationdMassChangeRate = get_dMassChangeDeformationGradientdMassChangeRate( );

            const floatVector *dMassChangeDeformationdDirectionVector = get_dMassChangeDeformationGradientdDirectionVector( );

            floatVector dRdAdditionalDOF( num_equations * num_additional_dof, 0 );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                dRdAdditionalDOF[ num_additional_dof * i + 0 ] = ( *dMassChangeDeformationdDensity )[ i ];

                dRdAdditionalDOF[ num_additional_dof * i + 1 ] = ( *dMassChangeDeformationdMassChangeRate )[ i ];

                for ( unsigned int j = 0; j < dim; j++ ){

                    dRdAdditionalDOF[ num_additional_dof * i + j + 2 ] = ( *dMassChangeDeformationdDirectionVector )[ dim * i + j ];

                }

            }

            setdRdAdditionalDOF( dRdAdditionalDOF );

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

        void residual::suggestInitialIterateValues( std::vector< unsigned int >   &indices,
                                                    std::vector< floatType > &values ){
            /*!
             * Suggest initial iterate values to try and improve convergence
             * 
             * \param &indices: The indices of the unknown vector to suggest initial values
             * \param &values: The values to suggest
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int configuration = *getMassChangeConfigurationIndex( );

            const floatVector *massChangeDeformationGradient = get_massChangeDeformationGradient( );

            indices = std::vector< unsigned int >( sot_dim, sot_dim * configuration );

            for ( unsigned int i = 0; i < sot_dim; i++ ){ indices[ i ] += i; }
            values = *massChangeDeformationGradient;

        }

    }

}
