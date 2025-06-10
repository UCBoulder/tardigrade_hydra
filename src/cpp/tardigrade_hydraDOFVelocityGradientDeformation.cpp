/**
  ******************************************************************************
  * \file tardigrade_hydraDOFVelocityGradientDeformation.cpp
  ******************************************************************************
  * An implementation of a deformation where the velocity gradient is in the
  * additionalDOF vector.
  ******************************************************************************
  */

#include<tardigrade_hydraDOFVelocityGradientDeformation.h>
#include<tardigrade_constitutive_tools.h>

namespace tardigradeHydra{

    namespace dofVelocityGradientDeformation{

        void residual::decomposeAdditionalDOF( ){
            /*!
             * Decompose the additional DOF vectors
             */

            constexpr unsigned int dim = 3;

            TARDIGRADE_ERROR_TOOLS_CHECK( hydra->getAdditionalDOF( )->size( ) >= getMassChangeVelocityGradientIndex( ) + dim * dim, "The additional DOF vector is of size " + std::to_string( hydra->getAdditionalDOF( )->size( ) ) + " which is less than the required size of " + std::to_string( getMassChangeVelocityGradientIndex( ) + dim * dim ) );

            TARDIGRADE_ERROR_TOOLS_CHECK( hydra->getPreviousAdditionalDOF( )->size( ) >= getMassChangeVelocityGradientIndex( ) + dim * dim, "The additional DOF vector is of size " + std::to_string( hydra->getPreviousAdditionalDOF( )->size( ) ) + " which is less than the required size of " + std::to_string( getMassChangeVelocityGradientIndex( ) + dim * dim ) );

            auto dofVelocityGradient         = get_setDataStorage_dofVelocityGradient( );

            auto previousMassChangeVelocityGradient = get_setDataStorage_previousMassChangeVelocityGradient( );

            *dofVelocityGradient.value = secondOrderTensor( hydra->getAdditionalDOF( )->begin( ) + getMassChangeVelocityGradientIndex( ),
                                                                   hydra->getAdditionalDOF( )->begin( ) + getMassChangeVelocityGradientIndex( ) + dim * dim );

            *previousMassChangeVelocityGradient.value = secondOrderTensor( hydra->getPreviousAdditionalDOF( )->begin( ) + getMassChangeVelocityGradientIndex( ),
                                                                           hydra->getPreviousAdditionalDOF( )->begin( ) + getMassChangeVelocityGradientIndex( ) + dim * dim );

        }

        void residual::setPrecedingDeformationGradient( const bool &isPrevious ){
            /*!
             * Set the preceding deformation gradient
             *
             * \param &isPrevious: Flag for whether to set the current (false) or previous (true) value
             */

            if ( isPrevious ){

                auto precedingDeformationGradient = get_setDataStorage_previousPrecedingDeformationGradient( );
                *precedingDeformationGradient.value = hydra->getPreviousPrecedingConfiguration( getMassChangeConfigurationIndex( ) );

            }
            else{

                auto precedingDeformationGradient = get_setDataStorage_precedingDeformationGradient( );
                *precedingDeformationGradient.value = hydra->getPrecedingConfiguration( getMassChangeConfigurationIndex( ) );

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

            const fourthOrderTensor *dF1dF;

            const floatVector *dF1dFn;

            floatVector dpFdFs;

            setDataStorageBase< secondOrderTensor > precedingDeformationGradient;

            setDataStorageBase< fourthOrderTensor > dpFdF;

            setDataStorageBase< floatVector > dpFdFn;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dF = hydra->get_previousdF1dF( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dFn = hydra->get_previousdF1dFn( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dpFdFs = hydra->getPreviousPrecedingConfigurationJacobian( getMassChangeConfigurationIndex( ) ) )

                auto precedingDeformationGradient = get_setDataStorage_previousPrecedingDeformationGradient( );
                *precedingDeformationGradient.value = hydra->getPreviousPrecedingConfiguration( getMassChangeConfigurationIndex( ) );

                dpFdF = get_setDataStorage_dPreviousPrecedingDeformationGradientdPreviousDeformationGradient( );

                dpFdFn = get_setDataStorage_dPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dF = hydra->get_dF1dF( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dFn = hydra->get_dF1dFn( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dpFdFs = hydra->getPrecedingConfigurationJacobian( getMassChangeConfigurationIndex( ) ) )

                auto precedingDeformationGradient = get_setDataStorage_precedingDeformationGradient( );
                *precedingDeformationGradient.value = hydra->getPrecedingConfiguration( getMassChangeConfigurationIndex( ) );

                dpFdF = get_setDataStorage_dPrecedingDeformationGradientdDeformationGradient( );

                dpFdFn = get_setDataStorage_dPrecedingDeformationGradientdSubDeformationGradients( );

            }

            dpFdF.zero( sot_dim * sot_dim );

            dpFdFn.zero( sot_dim * sot_dim * ( num_configs - 1 ) );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int j = 0; j < sot_dim; j++ ){

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        ( *dpFdF.value )[ sot_dim * i + k ] += dpFdFs[ num_configs * sot_dim * i + j ] * ( *dF1dF )[ sot_dim * j + k ];

                    }

                }

            }

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){

                    ( *dpFdFn.value )[ ( num_configs - 1 ) * sot_dim * i + j ] += dpFdFs[ num_configs * sot_dim * i + j + sot_dim ];

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        ( *dpFdFn.value )[ ( num_configs - 1 ) * sot_dim * i + j ] += dpFdFs[ num_configs * sot_dim * i + k ] * ( *dF1dFn )[ ( num_configs - 1 ) * sot_dim * k + j ];

                    }

                }

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

            const secondOrderTensor *velocityGradient;

            const secondOrderTensor *precedingDeformationGradient;

            setDataStorageBase< secondOrderTensor > intermediateVelocityGradient;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradient = get_previousMassChangeVelocityGradient( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingDeformationGradient = get_previousPrecedingDeformationGradient( ) )

                intermediateVelocityGradient = get_setDataStorage_previousMassChangeIntermediateVelocityGradient( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradient = get_dofVelocityGradient( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingDeformationGradient = get_precedingDeformationGradient( ) )

                intermediateVelocityGradient = get_setDataStorage_dofIntermediateVelocityGradient( );

            }

            tardigradeConstitutiveTools::pullBackVelocityGradient( *velocityGradient, *precedingDeformationGradient, *intermediateVelocityGradient.value );

        }

        void residual::setMassChangeIntermediateVelocityGradientDerivatives( const bool &isPrevious ){
            /*!
             * Set the derivatives of the velocity gradient in the intermediate configuration
             *
             * \param &isPrevious: Flag for whether this is being computed for the current or previous timestep
             */

            constexpr unsigned int dim = 3;

            auto sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const secondOrderTensor *velocityGradient;

            const secondOrderTensor *precedingDeformationGradient;

            const fourthOrderTensor *dPFdF;

            const floatVector *dPFdFn;

            setDataStorageBase< secondOrderTensor > intermediateVelocityGradient;

            setDataStorageBase< secondOrderTensor > dILdL;

            setDataStorageBase< fourthOrderTensor > dILdF;

            setDataStorageBase< fourthOrderTensor > dILdFn;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dPFdF = get_dPreviousPrecedingDeformationGradientdPreviousDeformationGradient( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dPFdFn = get_dPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradient = get_previousMassChangeVelocityGradient( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingDeformationGradient = get_previousPrecedingDeformationGradient( ) )

                intermediateVelocityGradient = get_setDataStorage_previousMassChangeIntermediateVelocityGradient( );

                dILdL                        = get_setDataStorage_dPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeVelocityGradient( );

                dILdF                        = get_setDataStorage_dPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient( );

                dILdFn                       = get_setDataStorage_dPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dPFdF = get_dPrecedingDeformationGradientdDeformationGradient( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( dPFdFn = get_dPrecedingDeformationGradientdSubDeformationGradients( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradient = get_dofVelocityGradient( ) )

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingDeformationGradient = get_precedingDeformationGradient( ) )

                intermediateVelocityGradient = get_setDataStorage_dofIntermediateVelocityGradient( );

                dILdL                        = get_setDataStorage_dMassChangeIntermediateVelocityGradientdMassChangeVelocityGradient( );

                dILdF                        = get_setDataStorage_dMassChangeIntermediateVelocityGradientdDeformationGradient( );

                dILdFn                       = get_setDataStorage_dMassChangeIntermediateVelocityGradientdSubDeformationGradients( );

            }

            fourthOrderTensor dILdPF;

            tardigradeConstitutiveTools::pullBackVelocityGradient( *velocityGradient, *precedingDeformationGradient, *intermediateVelocityGradient.value, *dILdL.value, dILdPF );

            dILdF.zero( sot_dim * sot_dim );

            dILdFn.zero( ( num_configs - 1 ) * sot_dim * sot_dim );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int j = 0; j < sot_dim; j++ ){

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        ( *dILdF.value )[ sot_dim * i + k ] += dILdPF[ sot_dim * i + j ] * ( *dPFdF )[ sot_dim * j + k ];

                    }

                    for ( unsigned int k = 0; k < ( num_configs - 1 ) * sot_dim; k++ ){

                        ( *dILdFn.value )[ ( num_configs - 1 ) * sot_dim * i + k ] += dILdPF[ sot_dim * i + j ] * ( *dPFdFn )[ ( num_configs - 1 ) * sot_dim * j + k ];

                    }

                }

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

        void residual::setdMassChangeIntermediateVelocityGradientdMassChangeVelocityGradient( ){
            /*!
             * Set the derivative of the current intermediate velocity gradient w.r.t. the mass change velocity gradient
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

        void residual::setdPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeVelocityGradient( ){
            /*!
             * Set the derivative of the previous intermediate velocity gradient w.r.t. the previous mass change velocity gradient
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

            const secondOrderTensor *intermediateVelocityGradient = get_dofIntermediateVelocityGradient( );

            const secondOrderTensor *previousIntermediateVelocityGradient = get_previousMassChangeIntermediateVelocityGradient( );

            const secondOrderTensor previousMassChangeDeformationGradient = hydra->getPreviousConfiguration( getMassChangeConfigurationIndex( ) );

            auto dofDeformationGradient = get_setDataStorage_dofDeformationGradient( );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::evolveFExponentialMap( *hydra->getDeltaTime( ), previousMassChangeDeformationGradient,
                                                                                              *previousIntermediateVelocityGradient, *intermediateVelocityGradient,
                                                                                              *dofDeformationGradient.value,
                                                                                              getIntegrationParameter( ) ) )

        }

        void residual::setMassChangeDeformationGradientDerivatives( const bool &computePrevious ){
            /*!
             * Compute the derivatives of the mass-change deformation gradient
             *
             * \param &computePrevious: Compute the gradients w.r.t. previous values
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const secondOrderTensor *intermediateVelocityGradient = get_dofIntermediateVelocityGradient( );

            const fourthOrderTensor *dILdL = get_dMassChangeIntermediateVelocityGradientdMassChangeVelocityGradient( );

            const fourthOrderTensor *dILdF = get_dMassChangeIntermediateVelocityGradientdDeformationGradient( );

            const floatVector *dILdFn      = get_dMassChangeIntermediateVelocityGradientdSubDeformationGradients( );

            const secondOrderTensor *previousIntermediateVelocityGradient = get_previousMassChangeIntermediateVelocityGradient( );

            const secondOrderTensor previousMassChangeDeformationGradient = hydra->getPreviousConfiguration( getMassChangeConfigurationIndex( ) );

            auto dofDeformationGradient = get_setDataStorage_dofDeformationGradient( );

            fourthOrderTensor dFmdIL;

            if ( computePrevious ){

                const fourthOrderTensor  *dILpdL = get_dPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeVelocityGradient( );

                const fourthOrderTensor *dILpdF  = get_dPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient( );

                const floatVector *dILpdFn       = get_dPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients( );

                fourthOrderTensor dFmdFp;

                fourthOrderTensor dFmdILp;

                TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::evolveFExponentialMap( *hydra->getDeltaTime( ), previousMassChangeDeformationGradient,
                                                                                                  *previousIntermediateVelocityGradient, *intermediateVelocityGradient,
                                                                                                  *dofDeformationGradient.value,
                                                                                                  dFmdIL, dFmdFp, dFmdILp,
                                                                                                  getIntegrationParameter( ) ) )

                auto dFmdPreviousL = get_setDataStorage_dMassChangeDeformationGradientdPreviousMassChangeVelocityGradient( );
                dFmdPreviousL.zero( sot_dim * sot_dim );

                auto dFmdPreviousF = get_setDataStorage_dMassChangeDeformationGradientdPreviousDeformationGradient( );
                dFmdPreviousF.zero( sot_dim * sot_dim );

                auto dFmdPreviousFn = get_setDataStorage_dMassChangeDeformationGradientdPreviousSubDeformationGradients( );
                dFmdPreviousFn.zero( sot_dim * sot_dim * ( num_configs - 1 ) );

                for ( unsigned int i = 0; i < sot_dim; i++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        ( *dFmdPreviousFn.value )[ ( num_configs - 1 ) * sot_dim * i + j + ( getMassChangeConfigurationIndex( ) - 1 ) * sot_dim ]
                            += dFmdFp[ sot_dim * i + j ];

                        for ( unsigned int k = 0; k < sot_dim; k++ ){

                            ( *dFmdPreviousL.value )[ sot_dim * i + k ] += dFmdILp[ sot_dim * i + j ] * ( *dILpdL )[ sot_dim * j + k ];

                            ( *dFmdPreviousF.value )[ sot_dim * i + k ] += dFmdILp[ sot_dim * i + j ] * ( *dILpdF )[ sot_dim * j + k ];

                        }

                        for ( unsigned int k = 0; k < ( num_configs - 1 ) * sot_dim; k++ ){

                            ( *dFmdPreviousFn.value )[ ( num_configs - 1 ) * sot_dim * i + k ] += dFmdILp[ sot_dim * i + j ] * ( *dILpdFn )[ ( num_configs - 1 ) * sot_dim * j + k ];

                        }

                    }

                }

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::evolveFExponentialMap( *hydra->getDeltaTime( ), previousMassChangeDeformationGradient,
                                                                                                  *previousIntermediateVelocityGradient, *intermediateVelocityGradient,
                                                                                                  *dofDeformationGradient.value,
                                                                                                  dFmdIL,
                                                                                                  getIntegrationParameter( ) ) )

            }

            auto dFmdL = get_setDataStorage_dMassChangeDeformationGradientdMassChangeVelocityGradient( );
            dFmdL.zero( sot_dim * sot_dim );

            auto dFmdF = get_setDataStorage_dMassChangeDeformationGradientdDeformationGradient( );
            dFmdF.zero( sot_dim * sot_dim );

            auto dFmdFn = get_setDataStorage_dMassChangeDeformationGradientdSubDeformationGradients( );
            dFmdFn.zero( sot_dim * sot_dim * ( num_configs - 1 ) );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int j = 0; j < sot_dim; j++ ){

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        ( *dFmdL.value )[ sot_dim * i + k ] += dFmdIL[ sot_dim * i + j ] * ( *dILdL )[ sot_dim * j + k ];

                        ( *dFmdF.value )[ sot_dim * i + k ] += dFmdIL[ sot_dim * i + j ] * ( *dILdF )[ sot_dim * j + k ];

                    }

                    for ( unsigned int k = 0; k < ( num_configs - 1 ) * sot_dim; k++ ){

                        ( *dFmdFn.value )[ ( num_configs - 1 ) * sot_dim * i + k ] += dFmdIL[ sot_dim * i + j ] * ( *dILdFn )[ ( num_configs - 1 ) * sot_dim * j + k ];

                    }

                }

            }

        }

        void residual::setdMassChangeDeformationGradientdMassChangeVelocityGradient( ){
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the mass change velocity gradient
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

        void residual::setdMassChangeDeformationGradientdPreviousMassChangeVelocityGradient( ){
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the previous mass change velocity gradient
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

            auto dofConfigurationIndex = getMassChangeConfigurationIndex( );

            auto residual = get_setDataStorage_residual( );

            *residual.value = *get_dofDeformationGradient( ) - secondOrderTensor( hydra->get_configurations( )->begin( ) +   dofConfigurationIndex * 9,
                                                                                         hydra->get_configurations( )->begin( ) + ( dofConfigurationIndex + 1 ) * 9 );

        }

        void residual::setJacobian( ){
            /*!
             * Set the values of the jacobian
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int num_unknowns = hydra->getNumUnknowns( );

            const unsigned int num_equations = *getNumEquations( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            auto jacobian = get_setDataStorage_jacobian( );
            jacobian.zero( num_equations * num_unknowns );

            const floatVector *dFmdFn = get_dMassChangeDeformationGradientdSubDeformationGradients( );

            for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                ( *jacobian.value )[ num_unknowns * i + sot_dim * getMassChangeConfigurationIndex( ) + i ] += -1;

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){

                    ( *jacobian.value )[ num_unknowns * i + j + sot_dim ] += ( *dFmdFn )[ ( num_configs - 1 ) * sot_dim * i + j ];

                }

            }

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            auto dRdT = get_setDataStorage_dRdT( );

            dRdT.zero( sot_dim );

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient
             */

            auto dRdF = get_setDataStorage_dRdF( );
            *dRdF.value = *get_dMassChangeDeformationGradientdDeformationGradient( );

        }

        void residual::setdRdAdditionalDOF( ){
            /*!
             * Set the additional derivatives
             */

            auto sot_dim = hydra->getSOTDimension( );

            auto num_equations = *getNumEquations( );

            auto num_additional_dof = hydra->getAdditionalDOF( )->size( );

            const fourthOrderTensor *dMassChangeDeformationdMassChangeVelocityGradient = get_dMassChangeDeformationGradientdMassChangeVelocityGradient( );

            auto dRdAdditionalDOF = get_setDataStorage_dRdAdditionalDOF( );
            dRdAdditionalDOF.zero( num_equations * num_additional_dof );

            auto offset = getMassChangeVelocityGradientIndex( );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int j = 0; j < sot_dim; j++ ){

                    ( *dRdAdditionalDOF.value )[ num_additional_dof * i + j + offset ] = ( *dMassChangeDeformationdMassChangeVelocityGradient )[ sot_dim * i + j ];

                }

            }

        }

        void residual::suggestInitialIterateValues( std::vector< unsigned int >   &indices,
                                                    std::vector< floatType > &values ){
            /*!
             * Suggest initial iterate values to try and improve convergence
             * 
             * \param &indices: The indices of the unknown vector to suggest initial values
             * \param &values: The values to suggest
             */

            auto sot_dim = hydra->getSOTDimension( );

            auto configuration = getMassChangeConfigurationIndex( );

            const secondOrderTensor *dofDeformationGradient = get_dofDeformationGradient( );

            indices = std::vector< unsigned int >( sot_dim, sot_dim * configuration );

            for ( unsigned int i = 0; i < sot_dim; i++ ){ indices[ i ] += i; }
            values = *dofDeformationGradient;

        }

    }

}
