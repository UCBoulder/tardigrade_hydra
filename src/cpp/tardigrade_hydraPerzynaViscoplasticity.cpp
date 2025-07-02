/**
  ******************************************************************************
  * \file tardigrade_hydraPerzynaViscoplasticity.h
  ******************************************************************************
  * An implementation of perzynaViscoplasticity using the hydra framework. Used
  * as an example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraPerzynaViscoplasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_stress_tools.h>

namespace tardigradeHydra{

    namespace perzynaViscoplasticity{

        void residual::setDrivingStress( ){
            /*!
             * Set the driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration.
             */

            setDrivingStress( false );

        }

        void residual::setdDrivingStressdCauchyStress( ){
            /*!
             * Set the derivative of the driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration with respect to the Cauchy Stress.
             */

            setdDrivingStressdCauchyStress( false );

        }

        void residual::setdDrivingStressdF( ){
            /*!
             * Set the derivative of the driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration with respect to the deformation
             * gradient
             */

            setdDrivingStressdF( false );

        }

        void residual::setdDrivingStressdSubFs( ){
            /*!
             * Set the derivative of the driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration with respect to the sub-deformation
             * gradients
             */

            setdDrivingStressdSubFs( false );

        }

        void residual::setPreviousDrivingStress( ){
            /*!
             * Set the driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration.
             */

            setDrivingStress( true );

        }

        void residual::setdPreviousDrivingStressdPreviousCauchyStress( ){
            /*!
             * Set the derivative of the previous driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration with respect to the previous Cauchy Stress.
             */

            setdDrivingStressdCauchyStress( true );

        }

        void residual::setdPreviousDrivingStressdPreviousF( ){
            /*!
             * Set the derivative of the previous driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration with respect to the previous deformation
             * gradient
             */

            setdDrivingStressdF( true );

        }

        void residual::setdPreviousDrivingStressdPreviousSubFs( ){
            /*!
             * Set the derivative of the previous driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration with respect to the previous sub-deformation
             * gradients
             */

            setdDrivingStressdSubFs( true );

        }

        void residual::setDrivingStress( const bool isPrevious ){
            /*!
             * Set the driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration.
             *
             * \param isPrevious: Flag for whether to compute this in the previous configuration
             */

            const floatVector *cauchyStress;

            floatVector precedingConfiguration;

            setDataStorageBase< secondOrderTensor > drivingStress;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfiguration = hydra->getPreviousPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getPreviousStress( ) );

                drivingStress = get_setDataStorage_previousDrivingStress( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfiguration = hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getStress( ) );

                drivingStress = get_setDataStorage_drivingStress( );

            }

            tardigradeConstitutiveTools::pullBackCauchyStress( *cauchyStress, precedingConfiguration, *drivingStress.value );

        }

        void residual::setDrivingStressDerivatives( const bool isPrevious ){
            /*!
             * Set the derivatives of the driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration.
             *
             * \param isPrevious: Flag for whether to compute this in the previous configuration
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            constexpr unsigned int tot_dim = sot_dim * dim;

            constexpr unsigned int fot_dim = tot_dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const secondOrderTensor *cauchyStress;

            secondOrderTensor precedingConfiguration;

            setDataStorageBase< secondOrderTensor > drivingStress;

            setDataStorageBase< fourthOrderTensor > dDrivingStressdCauchyStress;

            setDataStorageBase< fourthOrderTensor > dDrivingStressdF;

            setDataStorageBase< floatVector > dDrivingStressdSubFs;

            floatVector precedingConfigurationJacobian;

            const floatVector *dF1dF;

            const floatVector *dF1dSubFs;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfiguration = hydra->getPreviousPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfigurationJacobian = hydra->getPreviousPrecedingConfigurationJacobian( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getPreviousStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dF = hydra->get_previousdF1dF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dSubFs = hydra->get_previousdF1dFn( ) );

                drivingStress = get_setDataStorage_previousDrivingStress( );

                dDrivingStressdCauchyStress = get_setDataStorage_dPreviousDrivingStressdPreviousCauchyStress( );

                dDrivingStressdF            = get_setDataStorage_dPreviousDrivingStressdPreviousF( );

                dDrivingStressdSubFs        = get_setDataStorage_dPreviousDrivingStressdPreviousSubFs( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfiguration = hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfigurationJacobian = hydra->getPrecedingConfigurationJacobian( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dF = hydra->get_dF1dF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dSubFs = hydra->get_dF1dFn( ) );

                drivingStress = get_setDataStorage_drivingStress( );

                dDrivingStressdCauchyStress = get_setDataStorage_dDrivingStressdCauchyStress( );

                dDrivingStressdF            = get_setDataStorage_dDrivingStressdF( );

                dDrivingStressdSubFs        = get_setDataStorage_dDrivingStressdSubFs( );

            }

            floatVector dDrivingStressdPrecedingF;

            tardigradeConstitutiveTools::pullBackCauchyStress( *cauchyStress, precedingConfiguration, *drivingStress.value,
                                                               *dDrivingStressdCauchyStress.value, dDrivingStressdPrecedingF );

            auto map_dDrivingStressdPrecedingF      = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dDrivingStressdPrecedingF.data( ) );
            auto map_precedingConfigurationJacobian = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( precedingConfigurationJacobian.data( ), num_configs * sot_dim );

            floatVector dDrivingStressdFn( num_configs * fot_dim, 0 );
            auto map_dDrivingStressdFn = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dDrivingStressdFn.data( ), num_configs * sot_dim );

            map_dDrivingStressdFn = ( map_dDrivingStressdPrecedingF * map_precedingConfigurationJacobian ).eval( );

            dDrivingStressdF.zero( fot_dim );

            dDrivingStressdSubFs.zero( sot_dim * ( num_configs - 1 ) * sot_dim );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int k = 0; k < sot_dim; k++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        ( *dDrivingStressdF.value )[ sot_dim * i + j ] += dDrivingStressdFn[ num_configs * sot_dim * i + k ] * ( *dF1dF )[ sot_dim * k + j ];

                    }

                }

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){ //TODO: This order messes a bit with the cache. Investigate how to make it faster.

                    ( *dDrivingStressdSubFs.value )[ ( num_configs - 1 ) * sot_dim * i + j ] += dDrivingStressdFn[ num_configs * sot_dim * i + j + sot_dim ];

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        ( *dDrivingStressdSubFs.value )[ ( num_configs - 1 ) * sot_dim * i + j ] += dDrivingStressdFn[  num_configs * sot_dim * i + k ] * ( *dF1dSubFs )[ ( num_configs - 1 ) * sot_dim * k + j ];

                    }

                }

            }

        }

        void residual::setdDrivingStressdCauchyStress( const bool isPrevious ){
            /*!
             * Set the derivatives of the driving stress w.r.t. the Cauchy stress
             * 
             * \param isPrevious: Flag for whether to compute this in the previous configuration
             */

            TARDIGRADE_ERROR_TOOLS_CATCH( setDrivingStressDerivatives( isPrevious ) );

        }

        void residual::setdDrivingStressdF( const bool isPrevious ){
            /*!
             * Set the derivatives of the driving stress w.r.t. the deformation gradient
             * 
             * \param isPrevious: Flag for whether to compute this in the previous configuration
             */

            TARDIGRADE_ERROR_TOOLS_CATCH( setDrivingStressDerivatives( isPrevious ) );

        }

        void residual::setdDrivingStressdSubFs( const bool isPrevious ){
            /*!
             * Set the derivatives of the driving stress w.r.t. the sub-deformation gradients
             * 
             * \param isPrevious: Flag for whether to compute this in the previous configuration
             */

            TARDIGRADE_ERROR_TOOLS_CATCH( setDrivingStressDerivatives( isPrevious ) );

        }

        void residual::setFlowDirection( ){
            /*!
             * Set the flow direction in the current configuration of the
             * plastic configuration.
             * 
             * Current formulation is Drucker-Prager based on the potential
             * 
             * \f$ g = \hat{\sigma} + A \bar{\sigma}\f$
             * 
             * where
             * 
             * \f$ \hat{\sigma} = \sqrt{\frac{3}{2} \sigma_{ij}^{\text{dev}} \sigma_{ij}^{\text{dev}} }\f$
             * \f$ \bar{\sigma} = \sigma_{ii}\f$
             * 
             * \f$ \sigma_{ij}^{\text{dev}} = \sigma_{ij} - \frac{1}{3} \bar{\sigma} \delta_{ij}\f$
             */

            setFlowDirection( false );

        }

        void residual::setdFlowDirectiondCauchyStress( const bool isPrevious ){
            /*!
             * Set the derivative of the flow direction in the current configuration of the
             * plastic configuration w.r.t. the Cauchy stress
             * 
             * Current formulation is Drucker-Prager based on the potential
             * 
             * \f$ g = \hat{\sigma} + A \bar{\sigma}\f$
             * 
             * where
             * 
             * \f$ \hat{\sigma} = \sqrt{\frac{3}{2} \sigma_{ij}^{\text{dev}} \sigma_{ij}^{\text{dev}} }\f$
             * \f$ \bar{\sigma} = \sigma_{ii}\f$
             * 
             * \f$ \sigma_{ij}^{\text{dev}} = \sigma_{ij} - \frac{1}{3} \bar{\sigma} \delta_{ij}\f$
             * 
             * \param isPrevious: Flag for whether to compute the values in the previous configuration
             */

            setFlowDirectionDerivatives( isPrevious );

        }

        void residual::setdFlowDirectiondF( const bool isPrevious ){
            /*!
             * Set the derivative of the flow direction in the current configuration of the
             * plastic configuration w.r.t. the deformation gradient
             * 
             * Current formulation is Drucker-Prager based on the potential
             * 
             * \f$ g = \hat{\sigma} + A \bar{\sigma}\f$
             * 
             * where
             * 
             * \f$ \hat{\sigma} = \sqrt{\frac{3}{2} \sigma_{ij}^{\text{dev}} \sigma_{ij}^{\text{dev}} }\f$
             * \f$ \bar{\sigma} = \sigma_{ii}\f$
             * 
             * \f$ \sigma_{ij}^{\text{dev}} = \sigma_{ij} - \frac{1}{3} \bar{\sigma} \delta_{ij}\f$
             * 
             * \param isPrevious: Flag for whether to compute the values in the previous configuration
             */

            setFlowDirectionDerivatives( isPrevious );

        }

        void residual::setdFlowDirectiondSubFs( const bool isPrevious ){
            /*!
             * Set the derivative of the flow direction in the current configuration of the
             * plastic configuration w.r.t. the sub-deformation gradients
             * 
             * Current formulation is Drucker-Prager based on the potential
             * 
             * \f$ g = \hat{\sigma} + A \bar{\sigma}\f$
             * 
             * where
             * 
             * \f$ \hat{\sigma} = \sqrt{\frac{3}{2} \sigma_{ij}^{\text{dev}} \sigma_{ij}^{\text{dev}} }\f$
             * \f$ \bar{\sigma} = \sigma_{ii}\f$
             * 
             * \f$ \sigma_{ij}^{\text{dev}} = \sigma_{ij} - \frac{1}{3} \bar{\sigma} \delta_{ij}\f$
             * 
             * \param isPrevious: Flag for whether to compute the values in the previous configuration
             */

            setFlowDirectionDerivatives( isPrevious );

        }

        void residual::setdFlowDirectiondCauchyStress( ){
            /*!
             * Set the derivative of the flow direction in the current configuration of the
             * plastic configuration w.r.t. the Cauchy stress
             * 
             * Current formulation is Drucker-Prager based on the potential
             * 
             * \f$ g = \hat{\sigma} + A \bar{\sigma}\f$
             * 
             * where
             * 
             * \f$ \hat{\sigma} = \sqrt{\frac{3}{2} \sigma_{ij}^{\text{dev}} \sigma_{ij}^{\text{dev}} }\f$
             * \f$ \bar{\sigma} = \sigma_{ii}\f$
             * 
             * \f$ \sigma_{ij}^{\text{dev}} = \sigma_{ij} - \frac{1}{3} \bar{\sigma} \delta_{ij}\f$
             */

            setdFlowDirectiondCauchyStress( false );

        }

        void residual::setdFlowDirectiondF( ){
            /*!
             * Set the derivative of the flow direction in the current configuration of the
             * plastic configuration w.r.t. the deformation gradient
             * 
             * Current formulation is Drucker-Prager based on the potential
             * 
             * \f$ g = \hat{\sigma} + A \bar{\sigma}\f$
             * 
             * where
             * 
             * \f$ \hat{\sigma} = \sqrt{\frac{3}{2} \sigma_{ij}^{\text{dev}} \sigma_{ij}^{\text{dev}} }\f$
             * \f$ \bar{\sigma} = \sigma_{ii}\f$
             * 
             * \f$ \sigma_{ij}^{\text{dev}} = \sigma_{ij} - \frac{1}{3} \bar{\sigma} \delta_{ij}\f$
             */

            setdFlowDirectiondF( false );

        }

        void residual::setdFlowDirectiondSubFs( ){
            /*!
             * Set the derivative of the flow direction in the current configuration of the
             * plastic configuration w.r.t. the sub-deformation gradients
             * 
             * Current formulation is Drucker-Prager based on the potential
             * 
             * \f$ g = \hat{\sigma} + A \bar{\sigma}\f$
             * 
             * where
             * 
             * \f$ \hat{\sigma} = \sqrt{\frac{3}{2} \sigma_{ij}^{\text{dev}} \sigma_{ij}^{\text{dev}} }\f$
             * \f$ \bar{\sigma} = \sigma_{ii}\f$
             * 
             * \f$ \sigma_{ij}^{\text{dev}} = \sigma_{ij} - \frac{1}{3} \bar{\sigma} \delta_{ij}\f$
             */

            setdFlowDirectiondSubFs( false );

        }

        void residual::setPreviousFlowDirection( ){
            /*!
             * Set the flow direction in the current configuration of the
             * plastic configuration.
             * 
             * Current formulation is Drucker-Prager based on the potential
             * 
             * \f$ g = \hat{\sigma} + A \bar{\sigma}\f$
             * 
             * where
             * 
             * \f$ \hat{\sigma} = \sqrt{\frac{3}{2} \sigma_{ij}^{\text{dev}} \sigma_{ij}^{\text{dev}} }\f$
             * \f$ \bar{\sigma} = \sigma_{ii}\f$
             * 
             * \f$ \sigma_{ij}^{\text{dev}} = \sigma_{ij} - \frac{1}{3} \bar{\sigma} \delta_{ij}\f$
             */

            setFlowDirection( true );

        }

        void residual::setdPreviousFlowDirectiondPreviousCauchyStress( ){
            /*!
             * Set the derivative of the previous flow direction in the current configuration of the
             * plastic configuration w.r.t. the previous Cauchy stress
             * 
             * Current formulation is Drucker-Prager based on the potential
             * 
             * \f$ g = \hat{\sigma} + A \bar{\sigma}\f$
             * 
             * where
             * 
             * \f$ \hat{\sigma} = \sqrt{\frac{3}{2} \sigma_{ij}^{\text{dev}} \sigma_{ij}^{\text{dev}} }\f$
             * \f$ \bar{\sigma} = \sigma_{ii}\f$
             * 
             * \f$ \sigma_{ij}^{\text{dev}} = \sigma_{ij} - \frac{1}{3} \bar{\sigma} \delta_{ij}\f$
             */

            setdFlowDirectiondCauchyStress( true );

        }

        void residual::setdPreviousFlowDirectiondPreviousF( ){
            /*!
             * Set the derivative of the previous flow direction in the current configuration of the
             * plastic configuration w.r.t. the previous deformation gradient
             * 
             * Current formulation is Drucker-Prager based on the potential
             * 
             * \f$ g = \hat{\sigma} + A \bar{\sigma}\f$
             * 
             * where
             * 
             * \f$ \hat{\sigma} = \sqrt{\frac{3}{2} \sigma_{ij}^{\text{dev}} \sigma_{ij}^{\text{dev}} }\f$
             * \f$ \bar{\sigma} = \sigma_{ii}\f$
             * 
             * \f$ \sigma_{ij}^{\text{dev}} = \sigma_{ij} - \frac{1}{3} \bar{\sigma} \delta_{ij}\f$
             */

            setdFlowDirectiondF( true );

        }

        void residual::setdPreviousFlowDirectiondPreviousSubFs( ){
            /*!
             * Set the derivative of the previous flow direction in the current configuration of the
             * plastic configuration w.r.t. the previous sub-deformation gradients
             * 
             * Current formulation is Drucker-Prager based on the potential
             * 
             * \f$ g = \hat{\sigma} + A \bar{\sigma}\f$
             * 
             * where
             * 
             * \f$ \hat{\sigma} = \sqrt{\frac{3}{2} \sigma_{ij}^{\text{dev}} \sigma_{ij}^{\text{dev}} }\f$
             * \f$ \bar{\sigma} = \sigma_{ii}\f$
             * 
             * \f$ \sigma_{ij}^{\text{dev}} = \sigma_{ij} - \frac{1}{3} \bar{\sigma} \delta_{ij}\f$
             */

            setdFlowDirectiondSubFs( true );

        }

        void residual::setFlowDirection( const bool isPrevious ){
            /*!
             * Set the flow direction in the current configuration of the
             * plastic configuration.
             * 
             * Current formulation is Drucker-Prager based on the potential
             * 
             * \f$ g = \hat{\sigma} + A \bar{\sigma}\f$
             * 
             * where
             * 
             * \f$ \hat{\sigma} = \sqrt{\frac{3}{2} \sigma_{ij}^{\text{dev}} \sigma_{ij}^{\text{dev}} }\f$
             * \f$ \bar{\sigma} = \sigma_{ii}\f$
             * 
             * \f$ \sigma_{ij}^{\text{dev}} = \sigma_{ij} - \frac{1}{3} \bar{\sigma} \delta_{ij}\f$
             * 
             * \param isPrevious: Flag for whether to compute the values in the previous configuration
             */

            const floatVector *drivingStress;

            const floatVector *flowParameters;

            floatType g;

            setDataStorageBase< floatVector > flowDirection;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_previousDrivingStress( ) );

                flowDirection = get_setDataStorage_previousFlowDirection( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_drivingStress( ) );

                flowDirection = get_setDataStorage_flowDirection( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( flowParameters = get_flowParameters( ) );

            floatVector dgdDrivingStress( drivingStress->size( ), 0 );

            flowDirection.zero( drivingStress->size( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::druckerPragerSurface( *drivingStress, ( *flowParameters )[ 1 ], ( *flowParameters )[ 0 ],
                                          g, dgdDrivingStress, *flowDirection.value ) );

        }

        void residual::setFlowDirectionDerivatives( const bool isPrevious ){
            /*!
             * Set the derivatives of the flow direction in the current configuration of the
             * plastic configuration.
             * 
             * Current formulation is Drucker-Prager based on the potential
             * 
             * \f$ g = \hat{\sigma} + A \bar{\sigma}\f$
             * 
             * where
             * 
             * \f$ \hat{\sigma} = \sqrt{\frac{3}{2} \sigma_{ij}^{\text{dev}} \sigma_{ij}^{\text{dev}} }\f$
             * \f$ \bar{\sigma} = \sigma_{ii}\f$
             * 
             * \f$ \sigma_{ij}^{\text{dev}} = \sigma_{ij} - \frac{1}{3} \bar{\sigma} \delta_{ij}\f$
             * 
             * \param isPrevious: Flag for whether to compute the values in the previous configuration
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector *drivingStress;

            const floatVector *flowParameters;

            const floatVector *dDrivingStressdCauchyStress;

            const floatVector *dDrivingStressdF;

            const floatVector *dDrivingStressdSubFs;

            floatType g;

            setDataStorageBase< floatVector > flowDirection;

            setDataStorageBase< fourthOrderTensor > dFlowDirectiondCauchyStress;

            setDataStorageBase< fourthOrderTensor > dFlowDirectiondF;

            setDataStorageBase< floatVector > dFlowDirectiondSubFs;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = get_dPreviousDrivingStressdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = get_dPreviousDrivingStressdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = get_dPreviousDrivingStressdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_previousDrivingStress( ) );

                flowDirection = get_setDataStorage_previousFlowDirection( );

                dFlowDirectiondCauchyStress = get_setDataStorage_dPreviousFlowDirectiondPreviousCauchyStress( );

                dFlowDirectiondF            = get_setDataStorage_dPreviousFlowDirectiondPreviousF( );

                dFlowDirectiondSubFs        = get_setDataStorage_dPreviousFlowDirectiondPreviousSubFs( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = get_dDrivingStressdCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = get_dDrivingStressdF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = get_dDrivingStressdSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_drivingStress( ) );

                flowDirection = get_setDataStorage_flowDirection( );

                dFlowDirectiondCauchyStress = get_setDataStorage_dFlowDirectiondCauchyStress( );

                dFlowDirectiondF            = get_setDataStorage_dFlowDirectiondF( );

                dFlowDirectiondSubFs        = get_setDataStorage_dFlowDirectiondSubFs( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( flowParameters = get_flowParameters( ) );

            floatVector dgdDrivingStress( drivingStress->size( ), 0 );

            flowDirection.zero( sot_dim );

            floatMatrix _dFlowDirectiondDrivingStress;

            floatVector dFlowDirectiondDrivingStress;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::druckerPragerSurface( *drivingStress, ( *flowParameters )[ 1 ], ( *flowParameters )[ 0 ],
                                          g, dgdDrivingStress, *flowDirection.value, _dFlowDirectiondDrivingStress ) );

            dFlowDirectiondDrivingStress = tardigradeVectorTools::appendVectors( _dFlowDirectiondDrivingStress );

            auto map_dFlowDirectiondDrivingStress = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dFlowDirectiondDrivingStress.data( ) );
            auto map_dDrivingStressdCauchyStress  = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dDrivingStressdCauchyStress->data( ) );
            auto map_dDrivingStressdF             = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dDrivingStressdF->data( ) );

            auto map_dDrivingStressdSubFs         = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dDrivingStressdSubFs->data( ), ( num_configs - 1 ) * sot_dim );

            auto map_dFlowDirectiondCauchyStress  = dFlowDirectiondCauchyStress.zeroMap< floatType, sot_dim, sot_dim >( );
            auto map_dFlowDirectiondF             = dFlowDirectiondF.zeroMap< floatType, sot_dim, sot_dim >( );
            auto map_dFlowDirectiondSubFs         = dFlowDirectiondSubFs.zeroMap< floatType, sot_dim >( ( num_configs - 1 ) * sot_dim );

            map_dFlowDirectiondCauchyStress = ( map_dFlowDirectiondDrivingStress * map_dDrivingStressdCauchyStress ).eval( );

            map_dFlowDirectiondF = ( map_dFlowDirectiondDrivingStress * map_dDrivingStressdF ).eval( );

            map_dFlowDirectiondSubFs = ( map_dFlowDirectiondDrivingStress * map_dDrivingStressdSubFs ).eval( );

        }

        void residual::setYieldFunction( ){
            /*!
             * Set the value of the yield function
             */

            setYieldFunction( false );

        }

        void residual::setdYieldFunctiondCauchyStress( const bool isPrevious ){
            /*!
             * Set the value of the derivative of the yield function w.r.t. the Cauchy stress
             * 
             * \param isPrevious: Flag for if the calculation is for the previous value (true) or the current value (false)
             */

            setYieldFunctionDerivatives( isPrevious );

        }

        void residual::setdYieldFunctiondF( const bool isPrevious ){
            /*!
             * Set the value of the derivative of the yield function w.r.t. the deformation gradient
             * 
             * \param isPrevious: Flag for if the calculation is for the previous value (true) or the current value (false)
             */

            setYieldFunctionDerivatives( isPrevious );

        }

        void residual::setdYieldFunctiondSubFs( const bool isPrevious ){
            /*!
             * Set the value of the derivative of the yield function w.r.t. the sub-deformation gradients
             * 
             * \param isPrevious: Flag for if the calculation is for the previous value (true) or the current value (false)
             */

            setYieldFunctionDerivatives( isPrevious );

        }

        void residual::setdYieldFunctiondStateVariables( const bool isPrevious ){
            /*!
             * Set the value of the derivative of the yield function w.r.t. the state variables
             * 
             * \param isPrevious: Flag for if the calculation is for the previous value (true) or the current value (false)
             */

            setYieldFunctionDerivatives( isPrevious );

        }

        void residual::setdYieldFunctiondCauchyStress( ){
            /*!
             * Set the value of the derivative of the yield function w.r.t. the Cauchy stress
             */

            setdYieldFunctiondCauchyStress( false );

        }

        void residual::setdYieldFunctiondF( ){
            /*!
             * Set the value of the derivative of the yield function w.r.t. the deformation gradient
             */

            setdYieldFunctiondF( false );

        }

        void residual::setdYieldFunctiondSubFs( ){
            /*!
             * Set the value of the derivative of the yield function w.r.t. the sub-deformation gradients
             */

            setdYieldFunctiondSubFs( false );

        }

        void residual::setdYieldFunctiondStateVariables( ){
            /*!
             * Set the value of the derivative of the yield function w.r.t. the state variables
             */

            setdYieldFunctiondStateVariables( false );

        }

        void residual::setPreviousYieldFunction( ){
            /*!
             * Set the value of the yield function
             */

            setYieldFunction( true );

        }

        void residual::setdPreviousYieldFunctiondPreviousCauchyStress( ){
            /*!
             * Set the value of the derivative of the previous yield function w.r.t. the previous Cauchy stress
             */

            setdYieldFunctiondCauchyStress( true );

        }

        void residual::setdPreviousYieldFunctiondPreviousF( ){
            /*!
             * Set the value of the derivative of the previous yield function w.r.t. the previous deformation gradient
             */

            setdYieldFunctiondF( true );

        }

        void residual::setdPreviousYieldFunctiondPreviousSubFs( ){
            /*!
             * Set the value of the derivative of the previous yield function w.r.t. the previous sub-deformation gradients
             */

            setdYieldFunctiondSubFs( true );

        }

        void residual::setdPreviousYieldFunctiondPreviousStateVariables( ){
            /*!
             * Set the value of the derivative of the previous yield function w.r.t. the previous state variables
             */

            setdYieldFunctiondStateVariables( true );

        }

        void residual::setYieldFunction( const bool isPrevious ){
            /*!
             * Set the value of the yield function
             * 
             * \param isPrevious: Flag for whether this is the previous timestep
             */

            const floatVector* drivingStress;

            const floatVector* yieldParameters;

            setDataStorageBase< floatType > yieldFunction;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_previousDrivingStress( ) );

                yieldFunction = get_setDataStorage_previousYieldFunction( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_drivingStress( ) );

                yieldFunction = get_setDataStorage_yieldFunction( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( yieldParameters = get_yieldParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::druckerPragerSurface( *drivingStress, ( *yieldParameters )[ 1 ], ( *yieldParameters )[ 0 ], *yieldFunction.value ) );

        }

        void residual::setYieldFunctionDerivatives( const bool isPrevious ){
            /*!
             * Set the value of the yield function derivatives
             * 
             * \param isPrevious: Flag for whether this is the previous timestep
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector* drivingStress;

            const floatVector* dDrivingStressdCauchyStress;

            const floatVector* dDrivingStressdF;

            const floatVector* dDrivingStressdSubFs;

            const floatVector* yieldParameters;

            setDataStorageBase< floatType > yieldFunction;

            setDataStorageBase< secondOrderTensor > dYieldFunctiondCauchyStress;

            setDataStorageBase< secondOrderTensor > dYieldFunctiondF;

            setDataStorageBase< floatVector > dYieldFunctiondSubFs;

            setDataStorageBase< floatVector > dYieldFunctiondStateVariables;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = get_dPreviousDrivingStressdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = get_dPreviousDrivingStressdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = get_dPreviousDrivingStressdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_previousDrivingStress( ) );

                yieldFunction = get_setDataStorage_previousYieldFunction( );

                dYieldFunctiondCauchyStress   = get_setDataStorage_dPreviousYieldFunctiondPreviousCauchyStress( );

                dYieldFunctiondF              = get_setDataStorage_dPreviousYieldFunctiondPreviousF( );

                dYieldFunctiondSubFs          = get_setDataStorage_dPreviousYieldFunctiondPreviousSubFs( );

                dYieldFunctiondStateVariables = get_setDataStorage_dPreviousYieldFunctiondPreviousStateVariables( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = get_dDrivingStressdCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = get_dDrivingStressdF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = get_dDrivingStressdSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_drivingStress( ) );

                yieldFunction = get_setDataStorage_yieldFunction( );

                dYieldFunctiondCauchyStress   = get_setDataStorage_dYieldFunctiondCauchyStress( );

                dYieldFunctiondF              = get_setDataStorage_dYieldFunctiondF( );

                dYieldFunctiondSubFs          = get_setDataStorage_dYieldFunctiondSubFs( );

                dYieldFunctiondStateVariables = get_setDataStorage_dYieldFunctiondStateVariables( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( yieldParameters = get_yieldParameters( ) );

            floatVector dYieldFunctiondDrivingStress( drivingStress->size( ), 0 );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::druckerPragerSurface( *drivingStress, ( *yieldParameters )[ 1 ], ( *yieldParameters )[ 0 ],
                                          *yieldFunction.value, dYieldFunctiondDrivingStress ) );

            auto map_dYieldFunctiondDrivingStress = getFixedSizeMatrixMap< floatType,       1, sot_dim >( dYieldFunctiondDrivingStress.data( ) );
            auto map_dDrivingStressdCauchyStress  = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dDrivingStressdCauchyStress->data( ) );
            auto map_dDrivingStressdF             = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dDrivingStressdF->data( ) );

            auto map_dDrivingStressdSubFs         = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( dDrivingStressdSubFs->data( ), ( num_configs - 1 ) * sot_dim );

            auto map_dYieldFunctiondCauchyStress  = dYieldFunctiondCauchyStress.zeroMap< floatType, 1, sot_dim >( );
            auto map_dYieldFunctiondF             = dYieldFunctiondF.zeroMap< floatType, 1, sot_dim >( );
            auto map_dYieldFunctiondSubFs         = dYieldFunctiondSubFs.zeroMap< floatType, 1 >( ( num_configs - 1 ) * sot_dim );

            map_dYieldFunctiondCauchyStress = ( map_dYieldFunctiondDrivingStress * map_dDrivingStressdCauchyStress ).eval( );

            map_dYieldFunctiondF = ( map_dYieldFunctiondDrivingStress * map_dDrivingStressdF ).eval( );

            map_dYieldFunctiondSubFs = ( map_dYieldFunctiondDrivingStress * map_dDrivingStressdSubFs ).eval( );

            dYieldFunctiondStateVariables.zero( get_stateVariables( )->size( ) );

        }

        void residual::setPlasticThermalMultiplier( ){
            /*!
             * Set the plastic thermal multiplier
             */

            setPlasticThermalMultiplier( false );

        }

        void residual::setdPlasticThermalMultiplierdT( ){
            /*!
             * Set the derivative of the plastic thermal multiplier w.r.t. the temperature
             */

            setdPlasticThermalMultiplierdT( false );

        }

        void residual::setPreviousPlasticThermalMultiplier( ){
            /*!
             * Set the previous plastic thermal multiplier
             */

            setPlasticThermalMultiplier( true );

        }

        void residual::setdPreviousPlasticThermalMultiplierdPreviousT( ){
            /*!
             * Set the derivative of the previous plastic thermal multiplier w.r.t. the pervious temperature
             */

            setdPlasticThermalMultiplierdT( true );

        }

        void residual::setdPlasticThermalMultiplierdT( const bool isPrevious ){
            /*!
             * Set the derivative of the plastic thermal multiplier w.r.t. the temperature
             * 
             * \param isPrevious: A flag for if the derivative is to be computed for the previous (True) or current (False) plastic thermal multiplier
             */

            setPlasticThermalMultiplierDerivatives( isPrevious );

        }

        void residual::setPlasticThermalMultiplier( const bool isPrevious ){
            /*!
             * Set the plastic thermal multiplier
             * 
             * \param isPrevious: A flag for if the values are to be computed for the previous (True) or current (False) plastic thermal multiplier
             */

            const floatType *temperature;

            const floatVector *temperatureParameters;

            setDataStorageBase< floatType > plasticThermalMultiplier;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( temperature = hydra->getPreviousTemperature( ) );

                plasticThermalMultiplier    = get_setDataStorage_previousPlasticThermalMultiplier( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( temperature = hydra->getTemperature( ) );

                plasticThermalMultiplier    = get_setDataStorage_plasticThermalMultiplier( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( temperatureParameters = get_thermalParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::WLF( *temperature, { ( *temperatureParameters )[ 2 ], ( *temperatureParameters )[ 0 ], ( *temperatureParameters )[ 1 ] },
                                          *plasticThermalMultiplier.value ) ); 

        }

        void residual::setPlasticThermalMultiplierDerivatives( const bool isPrevious ){
            /*!
             * Set the derivatives of the plastic thermal multiplier
             * 
             * \param isPrevious: A flag for if the values are to be computed for the previous (True) or current (False) plastic thermal multiplier
             */

            const floatType *temperature;

            const floatVector *temperatureParameters;

            setDataStorageBase< floatType > plasticThermalMultiplier;

            setDataStorageBase< floatType > dPlasticThermalMultiplierdT;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( temperature = hydra->getPreviousTemperature( ) );

                plasticThermalMultiplier    = get_setDataStorage_previousPlasticThermalMultiplier( );

                dPlasticThermalMultiplierdT = get_setDataStorage_dPreviousPlasticThermalMultiplierdPreviousT( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( temperature = hydra->getTemperature( ) );

                plasticThermalMultiplier = get_setDataStorage_plasticThermalMultiplier( );

                dPlasticThermalMultiplierdT = get_setDataStorage_dPlasticThermalMultiplierdT( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( temperatureParameters = get_thermalParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::WLF( *temperature, { ( *temperatureParameters )[ 2 ], ( *temperatureParameters )[ 0 ], ( *temperatureParameters )[ 1 ] },
                                          *plasticThermalMultiplier.value, *dPlasticThermalMultiplierdT.value ) ); 

        }

        void residual::setDragStress( ){
            /*!
             * Set the value of the drag stress
             */

            setDragStress( false );

        }

        void residual::setdDragStressdStateVariables( ){
            /*!
             * Set the value of the derivative of the drag stress w.r.t. the state variables
             */

            setdDragStressdStateVariables( false );

        }

        void residual::setdDragStressdStateVariables( const bool isPrevious ){
            /*!
             * Set the value of the derivative of the drag stress w.r.t. the state variables
             * 
             * \param isPrevious: Flag for whether the derivative should be computed of the previous (true) or current (false) values
             */

            setDragStressDerivatives( isPrevious );

        }

        void residual::setPreviousDragStress( ){
            /*!
             * Set the value of the drag stress
             */

            setDragStress( true );

        }

        void residual::setdPreviousDragStressdPreviousStateVariables( ){
            /*!
             * Set the value of the derivative of the previous drag stress w.r.t. the previous state variables
             */

            setdDragStressdStateVariables( true );

        }

        void residual::setDragStress( const bool isPrevious ){
            /*!
             * Set the value of the drag stress
             * 
             * \param isPrevious: Flag for whether to compute the values for the
             *     previous timestep
             */

            const floatVector *stateVariables;

            const floatVector *dragStressParameters;

            setDataStorageBase< floatType > dragStress;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_previousStateVariables( ) );

                dragStress = get_setDataStorage_previousDragStress( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_stateVariables( ) );

                dragStress = get_setDataStorage_dragStress( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( dragStressParameters = get_dragStressParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::linearHardening( *stateVariables, floatVector( dragStressParameters->begin( ) + 1, dragStressParameters->end( ) ), ( *dragStressParameters )[ 0 ], *dragStress.value ) );

        }

        void residual::setDragStressDerivatives( const bool isPrevious ){
            /*!
             * Set the value of the derivatives of the drag stress
             * 
             * \param isPrevious: Flag for whether to compute the values for the
             *     previous timestep
             */

            const floatVector *stateVariables;

            const floatVector *dragStressParameters;

            setDataStorageBase< floatType > dragStress;

            setDataStorageBase< floatVector > dDragStressdStateVariables;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_previousStateVariables( ) );

                dragStress = get_setDataStorage_previousDragStress( );

                dDragStressdStateVariables = get_setDataStorage_dPreviousDragStressdPreviousStateVariables( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_stateVariables( ) );

                dragStress = get_setDataStorage_dragStress( );

                dDragStressdStateVariables = get_setDataStorage_dDragStressdStateVariables( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( dragStressParameters = get_dragStressParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::linearHardening( *stateVariables, floatVector( dragStressParameters->begin( ) + 1, dragStressParameters->end( ) ), ( *dragStressParameters )[ 0 ], *dragStress.value, *dDragStressdStateVariables.value ) );

        }

        void residual::setHardeningFunction( ){
            /*!
             * Set the value of the hardening function
             */

            setHardeningFunction( false );

        }

        void residual::setdHardeningFunctiondStateVariables( ){
            /*!
             * Set the value of the derivative of the hardening function w.r.t. the state variables
             */

            setdHardeningFunctiondStateVariables( false );

        }

        void residual::setdHardeningFunctiondStateVariables( const bool isPrevious ){
            /*!
             * Set the value of the derivative of the hardening function w.r.t. the state variables
             * 
             * \param isPrevious: Flag for if the derivative should be taken of the previous hardening function (true) or the current one (false)
             */

            setHardeningFunctionDerivatives( isPrevious );

        }

        void residual::setPreviousHardeningFunction( ){
            /*!
             * Set the value of the hardening function
             */

            setHardeningFunction( true );

        }

        void residual::setdPreviousHardeningFunctiondPreviousStateVariables( ){
            /*!
             * Set the value of the derivative of the previous hardening function w.r.t. the previous state variables
             */

            setdHardeningFunctiondStateVariables( true );

        }

        void residual::setHardeningFunction( const bool isPrevious ){
            /*!
             * Set the value of the hardening function
             * 
             * \param &isPrevious: Flag for whether to compute the values for the
             *     previous timestep
             */

            const floatVector *stateVariables;

            const floatVector *hardeningParameters;

            setDataStorageBase< floatVector > hardeningFunction;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_previousStateVariables( ) );

                hardeningFunction = get_setDataStorage_previousHardeningFunction( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_stateVariables( ) );

                hardeningFunction = get_setDataStorage_hardeningFunction( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( hardeningParameters = get_hardeningParameters( ) );

            floatType _hardeningFunction;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::linearHardening( *stateVariables, floatVector( hardeningParameters->begin( ) + 1, hardeningParameters->end( ) ), ( *hardeningParameters )[ 0 ], _hardeningFunction ) );

            *hardeningFunction.value = { _hardeningFunction };

        }

        void residual::setHardeningFunctionDerivatives( const bool isPrevious ){
            /*!
             * Set the value of the derivatives of the hardening function
             * 
             * \param &isPrevious: Flag for whether to compute the values for the
             *     previous timestep
             */

            const floatVector *stateVariables;

            const floatVector *hardeningParameters;

            setDataStorageBase< floatVector > hardeningFunction;

            setDataStorageBase< floatVector > dHardeningFunctiondStateVariables;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_previousStateVariables( ) );

                hardeningFunction = get_setDataStorage_previousHardeningFunction( );

                dHardeningFunctiondStateVariables = get_setDataStorage_dPreviousHardeningFunctiondPreviousStateVariables( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_stateVariables( ) );

                hardeningFunction = get_setDataStorage_hardeningFunction( );

                dHardeningFunctiondStateVariables = get_setDataStorage_dHardeningFunctiondStateVariables( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( hardeningParameters = get_hardeningParameters( ) );

            floatType _hardeningFunction;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::linearHardening( *stateVariables, floatVector( hardeningParameters->begin( ) + 1, hardeningParameters->end( ) ), ( *hardeningParameters )[ 0 ], _hardeningFunction, *dHardeningFunctiondStateVariables.value ) );

            *hardeningFunction.value = { _hardeningFunction };

        }

        void residual::setPlasticMultiplier( ){
            /*!
             * Set the plastic multiplier in the current configuration of the
             * plastic configuration
             */

            setPlasticMultiplier( false );

        }

        void residual::setdPlasticMultiplierdCauchyStress( ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the Cauchy stress
             */

            setdPlasticMultiplierdCauchyStress( false );

        }

        void residual::setdPlasticMultiplierdF( ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the deformation gradient
             */

            setdPlasticMultiplierdF( false );

        }

        void residual::setdPlasticMultiplierdSubFs( ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the sub-deformation gradients
             */

            setdPlasticMultiplierdSubFs( false );

        }

        void residual::setdPlasticMultiplierdT( ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the temperature
             */

            setdPlasticMultiplierdT( false );

        }

        void residual::setdPlasticMultiplierdStateVariables( ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the state variables
             */

            setdPlasticMultiplierdStateVariables( false );

        }

        void residual::setdPlasticMultiplierdCauchyStress( const bool isPrevious ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the Cauchy stress
             * 
             * \param isPrevious: Flag for whether to compute the derivatives of the current (false) or previous (true) value
             */

            setPlasticMultiplierDerivatives( isPrevious );

        }

        void residual::setdPlasticMultiplierdF( const bool isPrevious ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the deformation gradient
             * 
             * \param isPrevious: Flag for whether to compute the derivatives of the current (false) or previous (true) value
             */

            setPlasticMultiplierDerivatives( isPrevious );

        }

        void residual::setdPlasticMultiplierdSubFs( const bool isPrevious ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the sub-deformation gradients
             * 
             * \param isPrevious: Flag for whether to compute the derivatives of the current (false) or previous (true) value
             */

            setPlasticMultiplierDerivatives( isPrevious );

        }

        void residual::setdPlasticMultiplierdT( const bool isPrevious ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the temperature
             * 
             * \param isPrevious: Flag for whether to compute the derivatives of the current (false) or previous (true) value
             */

            setPlasticMultiplierDerivatives( isPrevious );

        }

        void residual::setdPlasticMultiplierdStateVariables( const bool isPrevious ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the state variables
             * 
             * \param isPrevious: Flag for whether to compute the derivatives of the current (false) or previous (true) value
             */

            setPlasticMultiplierDerivatives( isPrevious );

        }

        void residual::setPreviousPlasticMultiplier( ){
            /*!
             * Set the plastic multiplier in the current configuration of the
             * plastic configuration
             */

            setPlasticMultiplier( true );

        }

        void residual::setdPreviousPlasticMultiplierdPreviousCauchyStress( ){
            /*!
             * Set the derivative of the previous plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the previous Cauchy stress
             */

            setdPlasticMultiplierdCauchyStress( true );

        }

        void residual::setdPreviousPlasticMultiplierdPreviousF( ){
            /*!
             * Set the derivative of the previous plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the previous deformation gradient
             */

            setdPlasticMultiplierdF( true );

        }

        void residual::setdPreviousPlasticMultiplierdPreviousSubFs( ){
            /*!
             * Set the derivative of the previous plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the previous sub-deformation gradients
             */

            setdPlasticMultiplierdSubFs( true );

        }

        void residual::setdPreviousPlasticMultiplierdPreviousT( ){
            /*!
             * Set the derivative of the previous plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the previous temperature
             */

            setdPlasticMultiplierdT( true );

        }

        void residual::setdPreviousPlasticMultiplierdPreviousStateVariables( ){
            /*!
             * Set the derivative of the previous plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the previous state variables
             */

            setdPlasticMultiplierdStateVariables( true );

        }

        void residual::setPlasticMultiplier( const bool isPrevious ){
            /*!
             * Set the plastic multiplier in the current configuration of the
             * plastic configuration
             * 
             * \param &isPrevious: Flag for whether to compute the plastic multiplier
             *     in the previous timestep
             */

            const floatType *yieldFunction;

            const floatType *dragStress;

            const floatType *plasticThermalMultiplier;

            const floatVector *perzynaParameters;

            setDataStorageBase< floatType > plasticMultiplier;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( yieldFunction = get_previousYieldFunction( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dragStress = get_previousDragStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticThermalMultiplier = get_previousPlasticThermalMultiplier( ) );

                plasticMultiplier = get_setDataStorage_previousPlasticMultiplier( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( yieldFunction = get_yieldFunction( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dragStress = get_dragStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticThermalMultiplier = get_plasticThermalMultiplier( ) );

                plasticMultiplier = get_setDataStorage_plasticMultiplier( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( perzynaParameters = get_perzynaParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::perzynaModel( *yieldFunction, *dragStress, *plasticThermalMultiplier, ( *perzynaParameters )[ 0 ], *plasticMultiplier.value ) );

        }

        void residual::setPlasticMultiplierDerivatives( const bool isPrevious ){
            /*!
             * Set the derivatives of the plastic multiplier in the current configuration of the
             * plastic configuration
             * 
             * \param &isPrevious: Flag for whether to compute the plastic multiplier
             *     in the previous timestep
             */

            const floatType *yieldFunction;

            const floatType *dragStress;

            const floatType *plasticThermalMultiplier;

            const floatVector *perzynaParameters;

            const floatVector *dYieldFunctiondCauchyStress;

            const floatVector *dYieldFunctiondF;

            const floatVector *dYieldFunctiondSubFs;

            const floatVector *dYieldFunctiondStateVariables;

            const floatVector *dDragStressdStateVariables;

            const floatType   *dPlasticThermalMultiplierdT;

            setDataStorageBase< floatType > plasticMultiplier;

            setDataStorageBase< secondOrderTensor > dPlasticMultiplierdCauchyStress;

            setDataStorageBase< secondOrderTensor > dPlasticMultiplierdF;

            setDataStorageBase< secondOrderTensor > dPlasticMultiplierdSubFs;

            setDataStorageBase< floatType > dPlasticMultiplierdT;

            setDataStorageBase< floatVector > dPlasticMultiplierdStateVariables;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondCauchyStress = get_dPreviousYieldFunctiondPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondF = get_dPreviousYieldFunctiondPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondSubFs = get_dPreviousYieldFunctiondPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondStateVariables = get_dPreviousYieldFunctiondPreviousStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDragStressdStateVariables = get_dPreviousDragStressdPreviousStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticThermalMultiplierdT = get_dPreviousPlasticThermalMultiplierdPreviousT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( yieldFunction = get_previousYieldFunction( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dragStress = get_previousDragStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticThermalMultiplier = get_previousPlasticThermalMultiplier( ) );

                plasticMultiplier                 = get_setDataStorage_previousPlasticMultiplier( );

                dPlasticMultiplierdCauchyStress   = get_setDataStorage_dPreviousPlasticMultiplierdPreviousCauchyStress( );

                dPlasticMultiplierdF              = get_setDataStorage_dPreviousPlasticMultiplierdPreviousF( );

                dPlasticMultiplierdSubFs          = get_setDataStorage_dPreviousPlasticMultiplierdPreviousSubFs( );

                dPlasticMultiplierdT              = get_setDataStorage_dPreviousPlasticMultiplierdPreviousT( );

                dPlasticMultiplierdStateVariables = get_setDataStorage_dPreviousPlasticMultiplierdPreviousStateVariables( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondCauchyStress = get_dYieldFunctiondCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondF = get_dYieldFunctiondF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondSubFs = get_dYieldFunctiondSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondStateVariables = get_dYieldFunctiondStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDragStressdStateVariables = get_dDragStressdStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticThermalMultiplierdT = get_dPlasticThermalMultiplierdT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( yieldFunction = get_yieldFunction( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dragStress = get_dragStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticThermalMultiplier = get_plasticThermalMultiplier( ) );

                plasticMultiplier                 = get_setDataStorage_plasticMultiplier( );

                dPlasticMultiplierdCauchyStress   = get_setDataStorage_dPlasticMultiplierdCauchyStress( );

                dPlasticMultiplierdF              = get_setDataStorage_dPlasticMultiplierdF( );

                dPlasticMultiplierdSubFs          = get_setDataStorage_dPlasticMultiplierdSubFs( );

                dPlasticMultiplierdT              = get_setDataStorage_dPlasticMultiplierdT( );

                dPlasticMultiplierdStateVariables = get_setDataStorage_dPlasticMultiplierdStateVariables( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( perzynaParameters = get_perzynaParameters( ) );

            floatType dPlasticMultiplierdYieldFunction;

            floatType dPlasticMultiplierdDragStress;

            floatType dPlasticMultiplierdPlasticThermalMultiplier;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::perzynaModel( *yieldFunction, *dragStress, *plasticThermalMultiplier, ( *perzynaParameters )[ 0 ], *plasticMultiplier.value, dPlasticMultiplierdYieldFunction, dPlasticMultiplierdDragStress, dPlasticMultiplierdPlasticThermalMultiplier ) );

            *dPlasticMultiplierdCauchyStress.value = dPlasticMultiplierdYieldFunction * ( *dYieldFunctiondCauchyStress );

            *dPlasticMultiplierdF.value = dPlasticMultiplierdYieldFunction * ( *dYieldFunctiondF );

            *dPlasticMultiplierdSubFs.value = dPlasticMultiplierdYieldFunction * ( *dYieldFunctiondSubFs );

            *dPlasticMultiplierdT.value = dPlasticMultiplierdPlasticThermalMultiplier * ( *dPlasticThermalMultiplierdT );

            *dPlasticMultiplierdStateVariables.value = dPlasticMultiplierdYieldFunction * ( *dYieldFunctiondStateVariables )
                                                     + dPlasticMultiplierdDragStress    * ( *dDragStressdStateVariables );

        }

        void residual::setVelocityGradient( ){
            /*!
             * Set the velocity gradient in the current configuration of the plastic
             * configuration
             */

            setVelocityGradient( false );

        }

        void residual::setdVelocityGradientdCauchyStress( ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the Cauchy stress
             */

            setdVelocityGradientdCauchyStress( false );

        }

        void residual::setdVelocityGradientdF( ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the deformation gradient
             */

            setdVelocityGradientdF( false );

        }

        void residual::setdVelocityGradientdSubFs( ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the sub-deformation gradients
             */

            setdVelocityGradientdSubFs( false );

        }

        void residual::setdVelocityGradientdT( ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the temperature
             */

            setdVelocityGradientdT( false );

        }

        void residual::setdVelocityGradientdStateVariables( ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the state variables
             */

            setdVelocityGradientdStateVariables( false );

        }

        void residual::setdVelocityGradientdCauchyStress( const bool isPrevious ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the Cauchy stress
             * 
             * \param isPrevious: A flag for if the gradient is computed of the previous (true) or current (false) value
             */

            setVelocityGradientDerivatives( isPrevious );

        }

        void residual::setdVelocityGradientdF( const bool isPrevious ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the deformation gradient
             * 
             * \param isPrevious: A flag for if the gradient is computed of the previous (true) or current (false) value
             */

            setVelocityGradientDerivatives( isPrevious );

        }

        void residual::setdVelocityGradientdSubFs( const bool isPrevious ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the sub-deformation gradients
             * 
             * \param isPrevious: A flag for if the gradient is computed of the previous (true) or current (false) value
             */

            setVelocityGradientDerivatives( isPrevious );

        }

        void residual::setdVelocityGradientdT( const bool isPrevious ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the temperature
             * 
             * \param isPrevious: A flag for if the gradient is computed of the previous (true) or current (false) value
             */

            setVelocityGradientDerivatives( isPrevious );

        }

        void residual::setdVelocityGradientdStateVariables( const bool isPrevious ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the state variables
             * 
             * \param isPrevious: A flag for if the gradient is computed of the previous (true) or current (false) value
             */

            setVelocityGradientDerivatives( isPrevious );

        }

        void residual::setPreviousVelocityGradient( ){
            /*!
             * Set the velocity gradient in the current configuration of the plastic
             * configuration
             */

            setVelocityGradient( true );

        }

        void residual::setdPreviousVelocityGradientdPreviousCauchyStress( ){
            /*!
             * Set the derivative of the previous velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the previous Cauchy stress
             */

            setdVelocityGradientdCauchyStress( true );

        }

        void residual::setdPreviousVelocityGradientdPreviousF( ){
            /*!
             * Set the derivative of the previous velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the previous deformation gradient
             */

            setdVelocityGradientdF( true );

        }

        void residual::setdPreviousVelocityGradientdPreviousSubFs( ){
            /*!
             * Set the derivative of the previous velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the previous sub-deformation gradients
             */

            setdVelocityGradientdSubFs( true );

        }

        void residual::setdPreviousVelocityGradientdPreviousT( ){
            /*!
             * Set the derivative of the previous velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the previous temperature
             */

            setdVelocityGradientdT( true );

        }

        void residual::setdPreviousVelocityGradientdPreviousStateVariables( ){
            /*!
             * Set the derivative of the previous velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the previous state variables
             */

            setdVelocityGradientdStateVariables( true );

        }

        void residual::setVelocityGradient( const bool isPrevious ){
            /*!
             * Set the velocity gradient in the current configuration of the plastic
             * configuration
             * 
             * \param isPrevious: Flag for whether to compute the value at the previous
             *     timestep.
             */

            const floatType *plasticMultiplier;

            const floatVector *flowDirection;

            setDataStorageBase< secondOrderTensor > velocityGradient;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = get_previousPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( flowDirection = get_previousFlowDirection( ) );

                velocityGradient = get_setDataStorage_previousVelocityGradient( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = get_plasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( flowDirection = get_flowDirection( ) );

                velocityGradient = get_setDataStorage_velocityGradient( );

            }

            *velocityGradient.value = ( *plasticMultiplier ) * ( *flowDirection );

        }

        void residual::setVelocityGradientDerivatives( const bool isPrevious ){
            /*!
             * Set the derivatives of the velocity gradient in the current configuration of the plastic
             * configuration
             * 
             * \param isPrevious: Flag for whether to compute the value at the previous
             *     timestep.
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatType *plasticMultiplier;

            const floatVector *dPlasticMultiplierdCauchyStress;

            const floatVector *dPlasticMultiplierdF;

            const floatVector *dPlasticMultiplierdSubFs;

            const floatType   *dPlasticMultiplierdT;

            const floatVector *dPlasticMultiplierdStateVariables;

            const floatVector *flowDirection;

            const floatVector *dFlowDirectiondCauchyStress;

            const floatVector *dFlowDirectiondF;

            const floatVector *dFlowDirectiondSubFs;

            setDataStorageBase< secondOrderTensor > velocityGradient;

            setDataStorageBase< fourthOrderTensor > dVelocityGradientdCauchyStress;

            setDataStorageBase< fourthOrderTensor > dVelocityGradientdF;

            setDataStorageBase< floatVector > dVelocityGradientdSubFs;

            setDataStorageBase< secondOrderTensor > dVelocityGradientdT;

            setDataStorageBase< floatVector > dVelocityGradientdStateVariables;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdCauchyStress   = get_dPreviousPlasticMultiplierdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdF              = get_dPreviousPlasticMultiplierdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdSubFs          = get_dPreviousPlasticMultiplierdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdT              = get_dPreviousPlasticMultiplierdPreviousT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdStateVariables = get_dPreviousPlasticMultiplierdPreviousStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dFlowDirectiondCauchyStress       = get_dPreviousFlowDirectiondPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dFlowDirectiondF                  = get_dPreviousFlowDirectiondPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dFlowDirectiondSubFs              = get_dPreviousFlowDirectiondPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = get_previousPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( flowDirection = get_previousFlowDirection( ) );

                velocityGradient = get_setDataStorage_previousVelocityGradient( );

                dVelocityGradientdCauchyStress   = get_setDataStorage_dPreviousVelocityGradientdPreviousCauchyStress( );

                dVelocityGradientdF              = get_setDataStorage_dPreviousVelocityGradientdPreviousF( );

                dVelocityGradientdSubFs          = get_setDataStorage_dPreviousVelocityGradientdPreviousSubFs( );

                dVelocityGradientdT              = get_setDataStorage_dPreviousVelocityGradientdPreviousT( );

                dVelocityGradientdStateVariables = get_setDataStorage_dPreviousVelocityGradientdPreviousStateVariables( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdCauchyStress   = get_dPlasticMultiplierdCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdF              = get_dPlasticMultiplierdF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdSubFs          = get_dPlasticMultiplierdSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdT              = get_dPlasticMultiplierdT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdStateVariables = get_dPlasticMultiplierdStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dFlowDirectiondCauchyStress       = get_dFlowDirectiondCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dFlowDirectiondF                  = get_dFlowDirectiondF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dFlowDirectiondSubFs              = get_dFlowDirectiondSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = get_plasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( flowDirection = get_flowDirection( ) );

                velocityGradient = get_setDataStorage_velocityGradient( );

                dVelocityGradientdCauchyStress   = get_setDataStorage_dVelocityGradientdCauchyStress( );

                dVelocityGradientdF              = get_setDataStorage_dVelocityGradientdF( );

                dVelocityGradientdSubFs          = get_setDataStorage_dVelocityGradientdSubFs( );

                dVelocityGradientdT              = get_setDataStorage_dVelocityGradientdT( );

                dVelocityGradientdStateVariables = get_setDataStorage_dVelocityGradientdStateVariables( );
            }

            auto map_flowDirection                     = getFixedSizeVectorMap< floatType, sot_dim >( flowDirection->data( ) );
            auto map_dPlasticMultiplierdCauchyStress   = getFixedSizeMatrixMap< floatType, 1, sot_dim >( dPlasticMultiplierdCauchyStress->data( ) );
            auto map_dPlasticMultiplierdF              = getFixedSizeMatrixMap< floatType, 1, sot_dim >( dPlasticMultiplierdF->data( ) );

            auto map_dPlasticMultiplierdSubFs          = getDynamicColumnSizeMatrixMap< floatType, 1 >( dPlasticMultiplierdSubFs->data( ), ( num_configs - 1 ) * sot_dim );
            auto map_dPlasticMultiplierdStateVariables = getDynamicColumnSizeMatrixMap< floatType, 1 >( dPlasticMultiplierdStateVariables->data( ), dPlasticMultiplierdStateVariables->size( ) );

            auto map_dVelocityGradientdCauchyStress    = dVelocityGradientdCauchyStress.zeroMap< floatType, sot_dim, sot_dim >( );
            auto map_dVelocityGradientdF               = dVelocityGradientdF.zeroMap< floatType, sot_dim, sot_dim >( );
            auto map_dVelocityGradientdSubFs           = dVelocityGradientdSubFs.zeroMap< floatType, sot_dim >( ( num_configs - 1 ) * sot_dim );
            auto map_dVelocityGradientdStateVariables  = dVelocityGradientdStateVariables.zeroMap< floatType, sot_dim >( dPlasticMultiplierdStateVariables->size( ) );

            *velocityGradient.value = ( *plasticMultiplier ) * ( *flowDirection );

            map_dVelocityGradientdCauchyStress     = ( map_flowDirection * map_dPlasticMultiplierdCauchyStress ).eval( );
            *dVelocityGradientdCauchyStress.value += ( *plasticMultiplier ) * ( *dFlowDirectiondCauchyStress );

            map_dVelocityGradientdF     = ( map_flowDirection * map_dPlasticMultiplierdF ).eval( );
            *dVelocityGradientdF.value += ( *plasticMultiplier ) * ( *dFlowDirectiondF );

            map_dVelocityGradientdSubFs     = ( map_flowDirection * map_dPlasticMultiplierdSubFs ).eval( );
            *dVelocityGradientdSubFs.value += ( *plasticMultiplier ) * ( *dFlowDirectiondSubFs );

            *dVelocityGradientdT.value = ( *flowDirection ) * ( *dPlasticMultiplierdT );

            map_dVelocityGradientdStateVariables = ( map_flowDirection * map_dPlasticMultiplierdStateVariables ).eval( );

        }

        void residual::setStateVariableEvolutionRates( ){
            /*! 
             * Set the value of the state variable evolution rates
             */

            setStateVariableEvolutionRates( false );

        }

        void residual::setdStateVariableEvolutionRatesdCauchyStress( ){
            /*! 
             * Set the value of the derivative of the state variable evolution rates w.r.t. the Cauchy stress
             */

            setdStateVariableEvolutionRatesdCauchyStress( false );

        }

        void residual::setdStateVariableEvolutionRatesdF( ){
            /*! 
             * Set the value of the derivative of the state variable evolution rates w.r.t. the deformation gradient
             */

            setdStateVariableEvolutionRatesdF( false );

        }

        void residual::setdStateVariableEvolutionRatesdSubFs( ){
            /*! 
             * Set the value of the derivative of the state variable evolution rates w.r.t. the sub-deformation gradients
             */

            setdStateVariableEvolutionRatesdSubFs( false );

        }

        void residual::setdStateVariableEvolutionRatesdT( ){
            /*! 
             * Set the value of the derivative of the state variable evolution rates w.r.t. the temperature
             */

            setdStateVariableEvolutionRatesdT( false );

        }

        void residual::setdStateVariableEvolutionRatesdStateVariables( ){
            /*! 
             * Set the value of the derivative of the state variable evolution rates w.r.t. the state variables
             */

            setdStateVariableEvolutionRatesdStateVariables( false );

        }

        void residual::setdStateVariableEvolutionRatesdCauchyStress( const bool isPrevious ){
            /*! 
             * Set the value of the derivative of the state variable evolution rates w.r.t. the Cauchy stress
             * 
             * \param isPrevious: Flag for whether the gradients are of the current (false) or previous (true) value 
             */

            setStateVariableEvolutionRateDerivatives( isPrevious );

        }

        void residual::setdStateVariableEvolutionRatesdF( const bool isPrevious ){
            /*! 
             * Set the value of the derivative of the state variable evolution rates w.r.t. the deformation gradient
             * 
             * \param isPrevious: Flag for whether the gradients are of the current (false) or previous (true) value 
             */

            setStateVariableEvolutionRateDerivatives( isPrevious );

        }

        void residual::setdStateVariableEvolutionRatesdSubFs( const bool isPrevious ){
            /*! 
             * Set the value of the derivative of the state variable evolution rates w.r.t. the sub-deformation gradients
             * 
             * \param isPrevious: Flag for whether the gradients are of the current (false) or previous (true) value 
             */

            setStateVariableEvolutionRateDerivatives( isPrevious );

        }

        void residual::setdStateVariableEvolutionRatesdT( const bool isPrevious ){
            /*! 
             * Set the value of the derivative of the state variable evolution rates w.r.t. the temperature
             * 
             * \param isPrevious: Flag for whether the gradients are of the current (false) or previous (true) value 
             */

            setStateVariableEvolutionRateDerivatives( isPrevious );

        }

        void residual::setdStateVariableEvolutionRatesdStateVariables( const bool isPrevious ){
            /*! 
             * Set the value of the derivative of the state variable evolution rates w.r.t. the state variables
             * 
             * \param isPrevious: Flag for whether the gradients are of the current (false) or previous (true) value 
             */

            setStateVariableEvolutionRateDerivatives( isPrevious );

        }

        void residual::setPreviousStateVariableEvolutionRates( ){
            /*! 
             * Set the value of the state variable evolution rates
             */

            setStateVariableEvolutionRates( true );

        }

        void residual::setdPreviousStateVariableEvolutionRatesdPreviousCauchyStress( ){
            /*! 
             * Set the value of the derivative of the previous state variable evolution rates w.r.t. the previous Cauchy stress
             */

            setdStateVariableEvolutionRatesdCauchyStress( true );

        }

        void residual::setdPreviousStateVariableEvolutionRatesdPreviousF( ){
            /*! 
             * Set the value of the derivative of the previous state variable evolution rates w.r.t. the previous deformation gradient
             */

            setdStateVariableEvolutionRatesdF( true );

        }

        void residual::setdPreviousStateVariableEvolutionRatesdPreviousSubFs( ){
            /*! 
             * Set the value of the previous derivative of the state variable evolution rates w.r.t. the previous sub-deformation gradients
             */

            setdStateVariableEvolutionRatesdSubFs( true );

        }

        void residual::setdPreviousStateVariableEvolutionRatesdPreviousT( ){
            /*! 
             * Set the value of the derivative of the previous state variable evolution rates w.r.t. the previous temperature
             */

            setdStateVariableEvolutionRatesdT( true );

        }

        void residual::setdPreviousStateVariableEvolutionRatesdPreviousStateVariables( ){
            /*! 
             * Set the value of the derivative of the previous state variable evolution rates w.r.t. the previous state variables
             */

            setdStateVariableEvolutionRatesdStateVariables( true );

        }

        void residual::setStateVariableEvolutionRates( const bool isPrevious ){
            /*!
             * Set the value of the state variable evolution rates
             * 
             * \param isPrevious: A flag to indicate if the previous evolution rate
             *     should be computed.
             */

            const floatType *plasticMultiplier;

            const floatVector *hardeningFunction;

            setDataStorageBase< floatVector > stateVariableEvolutionRates;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = get_previousPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( hardeningFunction = get_previousHardeningFunction( ) );

                stateVariableEvolutionRates = get_setDataStorage_previousStateVariableEvolutionRates( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = get_plasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( hardeningFunction = get_hardeningFunction( ) );

                stateVariableEvolutionRates = get_setDataStorage_stateVariableEvolutionRates( );

            }

            *stateVariableEvolutionRates.value = ( *plasticMultiplier ) * ( *hardeningFunction );

        }

        void residual::setStateVariableEvolutionRateDerivatives( const bool isPrevious ){
            /*!
             * Set the value of the derivatives of the state variable evolution rates
             * 
             * \param isPrevious: A flag to indicate if the previous evolution rate
             *     should be computed.
             */

            const unsigned int dim = hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int num_stateVariables = get_stateVariables( )->size( );

            const floatType *plasticMultiplier;

            const floatVector *hardeningFunction;

            const floatVector *dPlasticMultiplierdCauchyStress;

            const floatVector *dPlasticMultiplierdF;

            const floatVector *dPlasticMultiplierdSubFs;

            const floatType   *dPlasticMultiplierdT;

            const floatVector *dPlasticMultiplierdStateVariables;

            const floatVector *dHardeningFunctiondStateVariables;

            setDataStorageBase< floatVector > stateVariableEvolutionRates;

            setDataStorageBase< floatVector > dStateVariableEvolutionRatesdCauchyStress;

            setDataStorageBase< floatVector > dStateVariableEvolutionRatesdF;

            setDataStorageBase< floatVector > dStateVariableEvolutionRatesdSubFs;

            setDataStorageBase< floatVector > dStateVariableEvolutionRatesdT;

            setDataStorageBase< floatVector > dStateVariableEvolutionRatesdStateVariables;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdCauchyStress = get_dPreviousPlasticMultiplierdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdF = get_dPreviousPlasticMultiplierdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdSubFs = get_dPreviousPlasticMultiplierdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdT = get_dPreviousPlasticMultiplierdPreviousT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdStateVariables = get_dPreviousPlasticMultiplierdPreviousStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dHardeningFunctiondStateVariables = get_dPreviousHardeningFunctiondPreviousStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = get_previousPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( hardeningFunction = get_previousHardeningFunction( ) );

                stateVariableEvolutionRates = get_setDataStorage_previousStateVariableEvolutionRates( );

                dStateVariableEvolutionRatesdCauchyStress   = get_setDataStorage_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( );

                dStateVariableEvolutionRatesdF              = get_setDataStorage_dPreviousStateVariableEvolutionRatesdPreviousF( );

                dStateVariableEvolutionRatesdSubFs          = get_setDataStorage_dPreviousStateVariableEvolutionRatesdPreviousSubFs( );

                dStateVariableEvolutionRatesdT              = get_setDataStorage_dPreviousStateVariableEvolutionRatesdPreviousT( );

                dStateVariableEvolutionRatesdStateVariables = get_setDataStorage_dPreviousStateVariableEvolutionRatesdPreviousStateVariables( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdCauchyStress = get_dPlasticMultiplierdCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdF = get_dPlasticMultiplierdF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdSubFs = get_dPlasticMultiplierdSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdT = get_dPlasticMultiplierdT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdStateVariables = get_dPlasticMultiplierdStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dHardeningFunctiondStateVariables = get_dHardeningFunctiondStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = get_plasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( hardeningFunction = get_hardeningFunction( ) );

                stateVariableEvolutionRates = get_setDataStorage_stateVariableEvolutionRates( );

                dStateVariableEvolutionRatesdCauchyStress   = get_setDataStorage_dStateVariableEvolutionRatesdCauchyStress( );

                dStateVariableEvolutionRatesdF              = get_setDataStorage_dStateVariableEvolutionRatesdF( );

                dStateVariableEvolutionRatesdSubFs          = get_setDataStorage_dStateVariableEvolutionRatesdSubFs( );

                dStateVariableEvolutionRatesdT              = get_setDataStorage_dStateVariableEvolutionRatesdT( );

                dStateVariableEvolutionRatesdStateVariables = get_setDataStorage_dStateVariableEvolutionRatesdStateVariables( );

            }

            *stateVariableEvolutionRates.value = ( *plasticMultiplier ) * ( *hardeningFunction );

            dStateVariableEvolutionRatesdCauchyStress.zero( num_stateVariables * sot_dim );

            dStateVariableEvolutionRatesdF.zero( num_stateVariables * sot_dim );

            dStateVariableEvolutionRatesdSubFs.zero( num_stateVariables * ( num_configs - 1 ) * sot_dim );

            dStateVariableEvolutionRatesdT.zero( num_stateVariables );

            dStateVariableEvolutionRatesdStateVariables.zero( num_stateVariables * num_stateVariables );

            const unsigned int iub = hardeningFunction->size( );
            for ( unsigned int i = 0; i < iub; i++ ){
                for ( unsigned int j = 0; j < sot_dim; j++ ){
                    ( *dStateVariableEvolutionRatesdCauchyStress.value )[ sot_dim * i + j ] += ( *hardeningFunction )[ i ] * ( *dPlasticMultiplierdCauchyStress )[ j ];
                    ( *dStateVariableEvolutionRatesdF.value )[ sot_dim * i + j ] += ( *hardeningFunction )[ i ] * ( *dPlasticMultiplierdF )[ j ];
                }
                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){
                    ( *dStateVariableEvolutionRatesdSubFs.value )[  ( num_configs - 1 ) * sot_dim * i + j ] += ( *dPlasticMultiplierdSubFs )[ j ] * ( *hardeningFunction )[ i ];
                }

                ( *dStateVariableEvolutionRatesdT.value )[ i ] = ( *dPlasticMultiplierdT ) * ( *hardeningFunction )[ i ];
            }

            for ( unsigned int i = 0; i < iub; i++ ){
                for ( unsigned int j = 0; j < num_stateVariables; j++ ){
                    ( *dStateVariableEvolutionRatesdStateVariables.value )[ num_stateVariables * i + j ] += ( *plasticMultiplier ) * ( *dHardeningFunctiondStateVariables )[ num_stateVariables * i + j ];
                    ( *dStateVariableEvolutionRatesdStateVariables.value )[ num_stateVariables * i + j ] += ( *hardeningFunction )[ i ] * ( *dPlasticMultiplierdStateVariables )[ j ];
                }
            }

        }

        void residual::setPlasticDeformationGradient( ){
            /*!
             * Set the plastic deformation gradient
             */

            const floatVector *velocityGradient;

            const floatVector *previousVelocityGradient;

            floatVector previousPlasticDeformationGradient;

            floatVector dFp;

            auto plasticDeformationGradient = get_setDataStorage_plasticDeformationGradient( );

            TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradient = get_velocityGradient( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousVelocityGradient = get_previousVelocityGradient( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousPlasticDeformationGradient = hydra->getPreviousConfiguration( *getPlasticConfigurationIndex( ) ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::evolveF( *hydra->getDeltaTime( ), previousPlasticDeformationGradient, *previousVelocityGradient, *velocityGradient, dFp,
                                          *plasticDeformationGradient.value, 1 - ( *getIntegrationParameter( ) ), 1 ) );

        }

        void residual::setdPlasticDeformationGradientdCauchyStress( ){
            /*!
             * set the derivative of the plastic deformation gradient w.r.t. the Cauchy stress
             */

            setPlasticDeformationGradientDerivatives( false );

        }

        void residual::setdPlasticDeformationGradientdF( ){
            /*!
             * set the derivative of the plastic deformation gradient w.r.t. the deformation gradient
             */

            setPlasticDeformationGradientDerivatives( false );

        }

        void residual::setdPlasticDeformationGradientdSubFs( ){
            /*!
             * set the derivative of the plastic deformation gradient w.r.t. the sub-deformation gradients
             */

            setPlasticDeformationGradientDerivatives( false );

        }

        void residual::setdPlasticDeformationGradientdT( ){
            /*!
             * set the derivative of the plastic deformation gradient w.r.t. the temperature
             */

            setPlasticDeformationGradientDerivatives( false );

        }

        void residual::setdPlasticDeformationGradientdStateVariables( ){
            /*!
             * set the derivative of the plastic deformation gradient w.r.t. the state variables
             */

            setPlasticDeformationGradientDerivatives( false );

        }

        void residual::setdPlasticDeformationGradientdPreviousCauchyStress( ){
            /*!
             * set the derivative of the plastic deformation gradient w.r.t. the previous Cauchy stress
             */

            setPlasticDeformationGradientDerivatives( true );

        }

        void residual::setdPlasticDeformationGradientdPreviousF( ){
            /*!
             * set the derivative of the plastic deformation gradient w.r.t. the previous deformation gradient
             */

            setPlasticDeformationGradientDerivatives( true );

        }

        void residual::setdPlasticDeformationGradientdPreviousSubFs( ){
            /*!
             * set the derivative of the plastic deformation gradient w.r.t. the previous sub-deformation gradients
             */

            setPlasticDeformationGradientDerivatives( true );

        }

        void residual::setdPlasticDeformationGradientdPreviousT( ){
            /*!
             * set the derivative of the plastic deformation gradient w.r.t. the previous temperature
             */

            setPlasticDeformationGradientDerivatives( true );

        }

        void residual::setdPlasticDeformationGradientdPreviousStateVariables( ){
            /*!
             * set the derivative of the plastic deformation gradient w.r.t. the previous state variables
             */

            setPlasticDeformationGradientDerivatives( true );

        }

        void residual::setPlasticDeformationGradientDerivatives( const bool setPreviousDerivatives ){
            /*!
             * Set the plastic deformation gradient derivatives
             * 
             * \param setPreviousDerivatives: Flag for if the previous derivatives should be set
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int plastic_config_index = *getPlasticConfigurationIndex( );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            const secondOrderTensor *velocityGradient;

            const secondOrderTensor *previousVelocityGradient;

            secondOrderTensor previousPlasticDeformationGradient;

            secondOrderTensor dFp;

            auto plasticDeformationGradient = get_setDataStorage_plasticDeformationGradient( );

            secondOrderTensor dFdL;

            TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradient = get_velocityGradient( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousVelocityGradient = get_previousVelocityGradient( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousPlasticDeformationGradient = hydra->getPreviousConfiguration( *getPlasticConfigurationIndex( ) ) );

            if ( setPreviousDerivatives ){

                fourthOrderTensor ddFdPreviousF;

                fourthOrderTensor dFdPreviousF;

                fourthOrderTensor dFdPreviousL;

                TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::evolveFFlatJ( *hydra->getDeltaTime( ), previousPlasticDeformationGradient, *previousVelocityGradient, *velocityGradient, dFp,
                                              *plasticDeformationGradient.value, dFdL, ddFdPreviousF, dFdPreviousF, dFdPreviousL, 1 - ( *getIntegrationParameter( ) ), 1 ) );

                auto map_dFdPreviousL = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dFdPreviousL.data( ) );

                auto map_dPreviousVelocityGradientdPreviousCauchyStress   = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_dPreviousVelocityGradientdPreviousCauchyStress( )->data( ) );
                auto map_dPreviousVelocityGradientdPreviousF              = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_dPreviousVelocityGradientdPreviousF( )->data( ) );

                auto map_dPreviousVelocityGradientdPreviousT              = getFixedSizeVectorMap< floatType, sot_dim >( get_dPreviousVelocityGradientdPreviousT( )->data( ) );

                auto map_dPreviousVelocityGradientdPreviousSubFs          = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( get_dPreviousVelocityGradientdPreviousSubFs( )->data( ), ( num_configs - 1 ) * sot_dim );
                auto map_dPreviousVelocityGradientdPreviousStateVariables = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( get_dPreviousVelocityGradientdPreviousStateVariables( )->data( ), num_isvs );

                auto dPlasticDeformationGradientdPreviousCauchyStress = get_setDataStorage_dPlasticDeformationGradientdPreviousCauchyStress( );
                auto map_dPlasticDeformationGradientdPreviousCauchyStress = dPlasticDeformationGradientdPreviousCauchyStress.zeroMap< floatType, sot_dim, sot_dim >( );

                auto dPlasticDeformationGradientdPreviousF = get_setDataStorage_dPlasticDeformationGradientdPreviousF( );
                auto map_dPlasticDeformationGradientdPreviousF = dPlasticDeformationGradientdPreviousF.zeroMap< floatType, sot_dim, sot_dim >( );

                auto dPlasticDeformationGradientdPreviousSubFs = get_setDataStorage_dPlasticDeformationGradientdPreviousSubFs( );
                auto map_dPlasticDeformationGradientdPreviousSubFs = dPlasticDeformationGradientdPreviousSubFs.zeroMap< floatType, sot_dim >( ( num_configs - 1 ) * sot_dim );

                auto dPlasticDeformationGradientdPreviousT = get_setDataStorage_dPlasticDeformationGradientdPreviousT( );
                auto map_dPlasticDeformationGradientdPreviousT = dPlasticDeformationGradientdPreviousT.zeroMap< floatType, sot_dim >( );

                auto dPlasticDeformationGradientdPreviousStateVariables = get_setDataStorage_dPlasticDeformationGradientdPreviousStateVariables( );
                auto map_dPlasticDeformationGradientdPreviousStateVariables = dPlasticDeformationGradientdPreviousStateVariables.zeroMap< floatType, sot_dim >( num_isvs );

                map_dPlasticDeformationGradientdPreviousCauchyStress = ( map_dFdPreviousL * map_dPreviousVelocityGradientdPreviousCauchyStress ).eval( );

                map_dPlasticDeformationGradientdPreviousF            = ( map_dFdPreviousL * map_dPreviousVelocityGradientdPreviousF ).eval( );

                map_dPlasticDeformationGradientdPreviousSubFs        = ( map_dFdPreviousL * map_dPreviousVelocityGradientdPreviousSubFs ).eval( );

                for ( unsigned int i = 0; i < sot_dim; i++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        ( *dPlasticDeformationGradientdPreviousSubFs.value )[ ( num_configs - 1 ) * sot_dim * i + j + ( plastic_config_index - 1 ) * sot_dim ] += dFdPreviousF[ sot_dim * i + j ];

                    }

                }

                map_dPlasticDeformationGradientdPreviousT = ( map_dFdPreviousL * map_dPreviousVelocityGradientdPreviousT ).eval( );

                map_dPlasticDeformationGradientdPreviousStateVariables = ( map_dFdPreviousL * map_dPreviousVelocityGradientdPreviousStateVariables ).eval( );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::evolveFFlatJ( *hydra->getDeltaTime( ), previousPlasticDeformationGradient, *previousVelocityGradient, *velocityGradient, dFp,
                                              *plasticDeformationGradient.value, dFdL, 1 - ( *getIntegrationParameter( ) ), 1 ) );

            }

            auto map_dFdL = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( dFdL.data( ) );

            auto map_dVelocityGradientdCauchyStress   = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_dVelocityGradientdCauchyStress( )->data( ) );
            auto map_dVelocityGradientdF              = getFixedSizeMatrixMap< floatType, sot_dim, sot_dim >( get_dVelocityGradientdF( )->data( ) );

            auto map_dVelocityGradientdT              = getFixedSizeVectorMap< floatType, sot_dim >( get_dVelocityGradientdT( )->data( ) );

            auto map_dVelocityGradientdSubFs          = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( get_dVelocityGradientdSubFs( )->data( ), ( num_configs - 1 ) * sot_dim );
            auto map_dVelocityGradientdStateVariables = getDynamicColumnSizeMatrixMap< floatType, sot_dim >( get_dVelocityGradientdStateVariables( )->data( ), num_isvs );

            auto dPlasticDeformationGradientdCauchyStress = get_setDataStorage_dPlasticDeformationGradientdCauchyStress( );
            auto map_dPlasticDeformationGradientdCauchyStress = dPlasticDeformationGradientdCauchyStress.zeroMap< floatType, sot_dim, sot_dim >( );

            auto dPlasticDeformationGradientdF = get_setDataStorage_dPlasticDeformationGradientdF( );
            auto map_dPlasticDeformationGradientdF = dPlasticDeformationGradientdF.zeroMap< floatType, sot_dim, sot_dim >( );

            auto dPlasticDeformationGradientdSubFs = get_setDataStorage_dPlasticDeformationGradientdSubFs( );
            auto map_dPlasticDeformationGradientdSubFs = dPlasticDeformationGradientdSubFs.zeroMap< floatType, sot_dim >( ( num_configs - 1 ) * sot_dim );

            auto dPlasticDeformationGradientdT = get_setDataStorage_dPlasticDeformationGradientdT( );
            auto map_dPlasticDeformationGradientdT = dPlasticDeformationGradientdT.zeroMap< floatType, sot_dim >( );

            auto dPlasticDeformationGradientdStateVariables = get_setDataStorage_dPlasticDeformationGradientdStateVariables( );
            auto map_dPlasticDeformationGradientdStateVariables = dPlasticDeformationGradientdStateVariables.zeroMap< floatType, sot_dim >( num_isvs );

            map_dPlasticDeformationGradientdCauchyStress   = ( map_dFdL * map_dVelocityGradientdCauchyStress ).eval( );

            map_dPlasticDeformationGradientdF              = ( map_dFdL * map_dVelocityGradientdF ).eval( );

            map_dPlasticDeformationGradientdSubFs          = ( map_dFdL * map_dVelocityGradientdSubFs ).eval( );

            map_dPlasticDeformationGradientdT              = ( map_dFdL * map_dVelocityGradientdT ).eval( );

            map_dPlasticDeformationGradientdStateVariables = ( map_dFdL * map_dVelocityGradientdStateVariables ).eval( );

        }

        void residual::setPlasticStateVariables( ){
            /*!
             * Set the plastic state variables
             */

            const floatVector *stateVariableEvolutionRates;

            const floatVector *previousStateVariableEvolutionRates;

            const floatVector *previousStateVariables;

            floatVector deltaPlasticStateVariables;

            auto plasticStateVariables = get_setDataStorage_plasticStateVariables( );

            TARDIGRADE_ERROR_TOOLS_CATCH( stateVariableEvolutionRates = get_stateVariableEvolutionRates( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousStateVariableEvolutionRates = get_previousStateVariableEvolutionRates( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousStateVariables = get_previousStateVariables( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::midpointEvolution( *hydra->getDeltaTime( ), *previousStateVariables, *previousStateVariableEvolutionRates, *stateVariableEvolutionRates, deltaPlasticStateVariables, *plasticStateVariables.value, ( 1 - *getIntegrationParameter( ) ) ) );

        }

        void residual::setdPlasticStateVariablesdCauchyStress( ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the Cauchy stress
             */

            setPlasticStateVariableDerivatives( false );

        }

        void residual::setdPlasticStateVariablesdF( ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the deformation gradient
             */

            setPlasticStateVariableDerivatives( false );

        }

        void residual::setdPlasticStateVariablesdSubFs( ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the sub-deformation gradients
             */

            setPlasticStateVariableDerivatives( false );

        }

        void residual::setdPlasticStateVariablesdT( ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the temperature
             */

            setPlasticStateVariableDerivatives( false );

        }

        void residual::setdPlasticStateVariablesdStateVariables( ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the state variables
             */

            setPlasticStateVariableDerivatives( false );

        }

        void residual::setdPlasticStateVariablesdPreviousCauchyStress( ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the previous Cauchy stress
             */

            setPlasticStateVariableDerivatives( true );

        }

        void residual::setdPlasticStateVariablesdPreviousF( ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the previous deformation gradient
             */

            setPlasticStateVariableDerivatives( true );

        }

        void residual::setdPlasticStateVariablesdPreviousSubFs( ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the previous sub-deformation gradients
             */

            setPlasticStateVariableDerivatives( true );

        }

        void residual::setdPlasticStateVariablesdPreviousT( ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the previous temperature
             */

            setPlasticStateVariableDerivatives( true );

        }

        void residual::setdPlasticStateVariablesdPreviousStateVariables( ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the previous state variables
             */

            setPlasticStateVariableDerivatives( true );

        }

        void residual::setPlasticStateVariableDerivatives( const bool setPreviousDerivatives ){
            /*!
             * Set the plastic state variables
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            const floatVector *stateVariableEvolutionRates;

            const floatVector *previousStateVariableEvolutionRates;

            const floatVector *previousStateVariables;

            const floatVector *dStateVariableEvolutionRatesdCauchyStress;

            const floatVector *dStateVariableEvolutionRatesdF;

            const floatVector *dStateVariableEvolutionRatesdSubFs;

            const floatVector *dStateVariableEvolutionRatesdT;

            const floatVector *dStateVariableEvolutionRatesdStateVariables;

            const floatVector *dPreviousStateVariableEvolutionRatesdPreviousCauchyStress = NULL;

            const floatVector *dPreviousStateVariableEvolutionRatesdPreviousF = NULL;

            const floatVector *dPreviousStateVariableEvolutionRatesdPreviousSubFs = NULL;

            const floatVector *dPreviousStateVariableEvolutionRatesdPreviousT = NULL;

            const floatVector *dPreviousStateVariableEvolutionRatesdPreviousStateVariables = NULL;

            floatVector deltaPlasticStateVariables;

            auto plasticStateVariables = get_setDataStorage_plasticStateVariables( );

            if ( setPreviousDerivatives ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dPreviousStateVariableEvolutionRatesdPreviousCauchyStress = get_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPreviousStateVariableEvolutionRatesdPreviousF = get_dPreviousStateVariableEvolutionRatesdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPreviousStateVariableEvolutionRatesdPreviousSubFs = get_dPreviousStateVariableEvolutionRatesdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPreviousStateVariableEvolutionRatesdPreviousT = get_dPreviousStateVariableEvolutionRatesdPreviousT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPreviousStateVariableEvolutionRatesdPreviousStateVariables = get_dPreviousStateVariableEvolutionRatesdPreviousStateVariables( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( dStateVariableEvolutionRatesdCauchyStress = get_dStateVariableEvolutionRatesdCauchyStress( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( dStateVariableEvolutionRatesdF = get_dStateVariableEvolutionRatesdF( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( dStateVariableEvolutionRatesdSubFs = get_dStateVariableEvolutionRatesdSubFs( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( dStateVariableEvolutionRatesdT = get_dStateVariableEvolutionRatesdT( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( dStateVariableEvolutionRatesdStateVariables = get_dStateVariableEvolutionRatesdStateVariables( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( stateVariableEvolutionRates = get_stateVariableEvolutionRates( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousStateVariableEvolutionRates = get_previousStateVariableEvolutionRates( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousStateVariables = get_previousStateVariables( ) );

            floatVector dXidXidot;

            if ( setPreviousDerivatives ){

                floatVector dXidXidotp;

                TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::midpointEvolutionFlatJ( *hydra->getDeltaTime( ), *previousStateVariables, *previousStateVariableEvolutionRates, *stateVariableEvolutionRates, deltaPlasticStateVariables, *plasticStateVariables.value, dXidXidot, dXidXidotp, ( 1 - *getIntegrationParameter( ) ) ) );

                auto map_dXidXidotp = getDynamicSizeMatrixMap< floatType >( dXidXidot.data( ), num_isvs, num_isvs );

                auto map_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress = getDynamicSizeMatrixMap< floatType >( dPreviousStateVariableEvolutionRatesdPreviousCauchyStress->data( ), num_isvs, sot_dim );

                auto map_dPreviousStateVariableEvolutionRatesdPreviousF = getDynamicSizeMatrixMap< floatType >( dPreviousStateVariableEvolutionRatesdPreviousF->data( ), num_isvs, sot_dim );

                auto map_dPreviousStateVariableEvolutionRatesdPreviousSubFs = getDynamicSizeMatrixMap< floatType >( dPreviousStateVariableEvolutionRatesdPreviousSubFs->data( ), num_isvs, ( num_configs - 1 ) * sot_dim );

                auto map_dPreviousStateVariableEvolutionRatesdPreviousT = getDynamicSizeVectorMap< floatType >( dPreviousStateVariableEvolutionRatesdPreviousT->data( ), num_isvs );

                auto map_dPreviousStateVariableEvolutionRatesdPreviousStateVariables = getDynamicSizeMatrixMap< floatType >( dPreviousStateVariableEvolutionRatesdPreviousStateVariables->data( ), num_isvs, num_isvs );

                auto dPlasticStateVariablesdPreviousCauchyStress = get_setDataStorage_dPlasticStateVariablesdPreviousCauchyStress( );
                auto map_dPlasticStateVariablesdPreviousCauchyStress = dPlasticStateVariablesdPreviousCauchyStress.zeroMap< floatType >( num_isvs, sot_dim );

                auto dPlasticStateVariablesdPreviousF = get_setDataStorage_dPlasticStateVariablesdPreviousF( );
                auto map_dPlasticStateVariablesdPreviousF = dPlasticStateVariablesdPreviousF.zeroMap< floatType >( num_isvs, sot_dim );

                auto dPlasticStateVariablesdPreviousSubFs = get_setDataStorage_dPlasticStateVariablesdPreviousSubFs( );
                auto map_dPlasticStateVariablesdPreviousSubFs = dPlasticStateVariablesdPreviousSubFs.zeroMap< floatType >( num_isvs, ( num_configs - 1 ) * sot_dim );

                auto dPlasticStateVariablesdPreviousT = get_setDataStorage_dPlasticStateVariablesdPreviousT( );
                auto map_dPlasticStateVariablesdPreviousT = dPlasticStateVariablesdPreviousT.zeroMap< floatType >( num_isvs );

                auto dPlasticStateVariablesdPreviousStateVariables = get_setDataStorage_dPlasticStateVariablesdPreviousStateVariables( );
                auto map_dPlasticStateVariablesdPreviousStateVariables = dPlasticStateVariablesdPreviousStateVariables.zeroMap< floatType >( num_isvs, num_isvs );

                map_dPlasticStateVariablesdPreviousCauchyStress   = ( map_dXidXidotp * map_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress ).eval( );

                map_dPlasticStateVariablesdPreviousF              = ( map_dXidXidotp * map_dPreviousStateVariableEvolutionRatesdPreviousF ).eval( );

                map_dPlasticStateVariablesdPreviousSubFs          = ( map_dXidXidotp * map_dPreviousStateVariableEvolutionRatesdPreviousSubFs ).eval( );

                map_dPlasticStateVariablesdPreviousT              = ( map_dXidXidotp * map_dPreviousStateVariableEvolutionRatesdPreviousT ).eval( );

                map_dPlasticStateVariablesdPreviousStateVariables = ( map_dXidXidotp * map_dPreviousStateVariableEvolutionRatesdPreviousStateVariables ).eval( );

                for ( unsigned int i = 0; i < num_isvs; i++ ){ ( *dPlasticStateVariablesdPreviousStateVariables.value )[ num_isvs * i + i ] += 1; }

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::midpointEvolutionFlatJ( *hydra->getDeltaTime( ), *previousStateVariables, *previousStateVariableEvolutionRates, *stateVariableEvolutionRates, deltaPlasticStateVariables, *plasticStateVariables.value, dXidXidot, ( 1 - *getIntegrationParameter( ) ) ) );

            }

            Eigen::Map< Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > map_dXidXidot( dXidXidot.data( ), num_isvs, num_isvs );

            auto map_dStateVariableEvolutionRatesdCauchyStress   = getDynamicSizeMatrixMap< floatType >( dStateVariableEvolutionRatesdCauchyStress->data( ), num_isvs, sot_dim );

            auto map_dStateVariableEvolutionRatesdF              = getDynamicSizeMatrixMap< floatType >( dStateVariableEvolutionRatesdF->data( ), num_isvs, sot_dim );

            auto map_dStateVariableEvolutionRatesdSubFs          = getDynamicSizeMatrixMap< floatType >( dStateVariableEvolutionRatesdSubFs->data( ), num_isvs, ( num_configs - 1 ) * sot_dim );

            auto map_dStateVariableEvolutionRatesdT              = getDynamicSizeVectorMap< floatType >( dStateVariableEvolutionRatesdT->data( ), num_isvs );

            auto map_dStateVariableEvolutionRatesdStateVariables = getDynamicSizeMatrixMap< floatType >( dStateVariableEvolutionRatesdStateVariables->data( ), num_isvs, num_isvs );

            auto dPlasticStateVariablesdCauchyStress = get_setDataStorage_dPlasticStateVariablesdCauchyStress( );
            auto map_dPlasticStateVariablesdCauchyStress = dPlasticStateVariablesdCauchyStress.zeroMap< floatType >( num_isvs, sot_dim );

            auto dPlasticStateVariablesdF = get_setDataStorage_dPlasticStateVariablesdF( );
            auto map_dPlasticStateVariablesdF = dPlasticStateVariablesdF.zeroMap< floatType >( num_isvs, sot_dim );

            auto dPlasticStateVariablesdSubFs = get_setDataStorage_dPlasticStateVariablesdSubFs( );
            auto map_dPlasticStateVariablesdSubFs = dPlasticStateVariablesdSubFs.zeroMap< floatType >( num_isvs, ( num_configs - 1 ) * sot_dim );

            auto dPlasticStateVariablesdT = get_setDataStorage_dPlasticStateVariablesdT( );
            auto map_dPlasticStateVariablesdT = dPlasticStateVariablesdT.zeroMap< floatType >( num_isvs );

            auto dPlasticStateVariablesdStateVariables = get_setDataStorage_dPlasticStateVariablesdStateVariables( );
            auto map_dPlasticStateVariablesdStateVariables = dPlasticStateVariablesdStateVariables.zeroMap< floatType >( num_isvs, num_isvs );

            map_dPlasticStateVariablesdCauchyStress   = ( map_dXidXidot * map_dStateVariableEvolutionRatesdCauchyStress ).eval( );

            map_dPlasticStateVariablesdF              = ( map_dXidXidot * map_dStateVariableEvolutionRatesdF ).eval( );

            map_dPlasticStateVariablesdSubFs          = ( map_dXidXidot * map_dStateVariableEvolutionRatesdSubFs ).eval( );

            map_dPlasticStateVariablesdT              = ( map_dXidXidot * map_dStateVariableEvolutionRatesdT ).eval( );

            map_dPlasticStateVariablesdStateVariables = ( map_dXidXidot * map_dStateVariableEvolutionRatesdStateVariables ).eval( );

        }

        void residual::setStateVariables( ){
            /*!
             * Set the state variables
             */

            setStateVariables( false );

        }

        void residual::setPreviousStateVariables( ){
            /*!
             * Set the state variables
             */

            setStateVariables( true );

        }

        void residual::setStateVariables( const bool isPrevious ){
            /*!
             * Set the state variables
             * 
             * \param isPrevious: Flag for whether the values are at the previous timestep.
             */

            const floatVector *allStateVariables;

            setDataStorageBase< floatVector > stateVariables;

            if ( isPrevious ){

                allStateVariables =  hydra->get_previousNonLinearSolveStateVariables( );

                stateVariables    = get_setDataStorage_previousStateVariables( );

            }
            else{

                allStateVariables =  hydra->get_nonLinearSolveStateVariables( );

                stateVariables    = get_setDataStorage_stateVariables( );

            }

            stateVariables.zero( getStateVariableIndices( )->size( ) );

            for ( auto index = getStateVariableIndices( )->begin( ); index != getStateVariableIndices( )->end( ); index++ ){

                if ( *index >= allStateVariables->size( ) ){

                    
                    std::string message = "The requested state variable is outside of the available range.\n";
                    message            += "  requested index: " + std::to_string( *index ) + "\n";
                    message            += "  total state variable number: " + std::to_string( allStateVariables->size( ) ) + "\n";

                    TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

                }

                ( *stateVariables.value )[ index - getStateVariableIndices( )->begin( ) ] = ( *allStateVariables )[ *index ];

            }

        }

        void residual::setResidual( ){
            /*!
             * Set the value of the residual
             * 
             * Defined as the residual's computed plastic deformation gradient minus the value stored in hydra's configurations
             * and the difference between the computed state variable's and hydra's stored values.
             */

            const unsigned int dim = hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int num_isvs = get_stateVariables( )->size( );

            auto residual = get_setDataStorage_residual( );
            residual.zero( *getNumEquations( ) );

            // Set the residual for the plastic deformation gradient
            for ( unsigned int i = 0; i < sot_dim; i++ ){

                ( *residual.value )[ i ] = ( *get_plasticDeformationGradient( ) )[ i ] - hydra->getConfiguration( *getPlasticConfigurationIndex( ) )[ i ];
    
            }

            // Set the residual for the plastic state variables
            for ( unsigned int i = 0; i < num_isvs; i++ ){

                ( *residual.value )[ i + sot_dim ] = ( *get_plasticStateVariables( ) )[ i ] - ( *get_stateVariables( ) )[ i ];
    
            }

        }

        void residual::setJacobian( ){
            /*!
             * Set the value of the Jacobian
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            const unsigned int num_unknowns = hydra->getNumUnknowns( );

            auto jacobian = get_setDataStorage_jacobian( );
            jacobian.zero( *getNumEquations( ) * num_unknowns );

            // Set the derivatives
            get_dPlasticDeformationGradientdCauchyStress( );

            get_dPlasticStateVariablesdCauchyStress( );

            // Get the Jacobians of the plastic deformation gradient
            for ( unsigned int i = 0; i < sot_dim; i++ ){
                unsigned int row = i;

                // Set the Jacobian with respect to the Cauchy stress
                for ( unsigned int j = 0; j < sot_dim; j++ ){
                    unsigned int col = j;

                    ( *jacobian.value )[ num_unknowns * row + col ] += ( *get_dPlasticDeformationGradientdCauchyStress( ) )[ sot_dim * i + j ];

                }

                // Set the Jacobian with respect to the sub-configurations
                ( *jacobian.value )[ num_unknowns * i + sot_dim + i ] -= 1;
                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){

                    unsigned int col = sot_dim + j;

                    ( *jacobian.value )[ num_unknowns * row + col ] += ( *get_dPlasticDeformationGradientdSubFs( ) )[ ( num_configs - 1 ) * sot_dim * i + j ];

                }

                // Set the Jacobian with respect to the state variables
                for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){
                    unsigned int col = sot_dim + ( num_configs - 1 ) * sot_dim + *ind;

                    ( *jacobian.value )[ num_unknowns * row + col ] += ( *get_dPlasticDeformationGradientdStateVariables( ) )[ num_isvs * i + ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ];

                }

            }

            // Get the Jacobians of the plastic state variables
            for ( unsigned int i = 0; i < num_isvs; i++ ){
                unsigned int row = sot_dim + i;

                // Set the Jacobian with respect to the Cauchy stress
                for ( unsigned int j = 0; j < sot_dim; j++ ){
                    unsigned int col = j;

                    ( *jacobian.value )[ num_unknowns * row + col ] += ( *get_dPlasticStateVariablesdCauchyStress( ) )[ sot_dim * i + j ];

                }

                // Set the Jacobian with respect to the other configurations
                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){
                    unsigned int col = sot_dim + j;

                    ( *jacobian.value )[ num_unknowns * row + col ] += ( *get_dPlasticStateVariablesdSubFs( ) )[ ( num_configs - 1 ) * sot_dim * i + j ];

                }

                // Set the Jacobian with respect to the state variables
                ( *jacobian.value )[ num_unknowns * row + sot_dim + ( num_configs - 1 ) * sot_dim + ( *getStateVariableIndices( ) )[ i ] ] -= 1;
                for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){
                    unsigned int col = sot_dim + ( num_configs - 1 ) * sot_dim + *ind;

                    ( *jacobian.value )[ num_unknowns * row + col ] += ( *get_dPlasticStateVariablesdStateVariables( ) )[ num_isvs * i + ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ];

                }

            }

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature.
             */

            auto dRdT = get_setDataStorage_dRdT( );
            dRdT.zero( *getNumEquations( ) );

            // Set the derivatives
            get_dPlasticDeformationGradientdT( );

            get_dPlasticStateVariablesdT( );

            // Get the Jacobians of the plastic deformation gradient
            for ( unsigned int i = 0; i < get_plasticDeformationGradient( )->size( ); i++ ){

                ( *dRdT.value )[ i ] = ( *get_dPlasticDeformationGradientdT( ) )[ i ];

            }

            // Get the Jacobians of the plastic state variables
            for ( unsigned int i = 0; i < get_plasticStateVariables( )->size( ); i++ ){

                ( *dRdT.value )[ i + get_plasticDeformationGradient( )->size( ) ] += ( *get_dPlasticStateVariablesdT( ) )[ i ];

            }

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient.
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            auto dRdF = get_setDataStorage_dRdF( );
            dRdF.zero( sot_dim * ( *getNumEquations( ) ) );

            // Set the derivatives
            get_dPlasticDeformationGradientdF( );

            get_dPlasticStateVariablesdF( );

            // Get the Jacobians of the plastic deformation gradient
            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int j = 0; j < sot_dim; j++ ){

                    ( *dRdF.value )[ sot_dim * i + j ] = ( *get_dPlasticDeformationGradientdF( ) )[ sot_dim * i + j ];

                }

            }

            // Get the Jacobians of the plastic state variables
            for ( unsigned int i = 0; i < num_isvs; i++ ){

                for ( unsigned int j = 0; j < sot_dim; j++ ){

                    ( *dRdF.value )[ sot_dim * ( i + sot_dim ) + j ] += ( *get_dPlasticStateVariablesdF( ) )[ sot_dim * i + j ];

                }

            }

        }

        const unsigned int* residual::getPlasticConfigurationIndex( ){
            /*!
             * Get plastic configuration index
             */

            return &_plasticConfigurationIndex;

        }

        const std::vector< unsigned int >* residual::getStateVariableIndices( ){
            /*!
             * Get state variable indices
             */

            return &_stateVariableIndices;

        }

        const floatType* residual::getIntegrationParameter( ){
            /*!
             * Get the integration parameter
             */

            return &_integrationParameter;

        }

        void residual::decomposeParameters( const floatVector &parameters ){
            /*!
             * Decompose the incoming parameter vector
             * 
             * \param &parameters: The incoming parameter vector. We assume a
             *     Perzyna viscoplastic driving stress with a Drucker-Prager
             *     yield and flow potential surface. The parameters are
             *     [ n, q0, q1, C1, C2, Tref, Y, A, B, h0, h1 ] where n is the
             *     Perzyna exponent, q0 is the initial drag stress, q1 is the
             *     drag modulus, C1, C2, and Tref are the temperature effect
             *     parameters, Y is the yield stress for the Drucker-Prager equation,
             *     A is the pressure term of the yield equation, B is the flow parameter,
             *     and h0 is the initial hardening modulus, and h1 is the linear hardening
             *     modulus.
             */

            constexpr unsigned int expectedSize = 11;

            TARDIGRADE_ERROR_TOOLS_CHECK( parameters.size( ) == expectedSize, "The parameters vector is not the correct length.\n  parameters: " + std::to_string( parameters.size( ) ) + "\n  required:   " + std::to_string( expectedSize ) + "\n" );

            set_perzynaParameters( { parameters[ 0 ] } );

            set_dragStressParameters( { parameters[ 1 ], parameters[ 2 ] } );

            set_thermalParameters( { parameters[ 3 ], parameters[ 4 ], parameters[ 5 ] } );

            set_yieldParameters( { parameters[ 6 ], parameters[ 7 ] } );

            set_flowParameters( { 0., parameters[ 8 ] } );

            set_hardeningParameters( { parameters[ 9 ], parameters[ 10 ] } );

        }

        void residual::addParameterizationInfo( std::string &parameterization_info ){
            /*!
             * Add parameterization information to the incoming string
             * 
             * \param &parameterization_info: The parameterization info string
             */

            std::stringstream ss;
            ss.precision(9);
            ss << std::scientific;

            ss << "class: tardigradeHydra::perzynaViscoplasticity::residual\n\n";
            ss << "name,                       description,       units, current value\n";
            ss << "   n,      the Perzyna exponential term,        none, " <<    ( *get_perzynaParameters( ) )[ 0 ] << "\n";
            ss << "  q0,           the initial drag stress,      stress, " << ( *get_dragStressParameters( ) )[ 0 ] << "\n";
            ss << "  Ep, the drag stress hardening modulus,      stress, " << ( *get_dragStressParameters( ) )[ 1 ] << "\n";
            ss << "  C1,              the WLF C1 parameter,        none, " <<    ( *get_thermalParameters( ) )[ 0 ] << "\n";
            ss << "  C2,              the WLF C2 parameter, temperature, " <<    ( *get_thermalParameters( ) )[ 1 ] << "\n";
            ss << "Tref,     the WLF reference temperature, temperature, " <<    ( *get_thermalParameters( ) )[ 2 ] << "\n";
            ss << "   Y,              initial yield stress,      stress, " <<      ( *get_yieldParameters( ) )[ 0 ] << "\n";
            ss << "   A,        yield pressure sensitivity,        none, " <<      ( *get_yieldParameters( ) )[ 1 ] << "\n";
            ss << "   B,         flow pressure sensitivity,        none, " <<       ( *get_flowParameters( ) )[ 1 ] << "\n";
            ss << " hi0,        isv initial evolution rate,        none, " <<  ( *get_hardeningParameters( ) )[ 0 ] << "\n";
            ss << " hi1,         isv linear evolution rate,        none, " <<  ( *get_hardeningParameters( ) )[ 1 ] << "\n";

            ss.unsetf(std::ios_base::floatfield);
            parameterization_info.append(ss.str());

        }

    }

}
