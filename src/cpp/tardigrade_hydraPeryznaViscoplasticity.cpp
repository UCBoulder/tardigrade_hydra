/**
  ******************************************************************************
  * \file tardigrade_hydraPeryznaViscoplasticity.h
  ******************************************************************************
  * An implementation of peryznaViscoplasticity using the hydra framework. Used
  * as an example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraPeryznaViscoplasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_stress_tools.h>

namespace tardigradeHydra{

    namespace peryznaViscoplasticity{

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

            floatVector drivingStress;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfiguration = hydra->getPreviousPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getPreviousStress( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfiguration = hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getStress( ) );

            }

            tardigradeConstitutiveTools::pullBackCauchyStress( *cauchyStress, precedingConfiguration, drivingStress );

            if ( isPrevious ){

                set_previousDrivingStress( drivingStress );

            }
            else{

                set_drivingStress( drivingStress );

            }

        }

        void residual::setDrivingStressDerivatives( const bool isPrevious ){
            /*!
             * Set the derivatives of the driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration.
             *
             * \param isPrevious: Flag for whether to compute this in the previous configuration
             */

            const unsigned int dim = *hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector *cauchyStress;

            floatVector precedingConfiguration;

            floatVector drivingStress;

            floatVector precedingConfigurationJacobian;

            const floatVector *dF1dF;

            const floatVector *dF1dSubFs;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfiguration = hydra->getPreviousPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfigurationJacobian = hydra->getPreviousPrecedingConfigurationJacobian( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getPreviousStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dF = hydra->get_previousdF1dF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dSubFs = hydra->get_previousdF1dFn( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfiguration = hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfigurationJacobian = hydra->getPrecedingConfigurationJacobian( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dF = hydra->get_dF1dF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dSubFs = hydra->get_dF1dFn( ) );

            }

            floatMatrix _dDrivingStressdCauchyStress;

            floatMatrix _dDrivingStressdPrecedingF;

            floatVector dDrivingStressdCauchyStress;

            floatVector dDrivingStressdPrecedingF;

            tardigradeConstitutiveTools::pullBackCauchyStress( *cauchyStress, precedingConfiguration, drivingStress,
                                                               _dDrivingStressdCauchyStress, _dDrivingStressdPrecedingF );

            dDrivingStressdCauchyStress = tardigradeVectorTools::appendVectors( _dDrivingStressdCauchyStress );

            dDrivingStressdPrecedingF = tardigradeVectorTools::appendVectors( _dDrivingStressdPrecedingF );

            floatVector dDrivingStressdFn = tardigradeVectorTools::matrixMultiply( dDrivingStressdPrecedingF, precedingConfigurationJacobian, sot_dim, sot_dim, sot_dim, num_configs * sot_dim );

            floatVector dDrivingStressdF( sot_dim * sot_dim, 0 );

            floatVector dDrivingStressdSubFs( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int k = 0; k < sot_dim; k++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        dDrivingStressdF[ sot_dim * i + j ] += dDrivingStressdFn[ num_configs * sot_dim * i + k ] * ( *dF1dF )[ sot_dim * k + j ];

                    }

                }

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){ //TODO: This order messes a bit with the cache. Investigate how to make it faster.

                    dDrivingStressdSubFs[ ( num_configs - 1 ) * sot_dim * i + j ] += dDrivingStressdFn[ num_configs * sot_dim * i + j + sot_dim ];

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        dDrivingStressdSubFs[ ( num_configs - 1 ) * sot_dim * i + j ] += dDrivingStressdFn[ num_configs * sot_dim * i + k ] * ( *dF1dSubFs )[ num_configs * sot_dim * k + j ];

                    }

                }

            }

            if ( isPrevious ){

                set_previousDrivingStress( drivingStress );

                set_dPreviousDrivingStressdPreviousCauchyStress( dDrivingStressdCauchyStress );

                set_dPreviousDrivingStressdPreviousF( dDrivingStressdF );

                set_dPreviousDrivingStressdPreviousSubFs( dDrivingStressdSubFs );

            }
            else{

                set_drivingStress( drivingStress );

                set_dDrivingStressdCauchyStress( dDrivingStressdCauchyStress );

                set_dDrivingStressdF( dDrivingStressdF );

                set_dDrivingStressdSubFs( dDrivingStressdSubFs );

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

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_previousDrivingStress( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_drivingStress( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( flowParameters = get_flowParameters( ) );

            floatVector dgdDrivingStress( drivingStress->size( ), 0 );

            floatVector flowDirection( drivingStress->size( ), 0 );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::druckerPragerSurface( *drivingStress, ( *flowParameters )[ 1 ], ( *flowParameters )[ 0 ], g, dgdDrivingStress, flowDirection ) );

            if ( isPrevious ){

                set_previousFlowDirection( flowDirection );

            }
            else{

                set_flowDirection( flowDirection );

            }

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

            const unsigned int dim = *hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector *drivingStress;

            const floatVector *flowParameters;

            const floatVector *dDrivingStressdCauchyStress;

            const floatVector *dDrivingStressdF;

            const floatVector *dDrivingStressdSubFs;

            floatType g;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = get_dPreviousDrivingStressdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = get_dPreviousDrivingStressdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = get_dPreviousDrivingStressdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_previousDrivingStress( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = get_dDrivingStressdCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = get_dDrivingStressdF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = get_dDrivingStressdSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_drivingStress( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( flowParameters = get_flowParameters( ) );

            floatVector dgdDrivingStress( drivingStress->size( ), 0 );

            floatVector flowDirection( drivingStress->size( ), 0 );

            floatMatrix _dFlowDirectiondDrivingStress;

            floatVector dFlowDirectiondDrivingStress;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::druckerPragerSurface( *drivingStress, ( *flowParameters )[ 1 ], ( *flowParameters )[ 0 ], g, dgdDrivingStress, flowDirection, _dFlowDirectiondDrivingStress ) );

            dFlowDirectiondDrivingStress = tardigradeVectorTools::appendVectors( _dFlowDirectiondDrivingStress );

            floatVector dFlowDirectiondCauchyStress = tardigradeVectorTools::matrixMultiply( dFlowDirectiondDrivingStress, *dDrivingStressdCauchyStress, sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dFlowDirectiondF            = tardigradeVectorTools::matrixMultiply( dFlowDirectiondDrivingStress, *dDrivingStressdF, sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dFlowDirectiondSubFs        = tardigradeVectorTools::matrixMultiply( dFlowDirectiondDrivingStress, *dDrivingStressdSubFs, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            if ( isPrevious ){

                set_previousFlowDirection( flowDirection );

                set_dPreviousFlowDirectiondPreviousCauchyStress( dFlowDirectiondCauchyStress );

                set_dPreviousFlowDirectiondPreviousF( dFlowDirectiondF );

                set_dPreviousFlowDirectiondPreviousSubFs( dFlowDirectiondSubFs );

            }
            else{

                set_flowDirection( flowDirection );

                set_dFlowDirectiondCauchyStress( dFlowDirectiondCauchyStress );

                set_dFlowDirectiondF( dFlowDirectiondF );

                set_dFlowDirectiondSubFs( dFlowDirectiondSubFs );

            }

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

        void residual::setYieldFunction( const bool isPrevious ){
            /*!
             * Set the value of the yield function
             * 
             * \param isPrevious: Flag for whether this is the previous timestep
             */

            const floatVector* drivingStress;

            const floatVector* yieldParameters;

            floatType yieldFunction;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_previousDrivingStress( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_drivingStress( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( yieldParameters = get_yieldParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::druckerPragerSurface( *drivingStress, ( *yieldParameters )[ 1 ], ( *yieldParameters )[ 0 ], yieldFunction ) );

            if ( isPrevious ){

                set_previousYieldFunction( yieldFunction );

            }
            else{

                set_yieldFunction( yieldFunction );

            }

        }

        void residual::setYieldFunctionDerivatives( const bool isPrevious ){
            /*!
             * Set the value of the yield function derivatives
             * 
             * \param isPrevious: Flag for whether this is the previous timestep
             */

            const unsigned int dim = *hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector* drivingStress;

            const floatVector* dDrivingStressdCauchyStress;

            const floatVector* dDrivingStressdF;

            const floatVector* dDrivingStressdSubFs;

            const floatVector* yieldParameters;

            floatType yieldFunction;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = get_dPreviousDrivingStressdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = get_dPreviousDrivingStressdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = get_dPreviousDrivingStressdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_previousDrivingStress( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = get_dDrivingStressdCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = get_dDrivingStressdF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = get_dDrivingStressdSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = get_drivingStress( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( yieldParameters = get_yieldParameters( ) );

            floatVector dYieldFunctiondDrivingStress( drivingStress->size( ), 0 );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::druckerPragerSurface( *drivingStress, ( *yieldParameters )[ 1 ], ( *yieldParameters )[ 0 ], yieldFunction, dYieldFunctiondDrivingStress ) );

            floatVector dYieldFunctiondCauchyStress = tardigradeVectorTools::matrixMultiply( dYieldFunctiondDrivingStress, *dDrivingStressdCauchyStress, 1, sot_dim, sot_dim, sot_dim );
 
            floatVector dYieldFunctiondF = tardigradeVectorTools::matrixMultiply( dYieldFunctiondDrivingStress, *dDrivingStressdF, 1, sot_dim, sot_dim, sot_dim );

            floatVector dYieldFunctiondSubFs = tardigradeVectorTools::matrixMultiply( dYieldFunctiondDrivingStress, *dDrivingStressdSubFs, 1, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            if ( isPrevious ){

                set_previousYieldFunction( yieldFunction );

                set_dPreviousYieldFunctiondPreviousCauchyStress( dYieldFunctiondCauchyStress );

                set_dPreviousYieldFunctiondPreviousF( dYieldFunctiondF );

                set_dPreviousYieldFunctiondPreviousSubFs( dYieldFunctiondSubFs );

            }
            else{

                set_yieldFunction( yieldFunction );

                set_dYieldFunctiondCauchyStress( dYieldFunctiondCauchyStress );

                set_dYieldFunctiondF( dYieldFunctiondF );

                set_dYieldFunctiondSubFs( dYieldFunctiondSubFs );

            }

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

            floatType plasticThermalMultiplier;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( temperature = hydra->getPreviousTemperature( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( temperature = hydra->getTemperature( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( temperatureParameters = get_thermalParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::WLF( *temperature, { ( *temperatureParameters )[ 2 ], ( *temperatureParameters )[ 0 ], ( *temperatureParameters )[ 1 ] }, plasticThermalMultiplier ) ); 

            if ( isPrevious ){

                set_previousPlasticThermalMultiplier( plasticThermalMultiplier );

            }
            else{

                set_plasticThermalMultiplier( plasticThermalMultiplier );

            }

        }

        void residual::setPlasticThermalMultiplierDerivatives( const bool isPrevious ){
            /*!
             * Set the derivatives of the plastic thermal multiplier
             * 
             * \param isPrevious: A flag for if the values are to be computed for the previous (True) or current (False) plastic thermal multiplier
             */

            const floatType *temperature;

            const floatVector *temperatureParameters;

            floatType plasticThermalMultiplier;

            floatType dPlasticThermalMultiplierdT;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( temperature = hydra->getPreviousTemperature( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( temperature = hydra->getTemperature( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( temperatureParameters = get_thermalParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::WLF( *temperature, { ( *temperatureParameters )[ 2 ], ( *temperatureParameters )[ 0 ], ( *temperatureParameters )[ 1 ] }, plasticThermalMultiplier, dPlasticThermalMultiplierdT ) ); 

            if ( isPrevious ){

                set_previousPlasticThermalMultiplier( plasticThermalMultiplier );

                set_dPreviousPlasticThermalMultiplierdPreviousT( dPlasticThermalMultiplierdT );

            }
            else{

                set_plasticThermalMultiplier( plasticThermalMultiplier );

                set_dPlasticThermalMultiplierdT( dPlasticThermalMultiplierdT );

            }

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

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_previousStateVariables( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_stateVariables( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( dragStressParameters = get_dragStressParameters( ) );

            floatType dragStress;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::linearHardening( *stateVariables, floatVector( dragStressParameters->begin( ) + 1, dragStressParameters->end( ) ), ( *dragStressParameters )[ 0 ], dragStress ) );

            if ( isPrevious ){

                set_previousDragStress( dragStress );

            }
            else{

                set_dragStress( dragStress );

            }

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

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_previousStateVariables( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_stateVariables( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( dragStressParameters = get_dragStressParameters( ) );

            floatType dragStress;

            floatVector dDragStressdStateVariables;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::linearHardening( *stateVariables, floatVector( dragStressParameters->begin( ) + 1, dragStressParameters->end( ) ), ( *dragStressParameters )[ 0 ], dragStress, dDragStressdStateVariables ) );

            if ( isPrevious ){

                set_previousDragStress( dragStress );

                set_dPreviousDragStressdPreviousStateVariables( dDragStressdStateVariables );

            }
            else{

                set_dragStress( dragStress );

                set_dDragStressdStateVariables( dDragStressdStateVariables );

            }

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

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_previousStateVariables( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_stateVariables( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( hardeningParameters = get_hardeningParameters( ) );

            floatType hardeningFunction;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::linearHardening( *stateVariables, floatVector( hardeningParameters->begin( ) + 1, hardeningParameters->end( ) ), ( *hardeningParameters )[ 0 ], hardeningFunction ) );

            if ( isPrevious ){

                set_previousHardeningFunction( hardeningFunction );

            }
            else{

                set_hardeningFunction( hardeningFunction );

            }

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

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_previousStateVariables( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = get_stateVariables( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( hardeningParameters = get_hardeningParameters( ) );

            floatType hardeningFunction;

            floatVector dHardeningFunctiondStateVariables;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::linearHardening( *stateVariables, floatVector( hardeningParameters->begin( ) + 1, hardeningParameters->end( ) ), ( *hardeningParameters )[ 0 ], hardeningFunction, dHardeningFunctiondStateVariables ) );

            if ( isPrevious ){

                set_previousHardeningFunction( hardeningFunction );

                set_dPreviousHardeningFunctiondPreviousStateVariables( dHardeningFunctiondStateVariables );

            }
            else{

                set_hardeningFunction( hardeningFunction );

                set_dHardeningFunctiondStateVariables( dHardeningFunctiondStateVariables );

            }

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

            const floatVector *peryznaParameters;

            floatType plasticMultiplier;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( yieldFunction = get_previousYieldFunction( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dragStress = get_previousDragStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticThermalMultiplier = get_previousPlasticThermalMultiplier( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( yieldFunction = get_yieldFunction( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dragStress = get_dragStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticThermalMultiplier = get_plasticThermalMultiplier( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( peryznaParameters = get_peryznaParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::peryznaModel( *yieldFunction, *dragStress, *plasticThermalMultiplier, ( *peryznaParameters )[ 0 ], plasticMultiplier ) );

            if ( isPrevious ){

                set_previousPlasticMultiplier( plasticMultiplier );

            }
            else{

                set_plasticMultiplier( plasticMultiplier );

            }

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

            const floatVector *peryznaParameters;

            floatType plasticMultiplier;

            const floatVector *dYieldFunctiondCauchyStress;

            const floatVector *dYieldFunctiondF;

            const floatVector *dYieldFunctiondSubFs;

            const floatVector *dDragStressdStateVariables;

            const floatType   *dPlasticThermalMultiplierdT;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondCauchyStress = get_dPreviousYieldFunctiondPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondF = get_dPreviousYieldFunctiondPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondSubFs = get_dPreviousYieldFunctiondPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDragStressdStateVariables = get_dPreviousDragStressdPreviousStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticThermalMultiplierdT = get_dPreviousPlasticThermalMultiplierdPreviousT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( yieldFunction = get_previousYieldFunction( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dragStress = get_previousDragStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticThermalMultiplier = get_previousPlasticThermalMultiplier( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondCauchyStress = get_dYieldFunctiondCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondF = get_dYieldFunctiondF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondSubFs = get_dYieldFunctiondSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDragStressdStateVariables = get_dDragStressdStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticThermalMultiplierdT = get_dPlasticThermalMultiplierdT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( yieldFunction = get_yieldFunction( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dragStress = get_dragStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticThermalMultiplier = get_plasticThermalMultiplier( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( peryznaParameters = get_peryznaParameters( ) );

            floatType dPlasticMultiplierdYieldFunction;

            floatType dPlasticMultiplierdDragStress;

            floatType dPlasticMultiplierdPlasticThermalMultiplier;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::peryznaModel( *yieldFunction, *dragStress, *plasticThermalMultiplier, ( *peryznaParameters )[ 0 ], plasticMultiplier, dPlasticMultiplierdYieldFunction, dPlasticMultiplierdDragStress, dPlasticMultiplierdPlasticThermalMultiplier ) );

            floatVector dPlasticMultiplierdCauchyStress = dPlasticMultiplierdYieldFunction * ( *dYieldFunctiondCauchyStress );

            floatVector dPlasticMultiplierdF = dPlasticMultiplierdYieldFunction * ( *dYieldFunctiondF );

            floatVector dPlasticMultiplierdSubFs = dPlasticMultiplierdYieldFunction * ( *dYieldFunctiondSubFs );

            floatType   dPlasticMultiplierdT = dPlasticMultiplierdPlasticThermalMultiplier * ( *dPlasticThermalMultiplierdT );

            floatVector dPlasticMultiplierdStateVariables = dPlasticMultiplierdDragStress * ( *dDragStressdStateVariables );

            if ( isPrevious ){

                set_previousPlasticMultiplier( plasticMultiplier );

                set_dPreviousPlasticMultiplierdPreviousCauchyStress( dPlasticMultiplierdCauchyStress );

                set_dPreviousPlasticMultiplierdPreviousF( dPlasticMultiplierdF );

                set_dPreviousPlasticMultiplierdPreviousSubFs( dPlasticMultiplierdSubFs );

                set_dPreviousPlasticMultiplierdPreviousT( dPlasticMultiplierdT );

                set_dPreviousPlasticMultiplierdPreviousStateVariables( dPlasticMultiplierdStateVariables );

            }
            else{

                set_plasticMultiplier( plasticMultiplier );

                set_dPlasticMultiplierdCauchyStress( dPlasticMultiplierdCauchyStress );

                set_dPlasticMultiplierdF( dPlasticMultiplierdF );

                set_dPlasticMultiplierdSubFs( dPlasticMultiplierdSubFs );

                set_dPlasticMultiplierdT( dPlasticMultiplierdT );

                set_dPlasticMultiplierdStateVariables( dPlasticMultiplierdStateVariables );

            }

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

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = get_previousPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( flowDirection = get_previousFlowDirection( ) );

                set_previousVelocityGradient( ( *plasticMultiplier ) * ( *flowDirection ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = get_plasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( flowDirection = get_flowDirection( ) );

                set_velocityGradient( ( *plasticMultiplier ) * ( *flowDirection ) );

            }

        }

        void residual::setVelocityGradientDerivatives( const bool isPrevious ){
            /*!
             * Set the derivatives of the velocity gradient in the current configuration of the plastic
             * configuration
             * 
             * \param isPrevious: Flag for whether to compute the value at the previous
             *     timestep.
             */

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

                set_previousVelocityGradient( ( *plasticMultiplier ) * ( *flowDirection ) );

                set_dPreviousVelocityGradientdPreviousCauchyStress( tardigradeVectorTools::appendVectors( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdCauchyStress ) ) + ( *plasticMultiplier ) * ( *dFlowDirectiondCauchyStress ) );

                set_dPreviousVelocityGradientdPreviousF( tardigradeVectorTools::appendVectors( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdF ) ) + ( *plasticMultiplier ) * ( *dFlowDirectiondF ) );

                set_dPreviousVelocityGradientdPreviousSubFs( tardigradeVectorTools::appendVectors( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdSubFs ) ) + ( *plasticMultiplier ) * ( *dFlowDirectiondSubFs ) );

                set_dPreviousVelocityGradientdPreviousT( ( *flowDirection ) * ( *dPlasticMultiplierdT ) );

                set_dPreviousVelocityGradientdPreviousStateVariables( tardigradeVectorTools::appendVectors( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdStateVariables ) ) );

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

                set_velocityGradient( ( *plasticMultiplier ) * ( *flowDirection ) );

                set_dVelocityGradientdCauchyStress( tardigradeVectorTools::appendVectors( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdCauchyStress ) ) + ( *plasticMultiplier ) * ( *dFlowDirectiondCauchyStress ) );

                set_dVelocityGradientdF( tardigradeVectorTools::appendVectors( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdF ) ) + ( *plasticMultiplier ) * ( *dFlowDirectiondF ) );

                set_dVelocityGradientdSubFs( tardigradeVectorTools::appendVectors( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdSubFs ) ) + ( *plasticMultiplier ) * ( *dFlowDirectiondSubFs ) );

                set_dVelocityGradientdT( ( *flowDirection ) * ( *dPlasticMultiplierdT ) );

                set_dVelocityGradientdStateVariables( tardigradeVectorTools::appendVectors( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdStateVariables ) ) );

            }

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

            const floatType *hardeningFunction;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = get_previousPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( hardeningFunction = get_previousHardeningFunction( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = get_plasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( hardeningFunction = get_hardeningFunction( ) );

            }

            floatVector stateVariableEvolutionRates = { ( *plasticMultiplier ) * ( *hardeningFunction ) };

            if ( isPrevious ){

                set_previousStateVariableEvolutionRates( stateVariableEvolutionRates );

            }
            else{

                set_stateVariableEvolutionRates( stateVariableEvolutionRates );

            }

        }

        void residual::setStateVariableEvolutionRateDerivatives( const bool isPrevious ){
            /*!
             * Set the value of the derivatives of the state variable evolution rates
             * 
             * \param isPrevious: A flag to indicate if the previous evolution rate
             *     should be computed.
             */

            const floatType *plasticMultiplier;

            const floatType *hardeningFunction;

            const floatVector *dPlasticMultiplierdCauchyStress;

            const floatVector *dPlasticMultiplierdF;

            const floatVector *dPlasticMultiplierdSubFs;

            const floatType   *dPlasticMultiplierdT;

            const floatVector *dPlasticMultiplierdStateVariables;

            const floatVector *dHardeningFunctiondStateVariables;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdCauchyStress = get_dPreviousPlasticMultiplierdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdF = get_dPreviousPlasticMultiplierdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdSubFs = get_dPreviousPlasticMultiplierdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdT = get_dPreviousPlasticMultiplierdPreviousT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdStateVariables = get_dPreviousPlasticMultiplierdPreviousStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dHardeningFunctiondStateVariables = get_dPreviousHardeningFunctiondPreviousStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = get_previousPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( hardeningFunction = get_previousHardeningFunction( ) );

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

            }

            floatVector stateVariableEvolutionRates = { ( *plasticMultiplier ) * ( *hardeningFunction ) };

            floatVector dStateVariableEvolutionRatesdCauchyStress = ( *dPlasticMultiplierdCauchyStress ) * ( *hardeningFunction );

            floatVector dStateVariableEvolutionRatesdF = ( *dPlasticMultiplierdF ) * ( *hardeningFunction );

            floatVector dStateVariableEvolutionRatesdSubFs = ( *dPlasticMultiplierdSubFs ) * ( *hardeningFunction );

            floatVector dStateVariableEvolutionRatesdT = { ( *dPlasticMultiplierdT ) * ( *hardeningFunction ) };

            floatVector dStateVariableEvolutionRatesdStateVariables = ( *dPlasticMultiplierdStateVariables ) * ( *hardeningFunction ) + ( *plasticMultiplier ) * ( *dHardeningFunctiondStateVariables );

            if ( isPrevious ){

                set_previousStateVariableEvolutionRates( stateVariableEvolutionRates );

                set_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( dStateVariableEvolutionRatesdCauchyStress );

                set_dPreviousStateVariableEvolutionRatesdPreviousF( dStateVariableEvolutionRatesdF );

                set_dPreviousStateVariableEvolutionRatesdPreviousSubFs( dStateVariableEvolutionRatesdSubFs );

                set_dPreviousStateVariableEvolutionRatesdPreviousT( dStateVariableEvolutionRatesdT );

                set_dPreviousStateVariableEvolutionRatesdPreviousStateVariables( dStateVariableEvolutionRatesdStateVariables );

            }
            else{

                set_stateVariableEvolutionRates( stateVariableEvolutionRates );

                set_dStateVariableEvolutionRatesdCauchyStress( dStateVariableEvolutionRatesdCauchyStress );

                set_dStateVariableEvolutionRatesdF( dStateVariableEvolutionRatesdF );

                set_dStateVariableEvolutionRatesdSubFs( dStateVariableEvolutionRatesdSubFs );

                set_dStateVariableEvolutionRatesdT( dStateVariableEvolutionRatesdT );

                set_dStateVariableEvolutionRatesdStateVariables( dStateVariableEvolutionRatesdStateVariables );

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

            floatVector plasticDeformationGradient;

            TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradient = get_velocityGradient( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousVelocityGradient = get_previousVelocityGradient( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousPlasticDeformationGradient = hydra->getPreviousConfiguration( *getPlasticConfigurationIndex( ) ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::evolveF( *hydra->getDeltaTime( ), previousPlasticDeformationGradient, *previousVelocityGradient, *velocityGradient, dFp, plasticDeformationGradient, 1 - ( *getIntegrationParameter( ) ), 1 ) );

            set_plasticDeformationGradient( plasticDeformationGradient );

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

            const unsigned int dim = *hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int plastic_config_index = *getPlasticConfigurationIndex( );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            const floatVector *velocityGradient;

            const floatVector *previousVelocityGradient;

            floatVector previousPlasticDeformationGradient;

            floatVector dFp;

            floatVector plasticDeformationGradient;

            floatMatrix _dFdL;

            floatVector dFdL;

            TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradient = get_velocityGradient( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousVelocityGradient = get_previousVelocityGradient( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousPlasticDeformationGradient = hydra->getPreviousConfiguration( *getPlasticConfigurationIndex( ) ) );

            if ( setPreviousDerivatives ){

                floatMatrix _ddFdPreviousF;

                floatMatrix _dFdPreviousF;

                floatMatrix _dFdPreviousL;

                floatVector ddFdPreviousF;

                floatVector dFdPreviousF;

                floatVector dFdPreviousL;

                TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::evolveF( *hydra->getDeltaTime( ), previousPlasticDeformationGradient, *previousVelocityGradient, *velocityGradient, dFp, plasticDeformationGradient, _dFdL, _ddFdPreviousF, _dFdPreviousF, _dFdPreviousL, 1 - ( *getIntegrationParameter( ) ), 1 ) );

                dFdL = tardigradeVectorTools::appendVectors( _dFdL );

                ddFdPreviousF = tardigradeVectorTools::appendVectors( _ddFdPreviousF );

                dFdPreviousF = tardigradeVectorTools::appendVectors( _dFdPreviousF );

                dFdPreviousL = tardigradeVectorTools::appendVectors( _dFdPreviousL );

                set_dPlasticDeformationGradientdPreviousCauchyStress( tardigradeVectorTools::matrixMultiply( dFdPreviousL, *get_dPreviousVelocityGradientdPreviousCauchyStress( ), sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dPlasticDeformationGradientdPreviousF( tardigradeVectorTools::matrixMultiply( dFdPreviousL, *get_dPreviousVelocityGradientdPreviousF( ), sot_dim, sot_dim, sot_dim, sot_dim ) );

                floatVector dPlasticDeformationGradientdPreviousSubFs = tardigradeVectorTools::matrixMultiply( dFdPreviousL, *get_dPreviousVelocityGradientdPreviousSubFs( ), sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

                for ( unsigned int i = 0; i < sot_dim; i++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        dPlasticDeformationGradientdPreviousSubFs[ ( num_configs - 1 ) * sot_dim * i + j + ( plastic_config_index - 1 ) * sot_dim ] += dFdPreviousF[ sot_dim * i + j ];

                    }

                }

                set_dPlasticDeformationGradientdPreviousSubFs( dPlasticDeformationGradientdPreviousSubFs );

                set_dPlasticDeformationGradientdPreviousT( tardigradeVectorTools::matrixMultiply( dFdPreviousL, *get_dPreviousVelocityGradientdPreviousT( ), sot_dim, sot_dim, sot_dim, 1 ) );

                set_dPlasticDeformationGradientdPreviousStateVariables( tardigradeVectorTools::matrixMultiply( dFdPreviousL, *get_dPreviousVelocityGradientdPreviousStateVariables( ), sot_dim, sot_dim, sot_dim, num_isvs ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::evolveF( *hydra->getDeltaTime( ), previousPlasticDeformationGradient, *previousVelocityGradient, *velocityGradient, dFp, plasticDeformationGradient, _dFdL, 1 - ( *getIntegrationParameter( ) ), 1 ) );

                dFdL = tardigradeVectorTools::appendVectors( _dFdL );

            }

            set_plasticDeformationGradient( plasticDeformationGradient );

            set_dPlasticDeformationGradientdCauchyStress( tardigradeVectorTools::matrixMultiply( dFdL, *get_dVelocityGradientdCauchyStress( ), sot_dim, sot_dim, sot_dim, sot_dim ) );

            set_dPlasticDeformationGradientdF( tardigradeVectorTools::matrixMultiply( dFdL, *get_dVelocityGradientdF( ), sot_dim, sot_dim, sot_dim, sot_dim ) );

            set_dPlasticDeformationGradientdSubFs( tardigradeVectorTools::matrixMultiply( dFdL, *get_dVelocityGradientdSubFs( ), sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

            set_dPlasticDeformationGradientdT( tardigradeVectorTools::matrixMultiply( dFdL, *get_dVelocityGradientdT( ), sot_dim, sot_dim, sot_dim, 1 ) );

            set_dPlasticDeformationGradientdStateVariables( tardigradeVectorTools::matrixMultiply( dFdL, *get_dVelocityGradientdStateVariables( ), sot_dim, sot_dim, sot_dim, num_isvs ) );

        }

        void residual::setPlasticStateVariables( ){
            /*!
             * Set the plastic state variables
             */

            const floatVector *stateVariableEvolutionRates;

            const floatVector *previousStateVariableEvolutionRates;

            const floatVector *previousStateVariables;

            floatVector deltaPlasticStateVariables;

            floatVector plasticStateVariables;

            TARDIGRADE_ERROR_TOOLS_CATCH( stateVariableEvolutionRates = get_stateVariableEvolutionRates( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousStateVariableEvolutionRates = get_previousStateVariableEvolutionRates( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousStateVariables = get_previousStateVariables( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolution( *hydra->getDeltaTime( ), *previousStateVariables, *previousStateVariableEvolutionRates, *stateVariableEvolutionRates, deltaPlasticStateVariables, plasticStateVariables, ( 1 - *getIntegrationParameter( ) ) ) );

            set_plasticStateVariables( plasticStateVariables );

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

            const unsigned int dim = *hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

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

            floatVector plasticStateVariables;

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

            floatMatrix _dXidXidot;

            floatVector dXidXidot;

            if ( setPreviousDerivatives ){

                floatMatrix _dXidXidotp;

                floatVector dXidXidotp;

                TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolution( *hydra->getDeltaTime( ), *previousStateVariables, *previousStateVariableEvolutionRates, *stateVariableEvolutionRates, deltaPlasticStateVariables, plasticStateVariables, _dXidXidot, _dXidXidotp, ( 1 - *getIntegrationParameter( ) ) ) );

                dXidXidot = tardigradeVectorTools::appendVectors( _dXidXidot );

                dXidXidotp = tardigradeVectorTools::appendVectors( _dXidXidotp );

                floatVector isv_eye( num_isvs * num_isvs, 0 );
                tardigradeVectorTools::eye( isv_eye );

                set_dPlasticStateVariablesdPreviousCauchyStress( tardigradeVectorTools::matrixMultiply( dXidXidotp, *dPreviousStateVariableEvolutionRatesdPreviousCauchyStress, num_isvs, num_isvs, num_isvs, sot_dim ) );

                set_dPlasticStateVariablesdPreviousF( tardigradeVectorTools::matrixMultiply( dXidXidotp, *dPreviousStateVariableEvolutionRatesdPreviousF, num_isvs, num_isvs, num_isvs, sot_dim ) );

                set_dPlasticStateVariablesdPreviousSubFs( tardigradeVectorTools::matrixMultiply( dXidXidotp, *dPreviousStateVariableEvolutionRatesdPreviousSubFs, num_isvs, num_isvs, num_isvs, ( num_configs - 1 ) * sot_dim ) );

                set_dPlasticStateVariablesdPreviousT( tardigradeVectorTools::matrixMultiply( dXidXidotp, *dPreviousStateVariableEvolutionRatesdPreviousT, num_isvs, num_isvs, num_isvs, 1 ) );

                set_dPlasticStateVariablesdPreviousStateVariables( tardigradeVectorTools::matrixMultiply( dXidXidotp, *dPreviousStateVariableEvolutionRatesdPreviousStateVariables, num_isvs, num_isvs, num_isvs, num_isvs ) + isv_eye );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolution( *hydra->getDeltaTime( ), *previousStateVariables, *previousStateVariableEvolutionRates, *stateVariableEvolutionRates, deltaPlasticStateVariables, plasticStateVariables, _dXidXidot, ( 1 - *getIntegrationParameter( ) ) ) );

                dXidXidot = tardigradeVectorTools::appendVectors( _dXidXidot );

            }

            set_dPlasticStateVariablesdCauchyStress( tardigradeVectorTools::matrixMultiply( dXidXidot, *dStateVariableEvolutionRatesdCauchyStress, num_isvs, num_isvs, num_isvs, sot_dim ) );

            set_dPlasticStateVariablesdF( tardigradeVectorTools::matrixMultiply( dXidXidot, *dStateVariableEvolutionRatesdF, num_isvs, num_isvs, num_isvs, sot_dim ) );

            set_dPlasticStateVariablesdSubFs( tardigradeVectorTools::matrixMultiply( dXidXidot, *dStateVariableEvolutionRatesdSubFs, num_isvs, num_isvs, num_isvs, ( num_configs - 1 ) * sot_dim ) );

            set_dPlasticStateVariablesdT( tardigradeVectorTools::matrixMultiply( dXidXidot, *dStateVariableEvolutionRatesdT, num_isvs, num_isvs, num_isvs, 1 ) );

            set_dPlasticStateVariablesdStateVariables( tardigradeVectorTools::matrixMultiply( dXidXidot, *dStateVariableEvolutionRatesdStateVariables, num_isvs, num_isvs, num_isvs, num_isvs ) );

            set_plasticStateVariables( plasticStateVariables );

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

            if ( isPrevious ){

                allStateVariables =  hydra->get_previousNonLinearSolveStateVariables( );

            }
            else{

                allStateVariables =  hydra->get_nonLinearSolveStateVariables( );

            }

            floatVector stateVariables( getStateVariableIndices( )->size( ), 0 );

            for ( auto index = getStateVariableIndices( )->begin( ); index != getStateVariableIndices( )->end( ); index++ ){

                if ( *index >= allStateVariables->size( ) ){

                    
                    std::string message = "The requested state variable is outside of the available range.\n";
                    message            += "  requested index: " + std::to_string( *index ) + "\n";
                    message            += "  total state variable number: " + std::to_string( allStateVariables->size( ) ) + "\n";

                    TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

                }

                stateVariables[ index - getStateVariableIndices( )->begin( ) ] = ( *allStateVariables )[ *index ];

            }

            if ( isPrevious ){

                set_previousStateVariables( stateVariables );

            }
            else{

                set_stateVariables( stateVariables );

            }

        }

        void residual::setResidual( ){
            /*!
             * Set the value of the residual
             * 
             * Defined as the residual's computed plastic deformation gradient minus the value stored in hydra's configurations
             * and the difference between the computed state variable's and hydra's stored values.
             */

            floatVector residual( *getNumEquations( ), 0 );

            // Set the residual for the plastic deformation gradient
            for ( unsigned int i = 0; i < get_plasticDeformationGradient( )->size( ); i++ ){

                residual[ i ] = ( *get_plasticDeformationGradient( ) )[ i ] - hydra->getConfiguration( *getPlasticConfigurationIndex( ) )[ i ];
    
            }

            // Set the residual for the plastic state variables
            for ( unsigned int i = 0; i < get_stateVariables( )->size( ); i++ ){

                residual[ i + get_plasticDeformationGradient( )->size( ) ] = ( *get_plasticStateVariables( ) )[ i ] - ( *get_stateVariables( ) )[ i ];
    
            }

            setResidual( residual );

        }

        void residual::setJacobian( ){
            /*!
             * Set the value of the Jacobian
             */

            const unsigned int dim = *hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            floatMatrix jacobian( *getNumEquations( ), floatVector( hydra->getUnknownVector( )->size( ), 0 ) );

            // Set the derivatives
            get_dPlasticDeformationGradientdCauchyStress( );

            get_dPlasticStateVariablesdCauchyStress( );

            // Get the Jacobians of the plastic deformation gradient
            for ( unsigned int i = 0; i < sot_dim; i++ ){
                unsigned int row = i;

                // Set the Jacobian with respect to the Cauchy stress
                for ( unsigned int j = 0; j < sot_dim; j++ ){
                    unsigned int col = j;

                    jacobian[ row ][ col ] += ( *get_dPlasticDeformationGradientdCauchyStress( ) )[ sot_dim * i + j ];

                }

                // Set the Jacobian with respect to the sub-configurations
                jacobian[ i ][ sot_dim + i ] -= 1;
                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){

                    unsigned int col = sot_dim + j;

                    jacobian[ row ][ col ] += ( *get_dPlasticDeformationGradientdSubFs( ) )[ ( num_configs - 1 ) * sot_dim * i + j ];

                }

                // Set the Jacobian with respect to the state variables
                for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){
                    unsigned int col = sot_dim + ( num_configs - 1 ) * sot_dim + *ind;

                    jacobian[ row ][ col ] += ( *get_dPlasticDeformationGradientdStateVariables( ) )[ num_isvs * i + ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ];

                }

            }

            // Get the Jacobians of the plastic state variables
            for ( unsigned int i = 0; i < num_isvs; i++ ){
                unsigned int row = sot_dim + i;

                // Set the Jacobian with respect to the Cauchy stress
                for ( unsigned int j = 0; j < sot_dim; j++ ){
                    unsigned int col = j;

                    jacobian[ row ][ col ] += ( *get_dPlasticStateVariablesdCauchyStress( ) )[ sot_dim * i + j ];

                }

                // Set the Jacobian with respect to the other configurations
                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){
                    unsigned int col = sot_dim + j;

                    jacobian[ row ][ col ] += ( *get_dPlasticStateVariablesdSubFs( ) )[ ( num_configs - 1 ) * sot_dim * i + j ];

                }

                // Set the Jacobian with respect to the state variables
                jacobian[ row ][ sot_dim + ( num_configs - 1 ) * sot_dim + ( *getStateVariableIndices( ) )[ i ] ] -= 1;
                for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){
                    unsigned int col = sot_dim + ( num_configs - 1 ) * sot_dim + *ind;

                    jacobian[ row ][ col ] += ( *get_dPlasticStateVariablesdStateVariables( ) )[ num_isvs * i + ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ];

                }

            }

            setJacobian( jacobian );

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature.
             */

            floatVector dRdT( *getNumEquations( ), 0 );

            // Set the derivatives
            get_dPlasticDeformationGradientdT( );

            get_dPlasticStateVariablesdT( );

            // Get the Jacobians of the plastic deformation gradient
            for ( unsigned int i = 0; i < get_plasticDeformationGradient( )->size( ); i++ ){

                dRdT[ i ] = ( *get_dPlasticDeformationGradientdT( ) )[ i ];

            }

            // Get the Jacobians of the plastic state variables
            for ( unsigned int i = 0; i < get_plasticStateVariables( )->size( ); i++ ){

                dRdT[ i + get_plasticDeformationGradient( )->size( ) ] += ( *get_dPlasticStateVariablesdT( ) )[ i ];

            }

            setdRdT( dRdT );

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient.
             */

            const unsigned int dim = *hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int num_isvs = get_plasticStateVariables( )->size( );

            floatMatrix dRdF( *getNumEquations( ), floatVector( sot_dim, 0 ) );

            // Set the derivatives
            get_dPlasticDeformationGradientdF( );

            get_dPlasticStateVariablesdF( );

            // Get the Jacobians of the plastic deformation gradient
            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int j = 0; j < sot_dim; j++ ){

                    dRdF[ i ][ j ] = ( *get_dPlasticDeformationGradientdF( ) )[ sot_dim * i + j ];

                }

            }

            // Get the Jacobians of the plastic state variables
            for ( unsigned int i = 0; i < num_isvs; i++ ){

                for ( unsigned int j = 0; j < sot_dim; j++ ){

                    dRdF[ i + sot_dim ][ j ] += ( *get_dPlasticStateVariablesdF( ) )[ sot_dim * i + j ];

                }

            }

            setdRdF( dRdF );

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
             *     Peryzna viscoplastic driving stress with a Drucker-Prager
             *     yield and flow potential surface. The parameters are
             *     [ n, q0, q1, C1, C2, Tref, Y, A, B, h0, h1 ] where n is the
             *     Peryzna exponent, q0 is the initial drag stress, q1 is the
             *     drag modulus, C1, C2, and Tref are the temperature effect
             *     parameters, Y is the yield stress for the Drucker-Prager equation,
             *     A is the pressure term of the yield equation, B is the flow parameter,
             *     and h0 is the initial hardening modulus, and h1 is the linear hardening
             *     modulus.
             */

            unsigned int expectedSize = 11;

            if ( parameters.size( ) != expectedSize ){

                std::string message = "The parameters vector is not the correct length.\n";
                message            += "  parameters: " + std::to_string( parameters.size( ) ) + "\n";
                message            += "  required:   " + std::to_string( expectedSize ) + "\n";

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

            }

            set_peryznaParameters( { parameters[ 0 ] } );

            set_dragStressParameters( { parameters[ 1 ], parameters[ 2 ] } );

            set_thermalParameters( { parameters[ 3 ], parameters[ 4 ], parameters[ 5 ] } );

            set_yieldParameters( { parameters[ 6 ], parameters[ 7 ] } );

            set_flowParameters( { 0., parameters[ 8 ] } );

            set_hardeningParameters( { parameters[ 9 ], parameters[ 10 ] } );

        }

    }

}
