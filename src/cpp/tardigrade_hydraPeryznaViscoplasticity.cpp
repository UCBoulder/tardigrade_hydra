/**
  ******************************************************************************
  * \file tardigrade-hydraPeryznaViscoplasticity.h
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

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getPreviousCauchyStress( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfiguration = hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getCauchyStress( ) );

            }

            tardigradeConstitutiveTools::pullBackCauchyStress( *cauchyStress, precedingConfiguration, drivingStress );

            if ( isPrevious ){

                setPreviousDrivingStress( drivingStress );

            }
            else{

                setDrivingStress( drivingStress );

            }

        }

        void residual::setDrivingStressDerivatives( const bool isPrevious ){
            /*!
             * Set the derivatives of the driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration.
             *
             * \param isPrevious: Flag for whether to compute this in the previous configuration
             */

            const floatVector *cauchyStress;

            floatVector precedingConfiguration;

            floatVector drivingStress;

            floatMatrix precedingConfigurationGradient;

            const floatMatrix *dF1dF;

            const floatMatrix *dF1dSubFs;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfiguration = hydra->getPreviousPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfigurationGradient = hydra->getPreviousPrecedingConfigurationGradient( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dF = hydra->getPreviousdF1dF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dSubFs = hydra->getPreviousdF1dFn( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfiguration = hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfigurationGradient = hydra->getPrecedingConfigurationGradient( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dF = hydra->getdF1dF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dSubFs = hydra->getdF1dFn( ) );

            }

            floatMatrix dDrivingStressdCauchyStress;

            floatMatrix dDrivingStressdPrecedingF;

            tardigradeConstitutiveTools::pullBackCauchyStress( *cauchyStress, precedingConfiguration, drivingStress,
                                                               dDrivingStressdCauchyStress, dDrivingStressdPrecedingF );

            floatMatrix dDrivingStressdFn = tardigradeVectorTools::dot( dDrivingStressdPrecedingF, precedingConfigurationGradient );

            floatMatrix dDrivingStressdF( drivingStress.size( ), floatVector( precedingConfiguration.size( ), 0 ) );

            floatMatrix dDrivingStressdSubFs( drivingStress.size( ), floatVector( ( *dF1dSubFs )[ 0 ].size( ), 0 ) );

            for ( unsigned int i = 0; i < ( *hydra->getDimension( ) ) * ( *hydra->getDimension( ) ); i++ ){

                for ( unsigned int j = 0; j < ( *hydra->getDimension( ) ) * ( *hydra->getDimension( ) ); j++ ){

                    for ( unsigned int k = 0; k < ( *hydra->getDimension( ) ) * ( *hydra->getDimension( ) ); k++ ){

                        dDrivingStressdF[ i ][ j ] += dDrivingStressdFn[ i ][ k ] * ( *dF1dF )[ k ][ j ];

                    }

                }

                for ( unsigned int j = 0; j < ( *dF1dSubFs )[ 0 ].size( ); j++ ){

                    dDrivingStressdSubFs[ i ][ j ] += dDrivingStressdFn[ i ][ j + precedingConfiguration.size( ) ];

                    for ( unsigned int k = 0; k < ( *hydra->getDimension( ) ) * ( *hydra->getDimension( ) ); k++ ){

                        dDrivingStressdSubFs[ i ][ j ] += dDrivingStressdFn[ i ][ k ] * ( *dF1dSubFs )[ k ][ j ];

                    }

                }

            }

            if ( isPrevious ){

                setPreviousDrivingStress( drivingStress );

                setdPreviousDrivingStressdPreviousCauchyStress( dDrivingStressdCauchyStress );

                setdPreviousDrivingStressdPreviousF( dDrivingStressdF );

                setdPreviousDrivingStressdPreviousSubFs( dDrivingStressdSubFs );

            }
            else{

                setDrivingStress( drivingStress );

                setdDrivingStressdCauchyStress( dDrivingStressdCauchyStress );

                setdDrivingStressdF( dDrivingStressdF );

                setdDrivingStressdSubFs( dDrivingStressdSubFs );

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

        void residual::setDrivingStress( const floatVector &drivingStress ){
            /*!
             * Set the driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration.
             * 
             * \param &drivingStress: The driving stress
             */

            _drivingStress.second = drivingStress;

            _drivingStress.first = true;

            addIterationData( &_drivingStress );

        }

        void residual::setdDrivingStressdCauchyStress( const floatMatrix &dDrivingStressdCauchyStress ){
            /*!
             * Set the derivative of the driving stress i.e. the Cauchy stress pulled back to the current
             * configuration of the plastic configuration with respect to the Cauchy stress
             * 
             * \param &dDrivingStressdCauchyStress: The derivative of the driving stress w.r.t. the Cauchy stress
             */

            _dDrivingStressdCauchyStress.second = dDrivingStressdCauchyStress;

            _dDrivingStressdCauchyStress.first = true;

            addIterationData( &_dDrivingStressdCauchyStress );

        }

        void residual::setdDrivingStressdF( const floatMatrix &dDrivingStressdF ){
            /*!
             * Set the derivative of the driving stress i.e. the Cauchy stress pulled back to the current
             * configuration of the plastic configuration with respect to the deformation gradient
             * 
             * \param &dDrivingStressdF: The derivative of the driving stress w.r.t. the deformation gradient
             */

            _dDrivingStressdF.second = dDrivingStressdF;

            _dDrivingStressdF.first = true;

            addIterationData( &_dDrivingStressdF );

        }

        void residual::setdDrivingStressdSubFs( const floatMatrix &dDrivingStressdSubFs ){
            /*!
             * Set the derivative of the driving stress i.e. the Cauchy stress pulled back to the current
             * configuration of the plastic configuration with respect to the sub-deformation gradients
             * 
             * \param &dDrivingStressdSubFs: The derivative of the driving stress w.r.t. the sub-deformation gradients
             */

            _dDrivingStressdSubFs.second = dDrivingStressdSubFs;

            _dDrivingStressdSubFs.first = true;

            addIterationData( &_dDrivingStressdSubFs );

        }

        void residual::setPreviousDrivingStress( const floatVector &previousDrivingStress ){
            /*!
             * Set the previous driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration.
             * 
             * \param &previousDrivingStress: The previous driving stress
             */

            _previousDrivingStress.second = previousDrivingStress;

            _previousDrivingStress.first = true;

        }

        void residual::setdPreviousDrivingStressdPreviousCauchyStress( const floatMatrix &dPreviousDrivingStressdPreviousCauchyStress ){
            /*!
             * Set the derivative of the previous driving stress i.e. the previous Cauchy stress pulled back to the previous current
             * configuration of the plastic configuration with respect to the previous Cauchy stress
             * 
             * \param &dPreviousDrivingStressdPreviousCauchyStress: The derivative of the previous driving stress w.r.t. the previous Cauchy stress
             */

            _dPreviousDrivingStressdPreviousCauchyStress.second = dPreviousDrivingStressdPreviousCauchyStress;

            _dPreviousDrivingStressdPreviousCauchyStress.first = true;

        }

        void residual::setdPreviousDrivingStressdPreviousF( const floatMatrix &dPreviousDrivingStressdPreviousF ){
            /*!
             * Set the derivative of the previous driving stress i.e. the Cauchy stress pulled back to the current
             * configuration of the plastic configuration with respect to the previous deformation gradient
             * 
             * \param &dPreviousDrivingStressdPreviousF: The derivative of the previous driving stress w.r.t. the previous deformation gradient
             */

            _dPreviousDrivingStressdPreviousF.second = dPreviousDrivingStressdPreviousF;

            _dPreviousDrivingStressdPreviousF.first = true;

        }

        void residual::setdPreviousDrivingStressdPreviousSubFs( const floatMatrix &dPreviousDrivingStressdPreviousSubFs ){
            /*!
             * Set the derivative of the previous driving stress i.e. the Cauchy stress pulled back to the current
             * configuration of the plastic configuration with respect to the previous sub-deformation gradients
             * 
             * \param &dPreviousDrivingStressdPreviousSubFs: The derivative of the previous driving stress w.r.t. the previous sub-deformation gradients
             */

            _dPreviousDrivingStressdPreviousSubFs.second = dPreviousDrivingStressdPreviousSubFs;

            _dPreviousDrivingStressdPreviousSubFs.first = true;

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

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = getPreviousDrivingStress( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = getDrivingStress( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( flowParameters = getFlowParameters( ) );

            floatVector dgdDrivingStress( drivingStress->size( ), 0 );

            floatVector flowDirection( drivingStress->size( ), 0 );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::druckerPragerSurface( *drivingStress, ( *flowParameters )[ 1 ], ( *flowParameters )[ 0 ], g, dgdDrivingStress, flowDirection ) );

            if ( isPrevious ){

                setPreviousFlowDirection( flowDirection );

            }
            else{

                setFlowDirection( flowDirection );

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

            const floatVector *drivingStress;

            const floatVector *flowParameters;

            const floatMatrix *dDrivingStressdCauchyStress;

            const floatMatrix *dDrivingStressdF;

            const floatMatrix *dDrivingStressdSubFs;

            floatType g;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = getdPreviousDrivingStressdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = getdPreviousDrivingStressdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = getdPreviousDrivingStressdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = getPreviousDrivingStress( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = getdDrivingStressdCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = getdDrivingStressdF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = getdDrivingStressdSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = getDrivingStress( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( flowParameters = getFlowParameters( ) );

            floatVector dgdDrivingStress( drivingStress->size( ), 0 );

            floatVector flowDirection( drivingStress->size( ), 0 );

            floatMatrix dFlowDirectiondDrivingStress;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::druckerPragerSurface( *drivingStress, ( *flowParameters )[ 1 ], ( *flowParameters )[ 0 ], g, dgdDrivingStress, flowDirection, dFlowDirectiondDrivingStress ) );

            floatMatrix dFlowDirectiondCauchyStress = tardigradeVectorTools::dot( dFlowDirectiondDrivingStress, *dDrivingStressdCauchyStress );

            floatMatrix dFlowDirectiondF            = tardigradeVectorTools::dot( dFlowDirectiondDrivingStress, *dDrivingStressdF );

            floatMatrix dFlowDirectiondSubFs        = tardigradeVectorTools::dot( dFlowDirectiondDrivingStress, *dDrivingStressdSubFs );

            if ( isPrevious ){

                setPreviousFlowDirection( flowDirection );

                setdPreviousFlowDirectiondPreviousCauchyStress( dFlowDirectiondCauchyStress );

                setdPreviousFlowDirectiondPreviousF( dFlowDirectiondF );

                setdPreviousFlowDirectiondPreviousSubFs( dFlowDirectiondSubFs );

            }
            else{

                setFlowDirection( flowDirection );

                setdFlowDirectiondCauchyStress( dFlowDirectiondCauchyStress );

                setdFlowDirectiondF( dFlowDirectiondF );

                setdFlowDirectiondSubFs( dFlowDirectiondSubFs );

            }

        }

        void residual::setFlowDirection( const floatVector &flowDirection ){
            /*!
             * Set the flow direction in the current configuration of the
             * plastic configuration.
             * 
             * \param &flowDirection: The flow direction in the current configuration
             *     of the plastic configuration.
             */

            _flowDirection.second = flowDirection;

            _flowDirection.first = true;

            addIterationData( &_flowDirection );

        }

        void residual::setdFlowDirectiondCauchyStress( const floatMatrix &dFlowDirectiondCauchyStress ){
            /*!
             * Set the derivative of the flow direction in the current configuration of the
             * plastic configuration w.r.t. the Cauchy stress
             * 
             * \param &dFlowDirectiondCauchyStress: The derivative of the flow direction in the current configuration
             *     of the plastic configuration with respect to the Cauchy stress
             */

            _dFlowDirectiondCauchyStress.second = dFlowDirectiondCauchyStress;

            _dFlowDirectiondCauchyStress.first = true;

            addIterationData( &_dFlowDirectiondCauchyStress );

        }

        void residual::setdFlowDirectiondF( const floatMatrix &dFlowDirectiondF ){
            /*!
             * Set the derivative of the flow direction in the current configuration of the
             * plastic configuration w.r.t. the deformation gradient
             * 
             * \param &dFlowDirectiondF: The derivative of the flow direction in the current configuration
             *     of the plastic configuration with respect to the deformation gradient
             */

            _dFlowDirectiondF.second = dFlowDirectiondF;

            _dFlowDirectiondF.first = true;

            addIterationData( &_dFlowDirectiondF );

        }

        void residual::setdFlowDirectiondSubFs( const floatMatrix &dFlowDirectiondSubFs ){
            /*!
             * Set the derivative of the flow direction in the current configuration of the
             * plastic configuration w.r.t. the sub-deformation gradients
             * 
             * \param &dFlowDirectiondSubFs: The derivative of the flow direction in the current configuration
             *     of the plastic configuration with respect to the sub-deformation gradients
             */

            _dFlowDirectiondSubFs.second = dFlowDirectiondSubFs;

            _dFlowDirectiondSubFs.first = true;

            addIterationData( &_dFlowDirectiondSubFs );

        }

        void residual::setPreviousFlowDirection( const floatVector &previousFlowDirection ){
            /*!
             * Set the previous flow direction in the current configuration of the
             * plastic configuration.
             * 
             * \param &previousFlowDirection: The previous flow direction in the current
             *     configuration of the plastic configuration.
             */

            _previousFlowDirection.second = previousFlowDirection;

            _previousFlowDirection.first = true;

        }

        void residual::setdPreviousFlowDirectiondPreviousCauchyStress( const floatMatrix &dPreviousFlowDirectiondPreviousCauchyStress ){
            /*!
             * Set the derivative of the previous flow direction in the current configuration of the
             * plastic configuration w.r.t. the previous Cauchy stress
             * 
             * \param &dPreviousFlowDirectiondPreviousCauchyStress: The derivative of the previous flow direction in the current configuration
             *     of the plastic configuration with respect to the previous Cauchy stress
             */

            _dPreviousFlowDirectiondPreviousCauchyStress.second = dPreviousFlowDirectiondPreviousCauchyStress;

            _dPreviousFlowDirectiondPreviousCauchyStress.first = true;

        }

        void residual::setdPreviousFlowDirectiondPreviousF( const floatMatrix &dPreviousFlowDirectiondPreviousF ){
            /*!
             * Set the derivative of the previous flow direction in the current configuration of the
             * plastic configuration w.r.t. the previous deformation gradient
             * 
             * \param &dPreviousFlowDirectiondPreviousF: The derivative of the previous flow direction in the current configuration
             *     of the plastic configuration with respect to the previous deformation gradient
             */

            _dPreviousFlowDirectiondPreviousF.second = dPreviousFlowDirectiondPreviousF;

            _dPreviousFlowDirectiondPreviousF.first = true;

        }

        void residual::setdPreviousFlowDirectiondPreviousSubFs( const floatMatrix &dPreviousFlowDirectiondPreviousSubFs ){
            /*!
             * Set the derivative of the previous flow direction in the current configuration of the
             * plastic configuration w.r.t. the previous sub-deformation gradients
             * 
             * \param &dPreviousFlowDirectiondPreviousSubFs: The derivative of the previous flow direction in the current configuration
             *     of the plastic configuration with respect to the previous sub-deformation gradients
             */

            _dPreviousFlowDirectiondPreviousSubFs.second = dPreviousFlowDirectiondPreviousSubFs;

            _dPreviousFlowDirectiondPreviousSubFs.first = true;

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

            setYieldFunctionDerivatives( false );

        }

        void residual::setdYieldFunctiondF( ){
            /*!
             * Set the value of the derivative of the yield function w.r.t. the deformation gradient
             */

            setYieldFunctionDerivatives( false );

        }

        void residual::setdYieldFunctiondSubFs( ){
            /*!
             * Set the value of the derivative of the yield function w.r.t. the sub-deformation gradients
             */

            setYieldFunctionDerivatives( false );

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

            setYieldFunctionDerivatives( true );

        }

        void residual::setdPreviousYieldFunctiondPreviousF( ){
            /*!
             * Set the value of the derivative of the previous yield function w.r.t. the previous deformation gradient
             */

            setYieldFunctionDerivatives( true );

        }

        void residual::setdPreviousYieldFunctiondPreviousSubFs( ){
            /*!
             * Set the value of the derivative of the previous yield function w.r.t. the previous sub-deformation gradients
             */

            setYieldFunctionDerivatives( true );

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

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = getPreviousDrivingStress( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = getDrivingStress( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( yieldParameters = getYieldParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::druckerPragerSurface( *drivingStress, ( *yieldParameters )[ 1 ], ( *yieldParameters )[ 0 ], yieldFunction ) );

            if ( isPrevious ){

                setPreviousYieldFunction( yieldFunction );

            }
            else{

                setYieldFunction( yieldFunction );

            }

        }

        void residual::setYieldFunctionDerivatives( const bool isPrevious ){
            /*!
             * Set the value of the yield function derivatives
             * 
             * \param isPrevious: Flag for whether this is the previous timestep
             */

            const floatVector* drivingStress;

            const floatMatrix* dDrivingStressdCauchyStress;

            const floatMatrix* dDrivingStressdF;

            const floatMatrix* dDrivingStressdSubFs;

            const floatVector* yieldParameters;

            floatType yieldFunction;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = getdPreviousDrivingStressdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = getdPreviousDrivingStressdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = getdPreviousDrivingStressdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = getPreviousDrivingStress( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdCauchyStress = getdDrivingStressdCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdF = getdDrivingStressdF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDrivingStressdSubFs = getdDrivingStressdSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( drivingStress = getDrivingStress( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( yieldParameters = getYieldParameters( ) );

            floatVector dYieldFunctiondDrivingStress( drivingStress->size( ), 0 );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::druckerPragerSurface( *drivingStress, ( *yieldParameters )[ 1 ], ( *yieldParameters )[ 0 ], yieldFunction, dYieldFunctiondDrivingStress ) );

            floatVector dYieldFunctiondCauchyStress = tardigradeVectorTools::Tdot( *dDrivingStressdCauchyStress, dYieldFunctiondDrivingStress );
 
            floatVector dYieldFunctiondF = tardigradeVectorTools::Tdot( *dDrivingStressdF, dYieldFunctiondDrivingStress );

            floatVector dYieldFunctiondSubFs = tardigradeVectorTools::Tdot( *dDrivingStressdSubFs, dYieldFunctiondDrivingStress );

            if ( isPrevious ){

                setPreviousYieldFunction( yieldFunction );

                setdPreviousYieldFunctiondPreviousCauchyStress( dYieldFunctiondCauchyStress );

                setdPreviousYieldFunctiondPreviousF( dYieldFunctiondF );

                setdPreviousYieldFunctiondPreviousSubFs( dYieldFunctiondSubFs );

            }
            else{

                setYieldFunction( yieldFunction );

                setdYieldFunctiondCauchyStress( dYieldFunctiondCauchyStress );

                setdYieldFunctiondF( dYieldFunctiondF );

                setdYieldFunctiondSubFs( dYieldFunctiondSubFs );

            }

        }


        void residual::setYieldFunction( const floatType &yieldFunction ){
            /*!
             * Set the value of the yield function
             */

            _yieldFunction.second = yieldFunction;

            _yieldFunction.first = true;

            addIterationData( &_yieldFunction );

        }

        void residual::setdYieldFunctiondCauchyStress( const floatVector &dYieldFunctiondCauchyStress ){
            /*!
             * Set the value of the derivative of the yield function w.r.t. the Cauchy stress
             */

            _dYieldFunctiondCauchyStress.second = dYieldFunctiondCauchyStress;

            _dYieldFunctiondCauchyStress.first = true;

            addIterationData( &_dYieldFunctiondCauchyStress );

        }

        void residual::setdYieldFunctiondF( const floatVector &dYieldFunctiondF ){
            /*!
             * Set the value of the derivative of the yield function w.r.t. the deformation gradient
             */

            _dYieldFunctiondF.second = dYieldFunctiondF;

            _dYieldFunctiondF.first = true;

            addIterationData( &_dYieldFunctiondF );

        }

        void residual::setdYieldFunctiondSubFs( const floatVector &dYieldFunctiondSubFs ){
            /*!
             * Set the value of the derivative of the yield function w.r.t. the sub-deformation gradients
             */

            _dYieldFunctiondSubFs.second = dYieldFunctiondSubFs;

            _dYieldFunctiondSubFs.first = true;

            addIterationData( &_dYieldFunctiondSubFs );

        }

        void residual::setPreviousYieldFunction( const floatType &previousYieldFunction ){
            /*!
             * Set the value of the previous yield function
             */

            _previousYieldFunction.second = previousYieldFunction;

            _previousYieldFunction.first = true;

        }

        void residual::setdPreviousYieldFunctiondPreviousCauchyStress( const floatVector &dPreviousYieldFunctiondPreviousCauchyStress ){
            /*!
             * Set the value of the derivative of the previous yield function w.r.t. the previous Cauchy stress
             */

            _dPreviousYieldFunctiondPreviousCauchyStress.second = dPreviousYieldFunctiondPreviousCauchyStress;

            _dPreviousYieldFunctiondPreviousCauchyStress.first = true;

        }

        void residual::setdPreviousYieldFunctiondPreviousF( const floatVector &dPreviousYieldFunctiondPreviousF ){
            /*!
             * Set the value of the derivative of the previous yield function w.r.t. the previous deformation gradient
             */

            _dPreviousYieldFunctiondPreviousF.second = dPreviousYieldFunctiondPreviousF;

            _dPreviousYieldFunctiondPreviousF.first = true;

        }

        void residual::setdPreviousYieldFunctiondPreviousSubFs( const floatVector &dPreviousYieldFunctiondPreviousSubFs ){
            /*!
             * Set the value of the derivative of the previous yield function w.r.t. the previous sub-deformation gradients
             */

            _dPreviousYieldFunctiondPreviousSubFs.second = dPreviousYieldFunctiondPreviousSubFs;

            _dPreviousYieldFunctiondPreviousSubFs.first = true;

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

            TARDIGRADE_ERROR_TOOLS_CATCH( temperatureParameters = getThermalParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::WLF( *temperature, { ( *temperatureParameters )[ 2 ], ( *temperatureParameters )[ 0 ], ( *temperatureParameters )[ 1 ] }, plasticThermalMultiplier ) ); 

            if ( isPrevious ){

                setPreviousPlasticThermalMultiplier( plasticThermalMultiplier );

            }
            else{

                setPlasticThermalMultiplier( plasticThermalMultiplier );

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

            TARDIGRADE_ERROR_TOOLS_CATCH( temperatureParameters = getThermalParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::WLF( *temperature, { ( *temperatureParameters )[ 2 ], ( *temperatureParameters )[ 0 ], ( *temperatureParameters )[ 1 ] }, plasticThermalMultiplier, dPlasticThermalMultiplierdT ) ); 

            if ( isPrevious ){

                setPreviousPlasticThermalMultiplier( plasticThermalMultiplier );

                setdPreviousPlasticThermalMultiplierdPreviousT( dPlasticThermalMultiplierdT );

            }
            else{

                setPlasticThermalMultiplier( plasticThermalMultiplier );

                setdPlasticThermalMultiplierdT( dPlasticThermalMultiplierdT );

            }

        }


        void residual::setPlasticThermalMultiplier( const floatType &plasticThermalMultiplier ){
            /*!
             * Set the plastic thermal multiplier
             * 
             * \param &plasticThermalMultiplier: The multiplicative term which gives the evolution equation temperature dependence
             */

            _plasticThermalMultiplier.second = plasticThermalMultiplier;

            _plasticThermalMultiplier.first = true;

            addIterationData( &_plasticThermalMultiplier );

        }

        void residual::setdPlasticThermalMultiplierdT( const floatType &dPlasticThermalMultiplierdT ){
            /*!
             * Set the derivative of the plastic thermal multiplier w.r.t. the temperature
             * 
             * \param &dPlasticThermalMultiplierdT: The derivative of the plastic thermal multiplier w.r.t. temperature
             */

            _dPlasticThermalMultiplierdT.second = dPlasticThermalMultiplierdT;

            _dPlasticThermalMultiplierdT.first = true;

            addIterationData( &_dPlasticThermalMultiplierdT );

        }

        void residual::setPreviousPlasticThermalMultiplier( const floatType &previousPlasticThermalMultiplier ){
            /*!
             * Set the previous plastic thermal multiplier
             * 
             * \param &previousPlasticThermalMultiplier: The previous value of the multiplicative term which gives the evolution equation temperature dependence
             */

            _previousPlasticThermalMultiplier.second = previousPlasticThermalMultiplier;

            _previousPlasticThermalMultiplier.first = true;

        }

        void residual::setdPreviousPlasticThermalMultiplierdPreviousT( const floatType &dPreviousPlasticThermalMultiplierdPreviousT ){
            /*!
             * Set the previous plastic thermal multiplier
             * 
             * \param &dPreviousPlasticThermalMultiplierdPreviousT: The derivative of the previous plastic thermal multiplier w.r.t. the previous temperature
             */

            _dPreviousPlasticThermalMultiplierdPreviousT.second = dPreviousPlasticThermalMultiplierdPreviousT;

            _dPreviousPlasticThermalMultiplierdPreviousT.first = true;

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

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = getPreviousStateVariables( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = getStateVariables( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( dragStressParameters = getDragStressParameters( ) );

            floatType dragStress;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::linearHardening( *stateVariables, floatVector( dragStressParameters->begin( ) + 1, dragStressParameters->end( ) ), ( *dragStressParameters )[ 0 ], dragStress ) );

            if ( isPrevious ){

                setPreviousDragStress( dragStress );

            }
            else{

                setDragStress( dragStress );

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

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = getPreviousStateVariables( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = getStateVariables( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( dragStressParameters = getDragStressParameters( ) );

            floatType dragStress;

            floatVector dDragStressdStateVariables;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::linearHardening( *stateVariables, floatVector( dragStressParameters->begin( ) + 1, dragStressParameters->end( ) ), ( *dragStressParameters )[ 0 ], dragStress, dDragStressdStateVariables ) );

            if ( isPrevious ){

                setPreviousDragStress( dragStress );

                setdPreviousDragStressdPreviousStateVariables( dDragStressdStateVariables );

            }
            else{

                setDragStress( dragStress );

                setdDragStressdStateVariables( dDragStressdStateVariables );

            }

        }

        void residual::setDragStress( const floatType &dragStress ){
            /*!
             * Set the drag stress
             * 
             * \param &dragStress: The value of the drag stress
             */

            _dragStress.second = dragStress;

            _dragStress.first = true;

            addIterationData( &_dragStress );

        }

        void residual::setdDragStressdStateVariables( const floatVector &dDragStressdStateVariables ){
            /*!
             * Set the derivative of the drag stress w.r.t. the state variables
             * 
             * \param &dDragStressdStateVariables: The value of the derivative of the drag stress w.r.t. the state variables
             */

            _dDragStressdStateVariables.second = dDragStressdStateVariables;

            _dDragStressdStateVariables.first = true;

            addIterationData( &_dDragStressdStateVariables );

        }

        void residual::setPreviousDragStress( const floatType &previousDragStress ){
            /*!
             * Set the previous drag stress
             * 
             * \param &previousDragStress: The value of the drag stress
             */

            _previousDragStress.second = previousDragStress;

            _previousDragStress.first = true;

        }

        void residual::setdPreviousDragStressdPreviousStateVariables( const floatVector &dPreviousDragStressdPreviousStateVariables ){
            /*!
             * Set the derivative of the previous drag stress w.r.t. the previous state variables
             * 
             * \param &dPreviousDragStressdPreviousStateVariables: The value of the derivative of the drag stress w.r.t. the previous state variables
             */

            _dPreviousDragStressdPreviousStateVariables.second = dPreviousDragStressdPreviousStateVariables;

            _dPreviousDragStressdPreviousStateVariables.first = true;

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

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = getPreviousStateVariables( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = getStateVariables( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( hardeningParameters = getHardeningParameters( ) );

            floatType hardeningFunction;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::linearHardening( *stateVariables, floatVector( hardeningParameters->begin( ) + 1, hardeningParameters->end( ) ), ( *hardeningParameters )[ 0 ], hardeningFunction ) );

            if ( isPrevious ){

                setPreviousHardeningFunction( hardeningFunction );

            }
            else{

                setHardeningFunction( hardeningFunction );

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

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = getPreviousStateVariables( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( stateVariables = getStateVariables( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( hardeningParameters = getHardeningParameters( ) );

            floatType hardeningFunction;

            floatVector dHardeningFunctiondStateVariables;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::linearHardening( *stateVariables, floatVector( hardeningParameters->begin( ) + 1, hardeningParameters->end( ) ), ( *hardeningParameters )[ 0 ], hardeningFunction, dHardeningFunctiondStateVariables ) );

            if ( isPrevious ){

                setPreviousHardeningFunction( hardeningFunction );

                setdPreviousHardeningFunctiondPreviousStateVariables( dHardeningFunctiondStateVariables );

            }
            else{

                setHardeningFunction( hardeningFunction );

                setdHardeningFunctiondStateVariables( dHardeningFunctiondStateVariables );

            }

        }

        void residual::setHardeningFunction( const floatType &hardeningFunction ){
            /*!
             * Set the hardening function
             * 
             * \param &hardeningFunction: The value of the hardening function
             */

            _hardeningFunction.second = hardeningFunction;

            _hardeningFunction.first = true;

            addIterationData( &_hardeningFunction );

        }

        void residual::setdHardeningFunctiondStateVariables( const floatVector &dHardeningFunctiondStateVariables ){
            /*!
             * Set the derivative of the hardening function w.r.t. the state variables
             * 
             * \param &dHardeningFunctiondStateVariables: The derivative of the hardening function w.r.t. the state variables
             */

            _dHardeningFunctiondStateVariables.second = dHardeningFunctiondStateVariables;

            _dHardeningFunctiondStateVariables.first = true;

            addIterationData( &_dHardeningFunctiondStateVariables );

        }

        void residual::setPreviousHardeningFunction( const floatType &previousHardeningFunction ){
            /*!
             * Set the previous hardening function
             * 
             * \param &previousHardeningFunction: The previous value of the hardening function
             */

            _previousHardeningFunction.second = previousHardeningFunction;

            _previousHardeningFunction.first = true;

        }

        void residual::setdPreviousHardeningFunctiondPreviousStateVariables( const floatVector &dPreviousHardeningFunctiondPreviousStateVariables ){
            /*!
             * Set the previous hardening function
             * 
             * \param &dPreviousHardeningFunctiondPreviousStateVariables: The derivative of the previous value of the hardening function w.r.t. the previous state variables
             */

            _dPreviousHardeningFunctiondPreviousStateVariables.second = dPreviousHardeningFunctiondPreviousStateVariables;

            _dPreviousHardeningFunctiondPreviousStateVariables.first = true;

        }

        void residual::setPlasticMultiplier( ){
            /*!
             * Set the plastic multiplier in the current configuration of the
             * plastic configuration
             */

            setPlasticMultiplier( false );

        }

        void residual::setPreviousPlasticMultiplier( ){
            /*!
             * Set the plastic multiplier in the current configuration of the
             * plastic configuration
             */

            setPlasticMultiplier( true );

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

                TARDIGRADE_ERROR_TOOLS_CATCH( yieldFunction = getPreviousYieldFunction( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dragStress = getPreviousDragStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticThermalMultiplier = getPreviousPlasticThermalMultiplier( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( yieldFunction = getYieldFunction( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dragStress = getDragStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticThermalMultiplier = getPlasticThermalMultiplier( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( peryznaParameters = getPeryznaParameters( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::peryznaModel( *yieldFunction, *dragStress, *plasticThermalMultiplier, ( *peryznaParameters )[ 0 ], plasticMultiplier ) );

            if ( isPrevious ){

                setPreviousPlasticMultiplier( plasticMultiplier );

            }
            else{

                setPlasticMultiplier( plasticMultiplier );

            }

        }

        void residual::setPlasticMultiplier( const floatType &plasticMultiplier ){
            /*!
             * Set the plastic multiplier in the current configuration of the
             * plastic configuration
             * 
             * \param &plasticMultiplier: The plastic multiplier in the current
             *     configuration of the plastic configuration
             */

            _plasticMultiplier.second = plasticMultiplier;

            _plasticMultiplier.first = true;

            addIterationData( &_plasticMultiplier );

        }

        void residual::setPreviousPlasticMultiplier( const floatType &previousPlasticMultiplier ){
            /*!
             * Set the previous plastic multiplier in the current configuration of the
             * plastic configuration
             * 
             * \param &previousPlasticMultiplier: The previous plastic multiplier in
             *     the current configuration of the plastic configuration
             */

            _previousPlasticMultiplier.second = previousPlasticMultiplier;

            _previousPlasticMultiplier.first = true;

        }

        void residual::setVelocityGradient( ){
            /*!
             * Set the velocity gradient in the current configuration of the plastic
             * configuration
             */

            setVelocityGradient( false );

        }

        void residual::setPreviousVelocityGradient( ){
            /*!
             * Set the velocity gradient in the current configuration of the plastic
             * configuration
             */

            setVelocityGradient( true );

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

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = getPreviousPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( flowDirection = getPreviousFlowDirection( ) );

                setPreviousVelocityGradient( ( *plasticMultiplier ) * ( *flowDirection ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = getPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( flowDirection = getFlowDirection( ) );

                setVelocityGradient( ( *plasticMultiplier ) * ( *flowDirection ) );

            }

        }

        void residual::setVelocityGradient( const floatVector &velocityGradient ){
            /*!
             * Set the velocity gradient in the current configuration of the plastic
             * configuration
             * 
             * \param &velocityGradient: The velocity gradient in the current
             *     configuration of the plastic configuration
             */

            _velocityGradient.second = velocityGradient;

            _velocityGradient.first = true;

            addIterationData( &_velocityGradient );

        }

        void residual::setPreviousVelocityGradient( const floatVector &previousVelocityGradient ){
            /*!
             * Set the previous velocity gradient in the current configuration of the plastic
             * configuration
             * 
             * \param &previousVelocityGradient: The velocity gradient in the current
             *     configuration of the plastic configuration
             */

            _previousVelocityGradient.second = previousVelocityGradient;

            _previousVelocityGradient.first = true;

        }

        void residual::setStateVariableEvolutionRates( ){
            /*! 
             * Set the value of the state variable evolution rates
             */

            setStateVariableEvolutionRates( false );

        }

        void residual::setPreviousStateVariableEvolutionRates( ){
            /*! 
             * Set the value of the state variable evolution rates
             */

            setStateVariableEvolutionRates( true );

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

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = getPreviousPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( hardeningFunction = getPreviousHardeningFunction( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = getPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( hardeningFunction = getHardeningFunction( ) );

            }

            floatVector stateVariableEvolutionRates = { ( *plasticMultiplier ) * ( *hardeningFunction ) };

            if ( isPrevious ){

                setPreviousStateVariableEvolutionRates( stateVariableEvolutionRates );

            }
            else{

                setStateVariableEvolutionRates( stateVariableEvolutionRates );

            }

        }

        void residual::setStateVariableEvolutionRates( const floatVector &stateVariableEvolutionRates ){
            /*!
             * Set the state variable evolution rate
             * 
             * \param &stateVariableEvolutionRates: The current state variable evolution rate
             */

            _stateVariableEvolutionRates.second = stateVariableEvolutionRates;

            _stateVariableEvolutionRates.first = true;

            addIterationData( &_stateVariableEvolutionRates );

        }

        void residual::setPreviousStateVariableEvolutionRates( const floatVector &previousStateVariableEvolutionRates ){
            /*!
             * Set the previous state variable evolution rate
             * 
             * \param &previousStateVariableEvolutionRates: The previous state variable evolution rate
             */

            _previousStateVariableEvolutionRates.second = previousStateVariableEvolutionRates;

            _previousStateVariableEvolutionRates.first = true;

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

            TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradient = getVelocityGradient( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousVelocityGradient = getPreviousVelocityGradient( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousPlasticDeformationGradient = hydra->getPreviousConfiguration( *getPlasticConfigurationIndex( ) ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::evolveF( *hydra->getDeltaTime( ), previousPlasticDeformationGradient, *previousVelocityGradient, *velocityGradient, dFp, plasticDeformationGradient, 1 - ( *getIntegrationParameter( ) ), 1 ) );

            setPlasticDeformationGradient( plasticDeformationGradient );

        }

        void residual::setPlasticDeformationGradient( const floatVector &plasticDeformationGradient ){
            /*!
             * Set the plastic deformation gradient
             *
             * \param &plasticDeformationGradient: The plastic deformation gradient
             */

            _plasticDeformationGradient.second = plasticDeformationGradient;

            _plasticDeformationGradient.first = true;

            addIterationData( &_plasticDeformationGradient );

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

            TARDIGRADE_ERROR_TOOLS_CATCH( stateVariableEvolutionRates = getStateVariableEvolutionRates( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousStateVariableEvolutionRates = getPreviousStateVariableEvolutionRates( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousStateVariables = getPreviousStateVariables( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolution( *hydra->getDeltaTime( ), *previousStateVariables, *previousStateVariableEvolutionRates, *stateVariableEvolutionRates, deltaPlasticStateVariables, plasticStateVariables, ( 1 - *getIntegrationParameter( ) ) ) );

            setPlasticStateVariables( plasticStateVariables );

        }

        void residual::setPlasticStateVariables( const floatVector &plasticStateVariables ){
            /*!
             * Set the plastic state variables
             * 
             * \param &plasticStateVariables
             */

            _plasticStateVariables.second = plasticStateVariables;

            _plasticDeformationGradient.first = true;

            addIterationData( &_plasticStateVariables );

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

                allStateVariables =  hydra->getPreviousNonLinearSolveStateVariables( );

            }
            else{

                allStateVariables =  hydra->getNonLinearSolveStateVariables( );

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

                setPreviousStateVariables( stateVariables );

            }
            else{

                setStateVariables( stateVariables );

            }

        }

        void residual::setStateVariables( const floatVector &stateVariables ){
            /*!
             * Set the state variables
             * 
             * \param &stateVariables: The state variables
             */

            _stateVariables.second = stateVariables;

            _stateVariables.first = true;

            addIterationData( &_stateVariables );

        }

        void residual::setPreviousStateVariables( const floatVector &previousStateVariables ){
            /*!
             * Set the previous values of the state variables
             * 
             * \param &previousStateVariables: The previous state variables
             */

            _previousStateVariables.second = previousStateVariables;

            _previousStateVariables.first = true;

        }

        void residual::setPeryznaParameters( const floatVector &peryznaParameters ){
            /*!
             * Set the Peryzna parameters
             * 
             * \param &peryznaParameters: The Peryzna parameters
             */

            _peryznaParameters.second = peryznaParameters;

            _peryznaParameters.first = true;

        }

        void residual::setDragStressParameters( const floatVector &dragStressParameters ){
            /*!
             * Set the drag stress parameters
             * 
             * \param &dragStressParameters: The drag stress parameters
             */

            _dragStressParameters.second = dragStressParameters;

            _dragStressParameters.first = true;

        }

        void residual::setThermalParameters( const floatVector &thermalParameters ){
            /*!
             * Set the thermal parameters
             * 
             * \param &thermalParameters: The thermal parameters
             */

            _thermalParameters.second = thermalParameters;

            _thermalParameters.first = true;

        }

        void residual::setYieldParameters( const floatVector &yieldParameters ){
            /*!
             * Set the yield parameters
             * 
             * \param &yieldParameters: The yield parameters
             */

            _yieldParameters.second = yieldParameters;

            _yieldParameters.first = true;

        }

        void residual::setFlowParameters( const floatVector &flowParameters ){
            /*!
             * Set the flow parameters
             * 
             * \param &flowParameters: The flow parameters
             */

            _flowParameters.second = flowParameters;

            _flowParameters.first = true;

        }

        void residual::setHardeningParameters( const floatVector &hardeningParameters ){
            /*!
             * Set the hardening parameters
             * 
             * \param &hardeningParameters: The hardening parameters
             */

            _hardeningParameters.second = hardeningParameters;

            _hardeningParameters.first = true;

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
            for ( unsigned int i = 0; i < getPlasticDeformationGradient( )->size( ); i++ ){

                residual[ i ] = ( *getPlasticDeformationGradient( ) )[ i ] - hydra->getConfiguration( *getPlasticConfigurationIndex( ) )[ i ];
    
            }

            // Set the residual for the plastic state variables
            for ( unsigned int i = 0; i < getStateVariables( )->size( ); i++ ){

                residual[ i + getPlasticDeformationGradient( )->size( ) ] = ( *getPlasticStateVariables( ) )[ i ] - ( *getStateVariables( ) )[ i ];
    
            }

            setResidual( residual );

        }

        void residual::setJacobian( ){
            /*!
             * Set the value of the Jacobian
             */

            throw "not implemented";

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature.
             */

            throw "not implemented";

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient.
             */

            throw "not implemented";

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

        const floatVector* residual::getDrivingStress( ){
            /*!
             * Get the driving stress
             */

            if ( !_drivingStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setDrivingStress( ) );

            }

            return &_drivingStress.second;

        }

        const floatMatrix* residual::getdDrivingStressdCauchyStress( ){
            /*!
             * Get the derivative of the driving stress w.r.t. the Cauchy stress
             */

            if ( !_dDrivingStressdCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdDrivingStressdCauchyStress( ) );

            }

            return &_dDrivingStressdCauchyStress.second;

        }

        const floatMatrix* residual::getdDrivingStressdF( ){
            /*!
             * Get the derivative of the driving stress w.r.t. the deformation gradient
             */

            if ( !_dDrivingStressdF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdDrivingStressdF( ) );

            }

            return &_dDrivingStressdF.second;

        }

        const floatMatrix* residual::getdDrivingStressdSubFs( ){
            /*!
             * Get the derivative of the driving stress w.r.t. the sub-deformation gradients
             */

            if ( !_dDrivingStressdSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdDrivingStressdSubFs( ) );

            }

            return &_dDrivingStressdSubFs.second;

        }

        const floatVector* residual::getFlowDirection( ){
            /*!
             * Get the flow direction
             */

            if ( !_flowDirection.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setFlowDirection( ) );

            }

            return &_flowDirection.second;

        }

        const floatMatrix* residual::getdFlowDirectiondCauchyStress( ){
            /*!
             * Get the derivative of the flow direction w.r.t. the Cauchy stress
             */

            if ( !_dFlowDirectiondCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdFlowDirectiondCauchyStress( ) );

            }

            return &_dFlowDirectiondCauchyStress.second;

        }

        const floatMatrix* residual::getdFlowDirectiondF( ){
            /*!
             * Get the derivative of the flow direction w.r.t. the deformation gradient
             */

            if ( !_dFlowDirectiondF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdFlowDirectiondF( ) );

            }

            return &_dFlowDirectiondF.second;

        }

        const floatMatrix* residual::getdFlowDirectiondSubFs( ){
            /*!
             * Get the derivative of the flow direction w.r.t. the sub-deformation gradients
             */

            if ( !_dFlowDirectiondSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdFlowDirectiondSubFs( ) );

            }

            return &_dFlowDirectiondSubFs.second;

        }

        const floatType* residual::getYieldFunction( ){
            /*!
             * Get the value of the yield function
             */

            if ( !_yieldFunction.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setYieldFunction( ) );

            }

            return &_yieldFunction.second;

        }

        const floatVector* residual::getdYieldFunctiondCauchyStress( ){
            /*!
             * Get the value of the derivative of the yield function w.r.t. the Cauchy stress
             */

            if ( !_dYieldFunctiondCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdYieldFunctiondCauchyStress( ) );

            }

            return &_dYieldFunctiondCauchyStress.second;

        }

        const floatVector* residual::getdYieldFunctiondF( ){
            /*!
             * Get the value of the derivative of the yield function w.r.t. the deformation gradient
             */

            if ( !_dYieldFunctiondF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdYieldFunctiondF( ) );

            }

            return &_dYieldFunctiondF.second;

        }

        const floatVector* residual::getdYieldFunctiondSubFs( ){
            /*!
             * Get the value of the derivative of the yield function w.r.t. the sub-deformation gradients
             */

            if ( !_dYieldFunctiondSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdYieldFunctiondSubFs( ) );

            }

            return &_dYieldFunctiondSubFs.second;

        }

        const floatVector* residual::getdPreviousYieldFunctiondPreviousCauchyStress( ){
            /*!
             * Get the value of the derivative of the previous yield function w.r.t. the previous Cauchy stress
             */

            if ( !_dPreviousYieldFunctiondPreviousCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousYieldFunctiondPreviousCauchyStress( ) );

            }

            return &_dPreviousYieldFunctiondPreviousCauchyStress.second;

        }

        const floatVector* residual::getdPreviousYieldFunctiondPreviousF( ){
            /*!
             * Get the value of the derivative of the previous yield function w.r.t. the previous deformation gradient
             */

            if ( !_dPreviousYieldFunctiondPreviousF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousYieldFunctiondPreviousF( ) );

            }

            return &_dPreviousYieldFunctiondPreviousF.second;

        }

        const floatVector* residual::getdPreviousYieldFunctiondPreviousSubFs( ){
            /*!
             * Get the value of the derivative of the previous yield function w.r.t. the previous sub-deformation gradients
             */

            if ( !_dPreviousYieldFunctiondPreviousSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousYieldFunctiondPreviousSubFs( ) );

            }

            return &_dPreviousYieldFunctiondPreviousSubFs.second;

        }

        const floatType* residual::getPlasticThermalMultiplier( ){
            /*!
             * Get the value of the plastic thermal multiplier
             */

            if ( !_plasticThermalMultiplier.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPlasticThermalMultiplier( ) );

            }

            return &_plasticThermalMultiplier.second;

        }

        const floatType* residual::getdPlasticThermalMultiplierdT( ){
            /*!
             * Get the value of the derivative of the plastic thermal multiplier w.r.t. the temperature
             */

            if ( !_dPlasticThermalMultiplierdT.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticThermalMultiplierdT( ) );

            }

            return &_dPlasticThermalMultiplierdT.second;

        }

        const floatType* residual::getDragStress( ){
            /*!
             * Get the drag stress
             */

            if ( !_dragStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setDragStress( ) );

            }

            return &_dragStress.second;

        }

        const floatVector* residual::getdDragStressdStateVariables( ){
            /*!
             * Get the derivative drag stress w.r.t. the state variables
             */

            if ( !_dDragStressdStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdDragStressdStateVariables( ) );

            }

            return &_dDragStressdStateVariables.second;

        }

        const floatType* residual::getHardeningFunction( ){
            /*!
             * Get the value of the hardening function
             */

            if ( !_hardeningFunction.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setHardeningFunction( ) );

            }

            return &_hardeningFunction.second;

        }

        const floatVector* residual::getdHardeningFunctiondStateVariables( ){
            /*!
             * Get the value of the derivative of the hardening function w.r.t. the state variables
             */

            if ( !_dHardeningFunctiondStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdHardeningFunctiondStateVariables( ) );

            }

            return &_dHardeningFunctiondStateVariables.second;

        }

        const floatType* residual::getPlasticMultiplier( ){
            /*!
             * Get the plastic multiplier
             */

            if ( !_plasticMultiplier.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPlasticMultiplier( ) );

            }

            return &_plasticMultiplier.second;

        }

        const floatVector* residual::getVelocityGradient( ){
            /*!
             * Get the velocity gradient
             */

            if ( !_velocityGradient.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setVelocityGradient( ) );

            }

            return &_velocityGradient.second;

        }

        const floatVector* residual::getStateVariableEvolutionRates( ){
            /*!
             * Get the state variable evolution rate
             */

            if ( !_stateVariableEvolutionRates.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setStateVariableEvolutionRates( ) );

            }

            return &_stateVariableEvolutionRates.second;

        }

        const floatVector* residual::getPlasticDeformationGradient( ){
            /*!
             * Get the plastic deformation gradient
             */

            if ( !_plasticDeformationGradient.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPlasticDeformationGradient( ) );

            }

            return &_plasticDeformationGradient.second;

        }

        const floatVector* residual::getPlasticStateVariables( ){
            /*!
             * Get the plastic state variables
             */

            if ( !_plasticStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPlasticStateVariables( ) );

            }

            return &_plasticStateVariables.second;

        }

        const floatVector* residual::getStateVariables( ){
            /*!
             * Get the state variables
             */

            if ( !_stateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setStateVariables( ) );

            }

            return &_stateVariables.second;

        }

        const floatVector* residual::getPreviousDrivingStress( ){
            /*!
             * Get the previous driving stress
             */

            if ( !_previousDrivingStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousDrivingStress( ) );

            }

            return &_previousDrivingStress.second;

        }

        const floatMatrix* residual::getdPreviousDrivingStressdPreviousCauchyStress( ){
            /*!
             * Get the derivative of the previous driving stress with respect to the previous Cauchy stress
             */

            if ( !_dPreviousDrivingStressdPreviousCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousDrivingStressdPreviousCauchyStress( ) );

            }

            return &_dPreviousDrivingStressdPreviousCauchyStress.second;

        }

        const floatMatrix* residual::getdPreviousDrivingStressdPreviousF( ){
            /*!
             * Get the derivative of the previous driving stress with respect to the previous deformation gradient
             */

            if ( !_dPreviousDrivingStressdPreviousF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousDrivingStressdPreviousF( ) );

            }

            return &_dPreviousDrivingStressdPreviousF.second;

        }

        const floatMatrix* residual::getdPreviousDrivingStressdPreviousSubFs( ){
            /*!
             * Get the derivative of the previous driving stress with respect to the previous sub-deformation gradients
             */

            if ( !_dPreviousDrivingStressdPreviousSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousDrivingStressdPreviousSubFs( ) );

            }

            return &_dPreviousDrivingStressdPreviousSubFs.second;

        }

        const floatVector* residual::getPreviousFlowDirection( ){
            /*!
             * Get the previous flow direction
             */

            if ( !_previousFlowDirection.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousFlowDirection( ) );

            }

            return &_previousFlowDirection.second;

        }

        const floatMatrix* residual::getdPreviousFlowDirectiondPreviousCauchyStress( ){
            /*!
             * Get the derivative of the previous flow direction w.r.t. the Cauchy stress
             */

            if ( !_dPreviousFlowDirectiondPreviousCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousFlowDirectiondPreviousCauchyStress( ) );

            }

            return &_dPreviousFlowDirectiondPreviousCauchyStress.second;

        }

        const floatMatrix* residual::getdPreviousFlowDirectiondPreviousF( ){
            /*!
             * Get the derivative of the previous flow direction w.r.t. the deformation gradient
             */

            if ( !_dPreviousFlowDirectiondPreviousF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousFlowDirectiondPreviousF( ) );

            }

            return &_dPreviousFlowDirectiondPreviousF.second;

        }

        const floatMatrix* residual::getdPreviousFlowDirectiondPreviousSubFs( ){
            /*!
             * Get the derivative of the previous flow direction w.r.t. the sub-deformation gradients
             */

            if ( !_dPreviousFlowDirectiondPreviousSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousFlowDirectiondPreviousSubFs( ) );

            }

            return &_dPreviousFlowDirectiondPreviousSubFs.second;

        }

        const floatType* residual::getPreviousYieldFunction( ){
            /*!
             * Get the previous value of the yield function
             */

            if ( !_previousYieldFunction.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousYieldFunction( ) );

            }

            return &_previousYieldFunction.second;

        }

        const floatType* residual::getPreviousPlasticThermalMultiplier( ){
            /*!
             * Get the previous value of the plastic thermal multiplier
             */

            if ( !_previousPlasticThermalMultiplier.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousPlasticThermalMultiplier( ) );

            }

            return &_previousPlasticThermalMultiplier.second;

        }

        const floatType* residual::getdPreviousPlasticThermalMultiplierdPreviousT( ){
            /*!
             * Get the previous value of the derivative of the plastic thermal multiplier w.r.t. the previous temperature
             */

            if ( !_dPreviousPlasticThermalMultiplierdPreviousT.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousPlasticThermalMultiplierdPreviousT( ) );

            }

            return &_dPreviousPlasticThermalMultiplierdPreviousT.second;

        }

        const floatType* residual::getPreviousDragStress( ){
            /*!
             * Get the previous value of the drag stress
             */

            if ( !_previousDragStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousDragStress( ) );

            }

            return &_previousDragStress.second;

        }

        const floatVector* residual::getdPreviousDragStressdPreviousStateVariables( ){
            /*!
             * Get the derivative of the previous value of the drag stress w.r.t. the previous state variables
             */

            if ( !_dPreviousDragStressdPreviousStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousDragStressdPreviousStateVariables( ) );

            }

            return &_dPreviousDragStressdPreviousStateVariables.second;

        }

        const floatType* residual::getPreviousHardeningFunction( ){
            /*!
             * Get the previous value of the hardening function
             */

            if ( !_previousHardeningFunction.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousHardeningFunction( ) );

            }

            return &_previousHardeningFunction.second;

        }

        const floatVector* residual::getdPreviousHardeningFunctiondPreviousStateVariables( ){
            /*!
             * Get the derivative of the previous value of the hardening function w.r.t. the previous state variables
             */

            if ( !_dPreviousHardeningFunctiondPreviousStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousHardeningFunctiondPreviousStateVariables( ) );

            }

            return &_dPreviousHardeningFunctiondPreviousStateVariables.second;

        }

        const floatType* residual::getPreviousPlasticMultiplier( ){
            /*!
             * Get the previous plastic multiplier
             */

            if ( !_previousPlasticMultiplier.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousPlasticMultiplier( ) );

            }

            return &_previousPlasticMultiplier.second;

        }

        const floatVector* residual::getPreviousVelocityGradient( ){
            /*!
             * Get the previous velocity gradient
             */

            if ( !_previousVelocityGradient.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousVelocityGradient( ) );

            }

            return &_previousVelocityGradient.second;

        }

        const floatVector* residual::getPreviousStateVariableEvolutionRates( ){
            /*!
             * Get the previous state variable evolution rate
             */

            if ( !_previousStateVariableEvolutionRates.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousStateVariableEvolutionRates( ) );

            }

            return &_previousStateVariableEvolutionRates.second;

        }

        const floatVector* residual::getPreviousStateVariables( ){
            /*!
             * Get the previous state variables
             */

            if ( !_previousStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousStateVariables( ) );

            }

            return &_previousStateVariables.second;

        }

        const floatVector* residual::getPeryznaParameters( ){
            /*!
             * Get the Peryzna parameters
             */

            if ( !_peryznaParameters.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "Peryzna parameters not defined but required" ) );

            }

            return &_peryznaParameters.second;

        }

        const floatVector* residual::getDragStressParameters( ){
            /*!
             * Get the drag stress parameters
             */

            if ( !_dragStressParameters.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "Drag stress parameters not defined but required" ) );

            }

            return &_dragStressParameters.second;

        }

        const floatVector* residual::getThermalParameters( ){
            /*!
             * Get the thermal parameters
             */

            if ( !_thermalParameters.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "Thermal parameters not defined but required" ) );

            }

            return &_thermalParameters.second;

        }

        const floatVector* residual::getYieldParameters( ){
            /*!
             * Get the yield parameters
             */

            if ( !_yieldParameters.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "Yield parameters not defined but required" ) );

            }

            return &_yieldParameters.second;

        }

        const floatVector* residual::getFlowParameters( ){
            /*!
             * Get the flow parameters
             */

            if ( !_flowParameters.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "Flow parameters not defined but required" ) );

            }

            return &_flowParameters.second;

        }

        const floatVector* residual::getHardeningParameters( ){
            /*!
             * Get the hardening parameters
             */

            if ( !_hardeningParameters.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "Hardening parameters not defined but required" ) );

            }

            return &_hardeningParameters.second;

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

            setPeryznaParameters( { parameters[ 0 ] } );

            setDragStressParameters( { parameters[ 1 ], parameters[ 2 ] } );

            setThermalParameters( { parameters[ 3 ], parameters[ 4 ], parameters[ 5 ] } );

            setYieldParameters( { parameters[ 6 ], parameters[ 7 ] } );

            setFlowParameters( { 0., parameters[ 8 ] } );

            setHardeningParameters( { parameters[ 9 ], parameters[ 10 ] } );

        }

    }

}
