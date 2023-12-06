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

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getPreviousStress( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfiguration = hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getStress( ) );

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

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getPreviousStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dF = hydra->getPreviousdF1dF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dF1dSubFs = hydra->getPreviousdF1dFn( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfiguration = hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfigurationGradient = hydra->getPrecedingConfigurationGradient( *getPlasticConfigurationIndex( ) ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = hydra->getStress( ) );

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

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondCauchyStress = getdPreviousYieldFunctiondPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondF = getdPreviousYieldFunctiondPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondSubFs = getdPreviousYieldFunctiondPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDragStressdStateVariables = getdPreviousDragStressdPreviousStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticThermalMultiplierdT = getdPreviousPlasticThermalMultiplierdPreviousT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( yieldFunction = getPreviousYieldFunction( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dragStress = getPreviousDragStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticThermalMultiplier = getPreviousPlasticThermalMultiplier( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondCauchyStress = getdYieldFunctiondCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondF = getdYieldFunctiondF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dYieldFunctiondSubFs = getdYieldFunctiondSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dDragStressdStateVariables = getdDragStressdStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticThermalMultiplierdT = getdPlasticThermalMultiplierdT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( yieldFunction = getYieldFunction( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dragStress = getDragStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticThermalMultiplier = getPlasticThermalMultiplier( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( peryznaParameters = getPeryznaParameters( ) );

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

                setPreviousPlasticMultiplier( plasticMultiplier );

                setdPreviousPlasticMultiplierdPreviousCauchyStress( dPlasticMultiplierdCauchyStress );

                setdPreviousPlasticMultiplierdPreviousF( dPlasticMultiplierdF );

                setdPreviousPlasticMultiplierdPreviousSubFs( dPlasticMultiplierdSubFs );

                setdPreviousPlasticMultiplierdPreviousT( dPlasticMultiplierdT );

                setdPreviousPlasticMultiplierdPreviousStateVariables( dPlasticMultiplierdStateVariables );

            }
            else{

                setPlasticMultiplier( plasticMultiplier );

                setdPlasticMultiplierdCauchyStress( dPlasticMultiplierdCauchyStress );

                setdPlasticMultiplierdF( dPlasticMultiplierdF );

                setdPlasticMultiplierdSubFs( dPlasticMultiplierdSubFs );

                setdPlasticMultiplierdT( dPlasticMultiplierdT );

                setdPlasticMultiplierdStateVariables( dPlasticMultiplierdStateVariables );

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

        void residual::setdPlasticMultiplierdCauchyStress( const floatVector &dPlasticMultiplierdCauchyStress ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration with respect to the Cauchy stress
             * 
             * \param &dPlasticMultiplierdCauchyStress: The derivative of the plastic multiplier in the current
             *     configuration of the plastic configuration w.r.t. the Cauchy stress
             */

            _dPlasticMultiplierdCauchyStress.second = dPlasticMultiplierdCauchyStress;

            _dPlasticMultiplierdCauchyStress.first = true;

            addIterationData( &_dPlasticMultiplierdCauchyStress );

        }

        void residual::setdPlasticMultiplierdF( const floatVector &dPlasticMultiplierdF ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration with respect to the deformation gradient
             * 
             * \param &dPlasticMultiplierdF: The derivative of the plastic multiplier in the current
             *     configuration of the plastic configuration w.r.t. the deformation gradient
             */

            _dPlasticMultiplierdF.second = dPlasticMultiplierdF;

            _dPlasticMultiplierdF.first = true;

            addIterationData( &_dPlasticMultiplierdF );

        }

        void residual::setdPlasticMultiplierdSubFs( const floatVector &dPlasticMultiplierdSubFs ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration with respect to the sub-deformation gradients
             * 
             * \param &dPlasticMultiplierdSubFs: The derivative of the plastic multiplier in the current
             *     configuration of the plastic configuration w.r.t. the sub-deformation gradient
             */

            _dPlasticMultiplierdSubFs.second = dPlasticMultiplierdSubFs;

            _dPlasticMultiplierdSubFs.first = true;

            addIterationData( &_dPlasticMultiplierdSubFs );

        }

        void residual::setdPlasticMultiplierdT( const floatType &dPlasticMultiplierdT ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration with respect to the temperature
             * 
             * \param &dPlasticMultiplierdT: The derivative of the plastic multiplier in the current
             *     configuration of the plastic configuration w.r.t. the temperature
             */

            _dPlasticMultiplierdT.second = dPlasticMultiplierdT;

            _dPlasticMultiplierdT.first = true;

            addIterationData( &_dPlasticMultiplierdT );

        }

        void residual::setdPlasticMultiplierdStateVariables( const floatVector &dPlasticMultiplierdStateVariables ){
            /*!
             * Set the derivative of the plastic multiplier in the current configuration of the
             * plastic configuration with respect to the state variables
             * 
             * \param &dPlasticMultiplierdStateVariables: The derivative of the plastic multiplier in the current
             *     configuration of the plastic configuration w.r.t. the state variables
             */

            _dPlasticMultiplierdStateVariables.second = dPlasticMultiplierdStateVariables;

            _dPlasticMultiplierdStateVariables.first = true;

            addIterationData( &_dPlasticMultiplierdStateVariables );

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

        void residual::setdPreviousPlasticMultiplierdPreviousCauchyStress( const floatVector &dPreviousPlasticMultiplierdPreviousCauchyStress ){
            /*!
             * Set the derivative of the previous plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the previous Cauchy stress
             * 
             * \param &dPreviousPlasticMultiplierdPreviousCauchyStress: The derivative of the previous plastic multiplier in
             *     the current configuration of the plastic configuration w.r.t. the previous Cauchy stress
             */

            _dPreviousPlasticMultiplierdPreviousCauchyStress.second = dPreviousPlasticMultiplierdPreviousCauchyStress;

            _dPreviousPlasticMultiplierdPreviousCauchyStress.first = true;

        }

        void residual::setdPreviousPlasticMultiplierdPreviousF( const floatVector &dPreviousPlasticMultiplierdPreviousF ){
            /*!
             * Set the derivative of the previous plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the previous deformation gradient
             * 
             * \param &dPreviousPlasticMultiplierdPreviousF: The derivative of the previous plastic multiplier in
             *     the current configuration of the plastic configuration w.r.t. the previous deformation gradient
             */

            _dPreviousPlasticMultiplierdPreviousF.second = dPreviousPlasticMultiplierdPreviousF;

            _dPreviousPlasticMultiplierdPreviousF.first = true;

        }

        void residual::setdPreviousPlasticMultiplierdPreviousSubFs( const floatVector &dPreviousPlasticMultiplierdPreviousSubFs ){
            /*!
             * Set the derivative of the previous plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the previous sub-deformation gradients
             * 
             * \param &dPreviousPlasticMultiplierdPreviousSubFs: The derivative of the previous plastic multiplier in
             *     the current configuration of the plastic configuration w.r.t. the previous sub-deformation gradients
             */

            _dPreviousPlasticMultiplierdPreviousSubFs.second = dPreviousPlasticMultiplierdPreviousSubFs;

            _dPreviousPlasticMultiplierdPreviousSubFs.first = true;

        }

        void residual::setdPreviousPlasticMultiplierdPreviousT( const floatType &dPreviousPlasticMultiplierdPreviousT ){
            /*!
             * Set the derivative of the previous plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the previous temperature
             * 
             * \param &dPreviousPlasticMultiplierdPreviousT: The derivative of the previous plastic multiplier in
             *     the current configuration of the plastic configuration w.r.t. the previous temperature
             */

            _dPreviousPlasticMultiplierdPreviousT.second = dPreviousPlasticMultiplierdPreviousT;

            _dPreviousPlasticMultiplierdPreviousT.first = true;

        }

        void residual::setdPreviousPlasticMultiplierdPreviousStateVariables( const floatVector &dPreviousPlasticMultiplierdPreviousStateVariables ){
            /*!
             * Set the derivative of the previous plastic multiplier in the current configuration of the
             * plastic configuration w.r.t. the previous state variables
             * 
             * \param &dPreviousPlasticMultiplierdPreviousStateVariables: The derivative of the previous plastic multiplier in
             *     the current configuration of the plastic configuration w.r.t. the previous state variables
             */

            _dPreviousPlasticMultiplierdPreviousStateVariables.second = dPreviousPlasticMultiplierdPreviousStateVariables;

            _dPreviousPlasticMultiplierdPreviousStateVariables.first = true;

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

            const floatMatrix *dFlowDirectiondCauchyStress;

            const floatMatrix *dFlowDirectiondF;

            const floatMatrix *dFlowDirectiondSubFs;

            if ( isPrevious ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdCauchyStress   = getdPreviousPlasticMultiplierdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdF              = getdPreviousPlasticMultiplierdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdSubFs          = getdPreviousPlasticMultiplierdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdT              = getdPreviousPlasticMultiplierdPreviousT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdStateVariables = getdPreviousPlasticMultiplierdPreviousStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dFlowDirectiondCauchyStress       = getdPreviousFlowDirectiondPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dFlowDirectiondF                  = getdPreviousFlowDirectiondPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dFlowDirectiondSubFs              = getdPreviousFlowDirectiondPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = getPreviousPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( flowDirection = getPreviousFlowDirection( ) );

                setPreviousVelocityGradient( ( *plasticMultiplier ) * ( *flowDirection ) );

                setdPreviousVelocityGradientdPreviousCauchyStress( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdCauchyStress ) + ( *plasticMultiplier ) * ( *dFlowDirectiondCauchyStress ) );

                setdPreviousVelocityGradientdPreviousF( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdF ) + ( *plasticMultiplier ) * ( *dFlowDirectiondF ) );

                setdPreviousVelocityGradientdPreviousSubFs( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdSubFs ) + ( *plasticMultiplier ) * ( *dFlowDirectiondSubFs ) );

                setdPreviousVelocityGradientdPreviousT( ( *flowDirection ) * ( *dPlasticMultiplierdT ) );

                setdPreviousVelocityGradientdPreviousStateVariables( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdStateVariables ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdCauchyStress   = getdPlasticMultiplierdCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdF              = getdPlasticMultiplierdF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdSubFs          = getdPlasticMultiplierdSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdT              = getdPlasticMultiplierdT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdStateVariables = getdPlasticMultiplierdStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dFlowDirectiondCauchyStress       = getdFlowDirectiondCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dFlowDirectiondF                  = getdFlowDirectiondF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dFlowDirectiondSubFs              = getdFlowDirectiondSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = getPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( flowDirection = getFlowDirection( ) );

                setVelocityGradient( ( *plasticMultiplier ) * ( *flowDirection ) );

                setdVelocityGradientdCauchyStress( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdCauchyStress ) + ( *plasticMultiplier ) * ( *dFlowDirectiondCauchyStress ) );

                setdVelocityGradientdF( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdF ) + ( *plasticMultiplier ) * ( *dFlowDirectiondF ) );

                setdVelocityGradientdSubFs( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdSubFs ) + ( *plasticMultiplier ) * ( *dFlowDirectiondSubFs ) );

                setdVelocityGradientdT( ( *flowDirection ) * ( *dPlasticMultiplierdT ) );

                setdVelocityGradientdStateVariables( tardigradeVectorTools::dyadic( *flowDirection, *dPlasticMultiplierdStateVariables ) );

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

        void residual::setdVelocityGradientdCauchyStress( const floatMatrix &dVelocityGradientdCauchyStress ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the Cauchy stress
             * 
             * \param &dVelocityGradientdCauchyStress: The derivative of the velocity gradient in the current
             *     configuration of the plastic configuration w.r.t. the Cauchy stress
             */

            _dVelocityGradientdCauchyStress.second = dVelocityGradientdCauchyStress;

            _dVelocityGradientdCauchyStress.first = true;

            addIterationData( &_dVelocityGradientdCauchyStress );

        }

        void residual::setdVelocityGradientdF( const floatMatrix &dVelocityGradientdF ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the deformation gradient
             * 
             * \param &dVelocityGradientdF: The derivative of the velocity gradient in the current
             *     configuration of the plastic configuration w.r.t. the deformation gradient
             */

            _dVelocityGradientdF.second = dVelocityGradientdF;

            _dVelocityGradientdF.first = true;

            addIterationData( &_dVelocityGradientdF );

        }

        void residual::setdVelocityGradientdSubFs( const floatMatrix &dVelocityGradientdSubFs ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the sub-deformation gradients
             * 
             * \param &dVelocityGradientdSubFs: The derivative of the velocity gradient in the current
             *     configuration of the plastic configuration w.r.t. the sub-deformation gradients
             */

            _dVelocityGradientdSubFs.second = dVelocityGradientdSubFs;

            _dVelocityGradientdSubFs.first = true;

            addIterationData( &_dVelocityGradientdSubFs );

        }

        void residual::setdVelocityGradientdT( const floatVector &dVelocityGradientdT ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the temperature
             * 
             * \param &dVelocityGradientdT: The derivative of the velocity gradient in the current
             *     configuration of the plastic configuration w.r.t. the temperature
             */

            _dVelocityGradientdT.second = dVelocityGradientdT;

            _dVelocityGradientdT.first = true;

            addIterationData( &_dVelocityGradientdT );

        }

        void residual::setdVelocityGradientdStateVariables( const floatMatrix &dVelocityGradientdStateVariables ){
            /*!
             * Set the derivative of the velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the state variables
             * 
             * \param &dVelocityGradientdStateVariables: The derivative of the velocity gradient in the current
             *     configuration of the plastic configuration w.r.t. the state variables
             */

            _dVelocityGradientdStateVariables.second = dVelocityGradientdStateVariables;

            _dVelocityGradientdStateVariables.first = true;

            addIterationData( &_dVelocityGradientdStateVariables );

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

        void residual::setdPreviousVelocityGradientdPreviousCauchyStress( const floatMatrix &dPreviousVelocityGradientdPreviousCauchyStress ){
            /*!
             * Set the derivative of the previous velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the previous Cauchy stress
             * 
             * \param &dPreviousVelocityGradientdPreviousCauchyStress: The derivative of the velocity gradient in the current
             *     w.r.t. the previous Cauchy stress
             */

            _dPreviousVelocityGradientdPreviousCauchyStress.second = dPreviousVelocityGradientdPreviousCauchyStress;

            _dPreviousVelocityGradientdPreviousCauchyStress.first = true;

        }

        void residual::setdPreviousVelocityGradientdPreviousF( const floatMatrix &dPreviousVelocityGradientdPreviousF ){
            /*!
             * Set the derivative of the previous velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the previous deformation gradient
             * 
             * \param &dPreviousVelocityGradientdPreviousF: The derivative of the velocity gradient in the current
             *     w.r.t. the previous deforamtion gradient
             */

            _dPreviousVelocityGradientdPreviousF.second = dPreviousVelocityGradientdPreviousF;

            _dPreviousVelocityGradientdPreviousF.first = true;

        }

        void residual::setdPreviousVelocityGradientdPreviousSubFs( const floatMatrix &dPreviousVelocityGradientdPreviousSubFs ){
            /*!
             * Set the derivative of the previous velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the previous sub-deformation gradients
             * 
             * \param &dPreviousVelocityGradientdPreviousSubFs: The derivative of the velocity gradient in the current
             *     w.r.t. the previous sub-deforamtion gradients
             */

            _dPreviousVelocityGradientdPreviousSubFs.second = dPreviousVelocityGradientdPreviousSubFs;

            _dPreviousVelocityGradientdPreviousSubFs.first = true;

        }

        void residual::setdPreviousVelocityGradientdPreviousT( const floatVector &dPreviousVelocityGradientdPreviousT ){
            /*!
             * Set the derivative of the previous velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the previous temperature
             * 
             * \param &dPreviousVelocityGradientdPreviousT: The derivative of the velocity gradient in the current
             *     w.r.t. the previous temperature
             */

            _dPreviousVelocityGradientdPreviousT.second = dPreviousVelocityGradientdPreviousT;

            _dPreviousVelocityGradientdPreviousT.first = true;

        }

        void residual::setdPreviousVelocityGradientdPreviousStateVariables( const floatMatrix &dPreviousVelocityGradientdPreviousStateVariables ){
            /*!
             * Set the derivative of the previous velocity gradient in the current configuration of the plastic
             * configuration w.r.t. the previous state variables
             * 
             * \param &dPreviousVelocityGradientdPreviousStateVariables: The derivative of the velocity gradient in the current
             *     w.r.t. the previous state variables
             */

            _dPreviousVelocityGradientdPreviousStateVariables.second = dPreviousVelocityGradientdPreviousStateVariables;

            _dPreviousVelocityGradientdPreviousStateVariables.first = true;

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

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdCauchyStress = getdPreviousPlasticMultiplierdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdF = getdPreviousPlasticMultiplierdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdSubFs = getdPreviousPlasticMultiplierdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdT = getdPreviousPlasticMultiplierdPreviousT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdStateVariables = getdPreviousPlasticMultiplierdPreviousStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dHardeningFunctiondStateVariables = getdPreviousHardeningFunctiondPreviousStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = getPreviousPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( hardeningFunction = getPreviousHardeningFunction( ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdCauchyStress = getdPlasticMultiplierdCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdF = getdPlasticMultiplierdF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdSubFs = getdPlasticMultiplierdSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdT = getdPlasticMultiplierdT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPlasticMultiplierdStateVariables = getdPlasticMultiplierdStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dHardeningFunctiondStateVariables = getdHardeningFunctiondStateVariables( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( plasticMultiplier = getPlasticMultiplier( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( hardeningFunction = getHardeningFunction( ) );

            }

            floatVector stateVariableEvolutionRates = { ( *plasticMultiplier ) * ( *hardeningFunction ) };

            floatMatrix dStateVariableEvolutionRatesdCauchyStress = { ( *dPlasticMultiplierdCauchyStress ) * ( *hardeningFunction ) };

            floatMatrix dStateVariableEvolutionRatesdF = { ( *dPlasticMultiplierdF ) * ( *hardeningFunction ) };

            floatMatrix dStateVariableEvolutionRatesdSubFs = { ( *dPlasticMultiplierdSubFs ) * ( *hardeningFunction ) };

            floatVector dStateVariableEvolutionRatesdT = { ( *dPlasticMultiplierdT ) * ( *hardeningFunction ) };

            floatMatrix dStateVariableEvolutionRatesdStateVariables = { ( *dPlasticMultiplierdStateVariables ) * ( *hardeningFunction ) + ( *plasticMultiplier ) * ( *dHardeningFunctiondStateVariables ) };

            if ( isPrevious ){

                setPreviousStateVariableEvolutionRates( stateVariableEvolutionRates );

                setdPreviousStateVariableEvolutionRatesdPreviousCauchyStress( dStateVariableEvolutionRatesdCauchyStress );

                setdPreviousStateVariableEvolutionRatesdPreviousF( dStateVariableEvolutionRatesdF );

                setdPreviousStateVariableEvolutionRatesdPreviousSubFs( dStateVariableEvolutionRatesdSubFs );

                setdPreviousStateVariableEvolutionRatesdPreviousT( dStateVariableEvolutionRatesdT );

                setdPreviousStateVariableEvolutionRatesdPreviousStateVariables( dStateVariableEvolutionRatesdStateVariables );

            }
            else{

                setStateVariableEvolutionRates( stateVariableEvolutionRates );

                setdStateVariableEvolutionRatesdCauchyStress( dStateVariableEvolutionRatesdCauchyStress );

                setdStateVariableEvolutionRatesdF( dStateVariableEvolutionRatesdF );

                setdStateVariableEvolutionRatesdSubFs( dStateVariableEvolutionRatesdSubFs );

                setdStateVariableEvolutionRatesdT( dStateVariableEvolutionRatesdT );

                setdStateVariableEvolutionRatesdStateVariables( dStateVariableEvolutionRatesdStateVariables );

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

        void residual::setdStateVariableEvolutionRatesdCauchyStress( const floatMatrix &dStateVariableEvolutionRatesdCauchyStress ){
            /*!
             * Set the derivative of the state variable evolution rates w.r.t. the Cauchy stress
             * 
             * \param &dStateVariableEvolutionRatesdCauchyStress: The value of the derivative
             */

            _dStateVariableEvolutionRatesdCauchyStress.second = dStateVariableEvolutionRatesdCauchyStress;

            _dStateVariableEvolutionRatesdCauchyStress.first = true;

            addIterationData( &_dStateVariableEvolutionRatesdCauchyStress );

        }

        void residual::setdStateVariableEvolutionRatesdF( const floatMatrix &dStateVariableEvolutionRatesdF ){
            /*!
             * Set the derivative of the state variable evolution rates w.r.t. the deformation gradient
             * 
             * \param &dStateVariableEvolutionRatesdF: The value of the derivative
             */

            _dStateVariableEvolutionRatesdF.second = dStateVariableEvolutionRatesdF;

            _dStateVariableEvolutionRatesdF.first = true;

            addIterationData( &_dStateVariableEvolutionRatesdF );

        }

        void residual::setdStateVariableEvolutionRatesdSubFs( const floatMatrix &dStateVariableEvolutionRatesdSubFs ){
            /*!
             * Set the derivative of the state variable evolution rates w.r.t. the sub-deformation gradients
             * 
             * \param &dStateVariableEvolutionRatesdSubFs: The value of the derivative
             */

            _dStateVariableEvolutionRatesdSubFs.second = dStateVariableEvolutionRatesdSubFs;

            _dStateVariableEvolutionRatesdSubFs.first = true;

            addIterationData( &_dStateVariableEvolutionRatesdSubFs );

        }

        void residual::setdStateVariableEvolutionRatesdT( const floatVector &dStateVariableEvolutionRatesdT ){
            /*!
             * Set the derivative of the state variable evolution rates w.r.t. the temperature
             * 
             * \param &dStateVariableEvolutionRatesdT: The value of the derivative
             */

            _dStateVariableEvolutionRatesdT.second = dStateVariableEvolutionRatesdT;

            _dStateVariableEvolutionRatesdT.first = true;

            addIterationData( &_dStateVariableEvolutionRatesdT );

        }

        void residual::setdStateVariableEvolutionRatesdStateVariables( const floatMatrix &dStateVariableEvolutionRatesdStateVariables ){
            /*!
             * Set the derivative of the state variable evolution rates w.r.t. the state variables
             * 
             * \param &dStateVariableEvolutionRatesdStateVariables: The value of the derivative
             */

            _dStateVariableEvolutionRatesdStateVariables.second = dStateVariableEvolutionRatesdStateVariables;

            _dStateVariableEvolutionRatesdStateVariables.first = true;

            addIterationData( &_dStateVariableEvolutionRatesdStateVariables );

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

        void residual::setdPreviousStateVariableEvolutionRatesdPreviousCauchyStress( const floatMatrix &dPreviousStateVariableEvolutionRatesdPreviousCauchyStress ){
            /*!
             * Set the derivative of the state variable evolution rates w.r.t. the Cauchy stress
             * 
             * \param &dPreviousStateVariableEvolutionRatesdPreviousCauchyStress: The value of the derivative
             */

            _dPreviousStateVariableEvolutionRatesdPreviousCauchyStress.second = dPreviousStateVariableEvolutionRatesdPreviousCauchyStress;

            _dPreviousStateVariableEvolutionRatesdPreviousCauchyStress.first = true;

        }

        void residual::setdPreviousStateVariableEvolutionRatesdPreviousF( const floatMatrix &dPreviousStateVariableEvolutionRatesdPreviousF ){
            /*!
             * Set the derivative of the state variable evolution rates w.r.t. the deformation gradient
             * 
             * \param &dPreviousStateVariableEvolutionRatesdPreviousF: The value of the derivative
             */

            _dPreviousStateVariableEvolutionRatesdPreviousF.second = dPreviousStateVariableEvolutionRatesdPreviousF;

            _dPreviousStateVariableEvolutionRatesdPreviousF.first = true;

        }

        void residual::setdPreviousStateVariableEvolutionRatesdPreviousSubFs( const floatMatrix &dPreviousStateVariableEvolutionRatesdPreviousSubFs ){
            /*!
             * Set the derivative of the state variable evolution rates w.r.t. the sub-deformation gradients
             * 
             * \param &dPreviousStateVariableEvolutionRatesdPreviousSubFs: The value of the derivative
             */

            _dPreviousStateVariableEvolutionRatesdPreviousSubFs.second = dPreviousStateVariableEvolutionRatesdPreviousSubFs;

            _dPreviousStateVariableEvolutionRatesdPreviousSubFs.first = true;

        }

        void residual::setdPreviousStateVariableEvolutionRatesdPreviousT( const floatVector &dPreviousStateVariableEvolutionRatesdPreviousT ){
            /*!
             * Set the derivative of the state variable evolution rates w.r.t. the temperature
             * 
             * \param &dPreviousStateVariableEvolutionRatesdPreviousT: The value of the derivative
             */

            _dPreviousStateVariableEvolutionRatesdPreviousT.second = dPreviousStateVariableEvolutionRatesdPreviousT;

            _dPreviousStateVariableEvolutionRatesdPreviousT.first = true;

        }

        void residual::setdPreviousStateVariableEvolutionRatesdPreviousStateVariables( const floatMatrix &dPreviousStateVariableEvolutionRatesdPreviousStateVariables ){
            /*!
             * Set the derivative of the state variable evolution rates w.r.t. the state variables
             * 
             * \param &dPreviousStateVariableEvolutionRatesdPreviousStateVariables: The value of the derivative
             */

            _dPreviousStateVariableEvolutionRatesdPreviousStateVariables.second = dPreviousStateVariableEvolutionRatesdPreviousStateVariables;

            _dPreviousStateVariableEvolutionRatesdPreviousStateVariables.first = true;

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

            const floatVector *velocityGradient;

            const floatVector *previousVelocityGradient;

            floatVector previousPlasticDeformationGradient;

            floatVector dFp;

            floatVector plasticDeformationGradient;

            floatMatrix dFdL;

            TARDIGRADE_ERROR_TOOLS_CATCH( velocityGradient = getVelocityGradient( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousVelocityGradient = getPreviousVelocityGradient( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousPlasticDeformationGradient = hydra->getPreviousConfiguration( *getPlasticConfigurationIndex( ) ) );

            if ( setPreviousDerivatives ){

                floatMatrix ddFdPreviousF;

                floatMatrix dFdPreviousF;

                floatMatrix dFdPreviousL;

                TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::evolveF( *hydra->getDeltaTime( ), previousPlasticDeformationGradient, *previousVelocityGradient, *velocityGradient, dFp, plasticDeformationGradient, dFdL, ddFdPreviousF, dFdPreviousF, dFdPreviousL, 1 - ( *getIntegrationParameter( ) ), 1 ) );

                setdPlasticDeformationGradientdPreviousCauchyStress( tardigradeVectorTools::dot( dFdPreviousL, *getdPreviousVelocityGradientdPreviousCauchyStress( ) ) );

                setdPlasticDeformationGradientdPreviousF( tardigradeVectorTools::dot( dFdPreviousL, *getdPreviousVelocityGradientdPreviousF( ) ) );

                floatMatrix dPlasticDeformationGradientdPreviousSubFs = tardigradeVectorTools::dot( dFdPreviousL, *getdPreviousVelocityGradientdPreviousSubFs( ) );

                for ( unsigned int i = 0; i < dFdPreviousF.size( ); i++ ){

                    for ( unsigned int j = 0; j < dFdPreviousF[ i ].size( ); j++ ){

                        dPlasticDeformationGradientdPreviousSubFs[ i ][ j + ( *getPlasticConfigurationIndex( ) - 1 ) * plasticDeformationGradient.size( ) ] += dFdPreviousF[ i ][ j ];

                    }

                }

                setdPlasticDeformationGradientdPreviousSubFs( dPlasticDeformationGradientdPreviousSubFs );

                setdPlasticDeformationGradientdPreviousT( tardigradeVectorTools::dot( dFdPreviousL, *getdPreviousVelocityGradientdPreviousT( ) ) );

                setdPlasticDeformationGradientdPreviousStateVariables( tardigradeVectorTools::dot( dFdPreviousL, *getdPreviousVelocityGradientdPreviousStateVariables( ) ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::evolveF( *hydra->getDeltaTime( ), previousPlasticDeformationGradient, *previousVelocityGradient, *velocityGradient, dFp, plasticDeformationGradient, dFdL, 1 - ( *getIntegrationParameter( ) ), 1 ) );
            }

            setPlasticDeformationGradient( plasticDeformationGradient );

            setdPlasticDeformationGradientdCauchyStress( tardigradeVectorTools::dot( dFdL, *getdVelocityGradientdCauchyStress( ) ) );

            setdPlasticDeformationGradientdF( tardigradeVectorTools::dot( dFdL, *getdVelocityGradientdF( ) ) );

            setdPlasticDeformationGradientdSubFs( tardigradeVectorTools::dot( dFdL, *getdVelocityGradientdSubFs( ) ) );

            setdPlasticDeformationGradientdT( tardigradeVectorTools::dot( dFdL, *getdVelocityGradientdT( ) ) );

            setdPlasticDeformationGradientdStateVariables( tardigradeVectorTools::dot( dFdL, *getdVelocityGradientdStateVariables( ) ) );

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

        void residual::setdPlasticDeformationGradientdCauchyStress( const floatMatrix &dPlasticDeformationGradientdCauchyStress ){
            /*!
             * Set the derivative of the plastic deformation gradient w.r.t. the Cauchy stress
             *
             * \param &dPlasticDeformationGradientdCauchyStress: The value of the derivative
             */

            _dPlasticDeformationGradientdCauchyStress.second = dPlasticDeformationGradientdCauchyStress;

            _dPlasticDeformationGradientdCauchyStress.first = true;

            addIterationData( &_dPlasticDeformationGradientdCauchyStress );

        }

        void residual::setdPlasticDeformationGradientdF( const floatMatrix &dPlasticDeformationGradientdF ){
            /*!
             * Set the derivative of the plastic deformation gradient w.r.t. the deformation gradient
             *
             * \param &dPlasticDeformationGradientdF: The value of the derivative
             */

            _dPlasticDeformationGradientdF.second = dPlasticDeformationGradientdF;

            _dPlasticDeformationGradientdF.first = true;

            addIterationData( &_dPlasticDeformationGradientdF );

        }

        void residual::setdPlasticDeformationGradientdSubFs( const floatMatrix &dPlasticDeformationGradientdSubFs ){
            /*!
             * Set the derivative of the plastic deformation gradient w.r.t. the sub-deformation gradients
             *
             * \param &dPlasticDeformationGradientdSubFs: The value of the derivative
             */

            _dPlasticDeformationGradientdSubFs.second = dPlasticDeformationGradientdSubFs;

            _dPlasticDeformationGradientdSubFs.first = true;

            addIterationData( &_dPlasticDeformationGradientdSubFs );

        }

        void residual::setdPlasticDeformationGradientdT( const floatVector &dPlasticDeformationGradientdT ){
            /*!
             * Set the derivative of the plastic deformation gradient w.r.t. the temperature
             *
             * \param &dPlasticDeformationGradientdT: The value of the derivative
             */

            _dPlasticDeformationGradientdT.second = dPlasticDeformationGradientdT;

            _dPlasticDeformationGradientdT.first = true;

            addIterationData( &_dPlasticDeformationGradientdT );

        }

        void residual::setdPlasticDeformationGradientdStateVariables( const floatMatrix &dPlasticDeformationGradientdStateVariables ){
            /*!
             * Set the derivative of the plastic deformation gradient w.r.t. the state variables
             *
             * \param &dPlasticDeformationGradientdStateVariables: The value of the derivative
             */

            _dPlasticDeformationGradientdStateVariables.second = dPlasticDeformationGradientdStateVariables;

            _dPlasticDeformationGradientdStateVariables.first = true;

            addIterationData( &_dPlasticDeformationGradientdStateVariables );

        }

        void residual::setdPlasticDeformationGradientdPreviousCauchyStress( const floatMatrix &dPlasticDeformationGradientdPreviousCauchyStress ){
            /*!
             * Set the derivative of the plastic deformation gradient w.r.t. the previous Cauchy stress
             *
             * \param &dPlasticDeformationGradientdPreviousCauchyStress: The value of the derivative
             */

            _dPlasticDeformationGradientdPreviousCauchyStress.second = dPlasticDeformationGradientdPreviousCauchyStress;

            _dPlasticDeformationGradientdPreviousCauchyStress.first = true;

            addIterationData( &_dPlasticDeformationGradientdPreviousCauchyStress );

        }

        void residual::setdPlasticDeformationGradientdPreviousF( const floatMatrix &dPlasticDeformationGradientdPreviousF ){
            /*!
             * Set the derivative of the plastic deformation gradient w.r.t. the previous deformation gradient
             *
             * \param &dPlasticDeformationGradientdPreviousF: The value of the derivative
             */

            _dPlasticDeformationGradientdPreviousF.second = dPlasticDeformationGradientdPreviousF;

            _dPlasticDeformationGradientdPreviousF.first = true;

            addIterationData( &_dPlasticDeformationGradientdPreviousF );

        }

        void residual::setdPlasticDeformationGradientdPreviousSubFs( const floatMatrix &dPlasticDeformationGradientdPreviousSubFs ){
            /*!
             * Set the derivative of the plastic deformation gradient w.r.t. the previous sub-deformation gradients
             *
             * \param &dPlasticDeformationGradientdPreviousSubFs: The value of the derivative
             */

            _dPlasticDeformationGradientdPreviousSubFs.second = dPlasticDeformationGradientdPreviousSubFs;

            _dPlasticDeformationGradientdPreviousSubFs.first = true;

            addIterationData( &_dPlasticDeformationGradientdPreviousSubFs );

        }

        void residual::setdPlasticDeformationGradientdPreviousT( const floatVector &dPlasticDeformationGradientdPreviousT ){
            /*!
             * Set the derivative of the plastic deformation gradient w.r.t. the previous temperature
             *
             * \param &dPlasticDeformationGradientdPreviousT: The value of the derivative
             */

            _dPlasticDeformationGradientdPreviousT.second = dPlasticDeformationGradientdPreviousT;

            _dPlasticDeformationGradientdPreviousT.first = true;

            addIterationData( &_dPlasticDeformationGradientdPreviousT );

        }

        void residual::setdPlasticDeformationGradientdPreviousStateVariables( const floatMatrix &dPlasticDeformationGradientdPreviousStateVariables ){
            /*!
             * Set the derivative of the plastic deformation gradient w.r.t. the previous state variables
             *
             * \param &dPlasticDeformationGradientdPreviousStateVariables: The value of the derivative
             */

            _dPlasticDeformationGradientdPreviousStateVariables.second = dPlasticDeformationGradientdPreviousStateVariables;

            _dPlasticDeformationGradientdPreviousStateVariables.first = true;

            addIterationData( &_dPlasticDeformationGradientdPreviousStateVariables );

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

            const floatVector *stateVariableEvolutionRates;

            const floatVector *previousStateVariableEvolutionRates;

            const floatVector *previousStateVariables;

            const floatMatrix *dStateVariableEvolutionRatesdCauchyStress;

            const floatMatrix *dStateVariableEvolutionRatesdF;

            const floatMatrix *dStateVariableEvolutionRatesdSubFs;

            const floatVector *dStateVariableEvolutionRatesdT;

            const floatMatrix *dStateVariableEvolutionRatesdStateVariables;

            const floatMatrix *dPreviousStateVariableEvolutionRatesdPreviousCauchyStress = NULL;

            const floatMatrix *dPreviousStateVariableEvolutionRatesdPreviousF = NULL;

            const floatMatrix *dPreviousStateVariableEvolutionRatesdPreviousSubFs = NULL;

            const floatVector *dPreviousStateVariableEvolutionRatesdPreviousT = NULL;

            const floatMatrix *dPreviousStateVariableEvolutionRatesdPreviousStateVariables = NULL;

            floatVector deltaPlasticStateVariables;

            floatVector plasticStateVariables;

            if ( setPreviousDerivatives ){

                TARDIGRADE_ERROR_TOOLS_CATCH( dPreviousStateVariableEvolutionRatesdPreviousCauchyStress = getdPreviousStateVariableEvolutionRatesdPreviousCauchyStress( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPreviousStateVariableEvolutionRatesdPreviousF = getdPreviousStateVariableEvolutionRatesdPreviousF( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPreviousStateVariableEvolutionRatesdPreviousSubFs = getdPreviousStateVariableEvolutionRatesdPreviousSubFs( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPreviousStateVariableEvolutionRatesdPreviousT = getdPreviousStateVariableEvolutionRatesdPreviousT( ) );

                TARDIGRADE_ERROR_TOOLS_CATCH( dPreviousStateVariableEvolutionRatesdPreviousStateVariables = getdPreviousStateVariableEvolutionRatesdPreviousStateVariables( ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( dStateVariableEvolutionRatesdCauchyStress = getdStateVariableEvolutionRatesdCauchyStress( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( dStateVariableEvolutionRatesdF = getdStateVariableEvolutionRatesdF( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( dStateVariableEvolutionRatesdSubFs = getdStateVariableEvolutionRatesdSubFs( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( dStateVariableEvolutionRatesdT = getdStateVariableEvolutionRatesdT( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( dStateVariableEvolutionRatesdStateVariables = getdStateVariableEvolutionRatesdStateVariables( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( stateVariableEvolutionRates = getStateVariableEvolutionRates( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousStateVariableEvolutionRates = getPreviousStateVariableEvolutionRates( ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( previousStateVariables = getPreviousStateVariables( ) );

            floatMatrix dXidXidot;

            if ( setPreviousDerivatives ){

                floatMatrix dXidXidotp;

                TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolution( *hydra->getDeltaTime( ), *previousStateVariables, *previousStateVariableEvolutionRates, *stateVariableEvolutionRates, deltaPlasticStateVariables, plasticStateVariables, dXidXidot, dXidXidotp, ( 1 - *getIntegrationParameter( ) ) ) );

                setdPlasticStateVariablesdPreviousCauchyStress( tardigradeVectorTools::dot( dXidXidotp, *dPreviousStateVariableEvolutionRatesdPreviousCauchyStress ) );

                setdPlasticStateVariablesdPreviousF( tardigradeVectorTools::dot( dXidXidotp, *dPreviousStateVariableEvolutionRatesdPreviousF ) );

                setdPlasticStateVariablesdPreviousSubFs( tardigradeVectorTools::dot( dXidXidotp, *dPreviousStateVariableEvolutionRatesdPreviousSubFs ) );

                setdPlasticStateVariablesdPreviousT( tardigradeVectorTools::dot( dXidXidotp, *dPreviousStateVariableEvolutionRatesdPreviousT ) );

                setdPlasticStateVariablesdPreviousStateVariables( tardigradeVectorTools::dot( dXidXidotp, *dPreviousStateVariableEvolutionRatesdPreviousStateVariables ) + tardigradeVectorTools::eye< floatType >( previousStateVariables->size( ) ) );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolution( *hydra->getDeltaTime( ), *previousStateVariables, *previousStateVariableEvolutionRates, *stateVariableEvolutionRates, deltaPlasticStateVariables, plasticStateVariables, dXidXidot, ( 1 - *getIntegrationParameter( ) ) ) );

            }

            setdPlasticStateVariablesdCauchyStress( tardigradeVectorTools::dot( dXidXidot, *dStateVariableEvolutionRatesdCauchyStress ) );

            setdPlasticStateVariablesdF( tardigradeVectorTools::dot( dXidXidot, *dStateVariableEvolutionRatesdF ) );

            setdPlasticStateVariablesdSubFs( tardigradeVectorTools::dot( dXidXidot, *dStateVariableEvolutionRatesdSubFs ) );

            setdPlasticStateVariablesdT( tardigradeVectorTools::dot( dXidXidot, *dStateVariableEvolutionRatesdT ) );

            setdPlasticStateVariablesdStateVariables( tardigradeVectorTools::dot( dXidXidot, *dStateVariableEvolutionRatesdStateVariables ) );

            setPlasticStateVariables( plasticStateVariables );

        }


        void residual::setPlasticStateVariables( const floatVector &plasticStateVariables ){
            /*!
             * Set the plastic state variables
             * 
             * \param &plasticStateVariables
             */

            _plasticStateVariables.second = plasticStateVariables;

            _plasticStateVariables.first = true;

            addIterationData( &_plasticStateVariables );

        }

        void residual::setdPlasticStateVariablesdCauchyStress( const floatMatrix &dPlasticStateVariablesdCauchyStress ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the Cauchy Stress
             * 
             * \param &dPlasticStateVariablesdCauchyStress: The derivative of the plastic state variables w.r.t. the Cauchy stress
             */

            _dPlasticStateVariablesdCauchyStress.second = dPlasticStateVariablesdCauchyStress;

            _dPlasticStateVariablesdCauchyStress.first = true;

            addIterationData( &_dPlasticStateVariablesdCauchyStress );

        }

        void residual::setdPlasticStateVariablesdF( const floatMatrix &dPlasticStateVariablesdF ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the deformation gradient
             * 
             * \param &dPlasticStateVariablesdF: The derivative of the plastic state variables w.r.t. the deformation gradient
             */

            _dPlasticStateVariablesdF.second = dPlasticStateVariablesdF;

            _dPlasticStateVariablesdF.first = true;

            addIterationData( &_dPlasticStateVariablesdF );

        }

        void residual::setdPlasticStateVariablesdSubFs( const floatMatrix &dPlasticStateVariablesdSubFs ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the sub-deformation gradients
             * 
             * \param &dPlasticStateVariablesdSubFs: The derivative of the plastic state variables w.r.t. the sub-deformation gradients
             */

            _dPlasticStateVariablesdSubFs.second = dPlasticStateVariablesdSubFs;

            _dPlasticStateVariablesdSubFs.first = true;

            addIterationData( &_dPlasticStateVariablesdSubFs );

        }

        void residual::setdPlasticStateVariablesdT( const floatVector &dPlasticStateVariablesdT ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the temperature
             * 
             * \param &dPlasticStateVariablesdT: The derivative of the plastic state variables w.r.t. the temperature
             */

            _dPlasticStateVariablesdT.second = dPlasticStateVariablesdT;

            _dPlasticStateVariablesdT.first = true;

            addIterationData( &_dPlasticStateVariablesdT );

        }

        void residual::setdPlasticStateVariablesdStateVariables( const floatMatrix &dPlasticStateVariablesdStateVariables ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the state variables
             * 
             * \param &dPlasticStateVariablesdStateVariables: The derivative of the plastic state variables w.r.t. the state variables
             */

            _dPlasticStateVariablesdStateVariables.second = dPlasticStateVariablesdStateVariables;

            _dPlasticStateVariablesdStateVariables.first = true;

            addIterationData( &_dPlasticStateVariablesdStateVariables );

        }

        void residual::setdPlasticStateVariablesdPreviousCauchyStress( const floatMatrix &dPlasticStateVariablesdPreviousCauchyStress ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the previous Cauchy Stress
             * 
             * \param &dPlasticStateVariablesdPreviousCauchyStress: The derivative of the plastic state variables w.r.t. the previous Cauchy stress
             */

            _dPlasticStateVariablesdPreviousCauchyStress.second = dPlasticStateVariablesdPreviousCauchyStress;

            _dPlasticStateVariablesdPreviousCauchyStress.first = true;

            addIterationData( &_dPlasticStateVariablesdPreviousCauchyStress );

        }

        void residual::setdPlasticStateVariablesdPreviousF( const floatMatrix &dPlasticStateVariablesdPreviousF ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the previous deformation gradient
             * 
             * \param &dPlasticStateVariablesdPreviousF: The derivative of the plastic state variables w.r.t. the previous deformation gradient
             */

            _dPlasticStateVariablesdPreviousF.second = dPlasticStateVariablesdPreviousF;

            _dPlasticStateVariablesdPreviousF.first = true;

            addIterationData( &_dPlasticStateVariablesdPreviousF );

        }

        void residual::setdPlasticStateVariablesdPreviousSubFs( const floatMatrix &dPlasticStateVariablesdPreviousSubFs ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the previous sub-deformation gradients
             * 
             * \param &dPlasticStateVariablesdPreviousSubFs: The derivative of the plastic state variables w.r.t. the previous sub-deformation gradients
             */

            _dPlasticStateVariablesdPreviousSubFs.second = dPlasticStateVariablesdPreviousSubFs;

            _dPlasticStateVariablesdPreviousSubFs.first = true;

            addIterationData( &_dPlasticStateVariablesdPreviousSubFs );

        }

        void residual::setdPlasticStateVariablesdPreviousT( const floatVector &dPlasticStateVariablesdPreviousT ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the previous temperature
             * 
             * \param &dPlasticStateVariablesdPreviousT: The derivative of the plastic state variables w.r.t. the previous temperature
             */

            _dPlasticStateVariablesdPreviousT.second = dPlasticStateVariablesdPreviousT;

            _dPlasticStateVariablesdPreviousT.first = true;

            addIterationData( &_dPlasticStateVariablesdPreviousT );

        }

        void residual::setdPlasticStateVariablesdPreviousStateVariables( const floatMatrix &dPlasticStateVariablesdPreviousStateVariables ){
            /*!
             * Set the derivative of the plastic state variables w.r.t. the previous state variables
             * 
             * \param &dPlasticStateVariablesdPreviousStateVariables: The derivative of the plastic state variables w.r.t. the previous state variables
             */

            _dPlasticStateVariablesdPreviousStateVariables.second = dPlasticStateVariablesdPreviousStateVariables;

            _dPlasticStateVariablesdPreviousStateVariables.first = true;

            addIterationData( &_dPlasticStateVariablesdPreviousStateVariables );

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

            floatMatrix jacobian( *getNumEquations( ), floatVector( hydra->getUnknownVector( )->size( ), 0 ) );

            // Set the derivatives
            getdPlasticDeformationGradientdCauchyStress( );

            getdPlasticStateVariablesdCauchyStress( );

            // Get the Jacobians of the plastic deformation gradient
            for ( unsigned int i = 0; i < getPlasticDeformationGradient( )->size( ); i++ ){
                unsigned int row = i;

                // Set the Jacobian with respect to the Cauchy stress
                for ( unsigned int j = 0; j < ( *getdPlasticDeformationGradientdCauchyStress( ) )[ i ].size( ); j++ ){
                    unsigned int col = j;

                    jacobian[ row ][ col ] += ( *getdPlasticDeformationGradientdCauchyStress( ) )[ i ][ j ];

                }

                // Set the Jacobian with respect to the sub-configurations
                jacobian[ i ][ ( *getdPlasticDeformationGradientdCauchyStress( ) )[ i ].size( ) + i ] -= 1;
                for ( unsigned int j = 0; j < ( *getdPlasticDeformationGradientdSubFs( ) )[ i ].size( ); j++ ){
                    unsigned int col = ( *getdPlasticDeformationGradientdCauchyStress( ) )[ i ].size( ) + j;

                    jacobian[ row ][ col ] += ( *getdPlasticDeformationGradientdSubFs( ) )[ i ][ j ];

                }

                // Set the Jacobian with respect to the state variables
                for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){
                    unsigned int col = ( *getdPlasticDeformationGradientdCauchyStress( ) )[ i ].size( ) + ( *getdPlasticDeformationGradientdSubFs( ) )[ i ].size( ) + *ind;

                    jacobian[ row ][ col ] += ( *getdPlasticDeformationGradientdStateVariables( ) )[ i ][ ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ];

                }

            }

            // Get the Jacobians of the plastic state variables
            for ( unsigned int i = 0; i < getPlasticStateVariables( )->size( ); i++ ){
                unsigned int row = getPlasticDeformationGradient( )->size( ) + i;

                // Set the Jacobian with respect to the Cauchy stress
                for ( unsigned int j = 0; j < ( *getdPlasticStateVariablesdCauchyStress( ) )[ i ].size( ); j++ ){
                    unsigned int col = j;

                    jacobian[ row ][ col ] += ( *getdPlasticStateVariablesdCauchyStress( ) )[ i ][ j ];

                }

                // Set the Jacobian with respect to the other configurations
                for ( unsigned int j = 0; j < ( *getdPlasticStateVariablesdSubFs( ) )[ i ].size( ); j++ ){
                    unsigned int col = ( *getdPlasticDeformationGradientdCauchyStress( ) )[ i ].size( ) + j;

                    jacobian[ row ][ col ] += ( *getdPlasticStateVariablesdSubFs( ) )[ i ][ j ];

                }

                // Set the Jacobian with respect to the state variables
                jacobian[ row ][ ( *getdPlasticStateVariablesdCauchyStress( ) )[ i ].size( ) + ( *getdPlasticStateVariablesdSubFs( ) )[ i ].size( ) + ( *getStateVariableIndices( ) )[ i ] ] -= 1;
                for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){
                    unsigned int col = ( *getdPlasticStateVariablesdCauchyStress( ) )[ i ].size( ) + ( *getdPlasticStateVariablesdSubFs( ) )[ i ].size( ) + *ind;

                    jacobian[ row ][ col ] += ( *getdPlasticStateVariablesdStateVariables( ) )[ i ][ ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ];

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
            getdPlasticDeformationGradientdT( );

            getdPlasticStateVariablesdT( );

            // Get the Jacobians of the plastic deformation gradient
            for ( unsigned int i = 0; i < getPlasticDeformationGradient( )->size( ); i++ ){

                dRdT[ i ] = ( *getdPlasticDeformationGradientdT( ) )[ i ];

            }

            // Get the Jacobians of the plastic state variables
            for ( unsigned int i = 0; i < getPlasticStateVariables( )->size( ); i++ ){

                dRdT[ i + getPlasticDeformationGradient( )->size( ) ] += ( *getdPlasticStateVariablesdT( ) )[ i ];

            }

            setdRdT( dRdT );

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient.
             */

            floatMatrix dRdF( *getNumEquations( ), floatVector( hydra->getDeformationGradient( )->size( ), 0 ) );

            // Set the derivatives
            getdPlasticDeformationGradientdF( );

            getdPlasticStateVariablesdF( );

            // Get the Jacobians of the plastic deformation gradient
            for ( unsigned int i = 0; i < getPlasticDeformationGradient( )->size( ); i++ ){

                for ( unsigned int j = 0; j < ( *getdPlasticDeformationGradientdF( ) )[ i ].size( ); j++ ){

                    dRdF[ i ][ j ] = ( *getdPlasticDeformationGradientdF( ) )[ i ][ j ];

                }

            }

            // Get the Jacobians of the plastic state variables
            for ( unsigned int i = 0; i < getPlasticStateVariables( )->size( ); i++ ){

                for ( unsigned int j = 0; j < ( *getdPlasticStateVariablesdF( ) )[ i ].size( ); j++ ){

                    dRdF[ i + getPlasticDeformationGradient( )->size( ) ][ j ] += ( *getdPlasticStateVariablesdF( ) )[ i ][ j ];

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

        const floatVector* residual::getdPlasticMultiplierdCauchyStress( ){
            /*!
             * Get the derivative of the plastic multiplier w.r.t. the Cauchy stress
             */

            if ( !_dPlasticMultiplierdCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticMultiplierdCauchyStress( ) );

            }

            return &_dPlasticMultiplierdCauchyStress.second;

        }

        const floatVector* residual::getdPlasticMultiplierdF( ){
            /*!
             * Get the derivative of the plastic multiplier w.r.t. the deformation gradient
             */

            if ( !_dPlasticMultiplierdF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticMultiplierdF( ) );

            }

            return &_dPlasticMultiplierdF.second;

        }

        const floatVector* residual::getdPlasticMultiplierdSubFs( ){
            /*!
             * Get the derivative of the plastic multiplier w.r.t. the sub-deformation gradients
             */

            if ( !_dPlasticMultiplierdSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticMultiplierdSubFs( ) );

            }

            return &_dPlasticMultiplierdSubFs.second;

        }

        const floatType* residual::getdPlasticMultiplierdT( ){
            /*!
             * Get the derivative of the plastic multiplier w.r.t. the temperature
             */

            if ( !_dPlasticMultiplierdT.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticMultiplierdT( ) );

            }

            return &_dPlasticMultiplierdT.second;

        }

        const floatVector* residual::getdPlasticMultiplierdStateVariables( ){
            /*!
             * Get the derivative of the plastic multiplier w.r.t. the state variables
             */

            if ( !_dPlasticMultiplierdStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticMultiplierdStateVariables( ) );

            }

            return &_dPlasticMultiplierdStateVariables.second;

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

        const floatMatrix* residual::getdVelocityGradientdCauchyStress( ){
            /*!
             * Get the derivative of the velocity gradient w.r.t. the Cauchy stress
             */

            if ( !_dVelocityGradientdCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdVelocityGradientdCauchyStress( ) );

            }

            return &_dVelocityGradientdCauchyStress.second;

        }

        const floatMatrix* residual::getdVelocityGradientdF( ){
            /*!
             * Get the derivative of the velocity gradient w.r.t. the deformation gradient
             */

            if ( !_dVelocityGradientdF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdVelocityGradientdF( ) );

            }

            return &_dVelocityGradientdF.second;

        }

        const floatMatrix* residual::getdVelocityGradientdSubFs( ){
            /*!
             * Get the derivative of the velocity gradient w.r.t. the sub-deformation gradients
             */

            if ( !_dVelocityGradientdSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdVelocityGradientdSubFs( ) );

            }

            return &_dVelocityGradientdSubFs.second;

        }

        const floatVector* residual::getdVelocityGradientdT( ){
            /*!
             * Get the derivative of the velocity gradient w.r.t. the temperature
             */

            if ( !_dVelocityGradientdT.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdVelocityGradientdT( ) );

            }

            return &_dVelocityGradientdT.second;

        }

        const floatMatrix* residual::getdVelocityGradientdStateVariables( ){
            /*!
             * Get the derivative of the velocity gradient w.r.t. the state variables
             */

            if ( !_dVelocityGradientdStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdVelocityGradientdStateVariables( ) );

            }

            return &_dVelocityGradientdStateVariables.second;

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

        const floatMatrix* residual::getdStateVariableEvolutionRatesdCauchyStress( ){
            /*!
             * Get the derivative of the state variable evolution rate w.r.t. the Cauchy stress
             */

            if ( !_dStateVariableEvolutionRatesdCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdStateVariableEvolutionRatesdCauchyStress( ) );

            }

            return &_dStateVariableEvolutionRatesdCauchyStress.second;

        }

        const floatMatrix* residual::getdStateVariableEvolutionRatesdF( ){
            /*!
             * Get the derivative of the state variable evolution rate w.r.t. the deformation gradient
             */

            if ( !_dStateVariableEvolutionRatesdF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdStateVariableEvolutionRatesdF( ) );

            }

            return &_dStateVariableEvolutionRatesdF.second;

        }

        const floatMatrix* residual::getdStateVariableEvolutionRatesdSubFs( ){
            /*!
             * Get the derivative of the state variable evolution rate w.r.t. the sub-deformation gradients
             */

            if ( !_dStateVariableEvolutionRatesdSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdStateVariableEvolutionRatesdSubFs( ) );

            }

            return &_dStateVariableEvolutionRatesdSubFs.second;

        }

        const floatVector* residual::getdStateVariableEvolutionRatesdT( ){
            /*!
             * Get the derivative of the state variable evolution rate w.r.t. the sub-deformation gradients
             */

            if ( !_dStateVariableEvolutionRatesdT.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdStateVariableEvolutionRatesdT( ) );

            }

            return &_dStateVariableEvolutionRatesdT.second;

        }

        const floatMatrix* residual::getdStateVariableEvolutionRatesdStateVariables( ){
            /*!
             * Get the derivative of the state variable evolution rate w.r.t. the state variables
             */

            if ( !_dStateVariableEvolutionRatesdStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdStateVariableEvolutionRatesdStateVariables( ) );

            }

            return &_dStateVariableEvolutionRatesdStateVariables.second;

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

        const floatMatrix* residual::getdPlasticDeformationGradientdCauchyStress( ){
            /*!
             * Get the derivative of the plastic deformation gradient w.r.t. the Cauchy stress
             */

            if ( !_dPlasticDeformationGradientdCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticDeformationGradientdCauchyStress( ) );

            }

            return &_dPlasticDeformationGradientdCauchyStress.second;

        }

        const floatMatrix* residual::getdPlasticDeformationGradientdF( ){
            /*!
             * Get the derivative of the plastic deformation gradient w.r.t. the deformation gradient
             */

            if ( !_dPlasticDeformationGradientdF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticDeformationGradientdF( ) );

            }

            return &_dPlasticDeformationGradientdF.second;

        }

        const floatMatrix* residual::getdPlasticDeformationGradientdSubFs( ){
            /*!
             * Get the derivative of the plastic deformation gradient w.r.t. the sub-deformation gradients
             */

            if ( !_dPlasticDeformationGradientdSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticDeformationGradientdSubFs( ) );

            }

            return &_dPlasticDeformationGradientdSubFs.second;

        }

        const floatVector* residual::getdPlasticDeformationGradientdT( ){
            /*!
             * Get the derivative of the plastic deformation gradient w.r.t. the sub-deformation gradients
             */

            if ( !_dPlasticDeformationGradientdT.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticDeformationGradientdT( ) );

            }

            return &_dPlasticDeformationGradientdT.second;

        }

        const floatMatrix* residual::getdPlasticDeformationGradientdStateVariables( ){
            /*!
             * Get the derivative of the plastic deformation gradient w.r.t. the state variables
             */

            if ( !_dPlasticDeformationGradientdStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticDeformationGradientdStateVariables( ) );

            }

            return &_dPlasticDeformationGradientdStateVariables.second;

        }

        const floatMatrix* residual::getdPlasticDeformationGradientdPreviousCauchyStress( ){
            /*!
             * Get the derivative of the plastic deformation gradient w.r.t. the previous Cauchy stress
             */

            if ( !_dPlasticDeformationGradientdPreviousCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticDeformationGradientdPreviousCauchyStress( ) );

            }

            return &_dPlasticDeformationGradientdPreviousCauchyStress.second;

        }

        const floatMatrix* residual::getdPlasticDeformationGradientdPreviousF( ){
            /*!
             * Get the derivative of the plastic deformation gradient w.r.t. the previous deformation gradient
             */

            if ( !_dPlasticDeformationGradientdPreviousF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticDeformationGradientdPreviousF( ) );

            }

            return &_dPlasticDeformationGradientdPreviousF.second;

        }

        const floatMatrix* residual::getdPlasticDeformationGradientdPreviousSubFs( ){
            /*!
             * Get the derivative of the plastic deformation gradient w.r.t. the previous sub-deformation gradients
             */

            if ( !_dPlasticDeformationGradientdPreviousSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticDeformationGradientdPreviousSubFs( ) );

            }

            return &_dPlasticDeformationGradientdPreviousSubFs.second;

        }

        const floatVector* residual::getdPlasticDeformationGradientdPreviousT( ){
            /*!
             * Get the derivative of the plastic deformation gradient w.r.t. the previous sub-deformation gradients
             */

            if ( !_dPlasticDeformationGradientdPreviousT.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticDeformationGradientdPreviousT( ) );

            }

            return &_dPlasticDeformationGradientdPreviousT.second;

        }

        const floatMatrix* residual::getdPlasticDeformationGradientdPreviousStateVariables( ){
            /*!
             * Get the derivative of the plastic deformation gradient w.r.t. the previous state variables
             */

            if ( !_dPlasticDeformationGradientdPreviousStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticDeformationGradientdPreviousStateVariables( ) );

            }

            return &_dPlasticDeformationGradientdPreviousStateVariables.second;

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

        const floatMatrix* residual::getdPlasticStateVariablesdCauchyStress( ){
            /*!
             * Get the derivative of the plastic state variables w.r.t. the Cauchy stress
             */

            if ( !_dPlasticStateVariablesdCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticStateVariablesdCauchyStress( ) );

            }

            return &_dPlasticStateVariablesdCauchyStress.second;

        }

        const floatMatrix* residual::getdPlasticStateVariablesdF( ){
            /*!
             * Get the derivative of the plastic state variables w.r.t. the deformation gradient
             */

            if ( !_dPlasticStateVariablesdF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticStateVariablesdF( ) );

            }

            return &_dPlasticStateVariablesdF.second;

        }

        const floatMatrix* residual::getdPlasticStateVariablesdSubFs( ){
            /*!
             * Get the derivative of the plastic state variables w.r.t. the sub-deformation gradients
             */

            if ( !_dPlasticStateVariablesdSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticStateVariablesdSubFs( ) );

            }

            return &_dPlasticStateVariablesdSubFs.second;

        }

        const floatVector* residual::getdPlasticStateVariablesdT( ){
            /*!
             * Get the derivative of the plastic state variables w.r.t. the temperature
             */

            if ( !_dPlasticStateVariablesdT.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticStateVariablesdT( ) );

            }

            return &_dPlasticStateVariablesdT.second;

        }

        const floatMatrix* residual::getdPlasticStateVariablesdStateVariables( ){
            /*!
             * Get the derivative of the plastic state variables w.r.t. the state variables
             */

            if ( !_dPlasticStateVariablesdStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticStateVariablesdStateVariables( ) );

            }

            return &_dPlasticStateVariablesdStateVariables.second;

        }

        const floatMatrix* residual::getdPlasticStateVariablesdPreviousCauchyStress( ){
            /*!
             * Get the derivative of the plastic state variables w.r.t. the previous Cauchy stress
             */

            if ( !_dPlasticStateVariablesdPreviousCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticStateVariablesdPreviousCauchyStress( ) );

            }

            return &_dPlasticStateVariablesdPreviousCauchyStress.second;

        }

        const floatMatrix* residual::getdPlasticStateVariablesdPreviousF( ){
            /*!
             * Get the derivative of the plastic state variables w.r.t. the previous deformation gradient
             */

            if ( !_dPlasticStateVariablesdPreviousF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticStateVariablesdPreviousF( ) );

            }

            return &_dPlasticStateVariablesdPreviousF.second;

        }

        const floatMatrix* residual::getdPlasticStateVariablesdPreviousSubFs( ){
            /*!
             * Get the derivative of the plastic state variables w.r.t. the previous sub-deformation gradients
             */

            if ( !_dPlasticStateVariablesdPreviousSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticStateVariablesdPreviousSubFs( ) );

            }

            return &_dPlasticStateVariablesdPreviousSubFs.second;

        }

        const floatVector* residual::getdPlasticStateVariablesdPreviousT( ){
            /*!
             * Get the derivative of the plastic state variables w.r.t. the previous temperature
             */

            if ( !_dPlasticStateVariablesdPreviousT.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticStateVariablesdPreviousT( ) );

            }

            return &_dPlasticStateVariablesdPreviousT.second;

        }

        const floatMatrix* residual::getdPlasticStateVariablesdPreviousStateVariables( ){
            /*!
             * Get the derivative of the plastic state variables w.r.t. the previous state variables
             */

            if ( !_dPlasticStateVariablesdPreviousStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPlasticStateVariablesdPreviousStateVariables( ) );

            }

            return &_dPlasticStateVariablesdPreviousStateVariables.second;

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

        const floatVector* residual::getdPreviousPlasticMultiplierdPreviousCauchyStress( ){
            /*!
             * Get the derivative of the previous plastic multiplier w.r.t. the previous Cauchy stress
             */

            if ( !_dPreviousPlasticMultiplierdPreviousCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousPlasticMultiplierdPreviousCauchyStress( ) );

            }

            return &_dPreviousPlasticMultiplierdPreviousCauchyStress.second;

        }

        const floatVector* residual::getdPreviousPlasticMultiplierdPreviousF( ){
            /*!
             * Get the derivative of the previous plastic multiplier w.r.t. the previous deformation gradient
             */

            if ( !_dPreviousPlasticMultiplierdPreviousF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousPlasticMultiplierdPreviousF( ) );

            }

            return &_dPreviousPlasticMultiplierdPreviousF.second;

        }

        const floatVector* residual::getdPreviousPlasticMultiplierdPreviousSubFs( ){
            /*!
             * Get the derivative of the previous plastic multiplier w.r.t. the previous sub-deformation gradients
             */

            if ( !_dPreviousPlasticMultiplierdPreviousSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousPlasticMultiplierdPreviousSubFs( ) );

            }

            return &_dPreviousPlasticMultiplierdPreviousSubFs.second;

        }

        const floatType* residual::getdPreviousPlasticMultiplierdPreviousT( ){
            /*!
             * Get the derivative of the previous plastic multiplier w.r.t. the previous temperature
             */

            if ( !_dPreviousPlasticMultiplierdPreviousT.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousPlasticMultiplierdPreviousT( ) );

            }

            return &_dPreviousPlasticMultiplierdPreviousT.second;

        }

        const floatVector* residual::getdPreviousPlasticMultiplierdPreviousStateVariables( ){
            /*!
             * Get the derivative of the previous plastic multiplier w.r.t. the previous state variables
             */

            if ( !_dPreviousPlasticMultiplierdPreviousStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousPlasticMultiplierdPreviousStateVariables( ) );

            }

            return &_dPreviousPlasticMultiplierdPreviousStateVariables.second;

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

        const floatMatrix* residual::getdPreviousVelocityGradientdPreviousCauchyStress( ){
            /*!
             * Get the derivative of the previous velocity gradient w.r.t.
             * the previous Cauchy stress.
             */

            if ( !_dPreviousVelocityGradientdPreviousCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousVelocityGradientdPreviousCauchyStress( ) );

            }

            return &_dPreviousVelocityGradientdPreviousCauchyStress.second;

        }

        const floatMatrix* residual::getdPreviousVelocityGradientdPreviousF( ){
            /*!
             * Get the derivative of the previous velocity gradient w.r.t.
             * the previous deformation gradient
             */

            if ( !_dPreviousVelocityGradientdPreviousF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousVelocityGradientdPreviousF( ) );

            }

            return &_dPreviousVelocityGradientdPreviousF.second;

        }

        const floatMatrix* residual::getdPreviousVelocityGradientdPreviousSubFs( ){
            /*!
             * Get the derivative of the previous velocity gradient w.r.t.
             * the previous sub-deformation gradients
             */

            if ( !_dPreviousVelocityGradientdPreviousSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousVelocityGradientdPreviousSubFs( ) );

            }

            return &_dPreviousVelocityGradientdPreviousSubFs.second;

        }

        const floatVector* residual::getdPreviousVelocityGradientdPreviousT( ){
            /*!
             * Get the derivative of the previous velocity gradient w.r.t.
             * the previous temperature
             */

            if ( !_dPreviousVelocityGradientdPreviousT.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousVelocityGradientdPreviousT( ) );

            }

            return &_dPreviousVelocityGradientdPreviousT.second;

        }

        const floatMatrix* residual::getdPreviousVelocityGradientdPreviousStateVariables( ){
            /*!
             * Get the derivative of the previous velocity gradient w.r.t.
             * the previous state variables
             */

            if ( !_dPreviousVelocityGradientdPreviousStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousVelocityGradientdPreviousStateVariables( ) );

            }

            return &_dPreviousVelocityGradientdPreviousStateVariables.second;

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

        const floatMatrix* residual::getdPreviousStateVariableEvolutionRatesdPreviousCauchyStress( ){
            /*!
             * Get the derivative of the previous state variable evolution rate w.r.t. the previous Cauchy stress
             */

            if ( !_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousStateVariableEvolutionRatesdPreviousCauchyStress( ) );

            }

            return &_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress.second;

        }

        const floatMatrix* residual::getdPreviousStateVariableEvolutionRatesdPreviousF( ){
            /*!
             * Get the derivative of the previous state variable evolution rate w.r.t. the previous deformation gradient
             */

            if ( !_dPreviousStateVariableEvolutionRatesdPreviousF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousStateVariableEvolutionRatesdPreviousF( ) );

            }

            return &_dPreviousStateVariableEvolutionRatesdPreviousF.second;

        }

        const floatMatrix* residual::getdPreviousStateVariableEvolutionRatesdPreviousSubFs( ){
            /*!
             * Get the derivative of the previous state variable evolution rate w.r.t. the previous sub-deformation gradients
             */

            if ( !_dPreviousStateVariableEvolutionRatesdPreviousSubFs.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousStateVariableEvolutionRatesdPreviousSubFs( ) );

            }

            return &_dPreviousStateVariableEvolutionRatesdPreviousSubFs.second;

        }

        const floatVector* residual::getdPreviousStateVariableEvolutionRatesdPreviousT( ){
            /*!
             * Get the derivative of the previous state variable evolution rate w.r.t. the previous sub-deformation gradients
             */

            if ( !_dPreviousStateVariableEvolutionRatesdPreviousT.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousStateVariableEvolutionRatesdPreviousT( ) );

            }

            return &_dPreviousStateVariableEvolutionRatesdPreviousT.second;

        }

        const floatMatrix* residual::getdPreviousStateVariableEvolutionRatesdPreviousStateVariables( ){
            /*!
             * Get the derivative of the previous state variable evolution rate w.r.t. the previous state variables
             */

            if ( !_dPreviousStateVariableEvolutionRatesdPreviousStateVariables.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPreviousStateVariableEvolutionRatesdPreviousStateVariables( ) );

            }

            return &_dPreviousStateVariableEvolutionRatesdPreviousStateVariables.second;

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
