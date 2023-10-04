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

namespace tardigradeHydra{

    namespace peryznaViscoplasticity{

        void residual::setDrivingStress( ){
            /*!
             * Set the driving stress i.e. the Cauchy stress pulled back to the
             * current configuration of the plastic configuration.
             */

            floatVector precedingConfiguration;

            TARDIGRADE_ERROR_TOOLS_CATCH( precedingConfiguration = hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

            floatVector drivingStress;

            tardigradeConstitutiveTools::pullBackCauchyStress( *hydra->getCauchyStress( ), precedingConfiguration, drivingStress );

            setDrivingStress( drivingStress );

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

            _yieldParameters.first = true;

        }

        void residual::setMixingParameters( const floatVector &mixingParameters ){
            /*!
             * Set the mixing parameters
             * 
             * \param &mixingParameters: The mixing parameters
             */

            _mixingParameters.second = mixingParameters;

            _mixingParameters.first = true;

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

        const floatVector* residual::getFlowDirection( ){
            /*!
             * Get the flow direction
             */

            if ( !_flowDirection.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setFlowDirection( ) );

            }

            return &_flowDirection.second;

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

        const floatVector* residual::getPlasticDeformationGradient( ){
            /*!
             * Get the plastic deformation gradient
             */

            if ( !_plasticDeformationGradient.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setPlasticDeformationGradient( ) );

            }

            return &_plasticDeformationGradient.second;

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

        const floatVector* residual::getMixingParameters( ){
            /*!
             * Get the mixing parameters
             */

            if ( !_mixingParameters.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "Mixing parameters not defined but required" ) );

            }

            return &_mixingParameters.second;

        }

        void residual::decomposeParameters( const floatVector &parameters ){
            /*!
             * Decompose the incoming parameter vector
             * 
             * \param &parameters: The incoming parameter vector. We assume a
             *     Peryzna viscoplastic driving stress with a Drucker-Prager
             *     yield and flow potential surface. The parameters are
             *     [ n, q0, q1, C1, C2, Tref, Y, A, B, H ] where n is the
             *     Peryzna exponent, q0 is the initial drag stress, q1 is the
             *     drag modulus, C1, C2, and Tref are the temperature effect
             *     parameters, Y is the yield stress for the Drucker-Prager equation,
             *     A is the pressure term of the yield equation, B is the flow parameter,
             *     and H is the mixing modulus.
             */

            unsigned int expectedSize = 10;

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

            setMixingParameters( { parameters[ 9 ] } );

        }

    }

}
