/**
  ******************************************************************************
  * \file tardigrade_hydraFourierHeatConduction.h
  ******************************************************************************
  * An implementation of Fourier heat conduction using the hydra framework.
  ******************************************************************************
  */

#include<tardigrade_hydraFourierHeatConduction.h>

namespace tardigradeHydra{

    namespace fourierHeatConduction{

        void residual::decomposeParameterVector( const floatVector &parameters ){
            /*!
             * Decompose the parameter vector
             * 
             * \param &parameters: The paramter vector. Assumed to be a vector of length _expected_parameter_size which defines the conductivity.
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(
                parameters.size( ) == get_expectedParameterVectorSize( ),
                "Parameter vector is expected to have a length of " + std::to_string( get_expectedParameterVectorSize( ) ) + " but has a length of " + std::to_string( parameters.size( ) )
            )

            setConductivityParameters( parameters );

        }

        void residual::setTemperatureGradient( ){
            /*!
             * Set the value of the temperature gradient. It's assumed this is stored in the additional DOF vector
             */

            auto temperatureGradientIndex = get_temperatureGradientIndex( );

            auto dim = hydra->getDimension( );

            auto temperature_gradient = get_setDataStorage_temperatureGradient( );

            temperature_gradient.zero( dim );

            std::copy(
                std::begin( *hydra->getAdditionalDOF( ) ) + temperatureGradientIndex,
                std::begin( *hydra->getAdditionalDOF( ) ) + temperatureGradientIndex + dim,
                temperature_gradient.begin( )
            );

        }

        void residual::setPreviousTemperatureGradient( ){
            /*!
             * Set the previous value of the temperature gradient. It's assumed this is stored in the additional DOF vector
             */

            auto temperatureGradientIndex = get_temperatureGradientIndex( );

            auto dim = hydra->getDimension( );

            auto temperature_gradient = get_setDataStorage_previousTemperatureGradient( );

            temperature_gradient.zero( dim );

            std::copy(
                std::begin( *hydra->getPreviousAdditionalDOF( ) ) + temperatureGradientIndex,
                std::begin( *hydra->getPreviousAdditionalDOF( ) ) + temperatureGradientIndex + dim,
                temperature_gradient.begin( )
            );

        }

        void residual::setConductivity( ){
            /*!
             * Set the conductivity
             */

            auto conductivity = get_setDataStorage_conductivity( );

            *conductivity.value = ( *get_conductivityParameters( ) )[ 0 ];

        }

        void residual::setPreviousConductivity( ){
            /*!
             * Set the previous conductivity
             */

            auto conductivity = get_setDataStorage_previousConductivity( );

            *conductivity.value = ( *get_conductivityParameters( ) )[ 0 ];

        }

        void residual::setHeatFlux( ){
            /*!
             * Set the heat flux
             */

            auto dim = hydra->getDimension( );

            auto heat_flux = get_setDataStorage_heatFlux( );

            heat_flux.zero( dim );

            std::transform(
                heat_flux.begin( ),
                heat_flux.end( ),
                std::begin( *get_temperatureGradient( ) ),
                heat_flux.begin( ),
                std::bind( std::multiplies<>( ), std::placeholders::_1, -( *get_conductivity( ) ) )
            );

        }

        void residual::setPreviousHeatFlux( ){
            /*!
             * Set the previous heat flux
             */

            auto dim = hydra->getDimension( );

            auto heat_flux = get_setDataStorage_previousHeatFlux( );

            heat_flux.zero( dim );

            std::transform(
                heat_flux.begin( ),
                heat_flux.end( ),
                std::begin( *get_previousTemperatureGradient( ) ),
                heat_flux.begin( ),
                std::bind( std::multiplies<>( ), std::placeholders::_1, -( *get_previousConductivity( ) ) )
            );

        }

        void residual::setdHeatFluxdGradT( ){
            /*!
             * Get the derivative of the heat flux w.r.t. the gradient of the temperature
             */
        }

        void residual::setResidual( ){
            /*!
             * Set the residual value
             */

            auto residual = get_setDataStorage_residual( );

            residual.zero( *getNumEquations( ) );

        }
    
        void residual::setJacobian( ){
            /*!
             * Set the Jacobian value
             */

            auto dim = hydra->getDimension( );

            auto sot_dim = dim * dim;

            const unsigned int num_unknowns = hydra->getNumUnknowns( );

            // Form the Jacobian
            auto jacobian = get_setDataStorage_jacobian( );

            jacobian.zero( sot_dim * num_unknowns );

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            auto dRdT = get_setDataStorage_dRdT( );

            dRdT.zero( *getNumEquations( ) );

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient
             */

            auto dRdF = get_setDataStorage_dRdF( );

            dRdF.zero( ( *getNumEquations( ) ) * hydra->getDimension( ) * hydra->getDimension( ) );

        }

    }

}
