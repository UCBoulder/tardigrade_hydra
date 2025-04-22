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
                std::begin( *get_temperatureGradient( ) ),
                std::end(   *get_temperatureGradient( ) ),
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
                std::begin( *get_previousTemperatureGradient( ) ),
                std::end(   *get_previousTemperatureGradient( ) ),
                heat_flux.begin( ),
                std::bind( std::multiplies<>( ), std::placeholders::_1, -( *get_previousConductivity( ) ) )
            );

        }

        void residual::setdHeatFluxdGradT( ){
            /*!
             * Get the derivative of the heat flux w.r.t. the gradient of the temperature
             */

            auto dim = hydra->getDimension( );

            auto sot_dim = dim * dim;

            auto dHeatFluxdGradT = get_setDataStorage_dHeatFluxdGradT( );

            dHeatFluxdGradT.zero( sot_dim );

            for ( unsigned int i = 0; i < dim; ++i ){
                ( *dHeatFluxdGradT.value )[ dim * i + i ] = -( *get_conductivity( ) );
            }

        }

        void residual::setResidual( ){
            /*!
             * Set the residual value
             */

            auto dim = hydra->getDimension( );

            auto residual = get_setDataStorage_residual( );

            residual.zero( *getNumEquations( ) );

            for ( unsigned int i = 0; i < dim; ++i ){

                ( *residual.value )[ i ] = ( *hydra->getUnknownVector( ) )[ get_heatFluxIndex( ) + i ] - ( *get_heatFlux( ) )[ i ];

            }

        }
    
        void residual::setJacobian( ){
            /*!
             * Set the Jacobian value
             */

            auto dim = hydra->getDimension( );

            auto num_unknowns = hydra->getNumUnknowns( );

            // Form the Jacobian
            auto jacobian = get_setDataStorage_jacobian( );

            jacobian.zero( num_unknowns * num_unknowns );

            for ( unsigned int i = 0; i < dim; ++i ){

                ( *jacobian.value )[ hydra->getNumUnknowns( ) * ( i ) + get_heatFluxIndex( ) + i ] = 1.;

            }

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

        void residual::setdRdAdditionalDOF( ){
            /*!
             * Set the derivative of the residual w.r.t. the additional DOF vector
             */

            auto dim = hydra->getDimension( );

            auto dRdAdditionalDOF = get_setDataStorage_dRdAdditionalDOF( );

            dRdAdditionalDOF.zero( ( *getNumEquations( ) ) * hydra->getNumAdditionalDOF( ) );

            for ( unsigned int i = 0; i < dim; ++i ){

                for ( unsigned int j = 0; j < dim; ++j ){

                    ( *dRdAdditionalDOF.value )[ hydra->getNumAdditionalDOF( ) * ( i ) + get_temperatureGradientIndex( ) + j ] -= ( *get_dHeatFluxdGradT( ) )[ dim * i + j ];

                }

            }

        }

    }

}
