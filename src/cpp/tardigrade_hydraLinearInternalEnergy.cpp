/**
  ******************************************************************************
  * \file tardigrade_hydraLinearInternalEnergy.h
  ******************************************************************************
  * An implementation of Fourier heat conduction using the hydra framework.
  ******************************************************************************
  */

#include<tardigrade_hydraLinearInternalEnergy.h>

namespace tardigradeHydra{

    namespace linearInternalEnergy{

        void residual::decomposeParameterVector( const floatVector &parameters ){
            /*!
             * Decompose the parameter vector
             * 
             * \param &parameters: The paramter vector. Assumed to be a vector of length _expected_parameter_size which defines the specific heat.
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(
                parameters.size( ) == get_expectedParameterVectorSize( ),
                "Parameter vector is expected to have a length of " + std::to_string( get_expectedParameterVectorSize( ) ) + " but has a length of " + std::to_string( parameters.size( ) )
            )

            setSpecificHeatParameters( parameters );

        }

        void residual::setSpecificHeat( ){
            /*!
             * Set the specific heat
             */

            auto specificHeat = get_SetDataStorage_specificHeat( );

            *specificHeat.value = ( *get_specificHeatParameters( ) )[ 0 ];

        }

        void residual::setPreviousSpecificHeat( ){
            /*!
             * Set the previous specific heat
             */

            auto specificHeat = get_SetDataStorage_previousSpecificHeat( );

            *specificHeat.value = ( *get_specificHeatParameters( ) )[ 0 ];

        }

        void residual::setInternalEnergy( ){
            /*!
             * Set the internal energy
             */

            auto internal_energy = get_SetDataStorage_internalEnergy( );

            ( *internal_energy.value ) = ( *get_specificHeat( ) ) * hydra->getTemperature( );

        }

        void residual::setPreviousInternalEnergy( ){
            /*!
             * Set the previous internal energy
             */

            auto internal_energy = get_SetDataStorage_previousInternalEnergy( );

            ( *internal_energy.value ) = ( *get_specificHeat( ) ) * hydra->getPreviousTemperature( );

        }

        void residual::setdInternalEnergydT( ){
            /*!
             * Get the derivative of the internal energy w.r.t. the temperature
             */

            auto dInternalEnergydT = get_SetDataStorage_dInternalEnergydT( );

            ( *dInternalEnergydT.value ) = *get_specificHeat( );

        }

        void residual::setResidual( ){
            /*!
             * Set the residual value
             */

            auto residual = get_SetDataStorage_residual( );

            residual.zero( getNumEquations( ) );

            ( *residual.value )[ 0 ] = ( *hydra->getUnknownVector( ) )[ get_internalEnergyIndex( ) ] - ( *get_internalEnergy( ) );

        }
    
        void residual::setJacobian( ){
            /*!
             * Set the Jacobian value
             */

            auto num_unknowns = hydra->getNumUnknowns( );

            // Form the Jacobian
            auto jacobian = get_SetDataStorage_jacobian( );

            jacobian.zero( getNumEquations( ) * num_unknowns );

            ( *jacobian.value )[ 0 + get_internalEnergyIndex( ) ] = 1.;

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            auto dRdT = get_SetDataStorage_dRdT( );

            dRdT.zero( getNumEquations( ) );

            ( *dRdT.value )[ 0 ] -= ( *get_dInternalEnergydT( ) );

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient
             */

            auto dRdF = get_SetDataStorage_dRdF( );

            dRdF.zero( getNumEquations( ) * hydra->getDimension( ) * hydra->getDimension( ) );

        }

        void residual::setdRdAdditionalDOF( ){
            /*!
             * Set the derivative of the residual w.r.t. the additional DOF vector
             */

            auto dRdAdditionalDOF = get_SetDataStorage_dRdAdditionalDOF( );

            dRdAdditionalDOF.zero( getNumEquations( ) * hydra->getNumAdditionalDOF( ) );

        }

    }

}
