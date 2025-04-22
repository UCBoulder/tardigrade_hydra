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

            auto specificHeat = get_setDataStorage_specificHeat( );

            *specificHeat.value = ( *get_specificHeatParameters( ) )[ 0 ];

        }

        void residual::setPreviousSpecificHeat( ){
            /*!
             * Set the previous specific heat
             */

            auto specificHeat = get_setDataStorage_previousSpecificHeat( );

            *specificHeat.value = ( *get_specificHeatParameters( ) )[ 0 ];

        }

        void residual::setInternalEnergy( ){
            /*!
             * Set the internal energy
             */

            auto internal_energy = get_setDataStorage_internalEnergy( );

            ( *internal_energy.value ) = ( *get_specificHeat( ) ) * ( *hydra->getTemperature( ) );

        }

        void residual::setPreviousInternalEnergy( ){
            /*!
             * Set the previous internal energy
             */

            auto internal_energy = get_setDataStorage_previousInternalEnergy( );

            ( *internal_energy.value ) = ( *get_specificHeat( ) ) * ( *hydra->getPreviousTemperature( ) );

        }

        void residual::setdInternalEnergydT( ){
            /*!
             * Get the derivative of the internal energy w.r.t. the temperature
             */

            auto dInternalEnergydT = get_setDataStorage_dInternalEnergydT( );

            ( *dInternalEnergydT.value ) = *get_specificHeat( );

        }

        void residual::setResidual( ){
            /*!
             * Set the residual value
             */

            auto residual = get_setDataStorage_residual( );

            residual.zero( *getNumEquations( ) );

            ( *residual.value )[ 0 ] = ( *hydra->getUnknownVector( ) )[ get_internalEnergyIndex( ) ] - ( *get_internalEnergy( ) );

        }
    
        void residual::setJacobian( ){
            /*!
             * Set the Jacobian value
             */

            auto num_unknowns = hydra->getNumUnknowns( );

            // Form the Jacobian
            auto jacobian = get_setDataStorage_jacobian( );

            jacobian.zero( num_unknowns * num_unknowns );

            ( *jacobian.value )[ 0 + get_internalEnergyIndex( ) ] = 1.;

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            auto dRdT = get_setDataStorage_dRdT( );

            dRdT.zero( *getNumEquations( ) );

            ( *dRdT.value )[ 0 ] -= ( *get_dInternalEnergydT( ) );

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

            auto dRdAdditionalDOF = get_setDataStorage_dRdAdditionalDOF( );

            dRdAdditionalDOF.zero( ( *getNumEquations( ) ) * hydra->getNumAdditionalDOF( ) );

        }

    }

}
