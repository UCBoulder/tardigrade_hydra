/**
  ******************************************************************************
  * \file tardigrade_hydraThermalExpansion.cpp
  ******************************************************************************
  * An implementation of thermal expansion using the hydra framework. Used as an
  * example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraThermalExpansion.h>
#include<tardigrade_constitutive_tools.h>

namespace tardigradeHydra{

    namespace thermalExpansion{

        void residual::setReferenceTemperature( const floatType &referenceTemperature ){
            /*!
             * Set the reference temperature\n";
             *
             * \param &referenceTemperature: The reference temperature
             */

            _referenceTemperature = referenceTemperature;

        }

        void residual::setLinearParameters( const floatVector &linearParameters ){
            /*!
             * Set the linear thermal expansion coefficients
             * 
             * \param &linearParameters: The expanded linear parameters in row-major format
             *     ex. for 3D [ c11, c12, c13, c21, c22, c23, c31, c32, c33 ]
             */

            _linearParameters = linearParameters;

        }

        void residual::setQuadraticParameters( const floatVector &quadraticParameters ){
            /*!
             * Set the quadratic thermal expansion coefficients
             * 
             * \param &quadraticParameters: The expanded quadratic parameters in row-major
             *     format ex. for 3D [ c11, c12, c13, c21, c22, c23, c31, c32, c33 ]
             */

            _quadraticParameters = quadraticParameters;

        }

        void residual::setThermalGreenLagrangeStrain( ){
            /*!
             * Set the thermal Green-Lagrange strain defined as
             * 
             * \f$ E^{\theta}_{IJ} = \frac{1}{2} \left( F_{iI}^{\theta} F_{iJ}^{\theta} - \delta_{IJ} \right) \f$
             */

            floatVector thermalGreenLagrangeStrain;

            floatVector dThermalGreenLagrangeStraindT;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::quadraticThermalExpansion( *hydra->getTemperature( ), *getReferenceTemperature( ), *getLinearParameters( ), *getQuadraticParameters( ), thermalGreenLagrangeStrain, dThermalGreenLagrangeStraindT ) );

            set_thermalGreenLagrangeStrain( thermalGreenLagrangeStrain );

            set_dThermalGreenLagrangeStraindT( dThermalGreenLagrangeStraindT );

        }

        void residual::setThermalDeformationGradient( ){
            /*!
             * Set the thermal deformation gradient
             */

            const unsigned int *dim = hydra->getDimension( );

            floatVector thermalDeformationGradient;

            floatMatrix dThermalGreenLagrangeStraindThermalDeformationGradient;

            floatVector eye( ( *hydra->getDimension( ) ) * ( *hydra->getDimension( ) ) );
            tardigradeVectorTools::eye( eye );

            TARDIGRADE_ERROR_TOOLS_CATCH( thermalDeformationGradient = tardigradeVectorTools::matrixSqrt( 2 * ( *get_thermalGreenLagrangeStrain( ) ) + eye, *hydra->getDimension( ), dThermalGreenLagrangeStraindThermalDeformationGradient ) );

            set_thermalDeformationGradient( thermalDeformationGradient );

            floatMatrix dThermalDeformationGradientdGreenLagrangeStrain;

            TARDIGRADE_ERROR_TOOLS_CATCH( dThermalDeformationGradientdGreenLagrangeStrain = tardigradeVectorTools::inflate( 2 * tardigradeVectorTools::inverse( tardigradeVectorTools::appendVectors( dThermalGreenLagrangeStraindThermalDeformationGradient ), ( *dim ) * ( *dim ), ( *dim ) * ( *dim ) ), ( *dim ) * ( *dim ), ( *dim ) * ( *dim ) ) );

            set_dThermalDeformationGradientdT( tardigradeVectorTools::dot( dThermalDeformationGradientdGreenLagrangeStrain, *get_dThermalGreenLagrangeStraindT( ) ) );

        }

        void residual::setdThermalGreenLagrangeStraindT( ){
            /*!
             * Set the derivative of the thermal Green-Lagrange strain w.r.t. the temperature
             */

            setThermalGreenLagrangeStrain( );

        }

        void residual::setdThermalDeformationGradientdT( ){
            /*!
             * Set the derivative of the thermal deformation gradient w.r.t. the temperature
             */

            setThermalDeformationGradient( );

        }

        void residual::setResidual( ){
            /*!
             * Set the value of the residual
             * 
             * Defined as the residual's computed thermal deformation gradient minus the value stored in hydra's configurations.
             */

            setResidual( *get_thermalDeformationGradient( ) - ( *hydra->get_configurations( ) )[ *getThermalConfigurationIndex( ) ] );

        }

        void residual::setJacobian( ){
            /*!
             * Set the values of the jacobian
             */

            const unsigned int *dim = hydra->getDimension( );

            floatMatrix jacobian( *getNumEquations( ), floatVector( hydra->getUnknownVector( )->size( ), 0 ) );

            for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                jacobian[ i ][ ( *dim ) * ( *dim ) * ( *getThermalConfigurationIndex( ) ) + i ] = -1;

            }

            setJacobian( jacobian );

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            setdRdT( *get_dThermalDeformationGradientdT( ) );

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient
             */

            setdRdF( floatMatrix( *getNumEquations( ), floatVector( hydra->getDeformationGradient( )->size( ), 0 ) ) );

        }

        void residual::decomposeParameters( const floatVector &parameters ){
            /*!
             * Decompose the parameter vector
             *
             * \param &parameters: The parameter vector. Assumed to be
             *     of the form ( referenceTemperature, linearParameters, quadraticParameters )
             *     where the parameters are defined as in
             *     tardigradeConstitutiveTools::quadraticThermalExpansion.
             */

            const unsigned int *dim = hydra->getDimension( );

            unsigned int parametersPerTerm = ( *dim ) + ( ( *dim ) * ( *dim ) - ( *dim ) ) / 2;

            unsigned int expectedSize = 1 + 2 * parametersPerTerm;

            if ( parameters.size( ) != expectedSize ){

                std::string message = "The parameters vector is not the expected length.\n";
                message            += "  expected:   " + std::to_string( expectedSize ) + "\n";
                message            += "  parameters: " + std::to_string( parameters.size( ) ) + "\n";
                message            += "Parameters should be organized as ( referenceTemperature, lP1, lp2..., qp1, qp2, ... )\n";
                message            += "where the number of parameters will be related to the dimension.\n";
                message            += "The thermal strain must be symmetric so only the upper triangular parameters are required.\n";
                message            += "e.x. for 3D 13 parameters should be given.\n";

            }

            setReferenceTemperature( parameters[ 0 ] );

            unsigned int index = 0;

            floatVector linearParameters( ( *dim ) * ( *dim ), 0 );

            floatVector quadraticParameters( ( *dim ) * ( *dim ), 0 );

            for ( unsigned int i = 0; i < ( *dim ); i++ ){

                for ( unsigned int j = i; j < ( *dim ); j++ ){

                    linearParameters[ ( *dim ) * i + j ] = linearParameters[ ( *dim ) * j + i ] = parameters[ 1 + index ];

                    quadraticParameters[ ( *dim ) * i + j ] = quadraticParameters[ ( *dim ) * j + i ] = parameters[ 1 + parametersPerTerm + index ];

                    index++;

                }

            }

            setLinearParameters( linearParameters );

            setQuadraticParameters( quadraticParameters );

        }

    }

}
