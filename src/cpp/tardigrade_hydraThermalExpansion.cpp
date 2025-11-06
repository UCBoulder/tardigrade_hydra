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

            auto thermalGreenLagrangeStrain = get_setDataStorage_thermalGreenLagrangeStrain( );

            auto dThermalGreenLagrangeStraindT = get_setDataStorage_dThermalGreenLagrangeStraindT( );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::quadraticThermalExpansion( *hydra->getTemperature( ), *getReferenceTemperature( ), *getLinearParameters( ), *getQuadraticParameters( ), *thermalGreenLagrangeStrain.value, *dThermalGreenLagrangeStraindT.value ) );

        }

        void residual::setThermalDeformationGradient( ){
            /*!
             * Set the thermal deformation gradient
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            constexpr unsigned int tot_dim = sot_dim * dim;

            constexpr unsigned int fot_dim = tot_dim * dim;

            auto thermalDeformationGradient = get_setDataStorage_thermalDeformationGradient( );

            auto dThermalDeformationGradientdT = get_setDataStorage_dThermalDeformationGradientdT( );

            floatVector dThermalGreenLagrangeStraindThermalDeformationGradient;

            floatMatrix _dThermalGreenLagrangeStraindThermalDeformationGradient;

            floatVector A = 2 * ( *get_thermalGreenLagrangeStrain( ) );
            for ( unsigned int i = 0; i < dim; i++ ){ A[ dim * i + i ] += 1; }

            TARDIGRADE_ERROR_TOOLS_CATCH( *thermalDeformationGradient.value = tardigradeVectorTools::matrixSqrt( A, dim, _dThermalGreenLagrangeStraindThermalDeformationGradient ) );

            dThermalGreenLagrangeStraindThermalDeformationGradient = tardigradeVectorTools::appendVectors( _dThermalGreenLagrangeStraindThermalDeformationGradient );

            floatVector dThermalDeformationGradientdGreenLagrangeStrain( fot_dim, 0 );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dThermalGreenLagrangeStraindThermalDeformationGradient( dThermalGreenLagrangeStraindThermalDeformationGradient.data( ), sot_dim, sot_dim );

            Eigen::Map< Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dThermalDeformationGradientdGreenLagrangeStrain( dThermalDeformationGradientdGreenLagrangeStrain.data( ), sot_dim, sot_dim );

            map_dThermalDeformationGradientdGreenLagrangeStrain = 2 * ( map_dThermalGreenLagrangeStraindThermalDeformationGradient.inverse( ) ).eval( );

            Eigen::Map< const Eigen::Vector< floatType, sot_dim > > map_dThermalGreenLagrangeStraindT( get_dThermalGreenLagrangeStraindT( )->data( ), sot_dim );

            dThermalDeformationGradientdT.zero( sot_dim );
            Eigen::Map< Eigen::Vector< floatType, sot_dim > > map_dThermalDeformationGradientdT( dThermalDeformationGradientdT.value->data( ), sot_dim );

            map_dThermalDeformationGradientdT = ( map_dThermalDeformationGradientdGreenLagrangeStrain * map_dThermalGreenLagrangeStraindT ).eval( );

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

            const unsigned int thermalConfigurationIndex = *getThermalConfigurationIndex( );

            auto residual = get_setDataStorage_residual( );
            *residual.value = *get_thermalDeformationGradient( ) - floatVector( hydra->get_configurations( )->begin( ) +   thermalConfigurationIndex * 9,
                                                                                hydra->get_configurations( )->begin( ) + ( thermalConfigurationIndex + 1 ) * 9 );

        }

        void residual::setJacobian( ){
            /*!
             * Set the values of the jacobian
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int num_unknowns = hydra->getNumUnknowns( );

            auto jacobian = get_setDataStorage_jacobian( );
            jacobian.zero( *getNumEquations( ) * num_unknowns );

            for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                ( *jacobian.value )[ num_unknowns * i + sot_dim * ( *getThermalConfigurationIndex( ) ) + i ] = -1;

            }

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            auto dRdT = get_setDataStorage_dRdT( );
            *dRdT.value = *get_dThermalDeformationGradientdT( );

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient
             */

            auto dRdF = get_setDataStorage_dRdF( );
            dRdF.zero( *getNumEquations( ) * hydra->getDeformationGradient( )->size( ) );

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

            const unsigned int dim = hydra->getDimension( );

            const unsigned int sot_dim = hydra->getSOTDimension( );

            unsigned int parametersPerTerm = dim + ( dim * dim - dim ) / 2;

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

            floatVector linearParameters( sot_dim, 0 );

            floatVector quadraticParameters( sot_dim, 0 );

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = i; j < dim; j++ ){

                    linearParameters[ dim * i + j ] = linearParameters[ dim * j + i ] = parameters[ 1 + index ];

                    quadraticParameters[ dim * i + j ] = quadraticParameters[ dim * j + i ] = parameters[ 1 + parametersPerTerm + index ];

                    index++;

                }

            }

            setLinearParameters( linearParameters );

            setQuadraticParameters( quadraticParameters );

        }

        void residual::addParameterizationInfo( std::string &parameterization_info ){
            /*!
             * Add the parameterization information to the incoming string
             * 
             * \param &parameterization_info: The incoming string
             */

            std::stringstream ss;
            ss.precision(9);
            ss << std::scientific;

            ss << "class: tardigradeHydra::thermalExpansion::residual\n\n";
            ss << "name,                                    description,       units, current value\n";
            ss << "Tref,                      The reference temperature, temperature, " << *getReferenceTemperature( ) << "\n";
            ss << "l_11,    The 11 linear thermal expansion coefficient,        none, " << ( *getLinearParameters( ) )[ 0 ] << "\n";
            ss << "l_12,    The 12 linear thermal expansion coefficient,        none, " << ( *getLinearParameters( ) )[ 1 ] << "\n";
            ss << "l_13,    The 13 linear thermal expansion coefficient,        none, " << ( *getLinearParameters( ) )[ 2 ] << "\n";
            ss << "l_21,    The 21 linear thermal expansion coefficient,        none, " << ( *getLinearParameters( ) )[ 3 ] << "\n";
            ss << "l_22,    The 22 linear thermal expansion coefficient,        none, " << ( *getLinearParameters( ) )[ 4 ] << "\n";
            ss << "l_23,    The 23 linear thermal expansion coefficient,        none, " << ( *getLinearParameters( ) )[ 5 ] << "\n";
            ss << "l_31,    The 31 linear thermal expansion coefficient,        none, " << ( *getLinearParameters( ) )[ 6 ] << "\n";
            ss << "l_32,    The 32 linear thermal expansion coefficient,        none, " << ( *getLinearParameters( ) )[ 7 ] << "\n";
            ss << "l_33,    The 33 linear thermal expansion coefficient,        none, " << ( *getLinearParameters( ) )[ 8 ] << "\n";
            ss << "q_11, The 11 quadratic thermal expansion coefficient,        none, " << ( *getQuadraticParameters( ) )[ 0 ] << "\n";
            ss << "q_12, The 12 quadratic thermal expansion coefficient,        none, " << ( *getQuadraticParameters( ) )[ 1 ] << "\n";
            ss << "q_13, The 13 quadratic thermal expansion coefficient,        none, " << ( *getQuadraticParameters( ) )[ 2 ] << "\n";
            ss << "q_21, The 21 quadratic thermal expansion coefficient,        none, " << ( *getQuadraticParameters( ) )[ 3 ] << "\n";
            ss << "q_22, The 22 quadratic thermal expansion coefficient,        none, " << ( *getQuadraticParameters( ) )[ 4 ] << "\n";
            ss << "q_23, The 23 quadratic thermal expansion coefficient,        none, " << ( *getQuadraticParameters( ) )[ 5 ] << "\n";
            ss << "q_31, The 31 quadratic thermal expansion coefficient,        none, " << ( *getQuadraticParameters( ) )[ 6 ] << "\n";
            ss << "q_32, The 32 quadratic thermal expansion coefficient,        none, " << ( *getQuadraticParameters( ) )[ 7 ] << "\n";
            ss << "q_33, The 33 quadratic thermal expansion coefficient,        none, " << ( *getQuadraticParameters( ) )[ 8 ] << "\n";

            ss.unsetf(std::ios_base::floatfield);
            parameterization_info.append(ss.str());

        }

    }

}
