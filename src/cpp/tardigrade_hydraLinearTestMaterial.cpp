/**
  ******************************************************************************
  * \file tardigrade_hydraLinearTestMaterial.cpp
  ******************************************************************************
  * An implementation of a linear test material using the hydra framework. Used
  * as an example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraLinearTestMaterial.h>
#include<tardigrade_constitutive_tools.h>

namespace tardigradeHydra{

    namespace linearTestMaterial{

        void residual::decomposeParameterVector( const floatVector &parameters ){
            /*!
             * Decompose the parameter vector
             * 
             * \param &parameters: The paramter vector. Assumed to be a vector of length 2 which defines lambda and mu.
             */

            TARDIGRADE_ERROR_TOOLS_EVAL(
                const unsigned int expected_parameter_size = ( *getNumEquations( ) ) * ( hydra->getNumAdditionalDOF( ) + 9 + 1 )
            )

            TARDIGRADE_ERROR_TOOLS_CHECK(
                parameters.size( ) == expected_parameter_size,
                "Parameter vector is expected to have a length of " + std::to_string( expected_parameter_size ) + " but has a length of " + std::to_string( parameters.size( ) )
            );

            // Decompose the parameter vector
            set_F_params( floatVector( std::begin( parameters ), std::begin( parameters ) + ( *getNumEquations( ) ) * ( 9 ) ) );
    
            set_T_params( floatVector( std::begin( parameters ) + ( *getNumEquations( ) ) * ( 9 ), std::begin( parameters ) + ( *getNumEquations( ) ) * ( 9 + 1 ) ) );
    
            set_add_dof_params( floatVector( std::begin( parameters ) + ( *getNumEquations( ) ) * ( 9 + 1 ), std::end( parameters ) ) );

        }

        void residual::setXPred( ){
            /*!
             * Compute all of the values of the X vector for the residual calculation
             */

            constexpr unsigned int dim = 3;

            auto F_params = get_F_params( );

            auto T_params = get_T_params( );

            auto add_dof_params = get_add_dof_params( );

            auto XPred = get_setDataStorage_XPred( );

            XPred.zero( *getNumEquations( ) );

            unsigned int i = 0;

            for
            (
                auto xi = XPred.begin( );
                xi != XPred.end( );
                ++i, ++xi
            )
            {

                unsigned int j = 0;

                // Add the contributions from the deformation gradient
                for
                (
                    auto v = std::begin( *hydra->getDeformationGradient( ) );
                    v != std::end( *hydra->getDeformationGradient( ) );
                    ++j, ++v
                )
                {
                
                    *xi += ( *F_params )[ dim * dim * i + j ] * ( *v );

                }

                // Add the contributions from the temperature
                *xi += ( *T_params )[ i ] * ( *( hydra->getTemperature( ) ) );

                // Add the contributions from the additional degrees of freedom
                j = 0;
                for
                (
                    auto v = std::begin( *hydra->getAdditionalDOF( ) );
                    v != std::end( *hydra->getAdditionalDOF( ) );
                    ++j, ++v
                )
                {

                    *xi += ( *add_dof_params )[ hydra->getNumAdditionalDOF( ) * i + j ] * ( *v );

                }

            }

        }

        void residual::setStress( ){
            /*!
             * Set the stress
             *
             * Currently uses the Cauchy stress
             */

            constexpr unsigned int dim = 3;

            auto stress = get_setDataStorage_stress( );

            stress.zero( dim * dim );

            std::copy(
                std::begin( *get_XPred( ) ),
                std::begin( *get_XPred( ) ) + dim * dim,
                std::begin( *stress.value )
            );

        }
    
        void residual::setPreviousStress( ){
            /*!
             * Set the previous stress
             *
             * Currently uses the Cauchy stress
             */

            TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "Not implemented" ) );

        }
    
        void residual::setResidual( ){
            /*!
             * Set the residual value
             */

            auto residual = get_setDataStorage_residual( );

            residual.zero( *getNumEquations( ) );

            std::transform(
                std::begin( *get_XPred( ) ),
                std::end( *get_XPred( ) ),
                std::begin( *( hydra->getUnknownVector( ) ) ),
                residual.begin( ),
                std::minus< floatType >{}
            );
    
        }
    
        void residual::setJacobian( ){
            /*!
             * Set the Jacobian value
             */

            auto jacobian = get_setDataStorage_jacobian( );

            jacobian.zero( ( *getNumEquations( ) ) * hydra->getNumUnknowns( ) );

            for
            (
                unsigned int i = 0;
                i < std::min( ( *getNumEquations( ) ), hydra->getNumUnknowns( ) );
                ++i
            )
            {
                ( *jacobian.value )[ ( *getNumEquations( ) ) * i + i ] = 1;
            }

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            auto dRdT = get_setDataStorage_dRdT( );

            dRdT.zero( ( *getNumEquations( ) ) );

            std::copy(
                std::begin( *get_T_params( ) ),
                std::end( *get_T_params( ) ),
                dRdT.begin( )
            );

        }

        void residual::setdRdF( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient
             */

            constexpr unsigned int dim = 3;

            auto dRdF = get_setDataStorage_dRdF( );

            dRdF.zero( ( *getNumEquations( ) ) * dim * dim );

            std::copy(
                std::begin( *get_F_params( ) ),
                std::end( *get_F_params( ) ),
                dRdF.begin( )
            );

        }

        void residual::setdRdAdditionalDOF( ){
            /*!
             * Set the derivative of the residual w.r.t. the additional DOF
             */

            auto dRdAdditionalDOF = get_setDataStorage_dRdAdditionalDOF( );

            dRdAdditionalDOF.zero( ( *getNumEquations( ) ) * ( hydra->getNumAdditionalDOF( ) ) );

            std::copy(
                std::begin( *get_add_dof_params( ) ),
                std::end( *get_add_dof_params( ) ),
                dRdAdditionalDOF.begin( )
            );

        }

    }

}
