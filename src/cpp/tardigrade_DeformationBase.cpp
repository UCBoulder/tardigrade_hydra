/**
  ******************************************************************************
  * \file tardigrade_DeformationBase.cpp
  ******************************************************************************
  * The base class for defining multiplicatively decomposed deformation
  ******************************************************************************
  */

#include "tardigrade_DeformationBase.h"
#include "tardigrade_error_tools.h"
#include "Eigen/Dense"

namespace tardigradeHydra{

    template<
        unsigned int rows,
        unsigned int inner,
        unsigned int columns,
        class A_iterator, class B_iterator, class C_iterator
    >
    void DeformationBase::_denseMatrixMultiply(
        const A_iterator &A_begin, const A_iterator &A_end,
        const B_iterator &B_begin, const B_iterator &B_end,
        C_iterator C_begin, C_iterator C_end
    ){
        /*!
         * Dense matrix multiplication of the form \f$ [A][B] = [C] \f$
         *
         * TODO This probably should be moved to tardigrade_vector_tools
         *
         * \param &A_begin: The starting iterator for the A matrix
         * \param &A_end: The stopping iterator for the A matrix
         * \param &B_begin: The starting iterator for the B matrix
         * \param &B_end: The stopping iterator for the B matrix
         * \param &C_begin: The starting iterator for the C matrix
         * \param &C_end: The stopping iterator for the C matrix
         */

        using C_type = typename std::iterator_traits<C_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( A_end - A_begin ) == rows * inner,
            "The size of matrix A is " + std::to_string( ( unsigned int )( A_end - A_begin ) ) + " but it should be " + std::to_string( rows * inner )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( B_end - B_begin ) == inner * columns,
            "The size of matrix B is " + std::to_string( ( unsigned int )( B_end - B_begin ) ) + " but it should be " + std::to_string( inner * columns )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( C_end - C_begin ) == rows * columns,
            "The size of matrix C is " + std::to_string( ( unsigned int )( C_end - C_begin ) ) + " but it should be " + std::to_string( rows * columns )
        );

        std::fill( C_begin, C_end, C_type( ) );

        for ( unsigned int i = 0; i < rows; ++i ){
            for ( unsigned int j = 0; j < inner; ++j ){
                for ( unsigned int k = 0; k < columns; ++k ){
                    *( C_begin + columns * i + k ) += ( *( A_begin + inner * i + j ) ) * ( *( B_begin + columns * j + k ) );
                }
            }
        }

    }

    template<
        unsigned int size,
        class configuration_iterator,
        class output_iterator
    >
    void DeformationBase::getNetConfiguration(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Construct a net configuration from the iterator using an assumed multiplicative decomposition
         * where each configuration is a square matrix of dimension size x size i.e.,
         *
         * \f$ A_{af} = A_{ab} A_{cd} \cdots A_{ef} \f$
         *
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<configuration_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == size * size,
            "The size of the output is " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " and should be " + std::to_string( size * size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( configurations_end - configurations_begin ) >= size * size,
            "The provided configurations have a size of " + std::to_string( ( unsigned int )( configurations_end - configurations_begin ) ) + " but must have a size of at least " + std::to_string( size * size )
        );


        if ( configurations_end != ( configurations_begin + size * size ) ){

            // Get preceeding net configuration
            std::array< output_type, size * size > Aminus;
            getNetConfiguration<size>( configurations_begin + size * size, configurations_end, std::begin( Aminus ), std::end( Aminus ) );

            // Update the output
            _denseMatrixMultiply<size,size,size>(
                configurations_begin, configurations_begin + size * size,
                std::begin( Aminus ), std::end( Aminus ),
                output_begin, output_end
            );

        }
        else{

            std::copy(
                configurations_begin, configurations_end, output_begin
            );

        }

    }

    template<
        unsigned int size,
        class configuration_iterator,
        class output_iterator
    >
    void DeformationBase::getLeadingNetConfigurationJacobian(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute the Jacobian of a net configuration with respect to the first configuration e.g.,
         *
         * Given \f$ [A] = [B] [C] [D] \f$, compute the derivative of \f$ [A] \f$ with respect to \f$ [B] \f$
         *
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<configuration_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * size * size ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should have a size of " + std::to_string( size * size * size * size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            configurations_end != configurations_begin,
            "The configurations vector has no size"
        );

        // Handle the case where the configuration array only contains one configuration
        std::fill( output_begin, output_end, output_type( ) );

        if ( output_end == ( output_begin + size * size ) ){

            for ( unsigned int i = 0; i < size * size; ++i ){ *( output_begin + size * i + i ) += 1; }

        }

        std::array< output_type, size * size > Aminus;

        getNetConfiguration<size>(
            configurations_begin + size * size, configurations_end,
            std::begin( Aminus ), std::end( Aminus )
        );

        for ( unsigned int i = 0; i < size; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int k = 0; k < size; ++k ){
                    *( output_begin + size * size * size * i + size * size * j + size * i + k ) += Aminus[ size * k + j ];
                }
            }
        }

    }

    template<
        unsigned int size,
        class configuration_iterator,
        class output_iterator
    >
    void DeformationBase::getTrailingNetConfigurationJacobian(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute the Jacobian of a net configuration with respect to the first configuration e.g.,
         *
         * Given \f$ [A] = [B] [C] [D] \f$, compute the derivative of \f$ [A] \f$ with respect to \f$ [D] \f$
         *
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<configuration_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * size * size ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should have a size of " + std::to_string( size * size * size * size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            configurations_end != configurations_begin,
            "The configurations vector has no size"
        );

        // Handle the case where the configuration array only contains one configuration
        std::fill( output_begin, output_end, output_type( ) );

        if ( output_end == ( output_begin + size * size ) ){

            for ( unsigned int i = 0; i < size * size; ++i ){ *( output_begin + size * i + i ) += 1; }

        }

        std::array< output_type, size * size > Aplus;

        getNetConfiguration<size>(
            configurations_begin, configurations_end - size * size,
            std::begin( Aplus ), std::end( Aplus )
        );

        for ( unsigned int i = 0; i < size; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int k = 0; k < size; ++k ){
                    *( output_begin + size * size * size * i + size * size * j + size * k + j ) += Aplus[ size * i + k ];
                }
            }
        }
    }

    template<
        unsigned int size,
        class configuration_iterator,
        class output_iterator
    >
    void DeformationBase::getNetConfigurationJacobian(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const unsigned int &configuration_index, output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute the Jacobian of a net configuration with respect to an internal configuration e.g.,
         *
         * Given \f$ [A] = [B] [C] [D] \f$, compute the derivative of \f$ [A] \f$ with respect to \f$ [C] \f$
         *
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_index: The index of the configuration to compute the Jacobian for
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<configuration_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * size * size ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should have a size of " + std::to_string( size * size * size * size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            configurations_end != configurations_begin,
            "The configurations vector has no size"
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( configurations_end - configurations_begin ) % ( size * size ) == 0,
            "The configurations iterator has a size of " + std::to_string( ( unsigned int )( configurations_end - configurations_begin ) ) + " which is not a multiple of " + std::to_string( size * size )
        )

        const unsigned int num_configurations = ( unsigned int )( configurations_end - configurations_begin ) / ( size * size );

        if ( configuration_index == 0 ){

            getLeadingNetConfigurationJacobian<size>(
                configurations_begin, configurations_end, output_begin, output_end
            );

        }
        else if ( ( configuration_index + 1 ) == num_configurations ){

            getTrailingNetConfigurationJacobian<size>(
                configurations_begin, configurations_end, output_begin, output_end
            );

        }
        else if ( ( 0 < configuration_index ) && ( configuration_index < ( num_configurations - 1 ) ) ){

            std::fill( output_begin, output_end, output_type( ) );
            
            // Get the prior and previous configurations
            std::array< output_type, size * size > Aplus, Aminus;
            getNetConfiguration<size>(
                configurations_begin, configurations_begin + size * size * configuration_index,
                std::begin( Aplus ), std::end( Aplus )
            );
            getNetConfiguration<size>(
                configurations_begin + size * size * ( configuration_index + 1 ), configurations_end,
                std::begin( Aminus ), std::end( Aminus )
            );
            
            // Assemble the Jacobian
            for ( unsigned int i = 0; i < size; ++i ){
                for ( unsigned int j = 0; j < size; ++j ){
                    for ( unsigned int k = 0; k < size; ++k ){
                        for ( unsigned int l = 0; l < size; ++l ){
                            *( output_begin + size * size * size * i + size * size * j + size * k + l ) += Aplus[ size * i + k ] * Aminus[ size * l + j ];
                        }
                    }
                }
            }

        }
        else{

            std::fill( output_begin, output_end, output_type( ) );

        }

    }

    template<
        unsigned int size,
        unsigned int dim,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class output_iterator
    >
    void DeformationBase::getNetConfigurationGradient(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute the gradient of a net configuration e.g., given a configuration
         *
         * \f$ [A] = [B][C][D] \f$
         *
         * compute
         *
         * \f$ \frac{\partial [A]}{\partial X} = \frac{\partial [B]}{\partial X} [C][D] + [B]\frac{\partial [C]}{\partial X} [D] + [B][C]\frac{\partial [D]}{\partial X} \f$
         *
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<configuration_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) == ( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) ),
            "The number of configurations from the sub configurations is " + std::to_string( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) + " but the number of configurations from the gradients is " + std::to_string( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * dim ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should be " + std::to_string( size * size * dim )
        )

        if ( configurations_end != ( configurations_begin + size * size ) ){

            std::fill( output_begin, output_end, output_type( ) );

            // Get the following configuration and its gradient
            std::array< output_type, size * size > Aminus;
            std::array< output_type, size * size * dim > dAminusdX;

            getNetConfiguration<size>(
                configurations_begin + size * size, configurations_end,
                std::begin( Aminus ), std::end( Aminus )
            );

            getNetConfigurationGradient<size,dim>(
                configurations_begin + size * size, configurations_end,
                configuration_gradients_begin + size * size * dim, configuration_gradients_end,
                std::begin( dAminusdX ), std::end( dAminusdX )
            );

            // Assemble the configuration gradient
            for ( unsigned int i = 0; i < size; ++i ){
                for ( unsigned int j = 0; j < size; ++j ){
                    for ( unsigned int a = 0; a < dim; ++a ){
                        for ( unsigned int l = 0; l < size; ++l ){
                            *( output_begin + size * dim * i + dim * j + a ) += ( *( configuration_gradients_begin + size * dim * i + dim * l + a ) ) * Aminus[ size * l + j ]
                                                                              + ( *( configurations_begin + size * i + l ) ) * dAminusdX[ size * dim * l + dim * j + a ];
                        }
                    }
                }
            }

        }
        else{

            std::copy(
                configuration_gradients_begin, configuration_gradients_end, output_begin
            );

        }

    }

    template<
        unsigned int size,
        unsigned int dim,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class output_iterator
    >
    void DeformationBase::getLeadingNetConfigurationGradientConfigurationJacobian(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute the Jacobian gradient of a net configuration with respect to the leading configuration e.g., given a configuration
         *
         * \f$ [A] = [B][C][D] \f$
         *
         * compute
         *
         * \f$ \frac{\partial^2 [A]}{\partial X \partial [B]} = \mathbbold{I} \left( \frac{\partial [C]}{\partial X} [D] + [C] \frac{\partial [D]}{\partial X}\right) \f$
         *
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<configuration_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) == ( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) ),
            "The number of configurations from the sub configurations is " + std::to_string( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) + " but the number of configurations from the gradients is " + std::to_string( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * dim * size * size ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should be " + std::to_string( size * size * dim * size * size )
        );

        // Initialize the output vector
        std::fill(
            output_begin, output_end, output_type( )
        );

        if ( ( unsigned int )( configurations_end - configurations_begin ) > ( size * size ) ){

            std::array< output_type, size * size * dim > dAminusdX;
            getNetConfigurationGradient<size,dim>(
                configurations_begin + size * size, configurations_end,
                configuration_gradients_begin + size * size * dim, configuration_gradients_end,
                std::begin( dAminusdX ), std::end( dAminusdX )
            );

            for ( unsigned int i = 0; i < size; ++i ){
                for ( unsigned int j = 0; j < size; ++j ){
                    for ( unsigned int a = 0; a < dim; ++a ){
                        for ( unsigned int b = 0; b < size; ++b ){
                            *( output_begin + size * dim * size * size * i + dim * size * size * j + size * size * a + size * i + b ) += dAminusdX[ size * dim * b + dim * j + a ];
                        }
                    }
                }
            }

        }

    }

    template<
        unsigned int size,
        unsigned int dim,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class output_iterator
    >
    void DeformationBase::getTrailingNetConfigurationGradientConfigurationJacobian(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute the Jacobian gradient of a net configuration with respect to the trailing configuration e.g., given a configuration
         *
         * \f$ [A] = [B][C][D] \f$
         *
         * compute
         *
         * \f$ \frac{\partial^2 [A]}{\partial X \partial [D]} = \left( \frac{\partial [B]}{\partial X} [C] + [B] \frac{\partial [C]}{\partial X}\right) \mathbbold{I} \f$
         *
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<configuration_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) == ( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) ),
            "The number of configurations from the sub configurations is " + std::to_string( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) + " but the number of configurations from the gradients is " + std::to_string( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * dim * size * size ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should be " + std::to_string( size * size * dim * size * size )
        );

        // Initialize the output vector
        std::fill(
            output_begin, output_end, output_type( )
        );

        if ( ( unsigned int )( configurations_end - configurations_begin ) > ( size * size ) ){

            std::array< output_type, size * size * dim > dAplusdX;
            getNetConfigurationGradient<size,dim>(
                configurations_begin, configurations_end - size * size,
                configuration_gradients_begin, configuration_gradients_end - size * size * dim,
                std::begin( dAplusdX ), std::end( dAplusdX )
            );

            for ( unsigned int i = 0; i < size; ++i ){
                for ( unsigned int j = 0; j < size; ++j ){
                    for ( unsigned int a = 0; a < dim; ++a ){
                        for ( unsigned int b = 0; b < size; ++b ){
                            *( output_begin + size * dim * size * size * i + dim * size * size * j + size * size * a + size * b + j ) += dAplusdX[ size * dim * i + dim * b + a ];
                        }
                    }
                }
            }

        }

    }

    template<
        unsigned int size,
        unsigned int dim,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class output_iterator
    >
    void DeformationBase::getNetConfigurationGradientConfigurationJacobian(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        const unsigned int &configuration_index,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute the Jacobian gradient of a net configuration with respect to an arbitrary configuration e.g., given a configuration
         *
         * \f$ [A] = [B][C][D] \f$
         *
         * compute
         *
         * \f$ \frac{\partial^2 [A]}{\partial X \partial [D]} = \left( \frac{\partial [B]}{\partial X} [C] + [B] \frac{\partial [C]}{\partial X}\right) \mathbbold{I} \f$
         *
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param &configuration_index: The index of the configuration to compute the Jacobian for
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<configuration_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) == ( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) ),
            "The number of configurations from the sub configurations is " + std::to_string( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) + " but the number of configurations from the gradients is " + std::to_string( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * dim * size * size ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should be " + std::to_string( size * size * dim * size * size )
        );

        const unsigned int num_configurations = ( unsigned int )( configurations_end - configurations_begin ) / ( size * size );

        if ( configuration_index == 0 ){

            getLeadingNetConfigurationGradientConfigurationJacobian<size,dim>(
                configurations_begin, configurations_end, configuration_gradients_begin, configuration_gradients_end, output_begin, output_end
            );

        }
        else if ( ( configuration_index + 1 ) == num_configurations ){

            getTrailingNetConfigurationGradientConfigurationJacobian<size,dim>(
                configurations_begin, configurations_end, configuration_gradients_begin, configuration_gradients_end, output_begin, output_end
            );

        }
        else if ( ( 0 < configuration_index ) && ( configuration_index < ( num_configurations - 1 ) ) ){

            std::fill( output_begin, output_end, output_type( ) );
            
            // Get the prior and previous configurations
            std::array< output_type, size * size > Aplus;
            std::array< output_type, size * size * dim > dAplusdX;
            getNetConfiguration<size>(
                configurations_begin, configurations_begin + size * size * configuration_index,
                std::begin( Aplus ), std::end( Aplus )
            );
            getNetConfigurationGradient<size,dim>(
                configurations_begin, configurations_begin + size * size * configuration_index,
                configuration_gradients_begin, configuration_gradients_begin + size * size * dim * configuration_index,
                std::begin( dAplusdX ), std::end( dAplusdX )
            );

            // Get the prior and previous configuration Jacobians
            std::array< output_type, size * size * size * size > J_Aminus;
            std::array< output_type, size * size * dim * size * size > J_dAminusdX;

            getLeadingNetConfigurationJacobian<size>(
                configurations_begin + size * size * configuration_index, configurations_end,
                std::begin( J_Aminus ), std::end( J_Aminus )
            );

            getLeadingNetConfigurationGradientConfigurationJacobian<size,dim>(
                configurations_begin + size * size * configuration_index, configurations_end,
                configuration_gradients_begin + size * size * dim * configuration_index, configuration_gradients_end,
                std::begin( J_dAminusdX ), std::end( J_dAminusdX )
            );

            // Assemble the Jacobian
            for ( unsigned int i = 0; i < size; ++i ){
                for ( unsigned int j = 0; j < size; ++j ){
                    for ( unsigned int a = 0; a < dim; ++a ){
                        for ( unsigned int bc = 0; bc < size * size; ++bc ){
                            for ( unsigned int l = 0; l < size; ++l ){
                                *( output_begin + size * dim * size * size * i + dim * size * size * j + size * size * a + bc )
                                    += dAplusdX[ size * dim * i + dim * l + a ] * J_Aminus[ size * size * size * l + size * size * j + bc ]
                                     + Aplus[ size * i + l ] * J_dAminusdX[ size * dim * size * size * l + dim * size * size * j + size * size * a + bc ];
                            }
                        }
                    }
                }
            }

        }
        else{

            std::fill( output_begin, output_end, output_type( ) );

        }

    }

    template<
        unsigned int size,
        unsigned int dim,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class output_iterator
    >
    void DeformationBase::getLeadingNetConfigurationGradientConfigurationGradientJacobian(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute the Jacobian gradient of a net configuration with respect to the leading configuration's gradient e.g., given a configuration
         *
         * \f$ [A] = [B][C][D] \f$
         *
         * compute
         *
         * \f$ \frac{\partial^2 [A]}{\partial X \partial \frac{\partial [B]}{\partial X} } = \mathbbold{I} [C] [D] \f$
         *
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<configuration_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) == ( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) ),
            "The number of configurations from the sub configurations is " + std::to_string( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) + " but the number of configurations from the gradients is " + std::to_string( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * dim * size * size * dim ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should be " + std::to_string( size * size * dim * size * size * dim )
        );

        // Initialize the output vector
        std::fill(
            output_begin, output_end, output_type( )
        );

        if ( ( unsigned int )( configurations_end - configurations_begin ) > ( size * size ) ){

            std::array< output_type, size * size > Aminus;
            getNetConfiguration<size>(
                configurations_begin + size * size, configurations_end,
                std::begin( Aminus ), std::end( Aminus )
            );

            for ( unsigned int i = 0; i < size; ++i ){
                for ( unsigned int j = 0; j < size; ++j ){
                    for ( unsigned int a = 0; a < dim; ++a ){
                        for ( unsigned int k = 0; k < size; ++k ){
                            *( output_begin + size * dim * size * size * dim * i + dim * size * size * dim * j + size * size * dim * a + size * dim * i + dim * k + a )
                                += Aminus[ size * k + j ];
                        }
                    }
                }
            }
        }
    }

    template<
        unsigned int size,
        unsigned int dim,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class output_iterator
    >
    void DeformationBase::getTrailingNetConfigurationGradientConfigurationGradientJacobian(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute the Jacobian gradient of a net configuration with respect to the trailing configuration's gradient e.g., given a configuration
         *
         * \f$ [A] = [B][C][D] \f$
         *
         * compute
         *
         * \f$ \frac{\partial^2 [A]}{\partial X \partial \frac{\partial [D]}{\partial X} } = [B][C]\mathbbold{I} \f$
         *
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<configuration_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) == ( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) ),
            "The number of configurations from the sub configurations is " + std::to_string( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) + " but the number of configurations from the gradients is " + std::to_string( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * dim * size * size * dim ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should be " + std::to_string( size * size * dim * size * size * dim )
        );

        // Initialize the output vector
        std::fill(
            output_begin, output_end, output_type( )
        );

        if ( ( unsigned int )( configurations_end - configurations_begin ) > ( size * size ) ){

            std::array< output_type, size * size > Aplus;
            getNetConfiguration<size>(
                configurations_begin, configurations_end - size * size,
                std::begin( Aplus ), std::end( Aplus )
            );

            for ( unsigned int i = 0; i < size; ++i ){
                for ( unsigned int ja = 0; ja < size * dim; ++ja ){
                    for ( unsigned int l = 0; l < size; ++l ){
                        *( output_begin + size * dim * size * size * dim * i + size * size * dim * ja + size * dim * l + ja ) += Aplus[ size * i + l ];
                    }
                }
            }

        }

    }

    template<
        unsigned int size,
        unsigned int dim,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class output_iterator
    >
    void DeformationBase::getNetConfigurationGradientConfigurationGradientJacobian(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        const unsigned int &configuration_index,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute the Jacobian gradient of a net configuration with respect to an arbitrary configuration gradient e.g., given a configuration
         *
         * \f$ [A] = [B][C][D] \f$
         *
         * compute
         *
         * \f$ \frac{\partial^2 [A]}{\partial X \partial \frac{\partial [C]}{\partial X} } = [b] \mathbbold{I} [D] \f$
         *
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param &configuration_index: The index of the configuration to compute the Jacobian for
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<configuration_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) == ( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) ),
            "The number of configurations from the sub configurations is " + std::to_string( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) + " but the number of configurations from the gradients is " + std::to_string( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * dim * size * size * dim ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should be " + std::to_string( size * size * dim * size * size * dim )
        );

        const unsigned int num_configurations = ( unsigned int )( configurations_end - configurations_begin ) / ( size * size );

        if ( configuration_index == 0 ){

            getLeadingNetConfigurationGradientConfigurationGradientJacobian<size,dim>(
                configurations_begin, configurations_end, configuration_gradients_begin, configuration_gradients_end, output_begin, output_end
            );

        }
        else if ( ( configuration_index + 1 ) == num_configurations ){

            getTrailingNetConfigurationGradientConfigurationGradientJacobian<size,dim>(
                configurations_begin, configurations_end, configuration_gradients_begin, configuration_gradients_end, output_begin, output_end
            );

        }
        else if ( ( 0 < configuration_index ) && ( configuration_index < ( num_configurations - 1 ) ) ){

            std::fill( output_begin, output_end, output_type( ) );
            
            // Get the prior and previous configurations
            std::array< output_type, size * size > Aplus, Aminus;
            getNetConfiguration<size>(
                configurations_begin, configurations_begin + size * size * configuration_index,
                std::begin( Aplus ), std::end( Aplus )
            );
            getNetConfiguration<size>(
                configurations_begin + size * size * ( configuration_index + 1 ), configurations_end,
                std::begin( Aminus ), std::end( Aminus )
            );

            // Assemble the Jacobian
            for ( unsigned int i = 0; i < size; ++i ){
                for ( unsigned int j = 0; j < size; ++j ){
                    for ( unsigned int a = 0; a < dim; ++a ){
                        for ( unsigned int k = 0; k < size; ++k ){
                            for ( unsigned int l = 0; l < size; ++l ){
                                *( output_begin + size * dim * size * size * dim * i + dim * size * size * dim * j + size * size * dim * a + size * dim * k + dim * l + a )
                                    += Aplus[ size * i + k ] * Aminus[ size * l + j ];
                            }
                        }
                    }
                }
            }

        }
        else{

            std::fill( output_begin, output_end, output_type( ) );

        }

    }

    template<
        unsigned int leading_rows,
        unsigned int size,
        class total_configuration_iterator,
        class configuration_iterator,
        class Aminus_inverse_iterator,
        class output_iterator
    >
    void DeformationBase::solveForLeadingConfiguration(
        const total_configuration_iterator   &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Solve for the leading configuration which would be required to achieve the total deformation i.e., if
         * the total deformation is \f$ [A] \f$ and we know the net deformation from the subsequent deformations
         * in the form of the configurations, then
         * 
         * \f$ [A] = [B] [A^{-}] \rightarrow [B] = [A] [A^{-}]^{-1} \f$
         *
         * \param &total_configuration_begin: The starting iterator of the total deformation
         * \param &total_configuration_end: The stopping iterator of the total deformation
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &Aminus_inverse_begin: The starting iterator of the inverse of the total deformation represented by the configurations
         * \param &Aminus_inverse_end: The stopping iterator of the inverse of the total deformation represented by the configuraitons
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type         = typename std::iterator_traits<output_iterator>::value_type;
        using Aminus_inverse_type = typename std::iterator_traits<Aminus_inverse_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( total_configuration_end - total_configuration_begin ) == ( unsigned int )( output_end - output_begin ),
            "The total deformation has a size of " + std::to_string( ( unsigned int )( total_configuration_end - total_configuration_begin ) ) + " but the output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == leading_rows * size,
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but it needs a size of " + std::to_string( leading_rows * size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( Aminus_inverse_end - Aminus_inverse_begin ) == size * size,
            "The inverse of the total deformation of the configurations has a size of " + std::to_string( ( unsigned int )( Aminus_inverse_end - Aminus_inverse_begin ) ) + " but it needs a size of " + std::to_string( size * size )
        );

        if ( configurations_end == configurations_begin ){

            std::fill( Aminus_inverse_begin, Aminus_inverse_end, Aminus_inverse_type( ) );

            std::copy(
                total_configuration_begin, total_configuration_end, output_begin
            );

        }
        else{

            std::array< output_type, size * size > Aminus;

            getNetConfiguration<size>(
                configurations_begin, configurations_end,
                std::begin( Aminus ), std::end( Aminus )
            );

            // TODO: Generalize this to a matrix solve rather than computing an inverse
            Eigen::Map< Eigen::Matrix<output_type, size, size, Eigen::RowMajor> > _Aminus( Aminus.data( ), size, size );
            Eigen::Map< Eigen::Matrix<output_type, size, size, Eigen::RowMajor> > _Aminus_inverse( &(*Aminus_inverse_begin ), size, size );
            _Aminus_inverse = _Aminus.inverse( );

            std::fill(
                output_begin, output_end, output_type( )
            );

            _denseMatrixMultiply<
                leading_rows, size, size
            >
            (
                total_configuration_begin, total_configuration_end,
                Aminus_inverse_begin, Aminus_inverse_end,
                output_begin, output_end
            );
        }

    }

    template<
        unsigned int leading_rows,
        unsigned int size,
        class total_configuration_iterator,
        class configuration_iterator,
        class output_iterator
    >
    void DeformationBase::solveForLeadingConfiguration(
        const total_configuration_iterator   &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Solve for the leading configuration which would be required to achieve the total deformation i.e., if
         * the total deformation is \f$ [A] \f$ and we know the net deformation from the subsequent deformations
         * in the form of the configurations, then
         * 
         * \f$ [A] = [B] [A^{-}] \rightarrow [B] = [A] [A^{-}]^{-1} \f$
         *
         * \param &total_configuration_begin: The starting iterator of the total deformation
         * \param &total_configuration_end: The stopping iterator of the total deformation
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using configuration_type = typename std::iterator_traits<configuration_iterator>::value_type;

        std::array< configuration_type, size * size > Aminus_inverse;

        solveForLeadingConfiguration<leading_rows, size>(
            total_configuration_begin, total_configuration_end, configurations_begin, configurations_end,
            std::begin( Aminus_inverse ), std::end( Aminus_inverse ),
            output_begin, output_end
        );

    }

    template<
        unsigned int leading_rows,
        unsigned int size,
        class total_configuration_iterator,
        class configuration_iterator,
        class output_iterator
    >
    void DeformationBase::solveForLeadingConfigurationDeformationJacobian(
        const total_configuration_iterator   &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Solve for the Jacobian of the leading configuration with respect to the total deformation i.e., if
         * the total deformation is \f$ [A] \f$ and we know the net deformation from the subsequent deformations
         * in the form of the configurations, then
         * 
         * \f$ [A] = [B] [A^{-}] \rightarrow [B] = [A] [A^{-}]^{-1} \f$
         *
         * \param &total_configuration_begin: The starting iterator of the total deformation
         * \param &total_configuration_end: The stopping iterator of the total deformation
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<output_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( total_configuration_end - total_configuration_begin ) * ( unsigned int )( total_configuration_end - total_configuration_begin ) == ( unsigned int )( output_end - output_begin ),
            "The Jacobian should have a size of " + std::to_string( ( unsigned int )( total_configuration_end - total_configuration_begin ) * ( unsigned int )( total_configuration_end - total_configuration_begin ) ) + " but the output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == leading_rows * size * leading_rows * size,
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but it needs a size of at least " + std::to_string( leading_rows * size * leading_rows * size )
        );

        if ( configurations_end == ( configurations_begin + size * size ) ){

            std::fill( output_begin, output_end, output_type( ) );

            for ( unsigned int i = 0; i < leading_rows * size; ++i ){

                *( output_begin + leading_rows * size * i + i ) += 1;

            }

        }
        else{

            std::array< output_type, size * size > Aminus, Aminus_inverse;

            getNetConfiguration<size>(
                configurations_begin, configurations_end,
                std::begin( Aminus ), std::end( Aminus )
            );

            // TODO: Generalize this to a matrix solve rather than computing an inverse
            Eigen::Map< Eigen::Matrix<output_type, size, size, Eigen::RowMajor> > _Aminus( Aminus.data( ), size, size );
            Eigen::Map< Eigen::Matrix<output_type, size, size, Eigen::RowMajor> > _Aminus_inverse( Aminus_inverse.data( ), size, size );
            _Aminus_inverse = _Aminus.inverse( );

            std::fill(
                output_begin, output_end, output_type( )
            );

            for ( unsigned int i = 0; i < leading_rows; ++i ){
                for ( unsigned int j = 0; j < size; ++j ){
                    for ( unsigned int a = 0; a < size; ++a ){
                        *( output_begin + size * leading_rows * size * i + leading_rows * size * j + size * i + a ) += Aminus_inverse[ size * a + j ];
                    }
                }
            }
        }

    }

    template<
        unsigned int leading_rows,
        unsigned int size,
        class total_configuration_iterator,
        class configuration_iterator,
        class output_iterator
    >
    void DeformationBase::solveForLeadingConfigurationConfigurationJacobian(
        const total_configuration_iterator   &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const unsigned int &configuration_index,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Solve for the Jacobian of the leading configuration with respect to the specified configuration i.e., if
         * the total deformation is \f$ [A] \f$ and we know the net deformation from the subsequent deformations
         * in the form of the configurations, then
         * 
         * \f$ [A] = [B] [A^{-}] \rightarrow [B] = [A] [A^{-}]^{-1} \f$
         *
         * \param &total_configuration_begin: The starting iterator of the total deformation
         * \param &total_configuration_end: The stopping iterator of the total deformation
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_index: The index of the configuration in the configuration iterator to
         *     compute the Jacobian with respect to
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<output_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == leading_rows * size * size * size,
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but it needs a size of at least " + std::to_string( size * size )
        );

        if ( configurations_end == ( configurations_begin + size * size ) ){

            std::fill( output_begin, output_end, output_type( ) );

        }
        else{

            std::array< output_type, leading_rows * size > leadingConfiguration;
            std::array< output_type, size * size > Aminus_inverse;
            std::array< output_type, leading_rows * size * size * size > intermediate_term;
            std::array< output_type, size * size * size * size > Aminus_jacobian;

            solveForLeadingConfiguration<leading_rows,size>(
                total_configuration_begin, total_configuration_end, configurations_begin, configurations_end,
                std::begin( Aminus_inverse ), std::end( Aminus_inverse ),
                std::begin( leadingConfiguration ), std::end( leadingConfiguration )
            );

            getNetConfigurationJacobian<size>(
                configurations_begin, configurations_end,
                configuration_index,
                std::begin( Aminus_jacobian ), std::end( Aminus_jacobian )
            );

            std::fill(
                std::begin( intermediate_term ), std::end( intermediate_term ), output_type( )
            );

            for ( unsigned int i = 0; i < leading_rows; ++i ){
                for ( unsigned int j = 0; j < size; ++j ){
                    for ( unsigned int c = 0; c < size; ++c ){
                        for ( unsigned int d = 0; d < size; ++d ){
                            intermediate_term[ size * size * size * i + size * size * j + size * c + d ]
                                -= leadingConfiguration[ size * i + c ] * Aminus_inverse[ size * d + j ];
                        }
                    }
                }
            }

            _denseMatrixMultiply<
                leading_rows * size,
                size * size,
                size * size
            >(
                std::begin( intermediate_term ), std::end( intermediate_term ),
                std::begin( Aminus_jacobian ),   std::end( Aminus_jacobian ),
                output_begin, output_end
            );

        }

    }

    template<
        unsigned int leading_rows,
        unsigned int size,
        unsigned int dim,
        class total_configuration_gradient_iterator,
        class leading_configuration_iterator,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class Aminus_inverse_iterator,
        class output_iterator
    >
    void DeformationBase::solveForLeadingConfigurationGradient(
        const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
        const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Solve for the leading configuration gradient which would be required to achieve the total configuration gradient i.e., if
         * the total deformation is \f$ [A] \f$ and we know the net deformation from the subsequent deformations
         * in the form of the configurations, then
         *
         * \f$ \frac{\partial [B]}{\partial X} = \frac{\partial [A]}{\partial X} A^{-} + [A] \frac{\partial [A]^{-}}{\partial X} \f$
         *
         * which means we can solve for \f$ \frac{\partial [A]}{\partial X} \f$ via
         *
         * \f$ \frac{\partial [A]}{\partial X} = \left(\frac{\partial [B]}{\partial X} - [A] \frac{\partial [A]^{-}}{\partial X}\right) \left([A]^{-}\right)^{-1} \f$
         *
         * \param &total_configuration_gradient_begin: The starting iterator of the total deformation gradient.
         *     Note that this deformation gradient is the derivative of the deformation \f$ [B] \f$ with respect
         *     to \f$ X \f$ rather than the standard deformation gradient from continuum (i.e., \f$ \bf{F} \f$)
         * \param &total_configuration_gradient_end: The stopping iterator of the total deformation gradient.
         *     Note that this deformation gradient is the derivative of the deformation \f$ [B] \f$ with respect
         *     to \f$ X \f$ rather than the standard deformation gradient from continuum (i.e., \f$ \bf{F} \f$)
         * \param &leading_configuration_begin: The starting iterator of the leading configuration
         * \param &leading_configuration_end: The stopping iterator of the leading configuration
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param Aminus_inverse_begin: The starting iterator of the inverse of the total trailing configuration
         * \param Aminus_inverse_end: The stopping iterator of the inverse of the total trailing configuration
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using configuration_type = typename std::iterator_traits<configuration_iterator>::value_type;
        using configuration_gradient_type = typename std::iterator_traits<configuration_gradient_iterator>::value_type;
        using output_type = typename std::iterator_traits<output_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) == leading_rows * size * dim,
            "The total deformation gradient has a size of " + std::to_string( ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) ) + " but must have a size of " + std::to_string( leading_rows * size * dim )
        )

        std::array< configuration_type, size * size > Aminus;
        std::array< configuration_gradient_type, size * size * dim > dAMinusdX;
        std::array< output_type, leading_rows * size * dim > intermediate_term1, intermediate_term2;

        // Compute the trailing configuration and it's gradient
        getNetConfiguration<size>(
            configurations_begin, configurations_end,
            std::begin( Aminus ), std::end( Aminus )
        );

        getNetConfigurationGradient<size,dim>(
            configurations_begin, configurations_end,
            configuration_gradients_begin, configuration_gradients_end,
            std::begin( dAMinusdX ), std::end( dAMinusdX )
        );

        // Assemble the intermediate terms
        _denseMatrixMultiply<leading_rows,size,size*dim>(
            leading_configuration_begin, leading_configuration_end,
            std::begin( dAMinusdX ), std::end( dAMinusdX ),
            std::begin( intermediate_term1 ), std::end( intermediate_term1 )
        );

        std::transform(
            total_configuration_gradient_begin, total_configuration_gradient_end,
            std::begin( intermediate_term1 ), std::begin( intermediate_term2 ),
            std::minus<>()
        );

        // TODO: Generalize this to a matrix solve rather than computing an inverse
        Eigen::Map< Eigen::Matrix<output_type, size, size, Eigen::RowMajor> > _Aminus( Aminus.data( ), size, size );
        Eigen::Map< Eigen::Matrix<output_type, size, size, Eigen::RowMajor> > _Aminus_inverse( &(*Aminus_inverse_begin), size, size );
        _Aminus_inverse = _Aminus.inverse( );

        std::fill(
            output_begin, output_end, output_type( )
        );

        for ( unsigned int i = 0; i < leading_rows; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int a = 0; a < dim; ++a ){
                    for ( unsigned int k = 0; k < size; ++k ){
                        *( output_begin + size * dim * i + dim * j + a ) += intermediate_term2[ size * dim * i + dim * k + a ] * ( *( Aminus_inverse_begin + size * k + j ) );
                    }
                }
            }
        }

    }

    template<
        unsigned int leading_rows,
        unsigned int size,
        unsigned int dim,
        class total_configuration_gradient_iterator,
        class leading_configuration_iterator,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class output_iterator
    >
    void DeformationBase::solveForLeadingConfigurationGradient(
        const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
        const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Solve for the leading configuration gradient which would be required to achieve the total configuration gradient i.e., if
         * the total deformation is \f$ [A] \f$ and we know the net deformation from the subsequent deformations
         * in the form of the configurations, then
         *
         * \f$ \frac{\partial [B]}{\partial X} = \frac{\partial [A]}{\partial X} A^{-} + [A] \frac{\partial [A]^{-}}{\partial X} \f$
         *
         * which means we can solve for \f$ \frac{\partial [A]}{\partial X} \f$ via
         *
         * \f$ \frac{\partial [A]}{\partial X} = \left(\frac{\partial [B]}{\partial X} - [A] \frac{\partial [A]^{-}}{\partial X}\right) \left([A]^{-}\right)^{-1} \f$
         *
         * \param &total_configuration_gradient_begin: The starting iterator of the total deformation gradient.
         *     Note that this deformation gradient is the derivative of the deformation \f$ [B] \f$ with respect
         *     to \f$ X \f$ rather than the standard deformation gradient from continuum (i.e., \f$ \bf{F} \f$)
         * \param &total_configuration_gradient_end: The stopping iterator of the total deformation gradient.
         *     Note that this deformation gradient is the derivative of the deformation \f$ [B] \f$ with respect
         *     to \f$ X \f$ rather than the standard deformation gradient from continuum (i.e., \f$ \bf{F} \f$)
         * \param &leading_configuration_begin: The starting iterator of the leading configuration
         * \param &leading_configuration_end: The stopping iterator of the leading configuration
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using configuration_type = typename std::iterator_traits<configuration_iterator>::value_type;
        std::array< configuration_type, size * size > Aminus;

        solveForLeadingConfigurationGradient<leading_rows, size, dim>(
            total_configuration_gradient_begin, total_configuration_gradient_end,
            leading_configuration_begin, leading_configuration_end,
            configurations_begin, configurations_end,
            configuration_gradients_begin, configuration_gradients_end,
            std::begin( Aminus ), std::end( Aminus ),
            output_begin, output_end
        );

    }

    template<
        unsigned int leading_rows,
        unsigned int size,
        unsigned int dim,
        class total_configuration_gradient_iterator,
        class leading_configuration_iterator,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class output_iterator
    >
    void DeformationBase::solveForLeadingConfigurationGradientTotalConfigurationGradientJacobian(
        const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
        const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Solve for the jacobian with respect to the total deformation gradient of the leading configuration gradient
         * which would be required to achieve the total configuration gradient i.e., if the total deformation is
         * \f$ [A] \f$ and we know the net deformation from the subsequent deformations in the form of the configurations, then
         *
         * \f$ \frac{\partial [B]}{\partial X} = \frac{\partial [A]}{\partial X} A^{-} + [A] \frac{\partial [A]^{-}}{\partial X} \f$
         *
         * which means we can solve for \f$ \frac{\partial [A]}{\partial X} \f$ via
         *
         * \f$ \frac{\partial [A]}{\partial X} = \left(\frac{\partial [B]}{\partial X} - [A] \frac{\partial [A]^{-}}{\partial X}\right) \left([A]^{-}\right)^{-1} \f$
         *
         * \param &total_configuration_gradient_begin: The starting iterator of the total deformation gradient.
         *     Note that this deformation gradient is the derivative of the deformation \f$ [B] \f$ with respect
         *     to \f$ X \f$ rather than the standard deformation gradient from continuum (i.e., \f$ \bf{F} \f$)
         * \param &total_configuration_gradient_end: The stopping iterator of the total deformation gradient.
         *     Note that this deformation gradient is the derivative of the deformation \f$ [B] \f$ with respect
         *     to \f$ X \f$ rather than the standard deformation gradient from continuum (i.e., \f$ \bf{F} \f$)
         * \param &leading_configuration_begin: The starting iterator of the leading configuration
         * \param &leading_configuration_end: The stopping iterator of the leading configuration
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using configuration_type = typename std::iterator_traits<configuration_iterator>::value_type;
        using output_type = typename std::iterator_traits<output_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) == leading_rows * size * dim,
            "The total deformation gradient has a size of " + std::to_string( ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) ) + " but must have a size of " + std::to_string( leading_rows * size * dim )
        )

        std::array< configuration_type, size * size > Aminus, Aminus_inverse;

        // Compute the trailing configuration and it's gradient
        getNetConfiguration<size>(
            configurations_begin, configurations_end,
            std::begin( Aminus ), std::end( Aminus )
        );

        // TODO: Generalize this to a matrix solve rather than computing an inverse
        Eigen::Map< Eigen::Matrix<output_type, size, size, Eigen::RowMajor> > _Aminus( Aminus.data( ), size, size );
        Eigen::Map< Eigen::Matrix<output_type, size, size, Eigen::RowMajor> > _Aminus_inverse( Aminus_inverse.data( ), size, size );
        _Aminus_inverse = _Aminus.inverse( );

        for ( unsigned int i = 0; i < leading_rows; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int a = 0; a < dim; ++a ){
                    for ( unsigned int b = 0; b < size; ++b ){
                        *( output_begin + size * dim * leading_rows * size * dim * i + dim * leading_rows * size * dim * j + leading_rows * size * dim * a + size * dim * i + dim * b + a )
                            += Aminus_inverse[ size * b + j ];
                    }
                }
            }
        }

    }

    template<
        unsigned int leading_rows,
        unsigned int size,
        unsigned int dim,
        class total_configuration_gradient_iterator,
        class leading_configuration_iterator,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class output_iterator
    >
    void DeformationBase::solveForLeadingConfigurationGradientLeadingConfigurationJacobian(
        const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
        const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Solve for the Jacobian with respect to the leading configuration of the leading configuration gradient
         * which would be required to achieve the total configuration gradient i.e., if the total deformation is
         * \f$ [A] \f$ and we know the net deformation from the subsequent deformations in the form of the configurations, then
         *
         * \f$ \frac{\partial [B]}{\partial X} = \frac{\partial [A]}{\partial X} A^{-} + [A] \frac{\partial [A]^{-}}{\partial X} \f$
         *
         * which means we can solve for \f$ \frac{\partial [A]}{\partial X} \f$ via
         *
         * \f$ \frac{\partial [A]}{\partial X} = \left(\frac{\partial [B]}{\partial X} - [A] \frac{\partial [A]^{-}}{\partial X}\right) \left([A]^{-}\right)^{-1} \f$
         *
         * \param &total_configuration_gradient_begin: The starting iterator of the total deformation gradient.
         *     Note that this deformation gradient is the derivative of the deformation \f$ [B] \f$ with respect
         *     to \f$ X \f$ rather than the standard deformation gradient from continuum (i.e., \f$ \bf{F} \f$)
         * \param &total_configuration_gradient_end: The stopping iterator of the total deformation gradient.
         *     Note that this deformation gradient is the derivative of the deformation \f$ [B] \f$ with respect
         *     to \f$ X \f$ rather than the standard deformation gradient from continuum (i.e., \f$ \bf{F} \f$)
         * \param &leading_configuration_begin: The starting iterator of the leading configuration
         * \param &leading_configuration_end: The stopping iterator of the leading configuration
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using configuration_type = typename std::iterator_traits<configuration_iterator>::value_type;
        using configuration_gradient_type = typename std::iterator_traits<configuration_gradient_iterator>::value_type;
        using output_type = typename std::iterator_traits<output_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) == leading_rows * size * dim,
            "The total deformation gradient has a size of " + std::to_string( ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) ) + " but must have a size of " + std::to_string( leading_rows * size * dim )
        )

        std::array< configuration_type, size * size > Aminus, Aminus_inverse;
        std::array< configuration_gradient_type, size * size * dim > dAMinusdX;

        // Compute the trailing configuration and it's gradient
        getNetConfiguration<size>(
            configurations_begin, configurations_end,
            std::begin( Aminus ), std::end( Aminus )
        );

        getNetConfigurationGradient<size,dim>(
            configurations_begin, configurations_end,
            configuration_gradients_begin, configuration_gradients_end,
            std::begin( dAMinusdX ), std::end( dAMinusdX )
        );

        // TODO: Generalize this to a matrix solve rather than computing an inverse
        Eigen::Map< Eigen::Matrix<output_type, size, size, Eigen::RowMajor> > _Aminus( Aminus.data( ), size, size );
        Eigen::Map< Eigen::Matrix<output_type, size, size, Eigen::RowMajor> > _Aminus_inverse( Aminus_inverse.data( ), size, size );
        _Aminus_inverse = _Aminus.inverse( );

        std::fill(
            output_begin, output_end, output_type( )
        );

        //J_ija = -A_ik dAminusdX_kla Aminus_inv_lj
        //J_ija_ib = -dAminusdX_bla Aminus_inv_lj

        for ( unsigned int i = 0; i < leading_rows; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int a = 0; a < dim; ++a ){
                    for ( unsigned int b = 0; b < size; ++b ){
                        for ( unsigned int l = 0; l < size; ++l ){
                            *( output_begin + size * dim * leading_rows * size * i + dim * leading_rows * size * j + leading_rows * size * a + size * i + b )
                                -= dAMinusdX[ size * dim * b + dim * l + a ] * Aminus_inverse[ size * l + j ];
                        }
                    }
                }
            }
        }

    }

}
