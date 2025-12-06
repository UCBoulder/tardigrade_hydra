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
        C_iterator C_begin, C_iterator C_end,
        const unsigned int offset, const unsigned int stride
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
         * \param offset: The column-wise shift applied to the output
         * \param stride: The size of each row of output
         */

        using C_type = typename std::iterator_traits<C_iterator>::value_type;

        std::fill( C_begin, C_end, C_type( ) );

        _denseMatrixMultiplyAccumulate<rows,inner,columns>(
            A_begin, A_end, B_begin, B_end, C_begin, C_end,
            offset, stride
        );

    }

    template<
        unsigned int rows,
        unsigned int inner,
        unsigned int columns,
        class A_iterator, class B_iterator, class C_iterator
    >
    void DeformationBase::_denseMatrixMultiplyAccumulate(
        const A_iterator &A_begin, const A_iterator &A_end,
        const B_iterator &B_begin, const B_iterator &B_end,
        C_iterator C_begin, C_iterator C_end,
        const unsigned int offset, const unsigned int stride
    ){
        /*!
         * Dense matrix multiplication of the form \f$ [A][B] = [C] \f$
         * where the values are accumulated into whatever is already in C
         *
         * TODO This probably should be moved to tardigrade_vector_tools
         *
         * \param &A_begin: The starting iterator for the A matrix
         * \param &A_end: The stopping iterator for the A matrix
         * \param &B_begin: The starting iterator for the B matrix
         * \param &B_end: The stopping iterator for the B matrix
         * \param &C_begin: The starting iterator for the C matrix
         * \param &C_end: The stopping iterator for the C matrix
         * \param offset: The column-wise shift applied to the output
         * \param stride: The size of each row of output
         */

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( A_end - A_begin ) == rows * inner,
            "The size of matrix A is " + std::to_string( ( unsigned int )( A_end - A_begin ) ) + " but it should be " + std::to_string( rows * inner )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( B_end - B_begin ) == inner * columns,
            "The size of matrix B is " + std::to_string( ( unsigned int )( B_end - B_begin ) ) + " but it should be " + std::to_string( inner * columns )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( C_end - C_begin ) == rows * stride,
            "The size of matrix C is " + std::to_string( ( unsigned int )( C_end - C_begin ) ) + " but it should be " + std::to_string( rows * stride )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            stride >= columns,
            "The stride is a size of " + std::to_string( stride ) + " and is less than the number of columns " + std::to_string( columns )
        );

        for ( unsigned int i = 0; i < rows; ++i ){
            for ( unsigned int j = 0; j < inner; ++j ){
                for ( unsigned int k = 0; k < columns; ++k ){
                    *( C_begin + stride * i + k + offset ) += ( *( A_begin + inner * i + j ) ) * ( *( B_begin + columns * j + k ) );
                }
            }
        }

    }

    template<
        unsigned int size,
        class A_inverse_iterator, class output_iterator
    >
    void DeformationBase::_assembledAinversedA(
        const A_inverse_iterator &A_inverse_begin, const A_inverse_iterator &A_inverse_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Construct the derivative of the inverse of a matrix w.r.t. the matrix
         *
         * \param &A_inverse_begin: The starting iterator of the inverse matrix
         * \param &A_inverse_end: The stopping iterator of the inverse matrix
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( A_inverse_end - A_inverse_begin ) == size * size,
            "The inverse matrix has a size of " + std::to_string( ( unsigned int )( A_inverse_end - A_inverse_begin ) ) + " but must have a size of " + std::to_string( size * size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == size * size * size * size,
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but must have a size of " + std::to_string( size * size * size * size )
        );

        for ( unsigned int i = 0; i < size; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int k = 0; k < size; ++k ){
                    for ( unsigned int l = 0; l < size; ++l ){
                        *( output_begin + size * size * size * i + size * size * j + size * k + l )
                            = -( *( A_inverse_begin + size * i + k ) ) * ( *( A_inverse_begin + size * l + j ) );
                    }
                }
            }
        }
    }

    template<
        unsigned int size,
        class A_iterator,
        class output_iterator
    >
    void DeformationBase::_compute_matrix_inverse(
        const A_iterator &A_begin, const A_iterator &A_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute the inverse of a matrix
         *
         * \param &A_begin: The starting iterator of the matrix
         * \param &A_end: The stopping iterator of the matrix
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<output_iterator>::value_type;
        Eigen::Map< Eigen::Matrix<output_type, size, size, Eigen::RowMajor> > _A( &(*A_begin), size, size );
        Eigen::Map< Eigen::Matrix<output_type, size, size, Eigen::RowMajor> > _A_inverse( &(*output_begin ), size, size );
        _A_inverse = _A.inverse( );

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

        using output_type = typename std::iterator_traits<output_iterator>::value_type;

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
        class Aminus_iterator,
        class output_iterator
    >
    void DeformationBase::_assemble_output_getLeadingNetConfigurationJacobian(
        const Aminus_iterator &Aminus_begin, const Aminus_iterator &Aminus_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Assemble the jacobian of getLeadingNetConfigurationJacobian
         *
         * Given \f$ [A] = [B] [C] [D] \f$, compute the derivative of \f$ [A] \f$ with respect to \f$ [B] \f$
         *
         * \param &Aminus_begin: The starting iterator of the trailing configuration
         * \param &Aminus_end: The stopping iterator of the trailing configuration
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<output_iterator>::value_type;

        std::fill( output_begin, output_end, output_type( ) );

        for ( unsigned int i = 0; i < size; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int k = 0; k < size; ++k ){
                    *( output_begin + size * size * size * i + size * size * j + size * i + k ) += ( *( Aminus_begin + size * k + j ) );
                }
            }
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

        using output_type = typename std::iterator_traits<output_iterator>::value_type;
        using configuration_type = typename std::iterator_traits<configuration_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * size * size ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should have a size of " + std::to_string( size * size * size * size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            configurations_end != configurations_begin,
            "The configurations vector has no size"
        );

        // Handle the case where the configuration array only contains one configuration

        if ( configurations_end == ( configurations_begin + size * size ) ){

            std::fill( output_begin, output_end, output_type( ) );
            for ( unsigned int i = 0; i < size * size; ++i ){ *( output_begin + size * size * i + i ) += 1; }
            return;

        }

        std::array< configuration_type, size * size > Aminus;

        getNetConfiguration<size>(
            configurations_begin + size * size, configurations_end,
            std::begin( Aminus ), std::end( Aminus )
        );

        _assemble_output_getLeadingNetConfigurationJacobian<size>(
            std::begin( Aminus ), std::end( Aminus ),
            output_begin, output_end
        );

    }

    template<
        unsigned int size,
        class Aplus_iterator,
        class output_iterator
    >
    void DeformationBase::_assemble_output_getTrailingNetConfigurationJacobian(
        const Aplus_iterator &Aplus_begin, const Aplus_iterator &Aplus_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Assemble the output for getTrailingNetConfigurationJacobian
         *
         * Given \f$ [A] = [B] [C] [D] \f$, compute the derivative of \f$ [A] \f$ with respect to \f$ [D] \f$
         *
         * \param &Aplus_begin: The starting iterator of the leading configurations
         * \param &Aplus_end: The stopping iterator of the leading configurations
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<output_iterator>::value_type;

        std::fill( output_begin, output_end, output_type( ) );

        for ( unsigned int i = 0; i < size; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int k = 0; k < size; ++k ){
                    *( output_begin + size * size * size * i + size * size * j + size * k + j ) += ( *( Aplus_begin + size * i + k ) );
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

        using output_type = typename std::iterator_traits<output_iterator>::value_type;
        using configuration_type = typename std::iterator_traits<configuration_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * size * size ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should have a size of " + std::to_string( size * size * size * size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            configurations_end != configurations_begin,
            "The configurations vector has no size"
        );

        // Handle the case where the configuration array only contains one configuration
        if ( configurations_end == ( configurations_begin + size * size ) ){

            std::fill( output_begin, output_end, output_type( ) );
            for ( unsigned int i = 0; i < size * size; ++i ){ *( output_begin + size * size * i + i ) += 1; }
            return;

        }

        std::array< configuration_type, size * size > Aplus;

        getNetConfiguration<size>(
            configurations_begin, configurations_end - size * size,
            std::begin( Aplus ), std::end( Aplus )
        );

        _assemble_output_getTrailingNetConfigurationJacobian<size>(
            std::begin( Aplus ), std::end( Aplus ),
            output_begin, output_end
        );

    }

    template<
        unsigned int size,
        class Aplus_iterator,
        class Aminus_iterator,
        class output_iterator
    >
    void DeformationBase::_assemble_output_getNetConfigurationJacobian(
        const Aplus_iterator &Aplus_begin, const Aplus_iterator &Aplus_end,
        const Aminus_iterator &Aminus_begin, const Aminus_iterator &Aminus_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Assemble the Jacobian for an arbitrarily located configuration
         *
         * \param Aplus_begin: The starting iterator of the leading configuration
         * \param Aplus_end: The stopping iterator of the leading configuration
         * \param Aminus_begin: The starting iterator of the trailing configuration
         * \param Aminus_end: The stopping iterator of the trailing configuration
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<output_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * size * size ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should have a size of " + std::to_string( size * size * size * size )
        );

        std::fill( output_begin, output_end, output_type( ) );

        // Assemble the Jacobian
        for ( unsigned int i = 0; i < size; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int k = 0; k < size; ++k ){
                    for ( unsigned int l = 0; l < size; ++l ){
                        *( output_begin + size * size * size * i + size * size * j + size * k + l ) += ( *( Aplus_begin + size * i + k ) ) * ( *( Aminus_begin + size * l + j ) );
                    }
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
        const unsigned int &configuration_index,
        output_iterator output_begin, output_iterator output_end
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

        using configuration_type = typename std::iterator_traits<configuration_iterator>::value_type;
        using output_type = typename std::iterator_traits<output_iterator>::value_type;

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
            std::array< configuration_type, size * size > Aplus, Aminus;

            getNetConfiguration<size>(
                configurations_begin, configurations_begin + size * size * configuration_index,
                std::begin( Aplus ), std::end( Aplus )
            );

            getNetConfiguration<size>(
                configurations_begin + size * size * ( configuration_index + 1 ), configurations_end,
                std::begin( Aminus ), std::end( Aminus )
            );

            _assemble_output_getNetConfigurationJacobian<size>(
                std::begin( Aplus ), std::end( Aplus ),
                std::begin( Aminus ), std::end( Aminus ),
                output_begin, output_end
            );

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
        class Aminus_iterator,
        class dAminusdX_iterator,
        class output_iterator
    >
    void DeformationBase::_assemble_dAdX_getNetConfigurationGradient(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        const Aminus_iterator &Aminus_begin, const Aminus_iterator &Aminus_end,
        const dAminusdX_iterator &dAminusdX_begin, const dAminusdX_iterator &dAminusdX_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param &Aminus_begin: The starting iterator of the net subsequent configuration
         * \param &Aminus_end: The stopping iterator of the net subsequent configuration
         * \param &dAminusdX_begin: The starting iterator of the net subsequent configuration gradient
         * \param &dAminusdX_end: The stopping iterator of the net subsequent configuration gradient
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<configuration_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * dim ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should be " + std::to_string( size * size * dim )
        )

        std::fill( output_begin, output_end, output_type( ) );

        // Assemble the configuration gradient
        for ( unsigned int i = 0; i < size; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int a = 0; a < dim; ++a ){
                    for ( unsigned int l = 0; l < size; ++l ){
                        *( output_begin + size * dim * i + dim * j + a ) += ( *( configuration_gradients_begin + size * dim * i + dim * l + a ) ) * ( *( Aminus_begin + size * l + j ) )
                                                                          + ( *( configurations_begin + size * i + l ) ) * ( *( dAminusdX_begin + size * dim * l + dim * j + a ) );
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

        using configuration_type = typename std::iterator_traits<configuration_iterator>::value_type;
        using configuration_gradient_type = typename std::iterator_traits<configuration_gradient_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) == ( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) ),
            "The number of configurations from the sub configurations is " + std::to_string( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) + " but the number of configurations from the gradients is " + std::to_string( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) )
        );

        if ( configurations_end != ( configurations_begin + size * size ) ){

            // Get the following configuration and its gradient
            std::array< configuration_type, size * size > Aminus;
            std::array< configuration_gradient_type, size * size * dim > dAminusdX;

            getNetConfiguration<size>(
                configurations_begin + size * size, configurations_end,
                std::begin( Aminus ), std::end( Aminus )
            );

            getNetConfigurationGradient<size,dim>(
                configurations_begin + size * size, configurations_end,
                configuration_gradients_begin + size * size * dim, configuration_gradients_end,
                std::begin( dAminusdX ), std::end( dAminusdX )
            );

            _assemble_dAdX_getNetConfigurationGradient<size,dim>(
                configurations_begin, configurations_end,
                configuration_gradients_begin, configuration_gradients_end,
                std::begin( Aminus ), std::end( Aminus ),
                std::begin( dAminusdX ), std::end( dAminusdX ),
                output_begin, output_end
            );

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
        class dAminusdX_iterator,
        class output_iterator
    >
    void DeformationBase::_assemble_output_getLeadingNetConfigurationGradientConfigurationJacobian(
        const dAminusdX_iterator &dAminusdX_begin, const dAminusdX_iterator &dAminusdX_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Assemble the output of getLeadingNetConfigurationGradientConfigurationJacobian
         * 
         * \param dAminusdX_begin: The starting iterator of the trailing configuration gradient
         * \param dAminusdX_end: The stopping iterator of the trailing configuration gradient
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<output_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * dim * size * size ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should be " + std::to_string( size * size * dim * size * size )
        );

        std::fill(
            output_begin, output_end, output_type( )
        );

        for ( unsigned int i = 0; i < size; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int a = 0; a < dim; ++a ){
                    for ( unsigned int b = 0; b < size; ++b ){
                        *( output_begin + size * dim * size * size * i + dim * size * size * j + size * size * a + size * i + b ) += ( *( dAminusdX_begin + size * dim * b + dim * j + a ) );
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

        using configuration_type = typename std::iterator_traits<configuration_iterator>::value_type;
        using output_type = typename std::iterator_traits<output_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) == ( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) ),
            "The number of configurations from the sub configurations is " + std::to_string( ( unsigned int )( configurations_end - configurations_begin ) / ( size * size ) ) + " but the number of configurations from the gradients is " + std::to_string( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) / ( size * size * dim ) )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == ( size * size * dim * size * size ),
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but should be " + std::to_string( size * size * dim * size * size )
        );

        if ( ( unsigned int )( configurations_end - configurations_begin ) > ( size * size ) ){

            std::fill(
                output_begin, output_end, output_type( )
            );
            std::array< configuration_type, size * size * dim > dAminusdX;
            getNetConfigurationGradient<size,dim>(
                configurations_begin + size * size, configurations_end,
                configuration_gradients_begin + size * size * dim, configuration_gradients_end,
                std::begin( dAminusdX ), std::end( dAminusdX )
            );

            _assemble_output_getLeadingNetConfigurationGradientConfigurationJacobian<size,dim>(
                std::begin( dAminusdX ), std::end( dAminusdX ),
                output_begin, output_end
            );

        }
        else{
            std::fill(
                output_begin, output_end, output_type( )
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

        using output_type = typename std::iterator_traits<output_iterator>::value_type;

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
        class Aminus_inverse_iterator,
        class output_iterator
    >
    void DeformationBase::_assemble_leading_configuration_solveForLeadingConfiguration(
        const total_configuration_iterator   &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
        Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Assemble the leading configuration
         *
         * \param &total_configuration_begin: The starting iterator of the total deformation
         * \param &total_configuration_end: The stopping iterator of the total deformation
         * \param &Aminus_inverse_begin: The starting iterator of the inverse of the total deformation represented by the configurations
         * \param &Aminus_inverse_end: The stopping iterator of the inverse of the total deformation represented by the configuraitons
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<output_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == leading_rows * size,
            "The output has a size of " + std::to_string( ( unsigned int )( output_end - output_begin ) ) + " but it needs a size of " + std::to_string( leading_rows * size )
        );

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
            _compute_matrix_inverse<size>(
                std::begin( Aminus ), std::end( Aminus ),
                Aminus_inverse_begin, Aminus_inverse_end
            );

            _assemble_leading_configuration_solveForLeadingConfiguration<leading_rows,size>(
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
    void DeformationBase::solveForLeadingConfigurationTotalConfigurationJacobian(
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
            _compute_matrix_inverse<size>(
                std::begin( Aminus ), std::end( Aminus ),
                std::begin( Aminus_inverse ), std::end( Aminus_inverse )
            );

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
        class leading_configuration_iterator,
        class Aminus_inverse_iterator,
        class output_iterator
    >
    void DeformationBase::_compute_intermediate_term_1_solveForLeadingConfigurationConfigurationJacobian(
        const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
        const Aminus_inverse_iterator &Aminus_inverse_begin, const Aminus_inverse_iterator &Aminus_inverse_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute intermediate term 1 for solveForLeadingConfigurationGradientConfigurationJacobian
         *
         * \param &leading_configuration_begin: The starting iterator of the leading configuration
         * \param &leading_configuration_end: The stopping iterator of the leading configuration
         * \param Aminus_inverse_begin: The starting iterator of the inverse of the total trailing configuration
         * \param Aminus_inverse_end: The stopping iterator of the inverse of the total trailing configuration
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using output_type = typename std::iterator_traits<output_iterator>::value_type;
        std::fill( output_begin, output_end, output_type( ) );

        for ( unsigned int i = 0; i < leading_rows; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int c = 0; c < size; ++c ){
                    for ( unsigned int d = 0; d < size; ++d ){
                        *( output_begin + size * size * size * i + size * size * j + size * c + d )
                            -= ( *( leading_configuration_begin + size * i + c ) ) * ( *( Aminus_inverse_begin + size * d + j ) );
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

        if ( configurations_end == configurations_begin ){

            std::fill( output_begin, output_end, output_type( ) );

        }
        else{

            std::array< output_type, leading_rows * size > leadingConfiguration;
            std::array< output_type, size * size > Aminus_inverse;
            std::array< output_type, leading_rows * size * size * size > intermediate_term = { output_type( ) };
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

            _compute_intermediate_term_1_solveForLeadingConfigurationConfigurationJacobian<leading_rows,size>(
                std::begin( leadingConfiguration ), std::end( leadingConfiguration ),
                std::begin( Aminus_inverse ), std::end( Aminus_inverse ),
                std::begin( intermediate_term ), std::end( intermediate_term )
            );

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
        class dAminusdX_iterator,
        class output_iterator
    >
    void DeformationBase::solveForLeadingConfigurationGradient(
        const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
        const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
        dAminusdX_iterator dAminusdX_begin, dAminusdX_iterator dAminusdX_end,
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
         * \param dAminusdX_begin: The starting iterator of the gradient of the total trailing configuration
         * \param dAminusdX_end: The stopping iterator of the gradient of the total trailing configuration
         * \param output_begin: The starting iterator of the output
         * \param output_end: The stopping iterator of the output
         */

        using configuration_type = typename std::iterator_traits<configuration_iterator>::value_type;
        using output_type = typename std::iterator_traits<output_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) == leading_rows * size * dim,
            "The total deformation gradient has a size of " + std::to_string( ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) ) + " but must have a size of " + std::to_string( leading_rows * size * dim )
        )

        std::array< configuration_type, size * size > Aminus;
        std::array< output_type, leading_rows * size * dim > intermediate_term1, intermediate_term2;

        // Compute the trailing configuration and it's gradient
        getNetConfiguration<size>(
            configurations_begin, configurations_end,
            std::begin( Aminus ), std::end( Aminus )
        );

        getNetConfigurationGradient<size,dim>(
            configurations_begin, configurations_end,
            configuration_gradients_begin, configuration_gradients_end,
            dAminusdX_begin, dAminusdX_end
        );

        // Assemble the intermediate terms
        _denseMatrixMultiply<leading_rows,size,size*dim>(
            leading_configuration_begin, leading_configuration_end,
            dAminusdX_begin, dAminusdX_end,
            std::begin( intermediate_term1 ), std::end( intermediate_term1 )
        );

        std::transform(
            total_configuration_gradient_begin, total_configuration_gradient_end,
            std::begin( intermediate_term1 ), std::begin( intermediate_term2 ),
            std::minus<>()
        );

        // TODO: Generalize this to a matrix solve rather than computing an inverse
        _compute_matrix_inverse<size>(
            std::begin( Aminus ), std::end( Aminus ),
            Aminus_inverse_begin, Aminus_inverse_end
        );

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
        using configuration_gradient_type = typename std::iterator_traits<configuration_gradient_iterator>::value_type;
        std::array< configuration_type, size * size > Aminus;
        std::array< configuration_gradient_type, size * size * dim > dAminusdx;

        solveForLeadingConfigurationGradient<leading_rows, size, dim>(
            total_configuration_gradient_begin, total_configuration_gradient_end,
            leading_configuration_begin, leading_configuration_end,
            configurations_begin, configurations_end,
            configuration_gradients_begin, configuration_gradients_end,
            std::begin( Aminus ), std::end( Aminus ),
            std::begin( dAminusdx ), std::end( dAminusdx ),
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

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) == leading_rows * size * dim,
            "The total deformation gradient has a size of " + std::to_string( ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) ) + " but must have a size of " + std::to_string( leading_rows * size * dim )
        )

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == leading_rows * size * dim * leading_rows * size * dim,
            "The jacobian has a size of " + std::to_string( ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) ) + " but must have a size of " + std::to_string( leading_rows * size * dim * leading_rows * size * dim )
        )

        std::array< configuration_type, size * size > Aminus, Aminus_inverse;

        // Compute the trailing configuration and it's gradient
        getNetConfiguration<size>(
            configurations_begin, configurations_end,
            std::begin( Aminus ), std::end( Aminus )
        );

        // TODO: Generalize this to a matrix solve rather than computing an inverse
        _compute_matrix_inverse<size>(
            std::begin( Aminus ), std::end( Aminus ),
            std::begin( Aminus_inverse ), std::end( Aminus_inverse )
        );

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

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == leading_rows * size * dim * leading_rows * size,
            "The jacobian has a size of " + std::to_string( ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) ) + " but must have a size of " + std::to_string( leading_rows * size * dim * leading_rows * size )
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
        _compute_matrix_inverse<size>(
            std::begin( Aminus ), std::end( Aminus ),
            std::begin( Aminus_inverse ), std::end( Aminus_inverse )
        );

        std::fill(
            output_begin, output_end, output_type( )
        );

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
    void DeformationBase::solveForLeadingConfigurationGradientConfigurationJacobian(
        const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
        const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        const unsigned int &configuration_index,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Solve for the Jacobian with respect to a configuration of the leading configuration gradient
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
         * \param &configuration_index: The index of the configuration in the configurations array to compute the Jacobian with respect to
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

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == leading_rows * size * dim * size * size,
            "The jacobian has a size of " + std::to_string( ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) ) + " but must have a size of " + std::to_string( leading_rows * size * dim * size * size )
        )

        std::array< configuration_gradient_type, leading_rows * size * dim > leadingConfigurationGradient;
        std::array< configuration_type, size * size * size * size > J_Aminus;
        std::array< configuration_type, size * size > Aminus_inverse;
        std::array< configuration_gradient_type, size * size * dim * size * size > J_dAminusdX;
        std::array< configuration_gradient_type, size * size * dim > dAminusdX;

        // Compute the Jacobian of the trailing configuration and it's gradient
        getNetConfigurationJacobian<size>(
            configurations_begin, configurations_end,
            configuration_index,
            std::begin( J_Aminus ), std::end( J_Aminus )
        );

        getNetConfigurationGradientConfigurationJacobian<size,dim>(
            configurations_begin, configurations_end,
            configuration_gradients_begin, configuration_gradients_end,
            configuration_index,
            std::begin( J_dAminusdX ), std::end( J_dAminusdX )
        );

        // Compute the leading configuration
        solveForLeadingConfigurationGradient<leading_rows, size, dim>(
            total_configuration_gradient_begin, total_configuration_gradient_end,
            leading_configuration_begin, leading_configuration_end,
            configurations_begin, configurations_end,
            configuration_gradients_begin, configuration_gradients_end,
            std::begin( Aminus_inverse ), std::end( Aminus_inverse ),
            std::begin( dAminusdX ), std::end( dAminusdX ),
            std::begin( leadingConfigurationGradient ), std::end( leadingConfigurationGradient )
        );

        // Assemble the Jacobian
        std::array< output_type, size * size * dim * size * size > intermediate_term1 = { output_type( ) };

        for ( unsigned int i = 0; i < size; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int k = 0; k < size; ++k ){
                    for ( unsigned int acd = 0; acd < dim * size * size; ++acd ){
                        intermediate_term1[ size * dim * size * size * i + dim * size * size * j + acd ]
                            -= J_dAminusdX[ size * dim * size * size * i + dim * size * size * k + acd ] * Aminus_inverse[ size * k + j ];
                    }
                }
            }
        }

        _denseMatrixMultiply<leading_rows,size,size*dim*size*size>(
            leading_configuration_begin, leading_configuration_end,
            std::begin( intermediate_term1 ), std::end( intermediate_term1 ),
            output_begin, output_end
        );

        std::array< output_type, leading_rows * size * dim * size * size > intermediate_term2 = { output_type( ) };

        for ( unsigned int i = 0; i < leading_rows; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int a = 0; a < dim; ++a ){
                    for ( unsigned int e = 0; e < size; ++e ){
                        for ( unsigned int f = 0; f < size; ++f ){
                            intermediate_term2[ size * dim * size * size * i + dim * size * size * j + size * size * a + size * e + f ]
                                -= leadingConfigurationGradient[ size * dim * i + dim * e + a ] * Aminus_inverse[ size * f + j ];
                        }
                    }
                }
            }
        }

        _denseMatrixMultiplyAccumulate<leading_rows * size * dim, size * size, size * size>(
            std::begin( intermediate_term2 ), std::end( intermediate_term2 ),
            std::begin( J_Aminus ), std::end( J_Aminus ),
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
    void DeformationBase::solveForLeadingConfigurationGradientConfigurationGradientJacobian(
        const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
        const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        const unsigned int &configuration_index,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Solve for the Jacobian with respect to a configuration gradient of the leading configuration gradient
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
         * \param &configuration_index: The index of the configuration in the configurations array to compute the Jacobian with respect to
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

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_end - output_begin ) == leading_rows * size * dim * size * size * dim,
            "The jacobian has a size of " + std::to_string( ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) ) + " but must have a size of " + std::to_string( leading_rows * size * dim * size * size * dim )
        )

        std::array< configuration_type, size * size > Aminus, Aminus_inverse;
        std::array< configuration_gradient_type, size * size * dim * size * size * dim > J_dAminusdX;

        getNetConfiguration<size>(
            configurations_begin, configurations_end,
            std::begin( Aminus ), std::end( Aminus )
        );

        // TODO: Generalize this to a matrix solve rather than computing an inverse
        _compute_matrix_inverse<size>(
            std::begin( Aminus ), std::end( Aminus ),
            std::begin( Aminus_inverse ), std::end( Aminus_inverse )
        );

        // Compute the Jacobian of the trailing configuration and it's gradient
        getNetConfigurationGradientConfigurationGradientJacobian<size,dim>(
            configurations_begin, configurations_end,
            configuration_gradients_begin, configuration_gradients_end,
            configuration_index,
            std::begin( J_dAminusdX ), std::end( J_dAminusdX )
        );

        std::array< output_type, size * size * dim * size * size * dim > intermediate_term = { output_type( ) };

        for ( unsigned int i = 0; i < size; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int k = 0; k < size; ++k ){
                    for ( unsigned int acde = 0; acde < dim * size * size * dim; ++acde ){
                        intermediate_term[ size * dim * size * size * dim * i + dim * size * size * dim * j + acde ]
                            -= J_dAminusdX[ size * dim * size * size * dim * i + dim * size * size * dim * k + acde ] * Aminus_inverse[ size * k + j ];
                    }
                }
            }
        }

        _denseMatrixMultiply<leading_rows,size,size*dim*size*size*dim>(
            leading_configuration_begin, leading_configuration_end,
            std::begin( intermediate_term ), std::end( intermediate_term ),
            output_begin, output_end
        );

    }

    template<
        unsigned int leading_rows,
        unsigned int size,
        unsigned int dim,
        class total_configuration_iterator,
        class total_configuration_gradient_iterator,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class output_leading_configuration_iterator,
        class output_leading_configuration_gradient_iterator
    >
    void DeformationBase::solveForAllLeading(
        const total_configuration_iterator &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
        const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        output_leading_configuration_iterator output_leading_configuration_begin, output_leading_configuration_iterator output_leading_configuration_end,
        output_leading_configuration_gradient_iterator output_leading_configuration_gradient_begin, output_leading_configuration_gradient_iterator output_leading_configuration_gradient_end
    ){
        /*!
         * Solve for the leading configuration and its gradient which would be required to achieve the total deformation i.e., if
         * the total deformation is \f$ [A] \f$ and we know the net deformation from the subsequent deformations
         * in the form of the configurations, then
         * 
         * \f$ [A] = [B] [A^{-}] \rightarrow [B] = [A] [A^{-}]^{-1} \f$
         *
         * \param &total_configuration_begin: The starting iterator of the total deformation
         * \param &total_configuration_end: The stopping iterator of the total deformation
         * \param &total_configuration_gradient_begin: The starting iterator of the total deformation gradient
         * \param &total_configuration_gradient_end: The stopping iterator of the total deformation gradient
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param &configurations_end: The stopping iterator of the configurations
         * \param output_leading_configuration_begin: The starting iterator of the leading configuration output
         * \param output_leading_configuration_end: The stopping iterator of the leading configuration output
         * \param output_leading_configuration_gradient_begin: The starting iterator of the leading configuration gradient output
         * \param output_leading_configuration_gradient_end: The stopping iterator of the leading configuration gradient output
         */

        using configuration_type = typename std::iterator_traits<configuration_iterator>::value_type;
        using configuration_gradient_type = typename std::iterator_traits<configuration_gradient_iterator>::value_type;

        std::array< configuration_type, size * size > Aminus_inverse;
        std::array< configuration_gradient_type, size * size * dim > dAminusdX;

        solveForAllLeading<leading_rows,size,dim>(
            total_configuration_begin, total_configuration_end,
            total_configuration_gradient_begin, total_configuration_gradient_end,
            configurations_begin, configurations_end,
            configuration_gradients_begin, configuration_gradients_end,
            std::begin( Aminus_inverse ), std::end( Aminus_inverse ),
            std::begin( dAminusdX ), std::end( dAminusdX ),
            output_leading_configuration_begin, output_leading_configuration_end,
            output_leading_configuration_gradient_begin, output_leading_configuration_gradient_end
        );
    }

    template<
        unsigned int leading_rows,
        unsigned int size,
        unsigned int dim,
        class total_configuration_iterator,
        class total_configuration_gradient_iterator,
        class Aminus_inverse_iterator,
        class dAminusdX_iterator,
        class output_leading_configuration_iterator,
        class output_leading_configuration_gradient_iterator
    >
    void DeformationBase::_sizeCheck_solveForAllLeading(
        const total_configuration_iterator &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
        const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
        Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
        dAminusdX_iterator dAminusdX_begin, dAminusdX_iterator dAminusdX_end,
        output_leading_configuration_iterator output_leading_configuration_begin, output_leading_configuration_iterator output_leading_configuration_end,
        output_leading_configuration_gradient_iterator output_leading_configuration_gradient_begin, output_leading_configuration_gradient_iterator output_leading_configuration_gradient_end
    ){
        /*!
         * Check the sizes for the solveForAllLeading function
         *
         *
         * \param &total_configuration_begin: The starting iterator of the total deformation
         * \param &total_configuration_end: The stopping iterator of the total deformation
         * \param &total_configuration_gradient_begin: The starting iterator of the total deformation gradient
         * \param &total_configuration_gradient_end: The stopping iterator of the total deformation gradient
         * \param Aminus_inverse_begin: The starting iterator for the net trailing configuration
         * \param Aminus_inverse_end: The stopping iterator for the net trailing configuration
         * \param dAminusdX_begin: The starting iterator for the net trailing configuration gradient
         * \param dAminusdX_end: The stopping iterator for the net trailing configuration gradient
         * \param output_leading_configuration_begin: The starting iterator of the leading configuration output
         * \param output_leading_configuration_end: The stopping iterator of the leading configuration output
         * \param output_leading_configuration_gradient_begin: The starting iterator of the leading configuration gradient output
         * \param output_leading_configuration_gradient_end: The stopping iterator of the leading configuration gradient output
         */

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( total_configuration_end - total_configuration_begin ) == ( leading_rows * size ),
            "The total configuration has a size of " + std::to_string( ( unsigned int )( total_configuration_end - total_configuration_begin ) ) + " but should have a size of " + std::to_string( leading_rows * size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) == ( leading_rows * size * dim ),
            "The total configuration gradient has a size of " + std::to_string( ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) ) + " but should have a size of " + std::to_string( leading_rows * size * dim )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( Aminus_inverse_end - Aminus_inverse_begin ) == size * size,
            "The inverse of the total deformation of the configurations has a size of " + std::to_string( ( unsigned int )( Aminus_inverse_end - Aminus_inverse_begin ) ) + " but it needs a size of " + std::to_string( size * size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( dAminusdX_end - dAminusdX_begin ) == size * size * dim,
            "The gradient of the total deformation of the configurations has a size of " + std::to_string( ( unsigned int )( dAminusdX_end - dAminusdX_begin ) ) + " but it needs a size of " + std::to_string( size * size * dim )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_leading_configuration_end - output_leading_configuration_begin ) == ( leading_rows * size ),
            "The leading configuration has a size of " + std::to_string( ( unsigned int )( output_leading_configuration_end - output_leading_configuration_begin ) ) + " but should have a size of " + std::to_string( leading_rows * size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_leading_configuration_gradient_end - output_leading_configuration_gradient_begin ) == ( leading_rows * size * dim ),
            "The leading configuration gradient has a size of " + std::to_string( ( unsigned int )( output_leading_configuration_gradient_end - output_leading_configuration_gradient_begin ) ) + " but should have a size of " + std::to_string( leading_rows * size * dim )
        );

    }

    template<
        unsigned int leading_rows,
        unsigned int size,
        unsigned int dim,
        class total_configuration_iterator,
        class total_configuration_gradient_iterator,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class Aminus_inverse_iterator,
        class dAminusdX_iterator,
        class output_leading_configuration_iterator,
        class output_leading_configuration_gradient_iterator
    >
    void DeformationBase::solveForAllLeading(
        const total_configuration_iterator &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
        const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
        dAminusdX_iterator dAminusdX_begin, dAminusdX_iterator dAminusdX_end,
        output_leading_configuration_iterator output_leading_configuration_begin, output_leading_configuration_iterator output_leading_configuration_end,
        output_leading_configuration_gradient_iterator output_leading_configuration_gradient_begin, output_leading_configuration_gradient_iterator output_leading_configuration_gradient_end
    ){
        /*!
         * Solve for the leading configuration and its gradient which would be required to achieve the total deformation i.e., if
         * the total deformation is \f$ [A] \f$ and we know the net deformation from the subsequent deformations
         * in the form of the configurations, then
         * 
         * \f$ [A] = [B] [A^{-}] \rightarrow [B] = [A] [A^{-}]^{-1} \f$
         *
         * \param &total_configuration_begin: The starting iterator of the total deformation
         * \param &total_configuration_end: The stopping iterator of the total deformation
         * \param &total_configuration_gradient_begin: The starting iterator of the total deformation gradient
         * \param &total_configuration_gradient_end: The stopping iterator of the total deformation gradient
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param Aminus_inverse_begin: The starting iterator for the net trailing configuration
         * \param Aminus_inverse_end: The stopping iterator for the net trailing configuration
         * \param dAminusdX_begin: The starting iterator for the net trailing configuration gradient
         * \param dAminusdX_end: The stopping iterator for the net trailing configuration gradient
         * \param output_leading_configuration_begin: The starting iterator of the leading configuration output
         * \param output_leading_configuration_end: The stopping iterator of the leading configuration output
         * \param output_leading_configuration_gradient_begin: The starting iterator of the leading configuration gradient output
         * \param output_leading_configuration_gradient_end: The stopping iterator of the leading configuration gradient output
         */

#ifndef TARDIGRADE_ERROR_TOOLS_OPT
        _sizeCheck_solveForAllLeading<leading_rows,size,dim>(
            total_configuration_begin, total_configuration_end,
            total_configuration_gradient_begin, total_configuration_gradient_end,
            Aminus_inverse_begin, Aminus_inverse_end,
            dAminusdX_begin, dAminusdX_end,
            output_leading_configuration_begin, output_leading_configuration_end,
            output_leading_configuration_gradient_begin, output_leading_configuration_gradient_end
        );
#endif

        using output_leading_configuration_gradient_type = typename std::iterator_traits<output_leading_configuration_gradient_iterator>::value_type;

        // Compute the leading configuration and its gradient

        solveForLeadingConfiguration<leading_rows,size>(
            total_configuration_begin, total_configuration_end,
            configurations_begin, configurations_end,
            Aminus_inverse_begin, Aminus_inverse_end,
            output_leading_configuration_begin, output_leading_configuration_end
        );

        std::array< output_leading_configuration_gradient_type, leading_rows * size * dim > intermediate_term1, intermediate_term2;

        // Compute the trailing configuration gradient
        getNetConfigurationGradient<size,dim>(
            configurations_begin, configurations_end,
            configuration_gradients_begin, configuration_gradients_end,
            dAminusdX_begin, dAminusdX_end
        );

        // Assemble the intermediate terms
        _denseMatrixMultiply<leading_rows,size,size*dim>(
            output_leading_configuration_begin, output_leading_configuration_end,
            dAminusdX_begin, dAminusdX_end,
            std::begin( intermediate_term1 ), std::end( intermediate_term1 )
        );

        std::transform(
            total_configuration_gradient_begin, total_configuration_gradient_end,
            std::begin( intermediate_term1 ), std::begin( intermediate_term2 ),
            std::minus<>()
        );

        std::fill(
            output_leading_configuration_gradient_begin, output_leading_configuration_gradient_end, output_leading_configuration_gradient_type( )
        );

        for ( unsigned int i = 0; i < leading_rows; ++i ){
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int a = 0; a < dim; ++a ){
                    for ( unsigned int k = 0; k < size; ++k ){
                        *( output_leading_configuration_gradient_begin + size * dim * i + dim * j + a ) += intermediate_term2[ size * dim * i + dim * k + a ] * ( *( Aminus_inverse_begin + size * k + j ) );
                    }
                }
            }
        }

    }

    template<
        unsigned int leading_rows,
        unsigned int size,
        unsigned int dim,
        class total_configuration_iterator,
        class total_configuration_gradient_iterator,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class output_leading_configuration_total_J_iterator,
        class output_leading_configuration_configurations_J_iterator,
        class output_leading_configuration_configuration_gradients_J_iterator,
        class output_leading_configuration_gradient_total_J_iterator,
        class output_leading_configuration_gradient_total_gradient_J_iterator,
        class output_leading_configuration_gradient_configurations_J_iterator,
        class output_leading_configuration_gradient_configuration_gradients_J_iterator
    >
    void DeformationBase::_sizeCheck_solveForAllLeadingJacobians(
        const total_configuration_iterator &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
        const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        output_leading_configuration_total_J_iterator output_leading_configuration_total_J_begin, output_leading_configuration_total_J_iterator output_leading_configuration_total_J_end,
        output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_begin, output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_end,
        output_leading_configuration_configuration_gradients_J_iterator output_leading_configuration_configuration_gradients_J_begin, output_leading_configuration_configuration_gradients_J_iterator output_leading_configuration_configuration_gradients_J_end,
        output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_begin, output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_end,
        output_leading_configuration_gradient_total_gradient_J_iterator output_leading_configuration_gradient_total_gradient_J_begin, output_leading_configuration_gradient_total_gradient_J_iterator output_leading_configuration_gradient_total_gradient_J_end,
        output_leading_configuration_gradient_configurations_J_iterator output_leading_configuration_gradient_configurations_J_begin, output_leading_configuration_gradient_configurations_J_iterator output_leading_configuration_gradient_configurations_J_end,
        output_leading_configuration_gradient_configuration_gradients_J_iterator output_leading_configuration_gradient_configuration_gradients_J_begin, output_leading_configuration_gradient_configuration_gradients_J_iterator output_leading_configuration_gradient_configuration_gradients_J_end
    ){
        /*!
         * Check the iterator sizes for solveForAllLeadingJacobians
         *
         * \param &total_configuration_begin: The starting iterator of the total deformation
         * \param &total_configuration_end: The stopping iterator of the total deformation
         * \param &total_configuration_gradient_begin: The starting iterator of the total deformation gradient
         * \param &total_configuration_gradient_end: The stopping iterator of the total deformation gradient
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param output_leading_configuration_total_J_begin: The starting iterator of the Jacobian of the leading configuration with respect to the total configuration output
         * \param output_leading_configuration_total_J_end: The stopping iterator of the Jacobian of the leading configuration with respect to the total configuration output
         * \param output_leading_configuration_configurations_J_begin: The starting iterator of the Jacobian of the leading configuration with respect to the configurations output
         * \param output_leading_configuration_configurations_J_end: The stopping iterator of the Jacobian of the leading configuration with respect to the configurations output
         * \param output_leading_configuration_configuration_gradients_J_begin: The starting iterator of the Jacobian of the leading configuration with respect to the configuration gradients output
         * \param output_leading_configuration_configuration_gradients_J_end: The stopping iterator of the Jacobian of the leading configuration with respect to the configuration gradients output
         * \param output_leading_configuration_gradient_total_J_begin: The starting iterator of the Jacobian of the leading configuration gradient with respect to the total configuration output
         * \param output_leading_configuration_gradient_total_J_end: The stopping iterator of the Jacobian of the leading configuration gradient with respect to the total configuration output
         * \param output_leading_configuration_gradient_total_gradient_J_begin: The starting iterator of the Jacobian of the leading configuration gradient with respect to the total configuration gradient output
         * \param output_leading_configuration_gradient_total_gradient_J_end: The stopping iterator of the Jacobian of the leading configuration gradient with respect to the total configuration gradient output
         * \param output_leading_configuration_gradient_configurations_J_begin: The starting iterator of the Jacobian of the leading configuration gradient with respect to the configurations output
         * \param output_leading_configuration_gradient_configurations_J_end: The stopping iterator of the Jacobian of the leading configuration gradient with respect to the configurations output
         * \param output_leading_configuration_gradient_configuration_gradients_J_begin: The starting iterator of the Jacobian of the leading configuration gradient with respect to the configuration gradients output
         * \param output_leading_configuration_gradient_configuration_gradients_J_end: The stopping iterator of the Jacobian of the leading configuration gradient with respect to the configuration gradients output
         */

        const unsigned int num_configs = ( configurations_end - configurations_begin ) / ( size * size );

        const unsigned int leading_configuration_size = leading_rows * size;
        const unsigned int leading_configuration_gradient_size = leading_rows * size * dim;
        const unsigned int configuration_size = size * size;
        const unsigned int configuration_gradient_size = size * size * dim;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( total_configuration_end - total_configuration_begin ) == leading_configuration_size,
            "The total configuration has a size of " + std::to_string( ( unsigned int )( total_configuration_end - total_configuration_begin ) ) + " but must have a size of " + std::to_string( leading_configuration_size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) == leading_configuration_gradient_size,
            "The total configuration gradient has a size of " + std::to_string( ( unsigned int )( total_configuration_gradient_end - total_configuration_gradient_begin ) ) + " but must have a size of " + std::to_string( leading_configuration_gradient_size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( configurations_end - configurations_begin ) == configuration_size * num_configs,
            "The configurations have a size of " + std::to_string( ( unsigned int )( configurations_end - configurations_begin ) ) + " but they must have a size of " + std::to_string( configuration_size * num_configs )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) == configuration_gradient_size * num_configs,
            "The configuration gradients have a size of " + std::to_string( ( unsigned int )( configuration_gradients_end - configuration_gradients_begin ) ) + " but they must have a size of " + std::to_string( configuration_gradient_size * num_configs )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_leading_configuration_total_J_end - output_leading_configuration_total_J_begin ) == leading_configuration_size * leading_configuration_size,
            "The jacobian of the leading configuration with respect to the total deformation has a size of " + std::to_string( ( unsigned int )( output_leading_configuration_total_J_end - output_leading_configuration_total_J_begin ) ) + " but must have a size of " + std::to_string( leading_configuration_size * leading_configuration_size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_leading_configuration_configurations_J_end - output_leading_configuration_configurations_J_begin ) == leading_configuration_size * configuration_size * num_configs,
            "The jacobian of the leading configuration with respect to the configurations has a size of " + std::to_string( ( unsigned int )( output_leading_configuration_configurations_J_end - output_leading_configuration_configurations_J_begin ) ) + " but must have a size of " + std::to_string( leading_configuration_size * configuration_size * num_configs )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_leading_configuration_configuration_gradients_J_end - output_leading_configuration_configuration_gradients_J_begin ) == leading_configuration_size * configuration_gradient_size * num_configs,
            "The jacobian of the leading configuration with respect to the configuration gradients has a size of " + std::to_string( ( unsigned int )( output_leading_configuration_configuration_gradients_J_end - output_leading_configuration_configuration_gradients_J_begin ) ) + " but must have a size of " + std::to_string( leading_configuration_size * configuration_gradient_size * num_configs )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_leading_configuration_gradient_total_J_end - output_leading_configuration_gradient_total_J_begin ) == leading_configuration_gradient_size * leading_configuration_size,
            "The jacobian of the leading configuration gradient with respect to the total deformation has a size of " + std::to_string( ( unsigned int )( output_leading_configuration_gradient_total_J_end - output_leading_configuration_gradient_total_J_begin ) ) + " but must have a size of " + std::to_string( leading_configuration_gradient_size * leading_configuration_size )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_leading_configuration_gradient_configurations_J_end - output_leading_configuration_gradient_configurations_J_begin ) == leading_configuration_gradient_size * configuration_size * num_configs,
            "The jacobian of the leading configuration gradient with respect to the configurations has a size of " + std::to_string( ( unsigned int )( output_leading_configuration_gradient_configurations_J_end - output_leading_configuration_gradient_configurations_J_begin ) ) + " but must have a size of " + std::to_string( leading_configuration_gradient_size * configuration_size * num_configs )
        );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            ( unsigned int )( output_leading_configuration_gradient_configuration_gradients_J_end - output_leading_configuration_gradient_configuration_gradients_J_begin ) == leading_configuration_gradient_size * configuration_gradient_size * num_configs,
            "The jacobian of the leading configuration gradient with respect to the configuration gradients has a size of " + std::to_string( ( unsigned int )( output_leading_configuration_gradient_configuration_gradients_J_end - output_leading_configuration_gradient_configuration_gradients_J_begin ) ) + " but must have a size of " + std::to_string( leading_configuration_gradient_size * configuration_gradient_size * num_configs )
        );

    }

    template<
        unsigned int leading_rows,
        unsigned int size,
        unsigned int dim,
        class total_configuration_iterator,
        class total_configuration_gradient_iterator,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class output_leading_configuration_total_J_iterator,
        class output_leading_configuration_configurations_J_iterator,
        class output_leading_configuration_configuration_gradients_J_iterator,
        class output_leading_configuration_gradient_total_J_iterator,
        class output_leading_configuration_gradient_total_gradient_J_iterator,
        class output_leading_configuration_gradient_configurations_J_iterator,
        class output_leading_configuration_gradient_configuration_gradients_J_iterator
    >
    void DeformationBase::_zeroOutputs_solveForAllLeadingJacobians(
        const total_configuration_iterator &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
        const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        output_leading_configuration_total_J_iterator output_leading_configuration_total_J_begin, output_leading_configuration_total_J_iterator output_leading_configuration_total_J_end,
        output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_begin, output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_end,
        output_leading_configuration_configuration_gradients_J_iterator output_leading_configuration_configuration_gradients_J_begin, output_leading_configuration_configuration_gradients_J_iterator output_leading_configuration_configuration_gradients_J_end,
        output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_begin, output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_end,
        output_leading_configuration_gradient_total_gradient_J_iterator output_leading_configuration_gradient_total_gradient_J_begin, output_leading_configuration_gradient_total_gradient_J_iterator output_leading_configuration_gradient_total_gradient_J_end,
        output_leading_configuration_gradient_configurations_J_iterator output_leading_configuration_gradient_configurations_J_begin, output_leading_configuration_gradient_configurations_J_iterator output_leading_configuration_gradient_configurations_J_end,
        output_leading_configuration_gradient_configuration_gradients_J_iterator output_leading_configuration_gradient_configuration_gradients_J_begin, output_leading_configuration_gradient_configuration_gradients_J_iterator output_leading_configuration_gradient_configuration_gradients_J_end
    ){
        /*!
         * Check the iterator sizes for solveForAllLeadingJacobians
         *
         * \param &total_configuration_begin: The starting iterator of the total deformation
         * \param &total_configuration_end: The stopping iterator of the total deformation
         * \param &total_configuration_gradient_begin: The starting iterator of the total deformation gradient
         * \param &total_configuration_gradient_end: The stopping iterator of the total deformation gradient
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param output_leading_configuration_total_J_begin: The starting iterator of the Jacobian of the leading configuration with respect to the total configuration output
         * \param output_leading_configuration_total_J_end: The stopping iterator of the Jacobian of the leading configuration with respect to the total configuration output
         * \param output_leading_configuration_configurations_J_begin: The starting iterator of the Jacobian of the leading configuration with respect to the configurations output
         * \param output_leading_configuration_configurations_J_end: The stopping iterator of the Jacobian of the leading configuration with respect to the configurations output
         * \param output_leading_configuration_configuration_gradients_J_begin: The starting iterator of the Jacobian of the leading configuration with respect to the configuration gradients output
         * \param output_leading_configuration_configuration_gradients_J_end: The stopping iterator of the Jacobian of the leading configuration with respect to the configuration gradients output
         * \param output_leading_configuration_gradient_total_J_begin: The starting iterator of the Jacobian of the leading configuration gradient with respect to the total configuration output
         * \param output_leading_configuration_gradient_total_J_end: The stopping iterator of the Jacobian of the leading configuration gradient with respect to the total configuration output
         * \param output_leading_configuration_gradient_total_gradient_J_begin: The starting iterator of the Jacobian of the leading configuration gradient with respect to the total configuration gradient output
         * \param output_leading_configuration_gradient_total_gradient_J_end: The stopping iterator of the Jacobian of the leading configuration gradient with respect to the total configuration gradient output
         * \param output_leading_configuration_gradient_configurations_J_begin: The starting iterator of the Jacobian of the leading configuration gradient with respect to the configurations output
         * \param output_leading_configuration_gradient_configurations_J_end: The stopping iterator of the Jacobian of the leading configuration gradient with respect to the configurations output
         * \param output_leading_configuration_gradient_configuration_gradients_J_begin: The starting iterator of the Jacobian of the leading configuration gradient with respect to the configuration gradients output
         * \param output_leading_configuration_gradient_configuration_gradients_J_end: The stopping iterator of the Jacobian of the leading configuration gradient with respect to the configuration gradients output
         */

        using output_lc_total_J_type = typename std::iterator_traits<output_leading_configuration_total_J_iterator>::value_type;
        using output_lc_configurations_J_type = typename std::iterator_traits<output_leading_configuration_configurations_J_iterator>::value_type;
        using output_lc_configuration_gradients_J_type = typename std::iterator_traits<output_leading_configuration_gradient_configuration_gradients_J_iterator>::value_type;
        using output_lcg_total_J_type = typename std::iterator_traits<output_leading_configuration_gradient_total_J_iterator>::value_type;
        using output_lcg_total_gradient_J_type = typename std::iterator_traits<output_leading_configuration_gradient_total_gradient_J_iterator>::value_type;
        using output_lcg_configurations_J_type = typename std::iterator_traits<output_leading_configuration_gradient_configurations_J_iterator>::value_type;
        using output_lcg_configuration_gradients_J_type = typename std::iterator_traits<output_leading_configuration_gradient_configuration_gradients_J_iterator>::value_type;

        std::fill(
            output_leading_configuration_total_J_begin, output_leading_configuration_total_J_end, output_lc_total_J_type( )
        );

        std::fill(
            output_leading_configuration_configurations_J_begin, output_leading_configuration_configurations_J_end, output_lc_configurations_J_type( )
        );

        std::fill(
            output_leading_configuration_configuration_gradients_J_begin, output_leading_configuration_configuration_gradients_J_end, output_lc_configuration_gradients_J_type( )
        );

        std::fill(
            output_leading_configuration_gradient_total_J_begin, output_leading_configuration_gradient_total_J_end, output_lcg_total_J_type( )
        );

        std::fill(
            output_leading_configuration_gradient_total_gradient_J_begin, output_leading_configuration_gradient_total_gradient_J_end, output_lcg_total_gradient_J_type( )
        );

        std::fill(
            output_leading_configuration_gradient_configurations_J_begin, output_leading_configuration_gradient_configurations_J_end, output_lcg_configurations_J_type( )
        );

        std::fill(
            output_leading_configuration_gradient_configuration_gradients_J_begin, output_leading_configuration_gradient_configuration_gradients_J_end, output_lcg_configuration_gradients_J_type( )
        );

    }

    template<
        unsigned int leading_rows,
        unsigned int size,
        unsigned int dim,
        class total_configuration_iterator,
        class total_configuration_gradient_iterator,
        class configuration_iterator,
        class configuration_gradient_iterator,
        class output_leading_configuration_total_J_iterator,
        class output_leading_configuration_configurations_J_iterator,
        class output_leading_configuration_configuration_gradients_J_iterator,
        class output_leading_configuration_gradient_total_J_iterator,
        class output_leading_configuration_gradient_total_gradient_J_iterator,
        class output_leading_configuration_gradient_configurations_J_iterator,
        class output_leading_configuration_gradient_configuration_gradients_J_iterator
    >
    void DeformationBase::solveForAllLeadingJacobians(
        const total_configuration_iterator &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
        const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
        output_leading_configuration_total_J_iterator output_leading_configuration_total_J_begin, output_leading_configuration_total_J_iterator output_leading_configuration_total_J_end,
        output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_begin, output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_end,
        output_leading_configuration_configuration_gradients_J_iterator output_leading_configuration_configuration_gradients_J_begin, output_leading_configuration_configuration_gradients_J_iterator output_leading_configuration_configuration_gradients_J_end,
        output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_begin, output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_end,
        output_leading_configuration_gradient_total_gradient_J_iterator output_leading_configuration_gradient_total_gradient_J_begin, output_leading_configuration_gradient_total_gradient_J_iterator output_leading_configuration_gradient_total_gradient_J_end,
        output_leading_configuration_gradient_configurations_J_iterator output_leading_configuration_gradient_configurations_J_begin, output_leading_configuration_gradient_configurations_J_iterator output_leading_configuration_gradient_configurations_J_end,
        output_leading_configuration_gradient_configuration_gradients_J_iterator output_leading_configuration_gradient_configuration_gradients_J_begin, output_leading_configuration_gradient_configuration_gradients_J_iterator output_leading_configuration_gradient_configuration_gradients_J_end
    ){
        /*!
         * Solve for all of the Jacobians of the leading configuration and its gradient which would be required to achieve the total deformation i.e., if
         * the total deformation is \f$ [A] \f$ and we know the net deformation from the subsequent deformations
         * in the form of the configurations, then
         * 
         * \f$ [A] = [B] [A^{-}] \rightarrow [B] = [A] [A^{-}]^{-1} \f$
         *
         * \f$ \frac{\partial [B]}{\partial X} = \left(\frac{\partial [A]}{\partial X} - [B] \frac{ \partial A^{-}}{\partial X}\right)[A^{-}]^{-1} \f$
         *
         * \param &total_configuration_begin: The starting iterator of the total deformation
         * \param &total_configuration_end: The stopping iterator of the total deformation
         * \param &total_configuration_gradient_begin: The starting iterator of the total deformation gradient
         * \param &total_configuration_gradient_end: The stopping iterator of the total deformation gradient
         * \param &configurations_begin: The starting iterator of the configurations
         * \param &configurations_end: The stopping iterator of the configurations
         * \param &configuration_gradients_end: The stopping iterator of the configuration gradients
         * \param &configuration_gradients_begin: The starting iterator of the configuration gradients
         * \param output_leading_configuration_total_J_begin: The starting iterator of the Jacobian of the leading configuration with respect to the total configuration output
         * \param output_leading_configuration_total_J_end: The stopping iterator of the Jacobian of the leading configuration with respect to the total configuration output
         * \param output_leading_configuration_configurations_J_begin: The starting iterator of the Jacobian of the leading configuration with respect to the configurations output
         * \param output_leading_configuration_configurations_J_end: The stopping iterator of the Jacobian of the leading configuration with respect to the configurations output
         * \param output_leading_configuration_configuration_gradients_J_begin: The starting iterator of the Jacobian of the leading configuration with respect to the configuration gradients output
         * \param output_leading_configuration_configuration_gradients_J_end: The stopping iterator of the Jacobian of the leading configuration with respect to the configuration gradients output
         * \param output_leading_configuration_gradient_total_J_begin: The starting iterator of the Jacobian of the leading configuration gradient with respect to the total configuration output
         * \param output_leading_configuration_gradient_total_J_end: The stopping iterator of the Jacobian of the leading configuration gradient with respect to the total configuration output
         * \param output_leading_configuration_gradient_total_gradient_J_begin: The starting iterator of the Jacobian of the leading configuration gradient with respect to the total configuration gradient output
         * \param output_leading_configuration_gradient_total_gradient_J_end: The stopping iterator of the Jacobian of the leading configuration gradient with respect to the total configuration gradient output
         * \param output_leading_configuration_gradient_configurations_J_begin: The starting iterator of the Jacobian of the leading configuration gradient with respect to the configurations output
         * \param output_leading_configuration_gradient_configurations_J_end: The stopping iterator of the Jacobian of the leading configuration gradient with respect to the configurations output
         * \param output_leading_configuration_gradient_configuration_gradients_J_begin: The starting iterator of the Jacobian of the leading configuration gradient with respect to the configuration gradients output
         * \param output_leading_configuration_gradient_configuration_gradients_J_end: The stopping iterator of the Jacobian of the leading configuration gradient with respect to the configuration gradients output
         */

        const unsigned int num_configs = ( configurations_end - configurations_begin ) / ( size * size );

        using configuration_type = typename std::iterator_traits<configuration_iterator>::value_type;
        using configuration_gradient_type = typename std::iterator_traits<configuration_gradient_iterator>::value_type;
        using output_lc_total_J_type = typename std::iterator_traits<output_leading_configuration_total_J_iterator>::value_type;
        using output_lc_configurations_J_type = typename std::iterator_traits<output_leading_configuration_configurations_J_iterator>::value_type;
        using output_lc_configuration_gradients_J_type = typename std::iterator_traits<output_leading_configuration_gradient_configuration_gradients_J_iterator>::value_type;
        using output_lcg_total_J_type = typename std::iterator_traits<output_leading_configuration_total_J_iterator>::value_type;
        using output_lcg_total_gradient_J_type = typename std::iterator_traits<output_leading_configuration_gradient_total_gradient_J_iterator>::value_type;
        using output_lcg_configurations_J_type = typename std::iterator_traits<output_leading_configuration_gradient_configurations_J_iterator>::value_type;
        using output_lcg_configuration_gradients_J_type = typename std::iterator_traits<output_leading_configuration_gradient_configuration_gradients_J_iterator>::value_type;

#ifndef TARDIGRADE_ERROR_TOOLS_OPT
        _sizeCheck_solveForAllLeadingJacobians<leading_rows,size,dim>(
            total_configuration_begin, total_configuration_end,
            total_configuration_gradient_begin, total_configuration_gradient_end,
            configurations_begin, configurations_end,
            configuration_gradients_begin, configuration_gradients_end,
            output_leading_configuration_total_J_begin, output_leading_configuration_total_J_end,
            output_leading_configuration_configurations_J_begin, output_leading_configuration_configurations_J_end,
            output_leading_configuration_configuration_gradients_J_begin, output_leading_configuration_configuration_gradients_J_end,
            output_leading_configuration_gradient_total_J_begin, output_leading_configuration_gradient_total_J_end,
            output_leading_configuration_gradient_total_gradient_J_begin, output_leading_configuration_gradient_total_gradient_J_end,
            output_leading_configuration_gradient_configurations_J_begin, output_leading_configuration_gradient_configurations_J_end,
            output_leading_configuration_gradient_configuration_gradients_J_begin, output_leading_configuration_gradient_configuration_gradients_J_end
        );
#endif
        // Initialize to zero
        _zeroOutputs_solveForAllLeadingJacobians<leading_rows,size,dim>(
            total_configuration_begin, total_configuration_end,
            total_configuration_gradient_begin, total_configuration_gradient_end,
            configurations_begin, configurations_end,
            configuration_gradients_begin, configuration_gradients_end,
            output_leading_configuration_total_J_begin, output_leading_configuration_total_J_end,
            output_leading_configuration_configurations_J_begin, output_leading_configuration_configurations_J_end,
            output_leading_configuration_configuration_gradients_J_begin, output_leading_configuration_configuration_gradients_J_end,
            output_leading_configuration_gradient_total_J_begin, output_leading_configuration_gradient_total_J_end,
            output_leading_configuration_gradient_total_gradient_J_begin, output_leading_configuration_gradient_total_gradient_J_end,
            output_leading_configuration_gradient_configurations_J_begin, output_leading_configuration_gradient_configurations_J_end,
            output_leading_configuration_gradient_configuration_gradients_J_begin, output_leading_configuration_gradient_configuration_gradients_J_end
        );

        if( num_configs == 0 ){

            // In this case the leading configuration is the total configuration
            for ( unsigned int i = 0; i < leading_rows * size; ++i ){

                *( output_leading_configuration_total_J_begin + leading_rows * size * i + i ) += 1;

            }

            for ( unsigned int i = 0; i < leading_rows * size * dim; ++i ){

                *( output_leading_configuration_gradient_total_gradient_J_begin + leading_rows * size * dim * i + i ) += 1;

            }

        }
        else{

            // Assemble the Jacobian Contributions
            std::array< output_lc_total_J_type, leading_rows * size > leading_configuration; //TODO: The type may not always be correct
            std::array< output_lcg_total_J_type, leading_rows * size * dim > leading_configuration_gradient; //TODO: The type may not always be correct
            std::array< configuration_type, size * size > Aminus_inverse;
            std::array< configuration_gradient_type, size * size * dim > dAminusdX;
            std::array< configuration_type, size * size * size * size > dAminusInversedA;

            // Compute the leading configuration and its gradients
            solveForAllLeading<leading_rows,size,dim>(
                total_configuration_begin, total_configuration_end,
                total_configuration_gradient_begin, total_configuration_gradient_end,
                configurations_begin, configurations_end,
                configuration_gradients_begin, configuration_gradients_end,
                std::begin( Aminus_inverse ), std::end( Aminus_inverse ),
                std::begin( dAminusdX ), std::end( dAminusdX ),
                std::begin( leading_configuration ), std::end( leading_configuration ),
                std::begin( leading_configuration_gradient ), std::end( leading_configuration_gradient )
            );

            _assembledAinversedA<size>(
                std::begin( Aminus_inverse ), std::end( Aminus_inverse ),
                std::begin( dAminusInversedA ), std::end( dAminusInversedA )
            );

            // JACOBIANS W.R.T. TOTAL CONFIGURATION
            // Assemble the Jacobians of the leading configuration
            for ( unsigned int i = 0; i < leading_rows; ++i ){
                for ( unsigned int j = 0; j < size; ++j ){
                    for ( unsigned int a = 0; a < size; ++a ){
                        *( output_leading_configuration_total_J_begin + size * leading_rows * size * i + leading_rows * size * j + size * i + a ) += Aminus_inverse[ size * a + j ];
                    }
                }
            }

            // Assemble the Jacobians of the leading configuration gradient
            for ( unsigned int i = 0; i < leading_rows; ++i ){
                for ( unsigned int l = 0; l < size; ++l ){
                    for ( unsigned int a = 0; a < dim; ++a ){
                        for ( unsigned int c = 0; c < size; ++c ){
                            for ( unsigned int kj = 0; kj < size * size; ++kj ){
                                *( output_leading_configuration_gradient_total_J_begin + size * dim * leading_rows * size * i + dim * leading_rows * size * l + leading_rows * size * a + size * i + c )
                                    += dAminusInversedA[ size * size * size * c + size * size * l + kj ] * dAminusdX[ dim * kj + a ];
                            }
                        }
                    }
                }
            }

            // JACOBIANS W.R.T. TOTAL CONFIGURATION GRADIENT
            // Assemble the Jacobians of the leading configuration gradient
            for ( unsigned int i = 0; i < leading_rows; ++i ){
                for ( unsigned int j = 0; j < size; ++j ){
                    for ( unsigned int a = 0; a < dim; ++a ){
                        for ( unsigned int b = 0; b < size; ++b ){
                            *( output_leading_configuration_gradient_total_gradient_J_begin + size * dim * leading_rows * size * dim * i + dim * leading_rows * size * dim * j + leading_rows * size * dim * a + size * dim * i + dim * b + a )
                                += Aminus_inverse[ size * b + j ];
                        }
                    }
                }
            }

            std::array< output_lc_configurations_J_type, leading_rows * size * size * size > intermediate_term1 = { output_lc_configurations_J_type( ) };
            _compute_intermediate_term_1_solveForLeadingConfigurationConfigurationJacobian<leading_rows,size>(
                std::begin( leading_configuration ), std::end( leading_configuration ),
                std::begin( Aminus_inverse ), std::end( Aminus_inverse ),
                std::begin( intermediate_term1 ), std::end( intermediate_term1 )
            );

            std::array< output_lcg_configurations_J_type, size * size * dim * size * size > intermediate_term2;
            std::array< output_lcg_configurations_J_type, leading_rows * size * dim * size * size > intermediate_term3;
            std::array< output_lc_configurations_J_type, size * size * size * size > Aminus_jacobian;
            std::array< output_lcg_configurations_J_type, size * size * dim * size * size > dAminusdX_jacobian;
            for ( unsigned int configuration_index = 0; configuration_index < num_configs; ++configuration_index ){

                // JACOBIANS W.R.T. CONFIGURATIONS
                // Assemble the Jacobians of the leading configuration
                getNetConfigurationJacobian<size>(
                    configurations_begin, configurations_end,
                    configuration_index,
                    std::begin( Aminus_jacobian ), std::end( Aminus_jacobian )
                );

                _denseMatrixMultiplyAccumulate<
                    leading_rows * size,
                    size * size,
                    size * size
                >(
                    std::begin( intermediate_term1 ), std::end( intermediate_term1 ),
                    std::begin( Aminus_jacobian ),   std::end( Aminus_jacobian ),
                    output_leading_configuration_configurations_J_begin, output_leading_configuration_configurations_J_end,
                    configuration_index * size * size, size * size * num_configs
                );

                // Assemble the Jacobians of the leading configuration gradient

//                std::fill(
//                    std::begin( intermediate_term2 ), std::end( intermediate_term2 ), output_lcg_configurations_J_type( )
//                );
//
//                for ( unsigned int i = 0; i < size; ++i ){
//                    for ( unsigned int j = 0; j < size; ++j ){
//                        for ( unsigned int k = 0; k < size; ++k ){
//                            for ( unsigned int acd = 0; acd < dim * size * size; ++acd ){
//                                intermediate_term2[ size * dim * size * size * i + dim * size * size * j + acd ]
//                                    -= dAminusdX_jacobian[ size * dim * size * size * i + dim * size * size * k + acd ] * Aminus_inverse[ size * k + j ];
//                            }
//                        }
//                    }
//                }
//
//                _denseMatrixMultiply<leading_rows,size,size*dim*size*size>(
//                    std::begin( leading_configuration ), std::begin( leading_configuration ),
//                    std::begin( intermediate_term2 ), std::end( intermediate_term2 ),
//                    output_leading_configuration_gradient_configurations_J_begin, output_leading_configuration_gradient_configurations_J_end
//                );
//
//                std::fill(
//                    std::begin( intermediate_term3 ), std::end( intermediate_term3 ), output_lcg_configurations_J_type( )
//                );
//
//                for ( unsigned int i = 0; i < leading_rows; ++i ){
//                    for ( unsigned int j = 0; j < size; ++j ){
//                        for ( unsigned int a = 0; a < dim; ++a ){
//                            for ( unsigned int e = 0; e < size; ++e ){
//                                for ( unsigned int f = 0; f < size; ++f ){
//                                    intermediate_term3[ size * dim * size * size * i + dim * size * size * j + size * size * a + size * e + f ]
//                                        -= leading_configuration_gradient[ size * dim * i + dim * e + a ] * Aminus_inverse[ size * f + j ];
//                                }
//                            }
//                        }
//                    }
//                }
//
//                _denseMatrixMultiplyAccumulate<leading_rows * size * dim, size * size, size * size>(
//                    std::begin( intermediate_term3 ), std::end( intermediate_term3 ),
//                    std::begin( Aminus_jacobian ), std::end( Aminus_jacobian ),
//                    output_leading_configuration_gradient_configurations_J_begin, output_leading_configuration_gradient_configurations_J_end
//                );

            }

        }

    }

}
