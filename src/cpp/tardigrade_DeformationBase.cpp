/**
  ******************************************************************************
  * \file tardigrade_DeformationBase.cpp
  ******************************************************************************
  * The base class for defining multiplicatively decomposed deformation
  ******************************************************************************
  */

#include "tardigrade_DeformationBase.h"
#include "tardigrade_error_tools.h"

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
    void DeformationBase::getSubConfiguration(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Construct a sub configuration from the iterator using an assumed multiplicative decomposition
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

            // Get preceeding sub-configuration
            std::array< output_type, size * size > Aminus;
            getSubConfiguration<size>( configurations_begin + size * size, configurations_end, std::begin( Aminus ), std::end( Aminus ) );

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
    void DeformationBase::getLeadingSubConfigurationJacobian(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute the Jacobian of a sub-configuration with respect to the first configuration e.g.,
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

        getSubConfiguration<size>(
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
    void DeformationBase::getTrailingSubConfigurationJacobian(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute the Jacobian of a sub-configuration with respect to the first configuration e.g.,
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

        getSubConfiguration<size>(
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
    void DeformationBase::getSubConfigurationJacobian(
        const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
        const unsigned int &configuration_index, output_iterator output_begin, output_iterator output_end
    ){
        /*!
         * Compute the Jacobian of a sub-configuration with respect to an internal configuration e.g.,
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

            getLeadingSubConfigurationJacobian<size>(
                configurations_begin, configurations_end, output_begin, output_end
            );

        }
        else if ( ( configuration_index + 1 ) == num_configurations ){

            getTrailingSubConfigurationJacobian<size>(
                configurations_begin, configurations_end, output_begin, output_end
            );

        }
        else if ( ( 0 < configuration_index ) && ( configuration_index < ( num_configurations - 1 ) ) ){

            std::fill( output_begin, output_end, output_type( ) );
            
            // Get the prior and previous configurations
            std::array< output_type, size * size > Aplus, Aminus;
            getSubConfiguration<size>(
                configurations_begin, configurations_begin + size * size * configuration_index,
                std::begin( Aplus ), std::end( Aplus )
            );
            getSubConfiguration<size>(
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

}
