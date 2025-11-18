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

            // Initialize the output
            std::fill(
                output_begin, output_end, output_type( )
            );

            // Get preceeding sub-configuration
            std::array< output_type, size * size > Aminus;
            getSubConfiguration<size>( configurations_begin + size * size, configurations_end, std::begin( Aminus ), std::end( Aminus ) );

            // Update the output
            for ( unsigned int i = 0; i < size; ++i ){
                for ( unsigned int j = 0; j < size; ++j ){
                    for ( unsigned int k = 0; k < size; ++k ){
                        *( output_begin + size * i + k ) += *( configurations_begin + size * i + j ) * Aminus[ size * j + k ];
                    }
                }
            }

        }
        else{

            std::copy(
                configurations_begin, configurations_end, output_begin
            );

        }

    }

}
