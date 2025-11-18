/**
  ******************************************************************************
  * \file tardigrade_DeformationBase.h
  ******************************************************************************
  * The base class for defining multiplicatively decomposed deformation
  ******************************************************************************
  */

#ifndef TARDIGRADE_DEFORMATIONBASE_H
#define TARDIGRADE_DEFORMATIONBASE_H

namespace tardigradeHydra{

    /*!
     * Base class for the decomposition of deformation
     */
    class DeformationBase{

        public:

            DeformationBase( ){
                /*!
                 * The base class for multiplicative deformation decomposition
                 *
                 * Provides utilities for decomposition which can then be used
                 * to create specific approaches (e.g., classical continuum,
                 * micromorphic continuum, etc.)
                 */
            }

            template<
                unsigned int size,
                class configuration_iterator,
                class output_iterator
            >
            void getSubConfiguration(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                output_iterator output_begin, output_iterator output_end
            );

    };

}

#include "tardigrade_DeformationBase.cpp"

#endif
