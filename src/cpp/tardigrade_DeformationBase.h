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

            template<
                unsigned int size,
                class configuration_iterator,
                class output_iterator
            >
            void getLeadingSubConfigurationJacobian(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int size,
                class configuration_iterator,
                class output_iterator
            >
            void getTrailingSubConfigurationJacobian(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int size,
                class configuration_iterator,
                class output_iterator
            >
            void getSubConfigurationJacobian(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const unsigned int &configuration_index, output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int size,
                unsigned int dim,
                class configuration_iterator,
                class configuration_gradient_iterator,
                class output_iterator
            >
            void getSubConfigurationGradient(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int size,
                unsigned int dim,
                class configuration_iterator,
                class configuration_gradient_iterator,
                class output_iterator
            >
            void getLeadingSubConfigurationGradientConfigurationJacobian(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int size,
                unsigned int dim,
                class configuration_iterator,
                class configuration_gradient_iterator,
                class output_iterator
            >
            void getTrailingSubConfigurationGradientConfigurationJacobian(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int size,
                unsigned int dim,
                class configuration_iterator,
                class configuration_gradient_iterator,
                class output_iterator
            >
            void getSubConfigurationGradientConfigurationJacobian(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                const unsigned int &configuration_index,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int size,
                unsigned int dim,
                class configuration_iterator,
                class configuration_gradient_iterator,
                class output_iterator
            >
            void getLeadingSubConfigurationGradientConfigurationGradientJacobian(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int size,
                unsigned int dim,
                class configuration_iterator,
                class configuration_gradient_iterator,
                class output_iterator
            >
            void getTrailingSubConfigurationGradientConfigurationGradientJacobian(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int size,
                unsigned int dim,
                class configuration_iterator,
                class configuration_gradient_iterator,
                class output_iterator
            >
            void getSubConfigurationGradientConfigurationGradientJacobian(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                const unsigned int &configuration_index,
                output_iterator output_begin, output_iterator output_end
            );

        protected:

            template<
                unsigned int rows,
                unsigned int columns,
                unsigned int inner,
                class A_iterator, class B_iterator, class C_iterator
            >
            void _denseMatrixMultiply(
                const A_iterator &A_begin, const A_iterator &A_end,
                const B_iterator &B_begin, const B_iterator &B_end,
                C_iterator C_begin, C_iterator C_end
            );

    };

}

#include "tardigrade_DeformationBase.cpp"

#endif
