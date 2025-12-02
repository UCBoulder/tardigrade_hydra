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
            void getNetConfiguration(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int size,
                class configuration_iterator,
                class output_iterator
            >
            void getLeadingNetConfigurationJacobian(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int size,
                class configuration_iterator,
                class output_iterator
            >
            void getTrailingNetConfigurationJacobian(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int size,
                class configuration_iterator,
                class output_iterator
            >
            void getNetConfigurationJacobian(
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
            void getNetConfigurationGradient(
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
            void getLeadingNetConfigurationGradientConfigurationJacobian(
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
            void getTrailingNetConfigurationGradientConfigurationJacobian(
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
            void getNetConfigurationGradientConfigurationJacobian(
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
            void getLeadingNetConfigurationGradientConfigurationGradientJacobian(
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
            void getTrailingNetConfigurationGradientConfigurationGradientJacobian(
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
            void getNetConfigurationGradientConfigurationGradientJacobian(
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                const unsigned int &configuration_index,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int deformation_rows,
                unsigned int size,
                class deformation_iterator,
                class configuration_iterator,
                class output_iterator
            >
            void solveForLeadingConfiguration(
                const deformation_iterator   &deformation_begin, const deformation_iterator &deformation_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int deformation_rows,
                unsigned int size,
                class deformation_iterator,
                class configuration_iterator,
                class Aminus_inverse_iterator,
                class output_iterator
            >
            void solveForLeadingConfiguration(
                const deformation_iterator   &deformation_begin, const deformation_iterator &deformation_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int deformation_rows,
                unsigned int size,
                class deformation_iterator,
                class configuration_iterator,
                class output_iterator
            >
            void solveForLeadingConfigurationDeformationJacobian(
                const deformation_iterator   &deformation_begin, const deformation_iterator &deformation_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int deformation_rows,
                unsigned int size,
                class deformation_iterator,
                class configuration_iterator,
                class output_iterator
            >
            void solveForLeadingConfigurationConfigurationJacobian(
                const deformation_iterator   &deformation_begin, const deformation_iterator &deformation_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const unsigned int &configuration_index,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int deformation_rows,
                unsigned int size,
                unsigned int dim,
                class deformation_gradient_iterator,
                class leading_configuration_iterator,
                class configuration_iterator,
                class configuration_gradient_iterator,
                class Aminus_inverse_iterator,
                class output_iterator
            >
            void solveForLeadingConfigurationGradient(
                const deformation_gradient_iterator &deformation_gradient_begin, const deformation_gradient_iterator &deformation_gradient_end,
                const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
                output_iterator output_begin, output_iterator output_end
            );

        protected:

            template<
                unsigned int rows,
                unsigned int inner,
                unsigned int columns,
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
