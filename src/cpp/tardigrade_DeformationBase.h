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
                unsigned int leading_rows,
                unsigned int size,
                class total_configuration_iterator,
                class configuration_iterator,
                class output_iterator
            >
            void solveForLeadingConfiguration(
                const total_configuration_iterator   &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int leading_rows,
                unsigned int size,
                class total_configuration_iterator,
                class configuration_iterator,
                class Aminus_inverse_iterator,
                class output_iterator
            >
            void solveForLeadingConfiguration(
                const total_configuration_iterator   &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int leading_rows,
                unsigned int size,
                class total_configuration_iterator,
                class configuration_iterator,
                class output_iterator
            >
            void solveForLeadingConfigurationTotalConfigurationJacobian(
                const total_configuration_iterator   &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                output_iterator output_begin, output_iterator output_end
            );

            template<
                unsigned int leading_rows,
                unsigned int size,
                class total_configuration_iterator,
                class configuration_iterator,
                class output_iterator
            >
            void solveForLeadingConfigurationConfigurationJacobian(
                const total_configuration_iterator   &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const unsigned int &configuration_index,
                output_iterator output_begin, output_iterator output_end
            );

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
            void solveForLeadingConfigurationGradient(
                const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
                const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
                dAminusdX_iterator dAminusdX_begin, dAminusdX_iterator dAminusdX_end,
                output_iterator output_begin, output_iterator output_end
            );

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
            void solveForLeadingConfigurationGradient(
                const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
                const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                output_iterator output_begin, output_iterator output_end
            );

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
            void solveForLeadingConfigurationGradientTotalConfigurationGradientJacobian(
                const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
                const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                output_iterator output_begin, output_iterator output_end
            );

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
            void solveForLeadingConfigurationGradientLeadingConfigurationJacobian(
                const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
                const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                output_iterator output_begin, output_iterator output_end
            );

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
            void solveForLeadingConfigurationGradientConfigurationJacobian(
                const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
                const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                const unsigned int &configuration_index,
                output_iterator output_begin, output_iterator output_end
            );

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
            void solveForLeadingConfigurationGradientConfigurationGradientJacobian(
                const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
                const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                const unsigned int &configuration_index,
                output_iterator output_begin, output_iterator output_end
            );

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
            void solveForAllLeading(
                const total_configuration_iterator &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
                const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
                dAminusdX_iterator dAminusdX_begin, dAminusdX_iterator dAminusdX_end,
                output_leading_configuration_iterator output_leading_configuration_begin, output_leading_configuration_iterator output_leading_configuration_end,
                output_leading_configuration_gradient_iterator output_leading_configuration_gradient_begin, output_leading_configuration_gradient_iterator output_leading_configuration_gradient_end
            );

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
            void solveForAllLeading(
                const total_configuration_iterator &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
                const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
                const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
                const configuration_gradient_iterator &configuration_gradients_begin, const configuration_gradient_iterator &configuration_gradients_end,
                output_leading_configuration_iterator output_leading_configuration_begin, output_leading_configuration_iterator output_leading_configuration_end,
                output_leading_configuration_gradient_iterator output_leading_configuration_gradient_begin, output_leading_configuration_gradient_iterator output_leading_configuration_gradient_end
            );

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
            void solveForAllLeadingJacobians(
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
                C_iterator C_begin, C_iterator C_end,
                const unsigned int offset = 0, const unsigned int stride = columns
            );

            template<
                unsigned int rows,
                unsigned int inner,
                unsigned int columns,
                class A_iterator, class B_iterator, class C_iterator
            >
            void _denseMatrixMultiplyAccumulate(
                const A_iterator &A_begin, const A_iterator &A_end,
                const B_iterator &B_begin, const B_iterator &B_end,
                C_iterator C_begin, C_iterator C_end,
                const unsigned int offset = 0, const unsigned int stride = columns
            );

            template<
                unsigned int size,
                class A_inverse_iterator, class output_iterator
            >
            void _assembledAinversedA(
                const A_inverse_iterator &A_inverse_begin, const A_inverse_iterator &A_inverse_end,
                output_iterator output_begin, output_iterator output_end
            );

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
            void _sizeCheck_solveForAllLeadingJacobians(
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
            );

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
            void _zeroOutputs_solveForAllLeadingJacobians(
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
            );

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
            void _sizeCheck_solveForAllLeading(
                const total_configuration_iterator &total_configuration_begin, const total_configuration_iterator &total_configuration_end,
                const total_configuration_gradient_iterator &total_configuration_gradient_begin, const total_configuration_gradient_iterator &total_configuration_gradient_end,
                Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
                dAminusdX_iterator dAminusdX_begin, dAminusdX_iterator dAminusdX_end,
                output_leading_configuration_iterator output_leading_configuration_begin, output_leading_configuration_iterator output_leading_configuration_end,
                output_leading_configuration_gradient_iterator output_leading_configuration_gradient_begin, output_leading_configuration_gradient_iterator output_leading_configuration_gradient_end
            );

            template<
                unsigned int leading_rows,
                unsigned int size,
                class leading_configuration_iterator,
                class Aminus_inverse_iterator,
                class output_iterator
            >
            void _compute_intermediate_term_1_solveForLeadingConfigurationConfigurationJacobian(
                const leading_configuration_iterator &leading_configuration_begin, const leading_configuration_iterator &leading_configuration_end,
                const Aminus_inverse_iterator &Aminus_inverse_begin, const Aminus_inverse_iterator &Aminus_inverse_end,
                output_iterator output_begin, output_iterator output_end
            );

    };

}

#include "tardigrade_DeformationBase.cpp"

#endif
