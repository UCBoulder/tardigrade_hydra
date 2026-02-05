/**
 ******************************************************************************
 * \file tardigrade_DeformationDecompositionBase.h
 ******************************************************************************
 * The base class for defining multiplicatively decomposed deformation
 ******************************************************************************
 */

#ifndef TARDIGRADE_DEFORMATIONDECOMPOSITIONBASE_H
#define TARDIGRADE_DEFORMATIONDECOMPOSITIONBASE_H

namespace tardigradeHydra {

    /*!
     * Base class for the decomposition of deformation
     */
    template <unsigned int _leading_rows, unsigned int _size, unsigned int _dim>
    class DeformationDecompositionBase {
       public:
        static constexpr unsigned int leading_rows =
            _leading_rows;                           //!< The number of rows in the leading configuration
        static constexpr unsigned int size = _size;  //!< The number of columns in the leading configuration and rows
                                                     //!< and columns for trailing configurations
        static constexpr unsigned int dim  = _dim;   //!< The dimension of the gradient

        DeformationDecompositionBase() {
            /*!
             * The base class for multiplicative deformation decomposition
             *
             * Provides utilities for decomposition which can then be used
             * to create specific approaches (e.g., classical continuum,
             * micromorphic continuum, etc.)
             *
             * The template parameters define the fundamental dimensions of
             * the deformation. If the deformation can be written as
             *
             * \f$ A_{ij} = B_{ik} C_{kj} \f$
             *
             * If the total configuration \f$ A \f$ has a dimension of
             * leading_rows x size then \f$ B \f$ will have a dimension
             * of leading_size x size and \f$ C \f$ will have a dimension
             * of size x size.
             *
             * The gradient can be written as
             *
             * \f$ A_{ij,a} = B_{ik,a} C_{kj} + B_{ik} C_{kj,a} \f$
             *
             * where the index \f$ a \f$ will have a dimension of dim.
             */
        }

        template <class configuration_iterator, class output_iterator>
        void getNetConfiguration(const configuration_iterator &configurations_begin,
                                 const configuration_iterator &configurations_end, output_iterator output_begin,
                                 output_iterator output_end);

        template <class configuration_iterator, class output_iterator>
        void accumulateLeadingNetConfigurationJacobian(const configuration_iterator &configurations_begin,
                                                       const configuration_iterator &configurations_end,
                                                       output_iterator output_begin, output_iterator output_end,
                                                       const unsigned int output_offset=0,
                                                       const unsigned int output_stride=size*size);

        template <class configuration_iterator, class output_iterator>
        void accumulateTrailingNetConfigurationJacobian(const configuration_iterator &configurations_begin,
                                                        const configuration_iterator &configurations_end,
                                                        output_iterator output_begin, output_iterator output_end,
                                                        const unsigned int output_offset=0,
                                                        const unsigned int output_stride=size*size);

        template <class configuration_iterator, class output_iterator>
        void accumulateNetConfigurationJacobian(const configuration_iterator &configurations_begin,
                                                const configuration_iterator &configurations_end,
                                                const unsigned int &configuration_index, output_iterator output_begin,
                                                output_iterator output_end, const unsigned int output_offset=0,
                                                const unsigned int output_stride = size * size);

        template <class configuration_iterator, class output_iterator>
        void getLeadingNetConfigurationJacobian(const configuration_iterator &configurations_begin,
                                                const configuration_iterator &configurations_end,
                                                output_iterator output_begin, output_iterator output_end);

        template <class configuration_iterator, class output_iterator>
        void getTrailingNetConfigurationJacobian(const configuration_iterator &configurations_begin,
                                                 const configuration_iterator &configurations_end,
                                                 output_iterator output_begin, output_iterator output_end);

        template <class configuration_iterator, class output_iterator>
        void getNetConfigurationJacobian(const configuration_iterator &configurations_begin,
                                         const configuration_iterator &configurations_end,
                                         const unsigned int &configuration_index, output_iterator output_begin,
                                         output_iterator output_end);

        template <class configuration_iterator, class output_iterator>
        void getNetConfigurationJacobian(const configuration_iterator &configurations_begin,
                                         const configuration_iterator &configurations_end,
                                         output_iterator output_begin,
                                         output_iterator output_end);

        template <class configuration_iterator, class configuration_gradient_iterator, class Aminus_iterator,
                  class dAminusdX_iterator, class output_iterator>
        void _assemble_dAdX_getNetConfigurationGradient(
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator &configuration_gradients_begin,
            const configuration_gradient_iterator &configuration_gradients_end, const Aminus_iterator &Aminus_begin,
            const Aminus_iterator &Aminus_end, const dAminusdX_iterator &dAminusdX_begin,
            const dAminusdX_iterator &dAminusdX_end, output_iterator output_begin, output_iterator output_end);

        template <class configuration_iterator, class configuration_gradient_iterator, class output_iterator>
        void getNetConfigurationGradient(const configuration_iterator          &configurations_begin,
                                         const configuration_iterator          &configurations_end,
                                         const configuration_gradient_iterator &configuration_gradients_begin,
                                         const configuration_gradient_iterator &configuration_gradients_end,
                                         output_iterator output_begin, output_iterator output_end);

        template <class configuration_iterator, class configuration_gradient_iterator, class output_iterator>
        void getLeadingNetConfigurationGradientConfigurationJacobian(
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator &configuration_gradients_begin,
            const configuration_gradient_iterator &configuration_gradients_end, output_iterator output_begin,
            output_iterator output_end);

        template <class configuration_iterator, class configuration_gradient_iterator, class output_iterator>
        void getTrailingNetConfigurationGradientConfigurationJacobian(
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator &configuration_gradients_begin,
            const configuration_gradient_iterator &configuration_gradients_end, output_iterator output_begin,
            output_iterator output_end);

        template <class configuration_iterator, class configuration_gradient_iterator, class output_iterator>
        void getNetConfigurationGradientConfigurationJacobian(
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator &configuration_gradients_begin,
            const configuration_gradient_iterator &configuration_gradients_end, const unsigned int &configuration_index,
            output_iterator output_begin, output_iterator output_end);

        template <class configuration_iterator, class configuration_gradient_iterator, class output_iterator>
        void getLeadingNetConfigurationGradientConfigurationGradientJacobian(
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator &configuration_gradients_begin,
            const configuration_gradient_iterator &configuration_gradients_end, output_iterator output_begin,
            output_iterator output_end);

        template <class configuration_iterator, class configuration_gradient_iterator, class output_iterator>
        void getTrailingNetConfigurationGradientConfigurationGradientJacobian(
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator &configuration_gradients_begin,
            const configuration_gradient_iterator &configuration_gradients_end, output_iterator output_begin,
            output_iterator output_end);

        template <class configuration_iterator, class configuration_gradient_iterator, class output_iterator>
        void getNetConfigurationGradientConfigurationGradientJacobian(
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator &configuration_gradients_begin,
            const configuration_gradient_iterator &configuration_gradients_end, const unsigned int &configuration_index,
            output_iterator output_begin, output_iterator output_end);

        template <class total_configuration_iterator, class configuration_iterator, class output_iterator>
        void solveForLeadingConfiguration(const total_configuration_iterator &total_configuration_begin,
                                          const total_configuration_iterator &total_configuration_end,
                                          const configuration_iterator       &configurations_begin,
                                          const configuration_iterator       &configurations_end,
                                          output_iterator output_begin, output_iterator output_end);

        template <class total_configuration_iterator, class configuration_iterator, class Aminus_inverse_iterator,
                  class output_iterator>
        void solveForLeadingConfiguration(const total_configuration_iterator &total_configuration_begin,
                                          const total_configuration_iterator &total_configuration_end,
                                          const configuration_iterator       &configurations_begin,
                                          const configuration_iterator       &configurations_end,
                                          Aminus_inverse_iterator             Aminus_inverse_begin,
                                          Aminus_inverse_iterator Aminus_inverse_end, output_iterator output_begin,
                                          output_iterator output_end);

        template <class total_configuration_iterator, class configuration_iterator, class output_iterator>
        void solveForLeadingConfigurationTotalConfigurationJacobian(
            const total_configuration_iterator &total_configuration_begin,
            const total_configuration_iterator &total_configuration_end,
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            output_iterator output_begin, output_iterator output_end);

        template <class total_configuration_iterator, class configuration_iterator, class output_iterator>
        void solveForLeadingConfigurationConfigurationJacobian(
            const total_configuration_iterator &total_configuration_begin,
            const total_configuration_iterator &total_configuration_end,
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const unsigned int &configuration_index, output_iterator output_begin, output_iterator output_end);

        template <class total_configuration_gradient_iterator, class leading_configuration_iterator,
                  class configuration_iterator, class configuration_gradient_iterator, class Aminus_inverse_iterator,
                  class dAminusdX_iterator, class output_iterator>
        void solveForLeadingConfigurationGradient(
            const total_configuration_gradient_iterator &total_configuration_gradient_begin,
            const total_configuration_gradient_iterator &total_configuration_gradient_end,
            const leading_configuration_iterator        &leading_configuration_begin,
            const leading_configuration_iterator        &leading_configuration_end,
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator &configuration_gradients_begin,
            const configuration_gradient_iterator &configuration_gradients_end,
            Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
            dAminusdX_iterator dAminusdX_begin, dAminusdX_iterator dAminusdX_end, output_iterator output_begin,
            output_iterator output_end);

        template <class total_configuration_gradient_iterator, class leading_configuration_iterator,
                  class configuration_iterator, class configuration_gradient_iterator, class output_iterator>
        void solveForLeadingConfigurationGradient(
            const total_configuration_gradient_iterator &total_configuration_gradient_begin,
            const total_configuration_gradient_iterator &total_configuration_gradient_end,
            const leading_configuration_iterator        &leading_configuration_begin,
            const leading_configuration_iterator        &leading_configuration_end,
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator &configuration_gradients_begin,
            const configuration_gradient_iterator &configuration_gradients_end, output_iterator output_begin,
            output_iterator output_end);

        template <class total_configuration_gradient_iterator, class leading_configuration_iterator,
                  class configuration_iterator, class configuration_gradient_iterator, class output_iterator>
        void solveForLeadingConfigurationGradientTotalConfigurationGradientJacobian(
            const total_configuration_gradient_iterator &total_configuration_gradient_begin,
            const total_configuration_gradient_iterator &total_configuration_gradient_end,
            const leading_configuration_iterator        &leading_configuration_begin,
            const leading_configuration_iterator        &leading_configuration_end,
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator &configuration_gradients_begin,
            const configuration_gradient_iterator &configuration_gradients_end, output_iterator output_begin,
            output_iterator output_end);

        template <class total_configuration_gradient_iterator, class leading_configuration_iterator,
                  class configuration_iterator, class configuration_gradient_iterator, class output_iterator>
        void solveForLeadingConfigurationGradientLeadingConfigurationJacobian(
            const total_configuration_gradient_iterator &total_configuration_gradient_begin,
            const total_configuration_gradient_iterator &total_configuration_gradient_end,
            const leading_configuration_iterator        &leading_configuration_begin,
            const leading_configuration_iterator        &leading_configuration_end,
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator &configuration_gradients_begin,
            const configuration_gradient_iterator &configuration_gradients_end, output_iterator output_begin,
            output_iterator output_end);

        template <class total_configuration_gradient_iterator, class leading_configuration_iterator,
                  class configuration_iterator, class configuration_gradient_iterator, class output_iterator>
        void solveForLeadingConfigurationGradientConfigurationJacobian(
            const total_configuration_gradient_iterator &total_configuration_gradient_begin,
            const total_configuration_gradient_iterator &total_configuration_gradient_end,
            const leading_configuration_iterator        &leading_configuration_begin,
            const leading_configuration_iterator        &leading_configuration_end,
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator &configuration_gradients_begin,
            const configuration_gradient_iterator &configuration_gradients_end, const unsigned int &configuration_index,
            output_iterator output_begin, output_iterator output_end);

        template <class total_configuration_gradient_iterator, class leading_configuration_iterator,
                  class configuration_iterator, class configuration_gradient_iterator, class output_iterator>
        void solveForLeadingConfigurationGradientConfigurationGradientJacobian(
            const total_configuration_gradient_iterator &total_configuration_gradient_begin,
            const total_configuration_gradient_iterator &total_configuration_gradient_end,
            const leading_configuration_iterator        &leading_configuration_begin,
            const leading_configuration_iterator        &leading_configuration_end,
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator &configuration_gradients_begin,
            const configuration_gradient_iterator &configuration_gradients_end, const unsigned int &configuration_index,
            output_iterator output_begin, output_iterator output_end);

        template <class total_configuration_iterator, class total_configuration_gradient_iterator,
                  class configuration_iterator, class configuration_gradient_iterator, class Aminus_inverse_iterator,
                  class dAminusdX_iterator, class output_leading_configuration_iterator,
                  class output_leading_configuration_gradient_iterator>
        void solveForAllLeading(
            const total_configuration_iterator          &total_configuration_begin,
            const total_configuration_iterator          &total_configuration_end,
            const total_configuration_gradient_iterator &total_configuration_gradient_begin,
            const total_configuration_gradient_iterator &total_configuration_gradient_end,
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator &configuration_gradients_begin,
            const configuration_gradient_iterator &configuration_gradients_end,
            Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
            dAminusdX_iterator dAminusdX_begin, dAminusdX_iterator dAminusdX_end,
            output_leading_configuration_iterator          output_leading_configuration_begin,
            output_leading_configuration_iterator          output_leading_configuration_end,
            output_leading_configuration_gradient_iterator output_leading_configuration_gradient_begin,
            output_leading_configuration_gradient_iterator output_leading_configuration_gradient_end);

        template <class total_configuration_iterator, class total_configuration_gradient_iterator,
                  class configuration_iterator, class configuration_gradient_iterator,
                  class output_leading_configuration_iterator, class output_leading_configuration_gradient_iterator>
        void solveForAllLeading(
            const total_configuration_iterator          &total_configuration_begin,
            const total_configuration_iterator          &total_configuration_end,
            const total_configuration_gradient_iterator &total_configuration_gradient_begin,
            const total_configuration_gradient_iterator &total_configuration_gradient_end,
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator         &configuration_gradients_begin,
            const configuration_gradient_iterator         &configuration_gradients_end,
            output_leading_configuration_iterator          output_leading_configuration_begin,
            output_leading_configuration_iterator          output_leading_configuration_end,
            output_leading_configuration_gradient_iterator output_leading_configuration_gradient_begin,
            output_leading_configuration_gradient_iterator output_leading_configuration_gradient_end);

        template <class total_configuration_iterator, class total_configuration_gradient_iterator,
                  class configuration_iterator, class configuration_gradient_iterator,
                  class output_leading_configuration_total_J_iterator,
                  class output_leading_configuration_configurations_J_iterator,
                  class output_leading_configuration_gradient_total_J_iterator,
                  class output_leading_configuration_gradient_total_gradient_J_iterator,
                  class output_leading_configuration_gradient_configurations_J_iterator,
                  class output_leading_configuration_gradient_configuration_gradients_J_iterator>
        void solveForAllLeadingJacobians(
            const total_configuration_iterator          &total_configuration_begin,
            const total_configuration_iterator          &total_configuration_end,
            const total_configuration_gradient_iterator &total_configuration_gradient_begin,
            const total_configuration_gradient_iterator &total_configuration_gradient_end,
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator                 &configuration_gradients_begin,
            const configuration_gradient_iterator                 &configuration_gradients_end,
            output_leading_configuration_total_J_iterator          output_leading_configuration_total_J_begin,
            output_leading_configuration_total_J_iterator          output_leading_configuration_total_J_end,
            output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_begin,
            output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_end,
            output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_begin,
            output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_end,
            output_leading_configuration_gradient_total_gradient_J_iterator
                output_leading_configuration_gradient_total_gradient_J_begin,
            output_leading_configuration_gradient_total_gradient_J_iterator
                output_leading_configuration_gradient_total_gradient_J_end,
            output_leading_configuration_gradient_configurations_J_iterator
                output_leading_configuration_gradient_configurations_J_begin,
            output_leading_configuration_gradient_configurations_J_iterator
                output_leading_configuration_gradient_configurations_J_end,
            output_leading_configuration_gradient_configuration_gradients_J_iterator
                output_leading_configuration_gradient_configuration_gradients_J_begin,
            output_leading_configuration_gradient_configuration_gradients_J_iterator
                output_leading_configuration_gradient_configuration_gradients_J_end);

       protected:
        template <unsigned int rows, unsigned int inner, unsigned int columns, class A_iterator, class B_iterator,
                  class C_iterator>
        void _denseMatrixMultiply(const A_iterator &A_begin, const A_iterator &A_end, const B_iterator &B_begin,
                                  const B_iterator &B_end, C_iterator C_begin, C_iterator C_end,
                                  const unsigned int A_offset = 0, const unsigned int A_stride = inner,
                                  const unsigned int B_offset = 0, const unsigned int B_stride = columns,
                                  const unsigned int output_offset = 0, const unsigned int output_stride = columns);

        template <unsigned int rows, unsigned int inner, unsigned int columns, class A_iterator, class B_iterator,
                  class C_iterator>
        void _denseMatrixMultiplyAccumulate(const A_iterator &A_begin, const A_iterator &A_end,
                                            const B_iterator &B_begin, const B_iterator &B_end, C_iterator C_begin,
                                            C_iterator C_end, const unsigned int A_offset = 0,
                                            const unsigned int A_stride = inner, const unsigned int B_offset = 0,
                                            const unsigned int B_stride = columns, const unsigned int output_offset = 0,
                                            const unsigned int output_stride = columns);

        template <unsigned int rows, unsigned int inner, unsigned int columns, unsigned int output_rows,
                  unsigned int output_columns, class A_iterator, class B_iterator, class C_iterator>
        void _denseMatrixMultiplyAccumulateReshape(
            const A_iterator &A_begin, const A_iterator &A_end, const B_iterator &B_begin, const B_iterator &B_end,
            C_iterator C_begin, C_iterator C_end, const unsigned int A_offset = 0, const unsigned int A_stride = inner,
            const unsigned int B_offset = 0, const unsigned int B_stride = columns,
            const unsigned int output_offset = 0, const unsigned int output_stride = output_columns);

        template <class A_inverse_iterator, class output_iterator>
        void _assembledAinversedA(const A_inverse_iterator &A_inverse_begin, const A_inverse_iterator &A_inverse_end,
                                  output_iterator output_begin, output_iterator output_end);

        template <unsigned int matrix_size, class A_iterator, class output_iterator>
        void _compute_matrix_inverse(const A_iterator &A_begin, const A_iterator &A_end, output_iterator output_begin,
                                     output_iterator output_end);

        template <class total_configuration_iterator, class Aminus_inverse_iterator, class output_iterator>
        void _assemble_leading_configuration_solveForLeadingConfiguration(
            const total_configuration_iterator &total_configuration_begin,
            const total_configuration_iterator &total_configuration_end, Aminus_inverse_iterator Aminus_inverse_begin,
            Aminus_inverse_iterator Aminus_inverse_end, output_iterator output_begin, output_iterator output_end);

        template <class total_configuration_iterator, class total_configuration_gradient_iterator,
                  class configuration_iterator, class configuration_gradient_iterator,
                  class output_leading_configuration_total_J_iterator,
                  class output_leading_configuration_configurations_J_iterator,
                  class output_leading_configuration_gradient_total_J_iterator,
                  class output_leading_configuration_gradient_total_gradient_J_iterator,
                  class output_leading_configuration_gradient_configurations_J_iterator,
                  class output_leading_configuration_gradient_configuration_gradients_J_iterator>
        void _sizeCheck_solveForAllLeadingJacobians(
            const total_configuration_iterator          &total_configuration_begin,
            const total_configuration_iterator          &total_configuration_end,
            const total_configuration_gradient_iterator &total_configuration_gradient_begin,
            const total_configuration_gradient_iterator &total_configuration_gradient_end,
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator                 &configuration_gradients_begin,
            const configuration_gradient_iterator                 &configuration_gradients_end,
            output_leading_configuration_total_J_iterator          output_leading_configuration_total_J_begin,
            output_leading_configuration_total_J_iterator          output_leading_configuration_total_J_end,
            output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_begin,
            output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_end,
            output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_begin,
            output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_end,
            output_leading_configuration_gradient_total_gradient_J_iterator
                output_leading_configuration_gradient_total_gradient_J_begin,
            output_leading_configuration_gradient_total_gradient_J_iterator
                output_leading_configuration_gradient_total_gradient_J_end,
            output_leading_configuration_gradient_configurations_J_iterator
                output_leading_configuration_gradient_configurations_J_begin,
            output_leading_configuration_gradient_configurations_J_iterator
                output_leading_configuration_gradient_configurations_J_end,
            output_leading_configuration_gradient_configuration_gradients_J_iterator
                output_leading_configuration_gradient_configuration_gradients_J_begin,
            output_leading_configuration_gradient_configuration_gradients_J_iterator
                output_leading_configuration_gradient_configuration_gradients_J_end);

        template <class output_leading_configuration_total_J_iterator,
                  class output_leading_configuration_configurations_J_iterator,
                  class output_leading_configuration_gradient_total_J_iterator,
                  class output_leading_configuration_gradient_total_gradient_J_iterator,
                  class output_leading_configuration_gradient_configurations_J_iterator,
                  class output_leading_configuration_gradient_configuration_gradients_J_iterator>
        void _zeroOutputs_solveForAllLeadingJacobians(
            output_leading_configuration_total_J_iterator          output_leading_configuration_total_J_begin,
            output_leading_configuration_total_J_iterator          output_leading_configuration_total_J_end,
            output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_begin,
            output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_end,
            output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_begin,
            output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_end,
            output_leading_configuration_gradient_total_gradient_J_iterator
                output_leading_configuration_gradient_total_gradient_J_begin,
            output_leading_configuration_gradient_total_gradient_J_iterator
                output_leading_configuration_gradient_total_gradient_J_end,
            output_leading_configuration_gradient_configurations_J_iterator
                output_leading_configuration_gradient_configurations_J_begin,
            output_leading_configuration_gradient_configurations_J_iterator
                output_leading_configuration_gradient_configurations_J_end,
            output_leading_configuration_gradient_configuration_gradients_J_iterator
                output_leading_configuration_gradient_configuration_gradients_J_begin,
            output_leading_configuration_gradient_configuration_gradients_J_iterator
                output_leading_configuration_gradient_configuration_gradients_J_end);

        template <class output_leading_configuration_total_J_iterator,
                  class output_leading_configuration_gradient_total_J_iterator,
                  class output_leading_configuration_gradient_total_gradient_J_iterator>
        void _zeroTotalOutputs_solveForAllLeadingJacobians(
            output_leading_configuration_total_J_iterator          output_leading_configuration_total_J_begin,
            output_leading_configuration_total_J_iterator          output_leading_configuration_total_J_end,
            output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_begin,
            output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_end,
            output_leading_configuration_gradient_total_gradient_J_iterator
                output_leading_configuration_gradient_total_gradient_J_begin,
            output_leading_configuration_gradient_total_gradient_J_iterator
                output_leading_configuration_gradient_total_gradient_J_end);

        template <class output_leading_configuration_configurations_J_iterator,
                  class output_leading_configuration_gradient_configurations_J_iterator,
                  class output_leading_configuration_gradient_configuration_gradients_J_iterator>
        void _zeroConfigurationOutputs_solveForAllLeadingJacobians(
            output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_begin,
            output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_end,
            output_leading_configuration_gradient_configurations_J_iterator
                output_leading_configuration_gradient_configurations_J_begin,
            output_leading_configuration_gradient_configurations_J_iterator
                output_leading_configuration_gradient_configurations_J_end,
            output_leading_configuration_gradient_configuration_gradients_J_iterator
                output_leading_configuration_gradient_configuration_gradients_J_begin,
            output_leading_configuration_gradient_configuration_gradients_J_iterator
                output_leading_configuration_gradient_configuration_gradients_J_end);

        template <class total_configuration_iterator, class total_configuration_gradient_iterator,
                  class Aminus_inverse_iterator, class dAminusdX_iterator, class output_leading_configuration_iterator,
                  class output_leading_configuration_gradient_iterator>
        void _sizeCheck_solveForAllLeading(
            const total_configuration_iterator          &total_configuration_begin,
            const total_configuration_iterator          &total_configuration_end,
            const total_configuration_gradient_iterator &total_configuration_gradient_begin,
            const total_configuration_gradient_iterator &total_configuration_gradient_end,
            Aminus_inverse_iterator Aminus_inverse_begin, Aminus_inverse_iterator Aminus_inverse_end,
            dAminusdX_iterator dAminusdX_begin, dAminusdX_iterator dAminusdX_end,
            output_leading_configuration_iterator          output_leading_configuration_begin,
            output_leading_configuration_iterator          output_leading_configuration_end,
            output_leading_configuration_gradient_iterator output_leading_configuration_gradient_begin,
            output_leading_configuration_gradient_iterator output_leading_configuration_gradient_end);

        template <class leading_configuration_iterator, class Aminus_inverse_iterator, class output_iterator>
        void _compute_intermediate_term_solveForLeadingConfigurationConfigurationJacobian(
            const leading_configuration_iterator &leading_configuration_begin,
            const leading_configuration_iterator &leading_configuration_end,
            const Aminus_inverse_iterator &Aminus_inverse_begin, const Aminus_inverse_iterator &Aminus_inverse_end,
            output_iterator output_begin, output_iterator output_end);

        template <class Aminus_iterator, class output_iterator>
        void _assemble_output_accumulateLeadingNetConfigurationJacobian(const Aminus_iterator &Aminus_begin,
                                                                        const Aminus_iterator &Aminus_end,
                                                                        output_iterator        output_begin,
                                                                        output_iterator        output_end,
                                                                        const unsigned int     output_offset = 0,
                                                                        const unsigned int     output_stride = size * size);

        template <class Aplus_iterator, class output_iterator>
        void _assemble_output_accumulateTrailingNetConfigurationJacobian(const Aplus_iterator &Aplus_begin,
                                                                         const Aplus_iterator &Aplus_end,
                                                                         output_iterator       output_begin,
                                                                         output_iterator       output_end,
                                                                         const unsigned int    output_offset = 0,
                                                                         const unsigned int    output_stride = size * size);

        template <class Aplus_iterator, class Aminus_iterator, class output_iterator>
        void _assemble_output_accumulateNetConfigurationJacobian(const Aplus_iterator  &Aplus_begin,
                                                                 const Aplus_iterator  &Aplus_end,
                                                                 const Aminus_iterator &Aminus_begin,
                                                                 const Aminus_iterator &Aminus_end,
                                                                 output_iterator output_begin, output_iterator output_end,
                                                                 const unsigned int     output_offset = 0,
                                                                 const unsigned int     output_stride = size * size);

        template <class dAminusdX_iterator, class output_iterator>
        void _assemble_output_getLeadingNetConfigurationGradientConfigurationJacobian(
            const dAminusdX_iterator &dAminusdX_begin, const dAminusdX_iterator &dAminusdX_end,
            output_iterator output_begin, output_iterator output_end);

        template <class dAplusdX_iterator, class output_iterator>
        void _assemble_output_getTrailingNetConfigurationGradientConfigurationJacobian(
            const dAplusdX_iterator &dAplusdX_begin, const dAplusdX_iterator &dAplusdX_end,
            output_iterator output_begin, output_iterator output_end);

        template <class Aplus_iterator, class dAplusdX_iterator, class Aminus_jacobian_iterator,
                  class dAminusdX_jacobian_iterator, class output_iterator>
        void _assemble_output_getNetConfigurationGradientConfigurationJacobian(
            const Aplus_iterator &Aplus_begin, const Aplus_iterator &Aplus_end, const dAplusdX_iterator &dAplusdX_begin,
            const dAplusdX_iterator &dAplusdX_end, const Aminus_jacobian_iterator &Aminus_jacobian_begin,
            const Aminus_jacobian_iterator    &Aminus_jacobian_end,
            const dAminusdX_jacobian_iterator &dAminusdX_jacobian_begin,
            const dAminusdX_jacobian_iterator &dAminusdX_jacobian_end, output_iterator output_begin,
            output_iterator output_end);

        template <class Aminus_iterator, class output_iterator>
        void _assemble_output_getLeadingNetConfigurationGradientConfigurationGradientJacobian(
            const Aminus_iterator &Aminus_begin, const Aminus_iterator &Aminus_end, output_iterator output_begin,
            output_iterator output_end);

        template <class Aplus_iterator, class output_iterator>
        void _assemble_output_getTrailingNetConfigurationGradientConfigurationGradientJacobian(
            const Aplus_iterator &Aplus_begin, const Aplus_iterator &Aplus_end, output_iterator output_begin,
            output_iterator output_end);

        template <class Aplus_iterator, class Aminus_iterator, class output_iterator>
        void _assemble_output_getNetConfigurationGradientConfigurationGradientJacobian(
            const Aplus_iterator &Aplus_begin, const Aplus_iterator &Aplus_end, const Aminus_iterator &Aminus_begin,
            const Aminus_iterator &Aminus_end, output_iterator output_begin, output_iterator output_end);

        template <class Aminus_inverse_iterator, class output_iterator>
        void _assemble_output_solveForLeadingConfigurationTotalConfigurationJacobian(
            const Aminus_inverse_iterator &Aminus_inverse_begin, const Aminus_inverse_iterator &Aminus_inverse_end,
            output_iterator output_begin, output_iterator output_end);

        template <class total_configuration_gradient_iterator, class leading_configuration_iterator,
                  class Aminus_inverse_iterator, class dAminusdX_iterator, class output_iterator>
        void _assemble_output_solveForLeadingConfigurationGradient(
            const total_configuration_gradient_iterator &total_configuration_gradient_begin,
            const total_configuration_gradient_iterator &total_configuration_gradient_end,
            const leading_configuration_iterator        &leading_configuration_begin,
            const leading_configuration_iterator        &leading_configuration_end,
            const Aminus_inverse_iterator &Aminus_inverse_begin, const Aminus_inverse_iterator &Aminus_inverse_end,
            const dAminusdX_iterator &dAminusdX_begin, const dAminusdX_iterator &dAminusdX_end,
            output_iterator output_begin, output_iterator output_end);

        template <class Aminus_inverse_iterator, class output_iterator>
        void _assemble_output_solveForLeadingConfigurationGradientTotalConfigurationGradientJacobian(
            const Aminus_inverse_iterator &Aminus_inverse_begin, const Aminus_inverse_iterator &Aminus_inverse_end,
            output_iterator output_begin, output_iterator output_end);

        template <class Aminus_inverse_iterator, class dAminusdX_iterator, class output_iterator>
        void _assemble_output_solveForLeadingConfigurationGradientLeadingConfigurationJacobian(
            const Aminus_inverse_iterator &Aminus_inverse_begin, const Aminus_inverse_iterator &Aminus_inverse_end,
            const dAminusdX_iterator &dAminusdX_begin, const dAminusdX_iterator &dAminusdX_end,
            output_iterator output_begin, output_iterator output_end);

        template <class leading_configuration_gradient_iterator, class Aminus_inverse_iterator, class output_iterator>
        void _assemble_intermediate_term_1_solveForLeadingConfigurationGradientConfigurationJacobian(
            const leading_configuration_gradient_iterator &leading_configuration_gradient_begin,
            const leading_configuration_gradient_iterator &leading_configuration_gradient_end,
            const Aminus_inverse_iterator &Aminus_inverse_begin, const Aminus_inverse_iterator &Aminus_inverse_end,
            output_iterator output_begin, output_iterator output_end);

        template <class leading_configuration_iterator, class Aminus_inverse_iterator, class output_iterator>
        void _assemble_intermediate_term_2_solveForLeadingConfigurationGradientConfigurationJacobian(
            const leading_configuration_iterator &leading_configuration_begin,
            const leading_configuration_iterator &leading_configuration_end,
            const Aminus_inverse_iterator &Aminus_inverse_begin, const Aminus_inverse_iterator &Aminus_inverse_end,
            output_iterator output_begin, output_iterator output_end);

        template <class Aminus_inverse_iterator, class dAminusdX_iterator, class output_iterator>
        void _assemble_leading_configuration_gradient_total_configuration_jacobian_solveforAllLeadingJacobians(
            const Aminus_inverse_iterator &Aminus_inverse_begin, const Aminus_inverse_iterator &Aminus_inverse_end,
            const dAminusdX_iterator &dAminusdX_begin, const dAminusdX_iterator &dAminusdX_end,
            output_iterator output_begin, output_iterator output_end);

        template <class Aminus_inverse_iterator, class dAminusdX_iterator,
                  class output_leading_configuration_total_J_iterator,
                  class output_leading_configuration_gradient_total_J_iterator,
                  class output_leading_configuration_gradient_total_gradient_J_iterator>
        void _assemble_total_jacobians_solveForAllLeadingJacobians(
            const Aminus_inverse_iterator &Aminus_inverse_begin, const Aminus_inverse_iterator &Aminus_inverse_end,
            const dAminusdX_iterator &dAminusdX_begin, const dAminusdX_iterator &dAminusdX_end,
            output_leading_configuration_total_J_iterator          output_leading_configuration_total_J_begin,
            output_leading_configuration_total_J_iterator          output_leading_configuration_total_J_end,
            output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_begin,
            output_leading_configuration_gradient_total_J_iterator output_leading_configuration_gradient_total_J_end,
            output_leading_configuration_gradient_total_gradient_J_iterator
                output_leading_configuration_gradient_total_gradient_J_begin,
            output_leading_configuration_gradient_total_gradient_J_iterator
                output_leading_configuration_gradient_total_gradient_J_end);

        template <class configuration_iterator, class configuration_gradient_iterator,
                  class intermediate_term1_iterator, class intermediate_term2_iterator,
                  class intermediate_term3_iterator, class intermediate_term4_iterator,
                  class output_leading_configuration_configurations_J_iterator,
                  class output_leading_configuration_gradient_configurations_J_iterator,
                  class output_leading_configuration_gradient_configuration_gradients_J_iterator>
        void _assemble_configuration_jacobians_solveForAllLeadingJacobians(
            const configuration_iterator &configurations_begin, const configuration_iterator &configurations_end,
            const configuration_gradient_iterator                 &configuration_gradients_begin,
            const configuration_gradient_iterator                 &configuration_gradients_end,
            const intermediate_term1_iterator                     &intermediate_term1_begin,
            const intermediate_term1_iterator                     &intermediate_term1_end,
            const intermediate_term2_iterator                     &intermediate_term2_begin,
            const intermediate_term2_iterator                     &intermediate_term2_end,
            const intermediate_term3_iterator                     &intermediate_term3_begin,
            const intermediate_term3_iterator                     &intermediate_term3_end,
            const intermediate_term4_iterator                     &intermediate_term4_begin,
            const intermediate_term4_iterator                     &intermediate_term4_end,
            output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_begin,
            output_leading_configuration_configurations_J_iterator output_leading_configuration_configurations_J_end,
            output_leading_configuration_gradient_configurations_J_iterator
                output_leading_configuration_gradient_configurations_J_begin,
            output_leading_configuration_gradient_configurations_J_iterator
                output_leading_configuration_gradient_configurations_J_end,
            output_leading_configuration_gradient_configuration_gradients_J_iterator
                output_leading_configuration_gradient_configuration_gradients_J_begin,
            output_leading_configuration_gradient_configuration_gradients_J_iterator
                output_leading_configuration_gradient_configuration_gradients_J_end);

        template <class intermediate_term2_iterator, class intermediate_term3_iterator,
                  class intermediate_term4_iterator, class leading_configuration_configurations_jacobian_iterator,
                  class Aminus_configuration_jacobian_iterator, class dAminusdX_configuration_jacobian_iterator,
                  class output_leading_configuration_gradient_configurations_jacobian_iterator>
        void _accumulate_output_leading_configuration_gradient_configurations_J_solveForAllLeadingJacobians(
            const unsigned int configuration_index, const unsigned int num_configs,
            const intermediate_term2_iterator &intermediate_term2_begin,
            const intermediate_term2_iterator &intermediate_term2_end,
            const intermediate_term3_iterator &intermediate_term3_begin,
            const intermediate_term3_iterator &intermediate_term3_end,
            const intermediate_term4_iterator &intermediate_term4_begin,
            const intermediate_term4_iterator &intermediate_term4_end,
            const leading_configuration_configurations_jacobian_iterator
                &leading_configuration_configurations_jacobian_begin,
            const leading_configuration_configurations_jacobian_iterator
                                                            &leading_configuration_configurations_jacobian_end,
            const Aminus_configuration_jacobian_iterator    &Aminus_configuration_jacobian_begin,
            const Aminus_configuration_jacobian_iterator    &Aminus_configuration_jacobian_end,
            const dAminusdX_configuration_jacobian_iterator &dAminusdX_configuration_jacobian_begin,
            const dAminusdX_configuration_jacobian_iterator &dAminusdX_configuration_jacobian_end,
            output_leading_configuration_gradient_configurations_jacobian_iterator
                output_leading_configuration_gradient_configurations_jacobian_begin,
            output_leading_configuration_gradient_configurations_jacobian_iterator
                output_leading_configuration_gradient_configurations_jacobian_end);

        template <class leading_configuration_iterator, class leading_configuration_gradient_iterator,
                  class Aminus_inverse_iterator, class dAminusdX_iterator, class output_intermediate_term1_iterator,
                  class output_intermediate_term2_iterator, class output_intermediate_term3_iterator,
                  class output_intermediate_term4_iterator>
        void _compute_intermediate_terms_solveForAllLeadingJacobians(
            const leading_configuration_iterator          &leading_configuration_begin,
            const leading_configuration_iterator          &leading_configuration_end,
            const leading_configuration_gradient_iterator &leading_configuration_gradient_begin,
            const leading_configuration_gradient_iterator &leading_configuration_gradient_end,
            const Aminus_inverse_iterator &Aminus_inverse_begin, const Aminus_inverse_iterator &Aminus_inverse_end,
            const dAminusdX_iterator &dAminusdX_begin, const dAminusdX_iterator &dAminusdX_end,
            output_intermediate_term1_iterator output_intermediate_term1_begin,
            output_intermediate_term1_iterator output_intermediate_term1_end,
            output_intermediate_term2_iterator output_intermediate_term2_begin,
            output_intermediate_term2_iterator output_intermediate_term2_end,
            output_intermediate_term3_iterator output_intermediate_term3_begin,
            output_intermediate_term3_iterator output_intermediate_term3_end,
            output_intermediate_term4_iterator output_intermediate_term4_begin,
            output_intermediate_term4_iterator output_intermediate_term4_end);
    };

}  // namespace tardigradeHydra

#include "tardigrade_DeformationDecompositionBase.cpp"

#endif
