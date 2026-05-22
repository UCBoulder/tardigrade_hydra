/**
 ******************************************************************************
 * \file tardigrade_DeformationEvolutionBase.tpp
 ******************************************************************************
 * The base class for defining an evolution equation for deformation
 ******************************************************************************
 */

#include "tardigrade_MatrixMap.h"
#include "tardigrade_vector_tools.h"

namespace tardigradeHydra {

    /*!
     * Form the LHS matrix for the deformation evolution
     *
     * \param &deltat: The change in time
     * \param &Ltp1_begin: The starting iterator of the velocity gradient
     * \param &Ltp1_end: The stopping iterator of the velocity gradient
     * \param LHS: The deformation LHS
     */
    template <class container, int size>
    template <typename deltat_type, class Ltp1_iterator>
    void DeformationEvolutionBase<container, size>::_formDeformationLHS(
        const deltat_type &deltat, const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
        std::array<typename std::iterator_traits<Ltp1_iterator>::value_type, size * size> &LHS) {
        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Ltp1_end - Ltp1_begin) == size * size,
                                     "The size of Ltp1 is " + std::to_string((unsigned int)(Ltp1_end - Ltp1_begin)) +
                                         " but it should be " + std::to_string(size * size))

        std::transform(Ltp1_begin, Ltp1_end, std::begin(LHS),
                       std::bind(std::multiplies<>(), std::placeholders::_1, -deltat * integration_parameter));

        for (unsigned int i = 0; i < size; ++i) {
            LHS[size * i + i] += 1;
        }
    }

    /*!
     * Form the linear solver for the deformation evolution
     *
     * \param &deltat: The change in time
     * \param &Ltp1_begin: The starting iterator of the velocity gradient
     * \param &Ltp1_end: The stopping iterator of the velocity gradient
     * \param &solver: The linear solver
     */
    template <class container, int size>
    template <typename deltat_type, class Ltp1_iterator, class solver_type>
    void DeformationEvolutionBase<container, size>::formDeformationSolver(const deltat_type   &deltat,
                                                                          const Ltp1_iterator &Ltp1_begin,
                                                                          const Ltp1_iterator &Ltp1_end,
                                                                          solver_type         &solver) {
        using Ltp1_type = typename std::iterator_traits<Ltp1_iterator>::value_type;

        std::array<Ltp1_type, size * size> LHS;

        _formDeformationLHS(deltat, Ltp1_begin, Ltp1_end, LHS);

        auto _LHS = getFixedSizeMatrixMap<Ltp1_type, size, size>(LHS.data());

        solver = tardigradeVectorTools::solverType<Ltp1_type, size, size>(_LHS);
    }

    /*!
     * Compute the deformation
     *
     * \param &deltat: The change in time
     * \param &Lt_begin: The starting iterator of the previous velocity gradient
     * \param &Lt_end: The stopping iterator of the previous velocity gradient
     * \param &Ltp1_begin: The starting iterator of the velocity gradient
     * \param &Ltp1_end: The stopping iterator of the velocity gradient
     * \param &Dt_begin: The starting iterator of the previous deformation
     * \param &Dt_end: The stopping iterator of the previous deformation
     * \param &Dtp1_begin: The starting iterator of the deformation
     * \param &Dtp1_end: The stopping iterator of the deformation
     */
    template <class container, int size>
    template <typename deltat_type, class Lt_iterator, class Ltp1_iterator, class Dt_iterator, class Dtp1_iterator>
    void DeformationEvolutionBase<container, size>::computeDeformation(
        const deltat_type &deltat, const Lt_iterator &Lt_begin, const Lt_iterator &Lt_end,
        const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end, const Dt_iterator &Dt_begin,
        const Dt_iterator &Dt_end, Dtp1_iterator Dtp1_begin, Dtp1_iterator Dtp1_end) {
        using Ltp1_type = typename std::iterator_traits<Ltp1_iterator>::value_type;
        using Lt_type   = typename std::iterator_traits<Lt_iterator>::value_type;
        using Dt_type   = typename std::iterator_traits<Dt_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Lt_end - Lt_begin) == size * size,
                                     "The size of Lt is " + std::to_string((unsigned int)(Lt_end - Lt_begin)) +
                                         " but it should be " + std::to_string(size * size))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Ltp1_end - Ltp1_begin) == size * size,
                                     "The size of Ltp1 is " + std::to_string((unsigned int)(Ltp1_end - Ltp1_begin)) +
                                         " but it should be " + std::to_string(size * size))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Dt_end - Dt_begin) == size * size,
                                     "The size of Dt is " + std::to_string((unsigned int)(Dt_end - Dt_begin)) +
                                         " but it should be " + std::to_string(size * size))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Dtp1_end - Dtp1_begin) == size * size,
                                     "The size of Dtp1 is " + std::to_string((unsigned int)(Dtp1_end - Dtp1_begin)) +
                                         " but it should be " + std::to_string(size * size))

        // Form the solver
        tardigradeVectorTools::solverType<Ltp1_type, size, size> linearSolver;
        formDeformationSolver(deltat, Ltp1_begin, Ltp1_end, linearSolver);

        // Form the right hand side
        std::array<Lt_type, size * size> previous_dF{};
        std::array<Dt_type, size * size> RHS{};

        std::transform(Lt_begin, Lt_end, std::begin(previous_dF),
                       std::bind(std::multiplies<>(), std::placeholders::_1, deltat * (1 - integration_parameter)));

        for (unsigned int i = 0; i < size; ++i) {
            previous_dF[size * i + i] += 1;
            for (unsigned int j = 0; j < size; ++j) {
                for (unsigned int I = 0; I < size; ++I) {
                    RHS[size * i + I] += previous_dF[dimension * i + j] * (*(Dt_begin + dimension * j + I));
                }
            }
        }

        auto _RHS  = getFixedSizeMatrixMap<Ltp1_type, size, size>(RHS.data());
        auto _Dtp1 = getFixedSizeMatrixMap<Ltp1_type, size, size>(Dtp1_begin);

        _Dtp1 = linearSolver.solve(_RHS);
    }

    /*!
     * Compute the derivative of the current deformation with respect to the current velocity gradient
     *
     * \param &deltat: The change in time
     * \param &Ltp1_begin: The starting iterator of the velocity gradient
     * \param &Ltp1_end: The stopping iterator of the velocity gradient
     * \param &Dtp1_begin: The starting iterator of the deformation
     * \param &Dtp1_end: The stopping iterator of the deformation
     * \param &dDtp1dLtp1_begin: The starting iterator of the derivative
     * \param &dDtp1dLtp1_end: The stopping iterator of the derivative
     */
    template <class container, int size>
    template <typename deltat_type, class Ltp1_iterator, class Dtp1_iterator, class dDtp1dLtp1_iterator>
    void DeformationEvolutionBase<container, size>::computeDeformation_dDtp1dLtp1(
        const deltat_type &deltat, const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
        const Dtp1_iterator &Dtp1_begin, const Dtp1_iterator &Dtp1_end, dDtp1dLtp1_iterator dDtp1dLtp1_begin,
        dDtp1dLtp1_iterator dDtp1dLtp1_end) {
        using dDtp1dLtp1_type = typename std::iterator_traits<dDtp1dLtp1_iterator>::value_type;
        using Ltp1_type       = typename std::iterator_traits<Ltp1_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dDtp1dLtp1_end - dDtp1dLtp1_begin) == size * size * size * size,
                                     "The size of dDtp1dLtp1 is " +
                                         std::to_string((unsigned int)(dDtp1dLtp1_end - dDtp1dLtp1_begin)) +
                                         " but should be " + std::to_string(size * size * size * size));

        std::array<dDtp1dLtp1_type, size * size * size * size> RHS{};

        for (unsigned int i = 0; i < size; ++i) {
            for (unsigned int I = 0; I < size; ++I) {
                for (unsigned int a = 0; a < size; ++a) {
                    RHS[size * size * size * i + size * size * I + size * i + a] +=
                        deltat * integration_parameter * (*(Dtp1_begin + size * a + I));
                }
            }
        }

        // Form the solver
        tardigradeVectorTools::solverType<Ltp1_type, size, size> linearSolver;
        formDeformationSolver(deltat, Ltp1_begin, Ltp1_end, linearSolver);

        auto _RHS        = getFixedSizeMatrixMap<dDtp1dLtp1_type, size, size * size * size>(RHS.data());
        auto _dDtp1dLtp1 = getFixedSizeMatrixMap<dDtp1dLtp1_type, size, size * size * size>(dDtp1dLtp1_begin);

        _dDtp1dLtp1 = linearSolver.solve(_RHS);
    }

    /*!
     * Compute the derivative of the current deformation with respect to the previous velocity gradient
     *
     * \param &deltat: The change in time
     * \param &Ltp1_begin: The starting iterator of the velocity gradient
     * \param &Ltp1_end: The stopping iterator of the velocity gradient
     * \param &Dt_begin: The starting iterator of the previous deformation
     * \param &Dt_end: The stopping iterator of the previous deformation
     * \param &dDtp1dLt_begin: The starting iterator of the derivative
     * \param &dDtp1dLt_end: The stopping iterator of the derivative
     */
    template <class container, int size>
    template <typename deltat_type, class Ltp1_iterator, class Dt_iterator, class dDtp1dLt_iterator>
    void DeformationEvolutionBase<container, size>::computeDeformation_dDtp1dLt(
        const deltat_type &deltat, const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
        const Dt_iterator &Dt_begin, const Dt_iterator &Dt_end, dDtp1dLt_iterator dDtp1dLt_begin,
        dDtp1dLt_iterator dDtp1dLt_end) {
        using dDtp1dLt_type = typename std::iterator_traits<dDtp1dLt_iterator>::value_type;
        using Ltp1_type     = typename std::iterator_traits<Ltp1_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dDtp1dLt_end - dDtp1dLt_begin) == size * size * size * size,
                                     "The size of dDtp1dLt is " +
                                         std::to_string((unsigned int)(dDtp1dLt_end - dDtp1dLt_begin)) +
                                         " but should be " + std::to_string(size * size * size * size));

        std::array<dDtp1dLt_type, size * size * size * size> RHS{};

        for (unsigned int i = 0; i < size; ++i) {
            for (unsigned int I = 0; I < size; ++I) {
                for (unsigned int a = 0; a < size; ++a) {
                    RHS[size * size * size * i + size * size * I + size * i + a] +=
                        deltat * (1. - integration_parameter) * (*(Dt_begin + size * a + I));
                }
            }
        }

        // Form the solver
        tardigradeVectorTools::solverType<Ltp1_type, size, size> linearSolver;
        formDeformationSolver(deltat, Ltp1_begin, Ltp1_end, linearSolver);

        auto _RHS      = getFixedSizeMatrixMap<dDtp1dLt_type, size, size * size * size>(RHS.data());
        auto _dDtp1dLt = getFixedSizeMatrixMap<dDtp1dLt_type, size, size * size * size>(dDtp1dLt_begin);

        _dDtp1dLt = linearSolver.solve(_RHS);
    }

    /*!
     * Compute the derivative of the current deformation with respect to the previous deformation
     *
     * \param &deltat: The change in time
     * \param &Ltp1_begin: The starting iterator of the velocity gradient
     * \param &Ltp1_end: The stopping iterator of the velocity gradient
     * \param &Lt_begin: The starting iterator of the previous velocity gradient
     * \param &Lt_end: The stopping iterator of the previous velocity gradient
     * \param &dDtp1dDt_begin: The starting iterator of the derivative
     * \param &dDtp1dDt_end: The stopping iterator of the derivative
     */
    template <class container, int size>
    template <typename deltat_type, class Ltp1_iterator, class Lt_iterator, class dDtp1dDt_iterator>
    void DeformationEvolutionBase<container, size>::computeDeformation_dDtp1dDt(
        const deltat_type &deltat, const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
        const Lt_iterator &Lt_begin, const Lt_iterator &Lt_end, dDtp1dDt_iterator dDtp1dDt_begin,
        dDtp1dDt_iterator dDtp1dDt_end) {
        using dDtp1dDt_type = typename std::iterator_traits<dDtp1dDt_iterator>::value_type;
        using Ltp1_type     = typename std::iterator_traits<Ltp1_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dDtp1dDt_end - dDtp1dDt_begin) == size * size * size * size,
                                     "The size of dDtp1dDt is " +
                                         std::to_string((unsigned int)(dDtp1dDt_end - dDtp1dDt_begin)) +
                                         " but should be " + std::to_string(size * size * size * size));

        std::array<dDtp1dDt_type, size * size * size * size> RHS{};

        for (unsigned int i = 0; i < size; ++i) {
            for (unsigned int I = 0; I < size; ++I) {
                RHS[size * size * size * i + size * size * I + size * i + I] += 1;
                for (unsigned int a = 0; a < size; ++a) {
                    RHS[size * size * size * i + size * size * I + size * a + I] +=
                        deltat * (1. - integration_parameter) * (*(Lt_begin + size * i + a));
                }
            }
        }

        // Form the solver
        tardigradeVectorTools::solverType<Ltp1_type, size, size> linearSolver;
        formDeformationSolver(deltat, Ltp1_begin, Ltp1_end, linearSolver);

        auto _RHS      = getFixedSizeMatrixMap<dDtp1dDt_type, size, size * size * size>(RHS.data());
        auto _dDtp1dDt = getFixedSizeMatrixMap<dDtp1dDt_type, size, size * size * size>(dDtp1dDt_begin);

        _dDtp1dDt = linearSolver.solve(_RHS);
    }

}  // namespace tardigradeHydra
