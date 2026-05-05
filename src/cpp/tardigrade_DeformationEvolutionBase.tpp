/**
 ******************************************************************************
 * \file tardigrade_DeformationEvolutionBase.tpp
 ******************************************************************************
 * The base class for defining an evolution equation for deformation
 ******************************************************************************
 */

#include "tardigrade_vector_tools.h"
#include "tardigrade_MatrixMap.h"

namespace tardigradeHydra {

    /*!
     * Form the LHS matrix for the deformation evolution
     *
     * \param &dt: The change in time
     * \param &Ltp1_begin: The starting iterator of the velocity gradient
     * \param &Ltp1_end: The stopping iterator of the velocity gradient
     * \param LHS: The deformation LHS
     */
    template<class container, int size>
    template<
    typename dt_type, class Ltp1_iterator
    >
    void DeformationEvolutionBase<container, size>::_formDeformationLHS(const dt_type &dt,
                             const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
                             std::array<typename std::iterator_traits<Ltp1_iterator>::value_type, size * size> &LHS){

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Ltp1_end - Ltp1_begin) == size*size,
                "The size of Ltp1 is " + std::to_string((unsigned int)(Ltp1_end - Ltp1_begin)) + " but it should be " + std::to_string(size*size))

        std::transform(Ltp1_begin, Ltp1_end, std::begin(LHS),
                std::bind(std::multiplies<>(), std::placeholders::_1, -dt * integration_parameter));

        for ( unsigned int i = 0; i < size; ++i ){
            LHS[size*i+i] += 1;
        }

    }

    /*!
     * Form the linear solver for the deformation evolution
     *
     * \param &dt: The change in time
     * \param &Ltp1_begin: The starting iterator of the velocity gradient
     * \param &Ltp1_end: The stopping iterator of the velocity gradient
     * \param &solver: The linear solver
     */
    template<class container, int size>
    template<
    typename dt_type, class Ltp1_iterator, class solver_type
    >
    void DeformationEvolutionBase<container,size>::formDeformationSolver(const dt_type &dt,
                                const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
                                solver_type &solver){

        using Ltp1_type = typename std::iterator_traits<Ltp1_iterator>::value_type;

        std::array<Ltp1_type, size * size> LHS;

        _formDeformationLHS(dt, Ltp1_begin, Ltp1_end, LHS);

        auto _LHS = getFixedSizeMatrixMap<Ltp1_type,size,size>(LHS.data());

        solver = tardigradeVectorTools::solverType<Ltp1_type, size, size>(_LHS);

    }

    /*!
     * Compute the deformation
     *
     * \param &dt: The change in time
     * \param &Lt_begin: The starting iterator of the previous velocity gradient
     * \param &Lt_end: The stopping iterator of the previous velocity gradient
     * \param &Ltp1_begin: The starting iterator of the velocity gradient
     * \param &Ltp1_end: The stopping iterator of the velocity gradient
     * \param &Ft_begin: The starting iterator of the previous deformation
     * \param &Ft_end: The stopping iterator of the previous deformation
     * \param &Ftp1_begin: The starting iterator of the deformation
     * \param &Ftp1_end: The stopping iterator of the deformation
     */
    template<class container, int size>
    template<
    typename dt_type, class Lt_iterator, class Ltp1_iterator, class Ft_iterator, class Ftp1_iterator
    >
    void DeformationEvolutionBase<container, size>::computeDeformation(
        const dt_type &dt,
        const Lt_iterator &Lt_begin, const Lt_iterator &Lt_end,
        const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
        const Ft_iterator &Ft_begin, const Ft_iterator &Ft_end,
        Ftp1_iterator Ftp1_begin, Ftp1_iterator Ftp1_end){

        using Ltp1_type = typename std::iterator_traits<Ltp1_iterator>::value_type;
        using Lt_type   = typename std::iterator_traits<Lt_iterator>::value_type;
        using Ft_type   = typename std::iterator_traits<Ft_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Lt_end - Lt_begin) == size*size,
                "The size of Lt is " + std::to_string((unsigned int)(Lt_end - Lt_begin)) + " but it should be " + std::to_string(size*size))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Ltp1_end - Ltp1_begin) == size*size,
                "The size of Ltp1 is " + std::to_string((unsigned int)(Ltp1_end - Ltp1_begin)) + " but it should be " + std::to_string(size*size))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Ft_end - Ft_begin) == size*size,
                "The size of Ft is " + std::to_string((unsigned int)(Ft_end - Ft_begin)) + " but it should be " + std::to_string(size*size))

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(Ftp1_end - Ftp1_begin) == size*size,
                "The size of Ftp1 is " + std::to_string((unsigned int)(Ftp1_end - Ftp1_begin)) + " but it should be " + std::to_string(size*size))

        //Form the solver
        tardigradeVectorTools::solverType<Ltp1_type, size, size> linearSolver;
        formDeformationSolver(dt, Ltp1_begin, Ltp1_end, linearSolver);

        //Form the right hand side
        std::array<Lt_type, size * size> previous_dF{};
        std::array<Ft_type, size * size> RHS{};

        std::transform(Lt_begin, Lt_end, std::begin(previous_dF),
                std::bind(std::multiplies<>(), std::placeholders::_1, dt * (1 - integration_parameter)));

        for ( unsigned int i = 0; i < size; ++i ){
            previous_dF[size*i+i] += 1;
            for ( unsigned int j = 0; j < size; ++j ){
                for ( unsigned int I = 0; I < size; ++I ){
                    RHS[size*i+I] += previous_dF[dimension * i + j] * (*(Ft_begin + dimension * j + I));
                }
            }
        }

        auto _RHS = getFixedSizeMatrixMap<Ltp1_type,size,size>(RHS.data());
        auto _Ftp1 = getFixedSizeMatrixMap<Ltp1_type,size,size>(Ftp1_begin);

        _Ftp1 = linearSolver.solve(_RHS);

    }

    /*!
     * Compute the derivative of the current deformation with respect to the current velocity gradient
     *
     * \param &dt: The change in time
     * \param &Ltp1_begin: The starting iterator of the velocity gradient
     * \param &Ltp1_end: The stopping iterator of the velocity gradient
     * \param &Ftp1_begin: The starting iterator of the deformation
     * \param &Ftp1_end: The stopping iterator of the deformation
     * \param &dFtp1dLtp1_begin: The starting iterator of the derivative
     * \param &dFtp1dLtp1_end: The stopping iterator of the derivative
     */
    template<class container, int size>
	template<
	typename dt_type, class Ltp1_iterator, class Ftp1_iterator, class dFtp1dLtp1_iterator
	>
    void DeformationEvolutionBase<container, size>::computeDeformation_dFtp1dLtp1(const dt_type &dt,
                                       const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
                                       const Ftp1_iterator &Ftp1_begin, const Ftp1_iterator &Ftp1_end,
                                       dFtp1dLtp1_iterator dFtp1dLtp1_begin, dFtp1dLtp1_iterator dFtp1dLtp1_end){

        using dFtp1dLtp1_type = typename std::iterator_traits<dFtp1dLtp1_iterator>::value_type;
        using Ltp1_type = typename std::iterator_traits<Ltp1_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dFtp1dLtp1_end - dFtp1dLtp1_begin) == size * size * size * size, "The size of dFtp1dLtp1 is " + std::to_string((unsigned int)(dFtp1dLtp1_end - dFtp1dLtp1_begin)) + " but should be " + std::to_string(size * size * size * size));

        std::array<dFtp1dLtp1_type, size * size * size * size> RHS{};

        for ( unsigned int i = 0; i < size; ++i ){
            for ( unsigned int I = 0; I < size; ++I ){
                for ( unsigned int a = 0; a < size; ++a ){

                    RHS[size * size * size * i + size * size * I + size * i + a] += dt * integration_parameter * ( *(Ftp1_begin + size * a + I) );

                }
            }
        }

        //Form the solver
        tardigradeVectorTools::solverType<Ltp1_type, size, size> linearSolver;
        formDeformationSolver(dt, Ltp1_begin, Ltp1_end, linearSolver);

        auto _RHS = getFixedSizeMatrixMap<dFtp1dLtp1_type, size, size * size * size>(RHS.data());
        auto _dFtp1dLtp1 = getFixedSizeMatrixMap<dFtp1dLtp1_type, size, size * size * size>(dFtp1dLtp1_begin);

        _dFtp1dLtp1 = linearSolver.solve(_RHS);

    }

    /*!
     * Compute the derivative of the current deformation with respect to the previous velocity gradient
     *
     * \param &dt: The change in time
     * \param &Ltp1_begin: The starting iterator of the velocity gradient
     * \param &Ltp1_end: The stopping iterator of the velocity gradient
     * \param &Ft_begin: The starting iterator of the previous deformation
     * \param &Ft_end: The stopping iterator of the previous deformation
     * \param &dFtp1dLt_begin: The starting iterator of the derivative
     * \param &dFtp1dLt_end: The stopping iterator of the derivative
     */
    template<class container, int size>
	template<
	typename dt_type, class Ltp1_iterator, class Ft_iterator, class dFtp1dLt_iterator
	>
    void DeformationEvolutionBase<container,size>::computeDeformation_dFtp1dLt(const dt_type &dt,
                                     const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
                                     const Ft_iterator &Ft_begin, const Ft_iterator &Ft_end,
                                     dFtp1dLt_iterator dFtp1dLt_begin, dFtp1dLt_iterator dFtp1dLt_end){

        using dFtp1dLt_type = typename std::iterator_traits<dFtp1dLt_iterator>::value_type;
        using Ltp1_type = typename std::iterator_traits<Ltp1_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dFtp1dLt_end - dFtp1dLt_begin) == size * size * size * size, "The size of dFtp1dLt is " + std::to_string((unsigned int)(dFtp1dLt_end - dFtp1dLt_begin)) + " but should be " + std::to_string(size * size * size * size));

        std::array<dFtp1dLt_type, size * size * size * size> RHS{};

        for ( unsigned int i = 0; i < size; ++i ){
            for ( unsigned int I = 0; I < size; ++I ){
                for ( unsigned int a = 0; a < size; ++a ){

                    RHS[size * size * size * i + size * size * I + size * i + a] += dt * (1. - integration_parameter) * ( *(Ft_begin + size * a + I) );

                }
            }
        }

        //Form the solver
        tardigradeVectorTools::solverType<Ltp1_type, size, size> linearSolver;
        formDeformationSolver(dt, Ltp1_begin, Ltp1_end, linearSolver);

        auto _RHS = getFixedSizeMatrixMap<dFtp1dLt_type, size, size * size * size>(RHS.data());
        auto _dFtp1dLt = getFixedSizeMatrixMap<dFtp1dLt_type, size, size * size * size>(dFtp1dLt_begin);

        _dFtp1dLt = linearSolver.solve(_RHS);

    }

    /*!
     * Compute the derivative of the current deformation with respect to the previous deformation
     *
     * \param &dt: The change in time
     * \param &Ltp1_begin: The starting iterator of the velocity gradient
     * \param &Ltp1_end: The stopping iterator of the velocity gradient
     * \param &Lt_begin: The starting iterator of the previous velocity gradient
     * \param &Lt_end: The stopping iterator of the previous velocity gradient
     * \param &dFtp1dFt_begin: The starting iterator of the derivative
     * \param &dFtp1dFt_end: The stopping iterator of the derivative
     */
    template<class container, int size>
	template<
	typename dt_type, class Ltp1_iterator, class Lt_iterator, class dFtp1dFt_iterator
	>
    void DeformationEvolutionBase<container,size>::computeDeformation_dFtp1dFt(const dt_type &dt,
                                     const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
                                     const Lt_iterator &Lt_begin, const Lt_iterator &Lt_end,
                                     dFtp1dFt_iterator dFtp1dFt_begin, dFtp1dFt_iterator dFtp1dFt_end){

        using dFtp1dFt_type = typename std::iterator_traits<dFtp1dFt_iterator>::value_type;
        using Ltp1_type = typename std::iterator_traits<Ltp1_iterator>::value_type;

        TARDIGRADE_ERROR_TOOLS_CHECK((unsigned int)(dFtp1dFt_end - dFtp1dFt_begin) == size * size * size * size, "The size of dFtp1dFt is " + std::to_string((unsigned int)(dFtp1dFt_end - dFtp1dFt_begin)) + " but should be " + std::to_string(size * size * size * size));

        std::array<dFtp1dFt_type, size * size * size * size> RHS{};

        for ( unsigned int i = 0; i < size; ++i ){
            for ( unsigned int I = 0; I < size; ++I ){
                RHS[size * size * size * i + size * size * I + size * i + I] += 1;
                for ( unsigned int a = 0; a < size; ++a ){

                    RHS[size * size * size * i + size * size * I + size * a + I] += dt * (1. - integration_parameter) * ( *(Lt_begin + size * i + a) );

                }
            }
        }

        //Form the solver
        tardigradeVectorTools::solverType<Ltp1_type, size, size> linearSolver;
        formDeformationSolver(dt, Ltp1_begin, Ltp1_end, linearSolver);

        auto _RHS = getFixedSizeMatrixMap<dFtp1dFt_type, size, size * size * size>(RHS.data());
        auto _dFtp1dFt = getFixedSizeMatrixMap<dFtp1dFt_type, size, size * size * size>(dFtp1dFt_begin);

        _dFtp1dFt = linearSolver.solve(_RHS);

    }

}
