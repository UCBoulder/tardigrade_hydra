/**
 ******************************************************************************
 * \file tardigrade_DeformationEvolutionBase.h
 ******************************************************************************
 * The base class for defining an evolution equation for deformation
 ******************************************************************************
 */

#ifndef TARDIGRADE_DEFORMATIONEVOLUTIONBASE_H
#define TARDIGRADE_DEFORMATIONEVOLUTIONBASE_H

#include "tardigrade_ResidualBase.h"

namespace tardigradeHydra {

    /*!
     * Base class for defining the evolution of a deformation
     *
     * Assumes an evolution equation of the form \f$ \dot{F}_{iI} = L_{ij} F_{jI} \f$
     * where \f$ F_{iI} \f$ is the deformation gradient i.e.,
     * \f$\frac{\partial x_i}{\partial X_I}$ \f$ and \f$ L_{ij} \f$
     * is the velocity gradient in \f$ F_{iI} \f$'s current configuration i.e.,
     * \f$\frac{\partial \dot{x}_i}{\partial x_j} \f$.
     *
     * The equation is integrated with a generalized
     * trapezoidal rule which can be expressed as
     *
     * \f$ F_{iI}^{t+1} = F_{iI}^t + \Delta t \dot{F}_{iI}^{t+\alpha} \f$
     *
     * \f$ \dot{F}_{iI}^{t+\alpha} = \left(1-\alpha\right) \dot{F}_{iI}^t + \alpha \dot{F}_{iI}^{t+1} \f$.
     *
     * This means
     *
     * \f$ F_{iI}^{t+1} = F_{iI}^t + \Delta t \left(\left(1-\alpha\right) \dot{F}_{iI}^t + \alpha \dot{F}_{iI}^{t+1} \f$\right) \f$
     *
     * \f$ F_{iI}^{t+1} = F_{iI}^t + \Delta t \left( 1 - \alpha \right) L_{ij}^t F_{jI}^t + \Delta t \alpha L_{ij}^{t+1} F_{jI}^{t+1} \f$
     *
     * \f$ \left(\delta_{ij} - \Delta t \alpha L_{ij}^{t+1} \right) F_{jI}^{t+1} = \left(\delta_{ij} + \Delta t \left(1 - \alpha \right) L_{ij}^{t}\right) F_{jI}^t \f$
     *
     * which can be solved for \f$F_{jI}^{t+1}\f$.
     *
     * By default the integration parameter (available as integration_paramter) is
     * 0.5 which enables second order accuracy. If it is set to 0.0, then the
     * evolution equation is explicitly integrated, and 1.0 then the integration
     * is fully implicit.
     *
     * In some cases (e.g., if the rate is very large) it may be helpful to use the
     * fully implicit method. A value of 0.5 can be shown to be marginally stable
     * in that it may oscillate, but the oscillations do not grow.
     */
    template <class container, int size>
    class DeformationEvolutionBase : public ResidualBase<container> {
        public:
	    using tardigradeHydra::ResidualBase<container>::ResidualBase;

        //! The spatial dimension
        using tardigradeHydra::ResidualBase<container>::dimension;

	    double integration_parameter = 0.5;

        template<
        typename dt_type, class Ltp1_iterator
        >
        void _formDeformationLHS(const dt_type &dt,
                                 const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
                                 std::array<typename std::iterator_traits<Ltp1_iterator>::value_type, size * size> &LHS);

        template<
        typename dt_type, class Ltp1_iterator, class solver_type
        >
        void formDeformationSolver(const dt_type &dt,
                                   const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
                                   solver_type &solver);

	    template<
		typename dt_type, class Lt_iterator, class Ltp1_iterator, class Ft_iterator, class Ftp1_iterator
	    >
        void computeDeformation(const dt_type &dt,
                                const Lt_iterator &Lt_begin, const Lt_iterator &Lt_end,
                                const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
                                const Ft_iterator &Ft_begin, const Ft_iterator &Ft_end,
                                Ftp1_iterator Ftp1_begin, Ftp1_iterator Ftp1_end);

	    template<
		typename dt_type, class Ltp1_iterator, class Ftp1_iterator, class dFtp1dLtp1_iterator
	    >
        void computeDeformation_dFtp1dLtp1(const dt_type &dt,
                                           const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
                                           const Ftp1_iterator &Ftp1_begin, const Ftp1_iterator &Ftp1_end,
                                           dFtp1dLtp1_iterator dFtp1dLtp1_begin, dFtp1dLtp1_iterator dFtp1dLtp1_end);

        protected:
	private:
    };

}

#include "tardigrade_DeformationEvolutionBase.tpp"

#endif
