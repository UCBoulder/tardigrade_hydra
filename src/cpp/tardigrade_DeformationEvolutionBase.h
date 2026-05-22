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
     * Assumes an evolution equation of the form \f$ \dot{D}_{iI} = L_{ij} D_{jI} \f$
     * where \f$ D_{iI} \f$ is the deformation gradient i.e.,
     * \f$\frac{\partial x_i}{\partial X_I}$ \f$ and \f$ L_{ij} \f$
     * is the velocity gradient in \f$ D_{iI} \f$'s current configuration e.g.,
     * \f$\frac{\partial \dot{x}_i}{\partial x_j} \f$.
     *
     * The equation is integrated with a generalized
     * trapezoidal rule which can be expressed as
     *
     * \f$ D_{iI}^{t+1} = D_{iI}^t + \Delta t \dot{D}_{iI}^{t+\alpha} \f$
     *
     * \f$ \dot{D}_{iI}^{t+\alpha} = \left(1-\alpha\right) \dot{D}_{iI}^t + \alpha \dot{D}_{iI}^{t+1} \f$.
     *
     * This means
     *
     * \f$ D_{iI}^{t+1} = D_{iI}^t + \Delta t \left(\left(1-\alpha\right) \dot{D}_{iI}^t + \alpha \dot{D}_{iI}^{t+1} \right) \f$
     *
     * \f$ D_{iI}^{t+1} = D_{iI}^t + \Delta t \left( 1 - \alpha \right) L_{ij}^t D_{jI}^t + \Delta t \alpha L_{ij}^{t+1}
     * D_{jI}^{t+1} \f$
     *
     * \f$ \left(\delta_{ij} - \Delta t \alpha L_{ij}^{t+1} \right) D_{jI}^{t+1} = \left(\delta_{ij} + \Delta t \left(1
     * - \alpha \right) L_{ij}^{t}\right) D_{jI}^t \f$
     *
     * which can be solved for \f$D_{jI}^{t+1}\f$.
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

        //! The integration parameter where 0 is fully explicit and 1 is fully implicit
        double integration_parameter = 0.5;

        template <typename deltat_type, class Ltp1_iterator>
        void _formDeformationLHS(
            const deltat_type &deltat, const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
            std::array<typename std::iterator_traits<Ltp1_iterator>::value_type, size * size> &LHS);

        template <typename deltat_type, class Ltp1_iterator, class solver_type>
        void formDeformationSolver(const deltat_type &deltat, const Ltp1_iterator &Ltp1_begin,
                                   const Ltp1_iterator &Ltp1_end, solver_type &solver);

        template <typename deltat_type, class Lt_iterator, class Ltp1_iterator, class Dt_iterator, class Dtp1_iterator>
        void computeDeformation(const deltat_type &deltat, const Lt_iterator &Lt_begin, const Lt_iterator &Lt_end,
                                const Ltp1_iterator &Ltp1_begin, const Ltp1_iterator &Ltp1_end,
                                const Dt_iterator &Dt_begin, const Dt_iterator &Dt_end, Dtp1_iterator Dtp1_begin,
                                Dtp1_iterator Dtp1_end);

        template <typename deltat_type, class Ltp1_iterator, class Dtp1_iterator, class dDtp1dLtp1_iterator>
        void computeDeformation_dDtp1dLtp1(const deltat_type &deltat, const Ltp1_iterator &Ltp1_begin,
                                           const Ltp1_iterator &Ltp1_end, const Dtp1_iterator &Dtp1_begin,
                                           const Dtp1_iterator &Dtp1_end, dDtp1dLtp1_iterator dDtp1dLtp1_begin,
                                           dDtp1dLtp1_iterator dDtp1dLtp1_end);

        template <typename deltat_type, class Ltp1_iterator, class Dt_iterator, class dDtp1dLt_iterator>
        void computeDeformation_dDtp1dLt(const deltat_type &deltat, const Ltp1_iterator &Ltp1_begin,
                                         const Ltp1_iterator &Ltp1_end, const Dt_iterator &Dt_begin,
                                         const Dt_iterator &Dt_end, dDtp1dLt_iterator dDtp1dLt_begin,
                                         dDtp1dLt_iterator dDtp1dLt_end);

        template <typename deltat_type, class Ltp1_iterator, class Lt_iterator, class dDtp1dDt_iterator>
        void computeDeformation_dDtp1dDt(const deltat_type &deltat, const Ltp1_iterator &Ltp1_begin,
                                         const Ltp1_iterator &Ltp1_end, const Lt_iterator &Lt_begin,
                                         const Lt_iterator &Lt_end, dDtp1dDt_iterator dDtp1dDt_begin,
                                         dDtp1dDt_iterator dDtp1dDt_end);

       protected:
       private:
    };

}  // namespace tardigradeHydra

#include "tardigrade_DeformationEvolutionBase.tpp"

#endif
