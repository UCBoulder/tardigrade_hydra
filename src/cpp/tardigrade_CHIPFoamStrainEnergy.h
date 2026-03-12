/**
 ******************************************************************************
 * \file tardigrade_CHIPFoamStrainEnergy.h
 ******************************************************************************
 * The CHIPFoam strain-energy potential
 ******************************************************************************
 */

#ifndef TARDIGRADE_CHIPFOAMSTRAINENERGY
#define TARDIGRADE_CHIPFOAMSTRAINENERGY

#include "tardigrade_HyperelasticBase.h"

namespace tardigradeHydra {

    /*!
     * The CHIPFoam strain-energy potential
     */
    class CHIPFoamStrainEnergy : public HyperelasticBase {
       public:
        /*!
         * Default residual
         */
        CHIPFoamStrainEnergy() : HyperelasticBase(), _parameters({}) {};

        /*!
         * Main utilization constructor
         *
         * \param *_hydra: A pointer to a hydraBase object
         * \param &_numEquations: The number of equations defined by the residual
         * \param &parameters: The parameter vector organized as
         *    Khat, Ghat, Jb
         */
        CHIPFoamStrainEnergy(hydraBase *_hydra, const unsigned int &_numEquations, const floatVector &parameters)
            : HyperelasticBase(_hydra, _numEquations), _parameters(parameters) {
            TARDIGRADE_ERROR_TOOLS_CHECK(_parameters.size() == 8, "The parameters vector must have a size of 8");

            setInitialized();
        }

        const floatType get_Khat();

        const floatType get_Ghat();

        const floatType get_Jb();

        const floatType get_C10();

        const floatType get_phi0();

        const floatType get_K();

        const floatType get_p0();

        const floatType get_gamma();

        const floatType compute_f(const floatType &J);

        const floatType compute_dfdJ(const floatType &J);

        const floatType compute_d2fdJ2(const floatType &J);

        const floatType compute_Jg(const floatType &Jbar, const floatType &Je);

        const floatType compute_dJgdJbar(const floatType &Jbar, const floatType &Je);

        const floatType compute_dJgdJe(const floatType &Jbar, const floatType &Je);

        const floatType compute_d2JgdJedJbar(const floatType &Jbar, const floatType &Je);

        const floatType compute_d2JgdJbar2(const floatType &Jbar, const floatType &Je);

        const floatType compute_pg(const floatType &Jbar, const floatType &Je);

        const floatType compute_dpgdJbar(const floatType &Jbar, const floatType &Je);

        const floatType compute_dpgdJe(const floatType &Jbar, const floatType &Je);

        const floatType compute_d2pgdJe2(const floatType &Jbar, const floatType &Je);

        const floatType compute_d2pgdJedJbar(const floatType &Jbar, const floatType &Je);

        const floatType compute_d2pgdJbar2(const floatType &Jbar, const floatType &Je);

        const floatType compute_ptilde(const floatType &Jbar, const floatType &Je);

        const floatType compute_dptildedJbar(const floatType &Jbar, const floatType &Je);

        const floatType compute_dptildedJe(const floatType &Jbar, const floatType &Je);

        const floatType compute_d2ptildedJe2(const floatType &Jbar, const floatType &Je);

        const floatType compute_d2ptildedJedJbar(const floatType &Jbar, const floatType &Je);

        const floatType compute_d2ptildedJbar2(const floatType &Jbar, const floatType &Je);

        const floatType compute_Jm(const floatType &Jbar, const floatType &Je);

        const floatType compute_dJmdJbar(const floatType &Jbar, const floatType &Je);

        const floatType compute_dJmdJe(const floatType &Jbar, const floatType &Je);

        const floatType compute_d2JmdJe2(const floatType &Jbar, const floatType &Je);

        const floatType compute_d2JmdJedJbar(const floatType &Jbar, const floatType &Je);

        const floatType compute_d2JmdJbar2(const floatType &Jbar, const floatType &Je);

        const floatType compute_Jbar_residual(const floatType &Jbar, const floatType &Je);

        const floatType compute_Jbar_dRdJbar(const floatType &Jbar, const floatType &Je);

        const floatType compute_Jbar_dRdJe(const floatType &Jbar, const floatType &Je);

        const floatType compute_Jbar_d2RdJe2(const floatType &Jbar, const floatType &Je);

        const floatType compute_Jbar_d2RdJedJbar(const floatType &Jbar, const floatType &Je);

        const floatType compute_Jbar_d2RdJbar2(const floatType &Jbar, const floatType &Je);

        const floatType Jbar_bisection(const floatType &Je, const floatType &lb, const floatType &ub, floatType tol_R = -1, floatType tol_dx = -1);

        const floatType Jbar_newton(const floatType &Je);

        //! The bisection method's relative tolerance
        floatType bisection_tolr = 1e-3;

        //! The bisection method's absolute tolerance
        floatType bisection_tola = 1e-3;

        //! The Newton method's relative tolerance
        floatType newton_tolr = 1e-9;

        //! The Newton method's absolute tolerance
        floatType newton_tola = 1e-9;

        //! The Newton method's maximum number of iterations
        floatType newton_maxiter = 10;

        //! The Newton method's line search alpha parameter
        floatType newton_lsalpha = 1e-4;

        //! The Newton method's maximum number of line-search iterations
        floatType newton_maxlsiter = 5;

       protected:
        //! The model parameters
        floatVector _parameters;

        virtual void setJe(bool isPrevious);

        virtual void setJe();

        virtual void setPreviousJe();

        virtual void setdJedFe(bool isPrevious);

        virtual void setdJedFe();

        virtual void setdPreviousJedPreviousFe();

        virtual void setd2JedFe2(bool isPrevious);

        virtual void setd2JedFe2();

        virtual void setd2PreviousJedPreviousFe2();

        virtual void setIbar1(bool isPrevious);

        virtual void setIbar1();

        virtual void setPreviousIbar1();

        virtual void setdIbar1dFe(bool isPrevious);

        virtual void setdIbar1dFe();

        virtual void setdPreviousIbar1dPreviousFe();

        virtual void setd2Ibar1dFe2(bool isPrevious);

        virtual void setd2Ibar1dFe2();

        virtual void setd2PreviousIbar1dPreviousFe2();

        virtual void setWLB(bool isPrevious);

        virtual void setWLB();

        virtual void setPreviousWLB();

        virtual void setWLBDerivatives(bool isPrevious);

        virtual void setWLBDerivatives();

        virtual void setPreviousWLBDerivatives();

        virtual void setWLBHessians(bool isPrevious);

        virtual void setWLBHessians();

        virtual void setPreviousWLBHessians();

        virtual void setJbar(bool isPrevious);

        virtual void setJbar();

        virtual void setPreviousJbar();

        virtual void setdJbardJe(bool isPrevious);

        virtual void setdJbardJe();

        virtual void setdPreviousJbardPreviousJe();

        virtual void setd2JbardJe2(bool isPrevious);

        virtual void setd2JbardJe2();

        virtual void setd2PreviousJbardPreviousJe2();

        virtual void setdJbardJe1(bool isPrevious);

        virtual void setdJbardJe1();

        virtual void setPreviousdJbardJe1();

        virtual void setddJbardJe1dJe(bool isPrevious);

        virtual void setddJbardJe1dJe();

        virtual void setdPreviousdJbardJe1dPreviousJe();

        virtual void setd2dJbardJe1dJe2(bool isPrevious);

        virtual void setd2dJbardJe1dJe2();

        virtual void setd2PreviousdJbardJe1dPreviousJe2();

        //! Check if the class has been initialized
        const bool isInitialized() { return is_initialized; }

        //! Set that the class has been initialized
        void setInitialized() { is_initialized = true; };

        //! Get the sign of a floating point number
        template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

       private:
        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, Je, floatType, setJe)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousJe, floatType, setPreviousJe)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dJedFe, secondOrderTensor, setdJedFe)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousJedPreviousFe, secondOrderTensor,
                                                  setdPreviousJedPreviousFe)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, d2JedFe2, fourthOrderTensor, setd2JedFe2)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2PreviousJedPreviousFe2, fourthOrderTensor,
                                                  setd2PreviousJedPreviousFe2)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, Ibar1, floatType, setIbar1)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousIbar1, floatType, setPreviousIbar1)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dIbar1dFe, secondOrderTensor, setdIbar1dFe)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousIbar1dPreviousFe, secondOrderTensor,
                                                  setdPreviousIbar1dPreviousFe)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, d2Ibar1dFe2, fourthOrderTensor, setd2Ibar1dFe2)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2PreviousIbar1dPreviousFe2, fourthOrderTensor,
                                                  setd2PreviousIbar1dPreviousFe2)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, WLB, floatType, setWLB);

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousWLB, floatType, setPreviousWLB);

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dWLBdD, floatVector, setWLBDerivatives);

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousWLBdPreviousD, floatVector,
                                                  setPreviousWLBDerivatives);

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, d2WLBdD2, floatVector, setWLBHessians);

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2PreviousWLBdPreviousD2, floatVector,
                                                  setPreviousWLBHessians);

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, Jbar, floatType, setJbar);

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousJbar, floatType, setPreviousJbar);

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dJbardJe, floatType, setdJbardJe);

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousJbardPreviousJe, floatType, setdPreviousJbardPreviousJe);

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, d2JbardJe2, floatType, setd2JbardJe2);

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2PreviousJbardPreviousJe2, floatType, setd2PreviousJbardPreviousJe2);

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dJbardJe1, floatType, setdJbardJe1);

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousdJbardJe1, floatType, setPreviousdJbardJe1);

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, ddJbardJe1dJe, floatType, setddJbardJe1dJe);

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousdJbardJe1dPreviousJe, floatType, setdPreviousdJbardJe1dPreviousJe);

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, d2dJbardJe1dJe2, floatType, setd2dJbardJe1dJe2);

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2PreviousdJbardJe1dPreviousJe2, floatType, setd2PreviousdJbardJe1dPreviousJe2);

        //! Whether the class has been initialized or not
        bool is_initialized = false;
    };

}  // namespace tardigradeHydra

#include "tardigrade_CHIPFoamStrainEnergy.tpp"

#endif
