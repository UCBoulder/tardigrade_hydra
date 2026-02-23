/**
 ******************************************************************************
 * \file tardigrade_HyperelasticBase.h
 ******************************************************************************
 * The base class for Hyperelastic materials
 ******************************************************************************
 */

#ifndef TARDIGRADE_HYPERELASTICBASE
#define TARDIGRADE_HYPERELASTICBASE

#include "tardigrade_ResidualBase.h"

namespace tardigradeHydra {

    /*!
     * A base class for Hyperelastic materials
     */
    class HyperelasticBase : public ResidualBase<> {

       public:

        using tardigradeHydra::ResidualBase<>::ResidualBase;
       
        protected:

            virtual void setFe(const bool isPrevious);

            virtual void setFe();

            virtual void setPreviousFe();

            virtual void setFeDerivatives(const bool isPrevious);

            virtual void setdFedF();

            virtual void setdFedFn();

            virtual void setPreviousdFedF();

            virtual void setPreviousdFedFn();

        private:

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, Fe, secondOrderTensor, setFe)

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dFedF, fourthOrderTensor, setdFedF)

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dFedFn, floatVector, setdFedFn)

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousFe, secondOrderTensor, setPreviousFe)

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousdFedF, fourthOrderTensor, setPreviousdFedF)

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousdFedFn, floatVector, setPreviousdFedFn)

    };

}

#include "tardigrade_ResidualBase.tpp"

#endif
