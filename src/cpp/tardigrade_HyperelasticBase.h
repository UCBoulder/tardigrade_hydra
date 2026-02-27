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

            virtual void setStrainEnergy(const bool isPrevious);

            virtual void setStrainEnergy();

            virtual void setPreviousStrainEnergy();

            virtual void setStrainEnergyJacobians(const bool isPrevious);

            virtual void setStrainEnergyJacobians();

            virtual void setPreviousStrainEnergyJacobians();

            virtual void setStrainEnergyHessians(const bool isPrevious);

            virtual void setStrainEnergyHessians();

            virtual void setPreviousStrainEnergyHessians();

            virtual void setCauchyStress(const bool isPrevious);

            virtual void setCauchyStress();

            virtual void setPreviousCauchyStress();

            virtual void setCauchyStressJacobians(const bool isPrevious);

            virtual void setdCauchyStressdF();

            virtual void setdCauchyStressdT();

            virtual void setdCauchyStressdFn();

            virtual void setdPreviousCauchyStressdPreviousF();

            virtual void setdPreviousCauchyStressdPreviousT();

            virtual void setdPreviousCauchyStressdPreviousFn();

            virtual void setStress() override;

            virtual void setPreviousStress() override;

            virtual void setResidual() override;

            virtual void setJacobian() override;

            virtual void setdRdT() override;

            virtual void setdRdF() override;

        private:

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, Fe, secondOrderTensor, setFe)

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dFedF, fourthOrderTensor, setdFedF)

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dFedFn, floatVector, setdFedFn)

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousFe, secondOrderTensor, setPreviousFe)

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousdFedF, fourthOrderTensor, setPreviousdFedF)

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousdFedFn, floatVector, setPreviousdFedFn)

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, strainEnergy, floatType, setStrainEnergy);

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousStrainEnergy, floatType, setPreviousStrainEnergy);

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dStrainEnergydFe, secondOrderTensor, setStrainEnergyJacobians);

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousStrainEnergydPreviousFe, secondOrderTensor, setPreviousStrainEnergyJacobians);

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2StrainEnergydFe2, fourthOrderTensor, setStrainEnergyHessians);

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2StrainEnergydFedT, secondOrderTensor, setStrainEnergyHessians);

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2PreviousStrainEnergydPreviousFe2, fourthOrderTensor, setPreviousStrainEnergyHessians);

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2PreviousStrainEnergydPreviousFedPreviousT, secondOrderTensor, setPreviousStrainEnergyHessians);

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, cauchyStress, secondOrderTensor, setCauchyStress);

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousCauchyStress, secondOrderTensor, setPreviousCauchyStress);

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dCauchyStressdF, fourthOrderTensor, setdCauchyStressdF);

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dCauchyStressdT, secondOrderTensor, setdCauchyStressdT);

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dCauchyStressdFn, floatVector, setdCauchyStressdFn);

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousCauchyStressdPreviousF, fourthOrderTensor, setdPreviousCauchyStressdPreviousF);

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousCauchyStressdPreviousT, secondOrderTensor, setdPreviousCauchyStressdPreviousT);

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousCauchyStressdPreviousFn, floatVector, setdPreviousCauchyStressdPreviousFn);

    };

}

#include "tardigrade_ResidualBase.tpp"

#endif
