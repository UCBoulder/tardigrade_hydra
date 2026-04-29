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

        virtual void setJe(bool isPrevious);

        virtual void setJe();

        virtual void setPreviousJe();

        virtual void setdJedFe(bool isPrevious);

        virtual void setdJedFe();

        virtual void setdPreviousJedPreviousFe();

        virtual void setd2JedFe2(bool isPrevious);

        virtual void setd2JedFe2();

        virtual void setd2PreviousJedPreviousFe2();

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

        template <typename T>
        inline T compute_I1(const bool isPrevious);

        template <class dI1dFe_iter>
        inline void compute_dI1dFe(const bool isPrevious, dI1dFe_iter dI1dFe_begin, dI1dFe_iter dI1dFe_end);

        template <class d2I1dFe2_iter>
        inline void compute_d2I1dFe2(d2I1dFe2_iter d2I1dFe2_begin, d2I1dFe2_iter d2I1dFe2_end);

        template <typename C_iter>
        void compute_right_cauchy_green_deformation_tensor(const bool isPrevious, C_iter C_begin, C_iter C_end);

        template <typename T>
        inline T compute_I2(const bool isPrevious);

        template <class dI2dFe_iter>
        inline void compute_dI2dFe(const bool isPrevious, dI2dFe_iter dI2dFe_begin, dI2dFe_iter dI2dFe_end);

        template <class d2I2dFe2_iter>
        inline void compute_d2I2dFe2(const bool isPrevious, d2I2dFe2_iter d2I2dFe2_begin, d2I2dFe2_iter d2I2dFe2_end);

        template <typename T>
        inline T compute_Ibar1(const bool isPrevious);

        template <class dIbar1dFe_iter>
        inline void compute_dIbar1dFe(const bool isPrevious, dIbar1dFe_iter dIbar1dFe_begin,
                                      dIbar1dFe_iter dIbar1dFe_end);

        template <class d2Ibar1dFe2_iter>
        inline void compute_d2Ibar1dFe2(const bool isPrevious, d2Ibar1dFe2_iter d2Ibar1dFe2_begin,
                                        d2Ibar1dFe2_iter d2Ibar1dFe2_end);

        template <typename T>
        inline T compute_Ibar2(const bool isPrevious);

        template <class dIbar2dFe_iter>
        inline void compute_dIbar2dFe(const bool isPrevious, dIbar2dFe_iter dIbar2dFe_begin,
                                      dIbar2dFe_iter dIbar2dFe_end);

        template <class d2Ibar2dFe2_iter>
        inline void compute_d2Ibar2dFe2(const bool isPrevious, d2Ibar2dFe2_iter d2Ibar2dFe2_begin,
                                        d2Ibar2dFe2_iter d2Ibar2dFe2_end);

       private:
        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, Fe, secondOrderTensor, setFe)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dFedF, fourthOrderTensor, setdFedF)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dFedFn, floatVector, setdFedFn)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousFe, secondOrderTensor, setPreviousFe)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousdFedF, fourthOrderTensor, setPreviousdFedF)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousdFedFn, floatVector, setPreviousdFedFn)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, Je, floatType, setJe)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousJe, floatType, setPreviousJe)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dJedFe, secondOrderTensor, setdJedFe)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousJedPreviousFe, secondOrderTensor,
                                                  setdPreviousJedPreviousFe)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, d2JedFe2, fourthOrderTensor, setd2JedFe2)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2PreviousJedPreviousFe2, fourthOrderTensor,
                                                  setd2PreviousJedPreviousFe2)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, strainEnergy, floatType, setStrainEnergy)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousStrainEnergy, floatType, setPreviousStrainEnergy)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dStrainEnergydFe, secondOrderTensor,
                                                  setStrainEnergyJacobians)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousStrainEnergydPreviousFe, secondOrderTensor,
                                                  setPreviousStrainEnergyJacobians)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousStrainEnergydPreviousT, floatType,
                                                  setPreviousStrainEnergyJacobians)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2StrainEnergydFe2, fourthOrderTensor,
                                                  setStrainEnergyHessians)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2StrainEnergydFedT, secondOrderTensor,
                                                  setStrainEnergyHessians)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2PreviousStrainEnergydPreviousFe2, fourthOrderTensor,
                                                  setPreviousStrainEnergyHessians)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, d2PreviousStrainEnergydPreviousFedPreviousT,
                                                  secondOrderTensor, setPreviousStrainEnergyHessians)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, cauchyStress, secondOrderTensor, setCauchyStress)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousCauchyStress, secondOrderTensor,
                                                  setPreviousCauchyStress)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dCauchyStressdF, fourthOrderTensor, setdCauchyStressdF)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dCauchyStressdT, secondOrderTensor, setdCauchyStressdT)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dCauchyStressdFn, floatVector, setdCauchyStressdFn)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousCauchyStressdPreviousF, fourthOrderTensor,
                                                  setdPreviousCauchyStressdPreviousF)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousCauchyStressdPreviousT, secondOrderTensor,
                                                  setdPreviousCauchyStressdPreviousT)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, dPreviousCauchyStressdPreviousFn, floatVector,
                                                  setdPreviousCauchyStressdPreviousFn)
    };

}  // namespace tardigradeHydra

#include "tardigrade_HyperelasticBase.tpp"

#endif
