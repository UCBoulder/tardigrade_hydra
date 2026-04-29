/**
 ******************************************************************************
 * \file tardigrade_NeoHookianStrainEnergy.cpp
 ******************************************************************************
 * A Neo-Hookian strain energy function
 ******************************************************************************
 */

#include "tardigrade_NeoHookianStrainEnergy.h"

namespace tardigradeHydra {

    /*!
     * Set the strain energy
     *
     * \param isPrevious: A flag for whether to set the current (false) or previous (true) strain energy
     */
    void NeoHookianStrainEnergy::setStrainEnergy(const bool isPrevious) {
        const floatType  Ibar1 = compute_Ibar1<floatType>(isPrevious);
        const floatType *Je;

        SetDataStorageBase<floatType> strainEnergy;

        if (isPrevious) {
            Je           = get_previousJe();
            strainEnergy = get_SetDataStorage_previousStrainEnergy();
        } else {
            Je           = get_Je();
            strainEnergy = get_SetDataStorage_strainEnergy();
        }

        *strainEnergy.value = _C10 * (Ibar1 - 3) + _D1 * (*Je - 1) * (*Je - 1);
    }

    /*!
     * Set the Jacobians of the strain energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) Jacobians of the strain energy
     */
    void NeoHookianStrainEnergy::setStrainEnergyJacobians(const bool isPrevious) {

        std::array<floatType, dimension * dimension> dIbar1dFe{};
        compute_dIbar1dFe(isPrevious, std::begin(dIbar1dFe), std::end(dIbar1dFe));
        const floatType         *Je;
        const secondOrderTensor *dJedFe;

        SetDataStorageBase<secondOrderTensor> dStrainEnergydFe;

        if (isPrevious) {
            Je               = get_previousJe();
            dJedFe           = get_dPreviousJedPreviousFe();
            dStrainEnergydFe = get_SetDataStorage_dPreviousStrainEnergydPreviousFe();
        } else {
            Je               = get_Je();
            dJedFe           = get_dJedFe();
            dStrainEnergydFe = get_SetDataStorage_dStrainEnergydFe();
        }

        dStrainEnergydFe.zero(dimension * dimension);

        for (unsigned int iI = 0; iI < dimension * dimension; ++iI) {
            (*dStrainEnergydFe.value)[iI] += _C10 * dIbar1dFe[iI] + 2 * _D1 * (*Je - 1) * (*dJedFe)[iI];
        }
    }

    /*!
     * Set the Hessians of the strain energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) Hessians of the strain energy
     */
    void NeoHookianStrainEnergy::setStrainEnergyHessians(const bool isPrevious) {

        std::array<floatType, dimension * dimension * dimension * dimension> d2Ibar1dFe2{};
        compute_d2Ibar1dFe2(isPrevious, std::begin(d2Ibar1dFe2), std::end(d2Ibar1dFe2));
        const floatType         *Je;
        const secondOrderTensor *dJedFe;
        const fourthOrderTensor *d2JedFe2;

        SetDataStorageBase<fourthOrderTensor> d2StrainEnergydFe2;
        SetDataStorageBase<secondOrderTensor> d2StrainEnergydFedT;

        if (isPrevious) {
            Je                  = get_previousJe();
            dJedFe              = get_dPreviousJedPreviousFe();
            d2JedFe2            = get_d2PreviousJedPreviousFe2();
            d2StrainEnergydFe2  = get_SetDataStorage_d2PreviousStrainEnergydPreviousFe2();
            d2StrainEnergydFedT = get_SetDataStorage_d2PreviousStrainEnergydPreviousFedPreviousT();
        } else {
            Je                  = get_Je();
            dJedFe              = get_dJedFe();
            d2JedFe2            = get_d2JedFe2();
            d2StrainEnergydFe2  = get_SetDataStorage_d2StrainEnergydFe2();
            d2StrainEnergydFedT = get_SetDataStorage_d2StrainEnergydFedT();
        }

        d2StrainEnergydFe2.zero(dimension * dimension * dimension * dimension);
        d2StrainEnergydFedT.zero(dimension * dimension);

        for (unsigned int iI = 0; iI < dimension * dimension; ++iI) {
            for (unsigned int jJ = 0; jJ < dimension * dimension; ++jJ) {
                (*d2StrainEnergydFe2.value)[dimension * dimension * iI + jJ] +=
                    _C10 * d2Ibar1dFe2[dimension * dimension * iI + jJ] + 2 * _D1 * (*dJedFe)[iI] * (*dJedFe)[jJ] +
                    2 * _D1 * (*Je - 1) * (*d2JedFe2)[dimension * dimension * iI + jJ];
            }
        }
    }
}  // namespace tardigradeHydra
