/**
 ******************************************************************************
 * \file tardigrade_MooneyRivlinStrainEnergy.cpp
 ******************************************************************************
 * A Mooney-Rivlin strain energy function
 ******************************************************************************
 */

#include "tardigrade_MooneyRivlinStrainEnergy.h"

namespace tardigradeHydra {

    /*!
     * Set the strain energy
     *
     * \param isPrevious: A flag for whether to set the current (false) or previous (true) strain energy
     */
    void MooneyRivlinStrainEnergy::setStrainEnergy(const bool isPrevious) {

        NeoHookianStrainEnergy::setStrainEnergy(isPrevious);

        const floatType  Ibar2 = compute_Ibar2<floatType>(isPrevious);

        SetDataStorageBase<floatType> strainEnergy;

        if (isPrevious) {
            strainEnergy = get_SetDataStorage_previousStrainEnergy();
        } else {
            strainEnergy = get_SetDataStorage_strainEnergy();
        }

        *strainEnergy.value += _C01 * (Ibar2 - 3);

    }

    /*!
     * Set the Jacobians of the strain energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) Jacobians of the strain energy
     */
    void MooneyRivlinStrainEnergy::setStrainEnergyJacobians(const bool isPrevious) {

        NeoHookianStrainEnergy::setStrainEnergyJacobians(isPrevious);

        constexpr unsigned int dim = 3;  // TODO: Replace with value from ResidualBase

        floatVector dIbar2dFe(dim * dim, 0);
        compute_dIbar2dFe(isPrevious, std::begin(dIbar2dFe), std::end(dIbar2dFe));

        SetDataStorageBase<secondOrderTensor> dStrainEnergydFe;

        if (isPrevious) {
            dStrainEnergydFe = get_SetDataStorage_dPreviousStrainEnergydPreviousFe();
        } else {
            dStrainEnergydFe = get_SetDataStorage_dStrainEnergydFe();
        }

        for (unsigned int iI = 0; iI < dim * dim; ++iI) {
            (*dStrainEnergydFe.value)[iI] += _C01 * dIbar2dFe[iI];
        }

    }

    /*!
     * Set the Hessians of the strain energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) Hessians of the strain energy
     */
    void MooneyRivlinStrainEnergy::setStrainEnergyHessians(const bool isPrevious) {

        NeoHookianStrainEnergy::setStrainEnergyHessians(isPrevious);

        constexpr unsigned int dim = 3;  // TODO: Replace with value from ResidualBase

        floatVector d2Ibar2dFe2(dim * dim * dim * dim, 0);
        compute_d2Ibar2dFe2(isPrevious, std::begin(d2Ibar2dFe2), std::end(d2Ibar2dFe2));

        SetDataStorageBase<fourthOrderTensor> d2StrainEnergydFe2;

        if (isPrevious) {
            d2StrainEnergydFe2  = get_SetDataStorage_d2PreviousStrainEnergydPreviousFe2();
        } else {
            d2StrainEnergydFe2  = get_SetDataStorage_d2StrainEnergydFe2();
        }

        for (unsigned int iIjJ = 0; iIjJ < dim * dim * dim * dim; ++iIjJ) {
            (*d2StrainEnergydFe2.value)[iIjJ] += _C01 * d2Ibar2dFe2[iIjJ];
        }

    }

}  // namespace tardigradeHydra
