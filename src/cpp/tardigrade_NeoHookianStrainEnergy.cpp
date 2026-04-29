/**
 ******************************************************************************
 * \file tardigrade_NeoHookianStrainEnergy.cpp
 ******************************************************************************
 * A Neo-Hookian strain energy function
 ******************************************************************************
 */

#include "tardigrade_NeoHookianStrainEnergy.h"

namespace tardigradeHydra{

    /*!
     * Set the strain energy
     *
     * \param isPrevious: A flag for whether to set the current (false) or previous (true) strain energy
     */
    void NeoHookianStrainEnergy::setStrainEnergy(const bool isPrevious) {

        const floatType C10 = _parameters[0];
        const floatType D1  = _parameters[1];

        const floatType Ibar1 = compute_Ibar1<floatType>(isPrevious);
        const floatType *Je;

        SetDataStorageBase<floatType> strainEnergy;

        if (isPrevious){
            Je = get_previousJe();
            strainEnergy = get_SetDataStorage_previousStrainEnergy();
        }
        else{
            Je = get_Je();
            strainEnergy = get_SetDataStorage_strainEnergy();
        }

        *strainEnergy.value = C10 * (Ibar1 - 3) + D1 * (*Je - 1) * (*Je - 1);
    }

    /*!
     * Set the Jacobians of the strain energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) Jacobians of the strain energy
     */
    void NeoHookianStrainEnergy::setStrainEnergyJacobians(const bool isPrevious) {

        constexpr unsigned int dim = 3; //TODO: Replace with value from ResidualBase

        const floatType C10 = _parameters[0];
        const floatType D1  = _parameters[1];

        floatVector dIbar1dFe(dim * dim, 0);
        compute_dIbar1dFe(isPrevious, std::begin(dIbar1dFe), std::end(dIbar1dFe));
        const floatType *Je;
        const secondOrderTensor *dJedFe;

        SetDataStorageBase<secondOrderTensor> dStrainEnergydFe;

        if (isPrevious){
            Je = get_previousJe();
            dJedFe = get_dPreviousJedPreviousFe();
            dStrainEnergydFe = get_SetDataStorage_dPreviousStrainEnergydPreviousFe();
        }
        else{
            Je = get_Je();
            dJedFe = get_dJedFe();
            dStrainEnergydFe = get_SetDataStorage_dStrainEnergydFe();
        }

        dStrainEnergydFe.zero(dim*dim);

        for ( unsigned int iI = 0; iI < dim * dim; ++iI ){
            (*dStrainEnergydFe.value)[iI] += C10 * dIbar1dFe[iI] + 2 * D1 * (*Je - 1) * (*dJedFe)[iI];
        }

    }

    /*!
     * Set the Hessians of the strain energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) Hessians of the strain energy
     */
    void NeoHookianStrainEnergy::setStrainEnergyHessians(const bool isPrevious) {
        constexpr unsigned int dim = 3; //TODO: Replace with value from ResidualBase

        const floatType C10 = _parameters[0];
        const floatType D1  = _parameters[1];

        floatVector d2Ibar1dFe2(dim * dim * dim * dim, 0);
        compute_d2Ibar1dFe2(isPrevious, std::begin(d2Ibar1dFe2), std::end(d2Ibar1dFe2));
        const floatType *Je;
        const secondOrderTensor *dJedFe;
        const fourthOrderTensor *d2JedFe2;

        SetDataStorageBase<fourthOrderTensor> d2StrainEnergydFe2;
        SetDataStorageBase<secondOrderTensor> d2StrainEnergydFedT;

        if (isPrevious){
            Je = get_previousJe();
            dJedFe = get_dPreviousJedPreviousFe();
            d2JedFe2 = get_d2PreviousJedPreviousFe2();
            d2StrainEnergydFe2 = get_SetDataStorage_d2PreviousStrainEnergydPreviousFe2();
            d2StrainEnergydFedT = get_SetDataStorage_d2PreviousStrainEnergydPreviousFedPreviousT();
        }
        else{
            Je = get_Je();
            dJedFe = get_dJedFe();
            d2JedFe2 = get_d2JedFe2();
            d2StrainEnergydFe2 = get_SetDataStorage_d2StrainEnergydFe2();
            d2StrainEnergydFedT = get_SetDataStorage_d2StrainEnergydFedT();
        }

        d2StrainEnergydFe2.zero(dim*dim*dim*dim);
        d2StrainEnergydFedT.zero(dim*dim);

        for ( unsigned int iI = 0; iI < dim * dim; ++iI ){
            for ( unsigned int jJ = 0; jJ < dim * dim; ++jJ ){
                (*d2StrainEnergydFe2.value)[dim*dim*iI+jJ] += C10 * d2Ibar1dFe2[dim*dim*iI+jJ]
                                                            + 2 * D1 * (*dJedFe)[iI] * (*dJedFe)[jJ]
                                                            + 2 * D1 * (*Je - 1) * (*d2JedFe2)[dim*dim*iI+jJ];
            }
        }

    }
}
