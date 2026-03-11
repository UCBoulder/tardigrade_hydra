/**
 ******************************************************************************
 * \file tardigrade_CHIPFoamStrainEnergy.tpp
 ******************************************************************************
 * The CHIPFoam strain-energy potential
 ******************************************************************************
 */

namespace tardigradeHydra {

    /*!
     * Compute the isochoric scale function
     */
    const floatType CHIPFoamStrainEnergy::compute_f(const floatType J){

        auto phi0 = get_phi0();

        floatType a = std::pow(phi0/(J-(1-phi0)),1./3);

        return (2*J - 1)*std::pow(J,-1./3) + (2-2*J-phi0) * a;

    }

    /*!
     * Compute the derivative of the isochoric scale function
     */
    const floatType CHIPFoamStrainEnergy::compute_dfdJ(const floatType J){

        auto phi0 = get_phi0();

        floatType a = std::pow(phi0/(J-(1-phi0)),1./3);
        floatType dadJ = -a/(3*(J - (1-phi0)));

        return 2 * std::pow(J,-1./3) - (2*J-1)/3*std::pow(J,-4./3) - 2 * a + (2-2*J-phi0) * dadJ;

    }

    /*!
     * Compute the second derivative of the isochoric scale function
     */
    const floatType CHIPFoamStrainEnergy::compute_d2fdJ2(const floatType J){

        auto phi0 = get_phi0();

        floatType a = std::pow(phi0/(J-(1-phi0)),1./3);
        floatType dadJ = -a/(3*(J - (1-phi0)));
        floatType d2adJ2 = -4 * dadJ/(3*(J - (1-phi0)));

        return -(4./3) * std::pow(J,-4./3) + (4./9) * (2 * J - 1) * std::pow(J,-7./3) - 4 * dadJ + (2 - 2*J-phi0) * d2adJ2;

    }

}

