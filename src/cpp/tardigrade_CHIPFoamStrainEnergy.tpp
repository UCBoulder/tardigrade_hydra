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
    const floatType CHIPFoamStrainEnergy::compute_f(const floatType &J){

        auto phi0 = get_phi0();

        floatType a = std::pow(phi0/(J-(1-phi0)),1./3);

        return (2*J - 1)*std::pow(J,-1./3) + (2-2*J-phi0) * a;

    }

    /*!
     * Compute the derivative of the isochoric scale function
     */
    const floatType CHIPFoamStrainEnergy::compute_dfdJ(const floatType &J){

        auto phi0 = get_phi0();

        floatType a = std::pow(phi0/(J-(1-phi0)),1./3);
        floatType dadJ = -a/(3*(J - (1-phi0)));

        return 2 * std::pow(J,-1./3) - (2*J-1)/3*std::pow(J,-4./3) - 2 * a + (2-2*J-phi0) * dadJ;

    }

    /*!
     * Compute the second derivative of the isochoric scale function
     */
    const floatType CHIPFoamStrainEnergy::compute_d2fdJ2(const floatType &J){

        auto phi0 = get_phi0();

        floatType a = std::pow(phi0/(J-(1-phi0)),1./3);
        floatType dadJ = -a/(3*(J - (1-phi0)));
        floatType d2adJ2 = -4 * dadJ/(3*(J - (1-phi0)));

        return -(4./3) * std::pow(J,-4./3) + (4./9) * (2 * J - 1) * std::pow(J,-7./3) - 4 * dadJ + (2 - 2*J-phi0) * d2adJ2;

    }

    /*!
     * Compute the gas relative volume
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_Jg(const floatType &Jbar, const floatType &Je){

        auto phi0 = get_phi0();

        return Je/Jbar * ( Jbar - 1 + phi0)/phi0;

    }

    /*!
     * Compute the derivative of the gas relative volume with respect to
     * the matrix volume-conserving compression
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_dJgdJbar(const floatType &Jbar, const floatType &Je){

        auto phi0 = get_phi0();

        auto Jg = compute_Jg(Jbar,Je);

        return -Jg/Jbar + Je/(Jbar * phi0);

    }

    /*!
     * Compute the derivative of the gas relative volume with respect to
     * the net elastic compression
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_dJgdJe(const floatType &Jbar, const floatType &Je){

        auto phi0 = get_phi0();

        return 1./Jbar * ( Jbar - 1 + phi0)/phi0;
    }

    /*!
     * Compute the second derivative of the gas relative volume with respect to
     * the matrix volume-conserving compression and the net elastic compression
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2JgdJedJbar(const floatType &Jbar, const floatType &Je){

        auto phi0 = get_phi0();

        return -compute_dJgdJe(Jbar,Je)/Jbar + 1./(Jbar * phi0);

    }

    /*!
     * Compute the second derivative of the gas relative volume with respect to
     * the matrix volume-conserving compression
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2JgdJbar2(const floatType &Jbar, const floatType &Je){

        auto phi0 = get_phi0();

        auto Jg = compute_Jg(Jbar, Je);

        auto dJgdJbar = compute_dJgdJbar(Jbar, Je);

        return -dJgdJbar/Jbar + Jg/(Jbar*Jbar) - Je/((Jbar * phi0)*(Jbar * phi0)) * phi0;

    }

    /*!
     * Compute the gas pressure
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_pg(const floatType &Jbar, const floatType &Je){

        auto p0 = get_p0();

        auto Jg = compute_Jg(Jbar, Je);

        auto gamma = get_gamma();

        return p0 * (std::pow(Jg, -gamma) - 1);

    }

    /*!
     * Compute the derivative of the gas pressure with respect to the volume-conserving compression
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_dpgdJbar(const floatType &Jbar, const floatType &Je){

        auto p0 = get_p0();

        auto Jg = compute_Jg(Jbar, Je);

        auto gamma = get_gamma();

        auto dJgdJbar = compute_dJgdJbar(Jbar, Je);

        return p0 * (-gamma * std::pow(Jg, -(gamma+1)) * dJgdJbar);

    }

    /*!
     * Compute the derivative of the gas pressure with respect to the elastic relative volume
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_dpgdJe(const floatType &Jbar, const floatType &Je){

        auto p0 = get_p0();

        auto Jg = compute_Jg(Jbar, Je);

        auto gamma = get_gamma();

        auto dJgdJe = compute_dJgdJe(Jbar, Je);

        return p0 * (-gamma * std::pow(Jg, -(gamma+1)) * dJgdJe);

    }

    /*!
     * Compute the second derivative of the gas pressure with respect to the net elastic relative volume
     * and the volume-conserving compression.
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2pgdJedJbar(const floatType &Jbar, const floatType &Je){

        auto p0 = get_p0();

        auto Jg = compute_Jg(Jbar, Je);

        auto gamma = get_gamma();

        auto dJgdJe = compute_dJgdJe(Jbar, Je);

        auto dJgdJbar = compute_dJgdJbar(Jbar, Je);

        auto d2JgdJedJbar = compute_d2JgdJedJbar(Jbar, Je);

        return p0 * (gamma * (gamma + 1) * std::pow(Jg, -(gamma+2)) * dJgdJe * dJgdJbar - gamma * std::pow(Jg, -(gamma+1)) * d2JgdJedJbar);

    }

    /*!
     * Compute the second derivative of the gas pressure with respect to the volume-conserving compression
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2pgdJbar2(const floatType &Jbar, const floatType &Je){

        auto p0 = get_p0();

        auto Jg = compute_Jg(Jbar, Je);

        auto gamma = get_gamma();

        auto dJgdJbar = compute_dJgdJbar(Jbar, Je);

        auto d2JgdJbar2 = compute_d2JgdJbar2(Jbar, Je);

        return p0 * (gamma * (gamma + 1) * std::pow(Jg, -(gamma+2)) * dJgdJbar * dJgdJbar - gamma * std::pow(Jg, -(gamma+1)) * d2JgdJbar2);

    }

    /*!
     * Compute the ptilde term
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_ptilde(const floatType &Jbar, const floatType &Je){

        auto pg = compute_pg(Jbar, Je);

        auto C10 = get_C10();

        auto phi0 = get_phi0();

        return pg + C10 * ( std::pow(phi0,1./3) * (4. * Jbar - 4. + 5. * phi0) * std::pow(Jbar - 1. + phi0,-4./3.) - (4. * Jbar - 1.)*(4.*Jbar + 1.) * std::pow(Jbar, -4./3.) / 3.);

    }

    /*!
     * Compute the derivative of the ptilde term with respect to the volume-conserving compression
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_dptildedJbar(const floatType &Jbar, const floatType &Je){

        auto dpgdJbar = compute_dpgdJbar(Jbar, Je);

        auto C10 = get_C10();

        auto phi0 = get_phi0();

        auto phi0_13 = std::pow(phi0,1./3);
        auto phi0_43 = std::pow(phi0,4./3);
        auto Jbar_13 = std::pow(Jbar,1./3);
        auto Jbar_73 = std::pow(Jbar,7./3);

        return dpgdJbar + ((-4*C10*Jbar*phi0_13 - 8*C10*phi0_43 + 4*C10*phi0_13) * std::pow((Jbar + phi0 - 1),-7./3.) / 3 - 32*C10/(9*Jbar_13) - 4*C10/(9*Jbar_73));

    }

    /*!
     * Compute the derivative of the ptilde term with respect to the elastic relative volume
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_dptildedJe(const floatType &Jbar, const floatType &Je){

        return compute_dpgdJe(Jbar, Je);

    }

    /*!
     * Compute the second derivative of the ptilde term with respect to the net elastic relative volume
     * and the volume-conserving compression.
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2ptildedJedJbar(const floatType &Jbar, const floatType &Je){

        return compute_d2pgdJedJbar(Jbar,Je);

    }

    /*!
     * Compute the second derivative of the ptilde term with respect to the volume-conserving compression
     *
     * \param Jbar: The matrix volume-conserving compression
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2ptildedJbar2(const floatType &Jbar, const floatType &Je){

        auto d2pgdJbar2 = compute_d2pgdJbar2(Jbar,Je);

        auto C10 = get_C10();

        auto phi0 = get_phi0();

        auto phi0_13 = std::pow(phi0,1./3);
        auto phi0_43 = std::pow(phi0,4./3);

        return d2pgdJbar2 + 16*C10*Jbar*phi0_13/(9*std::pow(Jbar + phi0 - 1,10./3)) + 44*C10*phi0_43/(9*std::pow(Jbar + phi0 - 1,10./3)) - 16*C10*phi0_13/(9*std::pow(Jbar + phi0 - 1,10./3)) + 32*C10/(27*std::pow(Jbar,4./3)) + 28*C10/(27*std::pow(Jbar,10./3));

    }

}
