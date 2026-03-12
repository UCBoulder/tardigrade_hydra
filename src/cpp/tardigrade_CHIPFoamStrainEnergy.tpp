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
    const floatType CHIPFoamStrainEnergy::compute_f(const floatType &J) {
        auto phi0 = get_phi0();

        floatType a = std::pow(phi0 / (J - (1 - phi0)), 1. / 3);

        return (2 * J - 1) * std::pow(J, -1. / 3) + (2 - 2 * J - phi0) * a;
    }

    /*!
     * Compute the derivative of the isochoric scale function
     */
    const floatType CHIPFoamStrainEnergy::compute_dfdJ(const floatType &J) {
        auto phi0 = get_phi0();

        floatType a    = std::pow(phi0 / (J - (1 - phi0)), 1. / 3);
        floatType dadJ = -a / (3 * (J - (1 - phi0)));

        return 2 * std::pow(J, -1. / 3) - (2 * J - 1) / 3 * std::pow(J, -4. / 3) - 2 * a + (2 - 2 * J - phi0) * dadJ;
    }

    /*!
     * Compute the second derivative of the isochoric scale function
     */
    const floatType CHIPFoamStrainEnergy::compute_d2fdJ2(const floatType &J) {
        auto phi0 = get_phi0();

        floatType a      = std::pow(phi0 / (J - (1 - phi0)), 1. / 3);
        floatType dadJ   = -a / (3 * (J - (1 - phi0)));
        floatType d2adJ2 = -4 * dadJ / (3 * (J - (1 - phi0)));

        return -(4. / 3) * std::pow(J, -4. / 3) + (4. / 9) * (2 * J - 1) * std::pow(J, -7. / 3) - 4 * dadJ +
               (2 - 2 * J - phi0) * d2adJ2;
    }

    /*!
     * Compute the gas relative volume
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_Jg(const floatType &Jbar, const floatType &Je) {
        auto phi0 = get_phi0();

        return Je / Jbar * (Jbar - 1 + phi0) / phi0;
    }

    /*!
     * Compute the derivative of the gas relative volume with respect to
     * the matrix volume-conserving relative volume
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_dJgdJbar(const floatType &Jbar, const floatType &Je) {
        auto phi0 = get_phi0();

        auto Jg = compute_Jg(Jbar, Je);

        return -Jg / Jbar + Je / (Jbar * phi0);
    }

    /*!
     * Compute the derivative of the gas relative volume with respect to
     * the net elastic relative volume
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_dJgdJe(const floatType &Jbar, const floatType &Je) {
        auto phi0 = get_phi0();

        return 1. / Jbar * (Jbar - 1 + phi0) / phi0;
    }

    /*!
     * Compute the second derivative of the gas relative volume with respect to
     * the matrix volume-conserving relative volume and the net elastic relative volume
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2JgdJedJbar(const floatType &Jbar, const floatType &Je) {
        auto phi0 = get_phi0();

        return -compute_dJgdJe(Jbar, Je) / Jbar + 1. / (Jbar * phi0);
    }

    /*!
     * Compute the second derivative of the gas relative volume with respect to
     * the matrix volume-conserving relative volume
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2JgdJbar2(const floatType &Jbar, const floatType &Je) {
        auto phi0 = get_phi0();

        auto Jg = compute_Jg(Jbar, Je);

        auto dJgdJbar = compute_dJgdJbar(Jbar, Je);

        return -dJgdJbar / Jbar + Jg / (Jbar * Jbar) - Je / ((Jbar * phi0) * (Jbar * phi0)) * phi0;
    }

    /*!
     * Compute the gas pressure
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_pg(const floatType &Jbar, const floatType &Je) {
        auto p0 = get_p0();

        auto Jg = compute_Jg(Jbar, Je);

        auto gamma = get_gamma();

        return p0 * (std::pow(Jg, -gamma) - 1);
    }

    /*!
     * Compute the derivative of the gas pressure with respect to the volume-conserving relative volume
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_dpgdJbar(const floatType &Jbar, const floatType &Je) {
        auto p0 = get_p0();

        auto Jg = compute_Jg(Jbar, Je);

        auto gamma = get_gamma();

        auto dJgdJbar = compute_dJgdJbar(Jbar, Je);

        return p0 * (-gamma * std::pow(Jg, -(gamma + 1)) * dJgdJbar);
    }

    /*!
     * Compute the derivative of the gas pressure with respect to the elastic relative volume
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_dpgdJe(const floatType &Jbar, const floatType &Je) {
        auto p0 = get_p0();

        auto Jg = compute_Jg(Jbar, Je);

        auto gamma = get_gamma();

        auto dJgdJe = compute_dJgdJe(Jbar, Je);

        return p0 * (-gamma * std::pow(Jg, -(gamma + 1)) * dJgdJe);
    }

    /*!
     * Compute the second derivative of the gas pressure with respect to the net elastic relative volume
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2pgdJe2(const floatType &Jbar, const floatType &Je) {
        auto p0 = get_p0();

        auto Jg = compute_Jg(Jbar, Je);

        auto gamma = get_gamma();

        auto dJgdJe = compute_dJgdJe(Jbar, Je);

        return p0 * (gamma * (gamma + 1) * std::pow(Jg, -(gamma + 2)) * dJgdJe * dJgdJe);
    }

    /*!
     * Compute the second derivative of the gas pressure with respect to the net elastic relative volume
     * and the volume-conserving relative volume.
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2pgdJedJbar(const floatType &Jbar, const floatType &Je) {
        auto p0 = get_p0();

        auto Jg = compute_Jg(Jbar, Je);

        auto gamma = get_gamma();

        auto dJgdJe = compute_dJgdJe(Jbar, Je);

        auto dJgdJbar = compute_dJgdJbar(Jbar, Je);

        auto d2JgdJedJbar = compute_d2JgdJedJbar(Jbar, Je);

        return p0 * (gamma * (gamma + 1) * std::pow(Jg, -(gamma + 2)) * dJgdJe * dJgdJbar -
                     gamma * std::pow(Jg, -(gamma + 1)) * d2JgdJedJbar);
    }

    /*!
     * Compute the second derivative of the gas pressure with respect to the volume-conserving relative volume
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2pgdJbar2(const floatType &Jbar, const floatType &Je) {
        auto p0 = get_p0();

        auto Jg = compute_Jg(Jbar, Je);

        auto gamma = get_gamma();

        auto dJgdJbar = compute_dJgdJbar(Jbar, Je);

        auto d2JgdJbar2 = compute_d2JgdJbar2(Jbar, Je);

        return p0 * (gamma * (gamma + 1) * std::pow(Jg, -(gamma + 2)) * dJgdJbar * dJgdJbar -
                     gamma * std::pow(Jg, -(gamma + 1)) * d2JgdJbar2);
    }

    /*!
     * Compute the ptilde term
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_ptilde(const floatType &Jbar, const floatType &Je) {
        auto pg = compute_pg(Jbar, Je);

        auto C10 = get_C10();

        auto phi0 = get_phi0();

        return pg +
               C10 * (std::pow(phi0, 1. / 3) * (4. * Jbar - 4. + 5. * phi0) * std::pow(Jbar - 1. + phi0, -4. / 3.) -
                      (4. * Jbar - 1.) * (4. * Jbar + 1.) * std::pow(Jbar, -4. / 3.) / 3.);
    }

    /*!
     * Compute the derivative of the ptilde term with respect to the volume-conserving relative volume
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_dptildedJbar(const floatType &Jbar, const floatType &Je) {
        auto dpgdJbar = compute_dpgdJbar(Jbar, Je);

        auto C10 = get_C10();

        auto phi0 = get_phi0();

        auto phi0_13 = std::pow(phi0, 1. / 3);
        auto phi0_43 = std::pow(phi0, 4. / 3);
        auto Jbar_13 = std::pow(Jbar, 1. / 3);
        auto Jbar_73 = std::pow(Jbar, 7. / 3);

        return dpgdJbar + ((-4 * C10 * Jbar * phi0_13 - 8 * C10 * phi0_43 + 4 * C10 * phi0_13) *
                               std::pow((Jbar + phi0 - 1), -7. / 3.) / 3 -
                           32 * C10 / (9 * Jbar_13) - 4 * C10 / (9 * Jbar_73));
    }

    /*!
     * Compute the derivative of the ptilde term with respect to the elastic relative volume
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_dptildedJe(const floatType &Jbar, const floatType &Je) {
        return compute_dpgdJe(Jbar, Je);
    }

    /*!
     * Compute the second derivative of the ptilde term with respect to the net elastic relative volume.
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2ptildedJe2(const floatType &Jbar, const floatType &Je) {
        return compute_d2pgdJe2(Jbar, Je);
    }

    /*!
     * Compute the second derivative of the ptilde term with respect to the net elastic relative volume
     * and the volume-conserving relative volume.
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2ptildedJedJbar(const floatType &Jbar, const floatType &Je) {
        return compute_d2pgdJedJbar(Jbar, Je);
    }

    /*!
     * Compute the second derivative of the ptilde term with respect to the volume-conserving relative volume
     *
     * \param Jbar: The matrix volume-conserving relative volume
     * \param Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2ptildedJbar2(const floatType &Jbar, const floatType &Je) {
        auto d2pgdJbar2 = compute_d2pgdJbar2(Jbar, Je);

        auto C10 = get_C10();

        auto phi0 = get_phi0();

        auto phi0_13  = std::pow(phi0, 1. / 3);
        auto phi0_43  = std::pow(phi0, 4. / 3);
        auto Jbar_43  = std::pow(Jbar, 4. / 3);
        auto Jbar_103 = std::pow(Jbar, 10. / 3);

        return d2pgdJbar2 +
               (16 * C10 * Jbar * phi0_13 + 44 * C10 * phi0_43 - 16 * C10 * phi0_13) /
                   (9 * std::pow(Jbar + phi0 - 1, 10. / 3)) +
               32 * C10 / (27 * Jbar_43) + 28 * C10 / (27 * Jbar_103);
    }

    /*!
     * Compute the relative volume of the parent material
     *
     * \param &Jbar: The matrix volume-conserving relative volume
     * \param &Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_Jm(const floatType &Jbar, const floatType &Je) {
        auto K = get_K();

        auto ptilde = compute_ptilde(Jbar, Je);

        return std::exp(-ptilde / K);
    }

    /*!
     * Compute the derivative of the relative volume of the parent material
     * with respect to the matrix volume-conserving relative volume
     *
     * \param &Jbar: The matrix volume-conserving relative volume
     * \param &Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_dJmdJbar(const floatType &Jbar, const floatType &Je) {
        auto K = get_K();

        auto Jm = compute_Jm(Jbar, Je);

        auto dptildedJbar = compute_dptildedJbar(Jbar, Je);

        return -Jm * dptildedJbar / K;
    }

    /*!
     * Compute the derivative of the relative volume of the parent material
     * with respect to the net elastic relative volume
     *
     * \param &Jbar: The matrix volume-conserving relative volume
     * \param &Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_dJmdJe(const floatType &Jbar, const floatType &Je) {
        auto K = get_K();

        auto Jm = compute_Jm(Jbar, Je);

        auto dptildedJe = compute_dptildedJe(Jbar, Je);

        return -Jm * dptildedJe / K;
    }

    /*!
     * Compute the second derivative of the relative volume of the parent material
     * with respect to the net elastic relative volume
     *
     * \param &Jbar: The matrix volume-conserving relative volume
     * \param &Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2JmdJe2(const floatType &Jbar, const floatType &Je) {
        auto K = get_K();

        auto Jm = compute_Jm(Jbar, Je);

        auto dJmdJe = compute_dJmdJe(Jbar, Je);

        auto dptildedJe = compute_dptildedJe(Jbar, Je);

        auto d2ptildedJe2 = compute_d2ptildedJe2(Jbar, Je);

        return -dJmdJe * dptildedJe / K - Jm * d2ptildedJe2 / K;
    }

    /*!
     * Compute the second derivative of the relative volume of the parent material
     * with respect to the net elastic relative volume and the matrix volume-conserving
     * relative volume
     *
     * \param &Jbar: The matrix volume-conserving relative volume
     * \param &Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2JmdJedJbar(const floatType &Jbar, const floatType &Je) {
        auto K = get_K();

        auto Jm = compute_Jm(Jbar, Je);

        auto dptildedJe = compute_dptildedJe(Jbar, Je);

        auto dJmdJbar = compute_dJmdJbar(Jbar, Je);

        auto d2ptildedJedJbar = compute_d2ptildedJedJbar(Jbar, Je);

        return -dJmdJbar * dptildedJe / K - Jm * d2ptildedJedJbar / K;
    }

    /*!
     * Compute the second derivative of the relative volume of the parent material
     * with respect to the matrix volume-conserving relative volume
     *
     * \param &Jbar: The matrix volume-conserving relative volume
     * \param &Je: The net elastic relative volume
     */
    const floatType CHIPFoamStrainEnergy::compute_d2JmdJbar2(const floatType &Jbar, const floatType &Je) {
        auto K = get_K();

        auto Jm = compute_Jm(Jbar, Je);

        auto dptildedJbar = compute_dptildedJbar(Jbar, Je);

        auto dJmdJbar = compute_dJmdJbar(Jbar, Je);

        auto d2ptildedJbar2 = compute_d2ptildedJbar2(Jbar, Je);

        return -dJmdJbar * dptildedJbar / K - Jm * d2ptildedJbar2 / K;
    }

    /*!
     * Compute the residual for the Jbar iteration
     *
     * \param &Jbar: The matrix volume-conserving relative volume
     * \param &Je: The net elastic relative volume
     */ 
    const floatType CHIPFoamStrainEnergy::compute_Jbar_residual(const floatType &Jbar, const floatType &Je) {

        auto Jm = compute_Jm(Jbar, Je);

        return Je / Jm - Jbar;
    }

    /*!
     * Compute the derivative of the residual for the Jbar iteration
     * with respect to the matrix volume-conserving relative volume
     *
     * \param &Jbar: The matrix volume-conserving relative volume
     * \param &Je: The net elastic relative volume
     */ 
    const floatType CHIPFoamStrainEnergy::compute_Jbar_dRdJbar(const floatType &Jbar, const floatType &Je) {

        auto Jm = compute_Jm(Jbar, Je);

        auto dJmdJbar = compute_dJmdJbar(Jbar, Je);

        return -Je / (Jm * Jm) * dJmdJbar - 1;

    }

    /*!
     * Compute the derivative of the residual for the Jbar iteration
     * with respect to the net elastic relative volume
     *
     * \param &Jbar: The matrix volume-conserving relative volume
     * \param &Je: The net elastic relative volume
     */ 
    const floatType CHIPFoamStrainEnergy::compute_Jbar_dRdJe(const floatType &Jbar, const floatType &Je) {

        auto Jm = compute_Jm(Jbar, Je);

        auto dJmdJe = compute_dJmdJe(Jbar, Je);

        return 1. / Jm - Je / (Jm * Jm) * dJmdJe;

    }

    /*!
     * Compute the second derivative of the residual for the Jbar iteration
     * with respect to the net elastic relative volume
     *
     * \param &Jbar: The matrix volume-conserving relative volume
     * \param &Je: The net elastic relative volume
     */ 
    const floatType CHIPFoamStrainEnergy::compute_Jbar_d2RdJe2(const floatType &Jbar, const floatType &Je) {

        auto Jm = compute_Jm(Jbar, Je);

        auto dJmdJe = compute_dJmdJe(Jbar, Je);

        auto d2JmdJe2 = compute_d2JmdJe2(Jbar, Je);

        return -1. / (Jm * Jm) * dJmdJe - 1 / (Jm * Jm) * dJmdJe + 2 * Je / (Jm * Jm * Jm) * dJmdJe * dJmdJe - Je / (Jm * Jm) * d2JmdJe2;

    }

    /*!
     * Compute the second derivative of the residual for the Jbar iteration
     * with respect to the net elastic relative volume and the matrix volume-conserving
     * relative volume
     *
     * \param &Jbar: The matrix volume-conserving relative volume
     * \param &Je: The net elastic relative volume
     */ 
    const floatType CHIPFoamStrainEnergy::compute_Jbar_d2RdJedJbar(const floatType &Jbar, const floatType &Je) {
        auto Jm = compute_Jm(Jbar, Je);

        auto dJmdJe = compute_dJmdJe(Jbar, Je);

        auto dJmdJbar = compute_dJmdJbar(Jbar, Je);

        auto d2JmdJedJbar = compute_d2JmdJedJbar(Jbar, Je);

        return -1. / (Jm * Jm) * dJmdJbar + 2 * Je / (Jm * Jm * Jm) * dJmdJe * dJmdJbar - Je / (Jm * Jm) * d2JmdJedJbar;
    }

    /*!
     * Compute the second derivative of the residual for the Jbar iteration
     * with respect to the the matrix volume-conserving relative volume
     *
     * \param &Jbar: The matrix volume-conserving relative volume
     * \param &Je: The net elastic relative volume
     */ 
    const floatType CHIPFoamStrainEnergy::compute_Jbar_d2RdJbar2(const floatType &Jbar, const floatType &Je) {

        auto Jm = compute_Jm(Jbar, Je);

        auto dJmdJbar = compute_dJmdJbar(Jbar, Je);

        auto d2JmdJbar2 = compute_d2JmdJbar2(Jbar, Je);

        return 2 * Je / (Jm * Jm * Jm) * dJmdJbar * dJmdJbar - Je / (Jm * Jm) * d2JmdJbar2;

    }

}  // namespace tardigradeHydra
