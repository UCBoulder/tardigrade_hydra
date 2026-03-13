/**
 ******************************************************************************
 * \file tardigrade_CHIPFoamStrainEnergy.cpp
 ******************************************************************************
 *
 ******************************************************************************
 */

#include "tardigrade_CHIPFoamStrainEnergy.h"

namespace tardigradeHydra {

    /*!
     * Get the shear modulus
     */
    const floatType CHIPFoamStrainEnergy::get_Ghat() {
        TARDIGRADE_ERROR_TOOLS_CHECK(isInitialized(), "The parameter vector has not been initialized");
        return _parameters[0];
    }

    /*!
     * Get the bulk modulus
     */
    const floatType CHIPFoamStrainEnergy::get_Khat() {
        TARDIGRADE_ERROR_TOOLS_CHECK(isInitialized(), "The parameter vector has not been initialized");
        return _parameters[1];
    }

    /*!
     * Get the buckling relative volume
     */
    const floatType CHIPFoamStrainEnergy::get_Jb() {
        TARDIGRADE_ERROR_TOOLS_CHECK(isInitialized(), "The parameter vector has not been initialized");
        return _parameters[2];
    }

    /*!
     * Get the effective Neo-Hookean modulus
     */
    const floatType CHIPFoamStrainEnergy::get_C10() {
        TARDIGRADE_ERROR_TOOLS_CHECK(isInitialized(), "The parameter vector has not been initialized");
        return _parameters[3];
    }

    /*!
     * Get the initial porosity
     */
    const floatType CHIPFoamStrainEnergy::get_phi0() {
        TARDIGRADE_ERROR_TOOLS_CHECK(isInitialized(), "The parameter vector has not been initialized");
        return _parameters[4];
    }

    /*!
     * Get the parent material effective bulk modulus
     */
    const floatType CHIPFoamStrainEnergy::get_K() {
        TARDIGRADE_ERROR_TOOLS_CHECK(isInitialized(), "The parameter vector has not been initialized");
        return _parameters[5];
    }

    /*!
     * Get the initial gas pressure
     */
    const floatType CHIPFoamStrainEnergy::get_p0() {
        TARDIGRADE_ERROR_TOOLS_CHECK(isInitialized(), "The parameter vector has not been initialized");
        return _parameters[6];
    }

    /*!
     * Get the ratio of gas specific heats
     */
    const floatType CHIPFoamStrainEnergy::get_gamma() {
        TARDIGRADE_ERROR_TOOLS_CHECK(isInitialized(), "The parameter vector has not been initialized");
        return _parameters[7];
    }

    /*!
     * Compute the Jacobian of the elastic deformation
     *
     * \param isPrevious: A flag for of the current (false) or previous (true) value should be computed
     */
    void CHIPFoamStrainEnergy::setJe(bool isPrevious) {
        constexpr unsigned int dim = 3;

        const secondOrderTensor *Fe;

        SetDataStorageBase<floatType> Je;

        if (isPrevious) {
            Fe = get_previousFe();

            Je = get_SetDataStorage_previousJe();

        } else {
            Fe = get_Fe();

            Je = get_SetDataStorage_Je();
        }

        Eigen::Map<const Eigen::Matrix<floatType, -1, -1, Eigen::RowMajor>> Femat(
            Fe->data(), dim, dim);  // TODO: Change this to a constant size when possible
        *Je.value = Femat.determinant();
    }

    /*!
     * Compute the Jacobian of the elastic deformation
     */
    void CHIPFoamStrainEnergy::setJe() { setJe(false); }

    /*!
     * Compute the previous Jacobian of the elastic deformation
     */
    void CHIPFoamStrainEnergy::setPreviousJe() { setJe(true); }

    /*!
     * Compute the derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true)
     *     value
     */
    void CHIPFoamStrainEnergy::setdJedFe(bool isPrevious) {
        constexpr unsigned int dim = 3;

        const secondOrderTensor *Fe;

        const floatType *Je;

        SetDataStorageBase<secondOrderTensor> dJedFe;

        if (isPrevious) {
            Fe = get_previousFe();

            Je = get_previousJe();

            dJedFe = get_SetDataStorage_dPreviousJedPreviousFe();

        } else {
            Fe = get_Fe();

            Je = get_Je();

            dJedFe = get_SetDataStorage_dJedFe();
        }

        secondOrderTensor                                                   invFe(dim * dim, 0);
        Eigen::Map<const Eigen::Matrix<floatType, -1, -1, Eigen::RowMajor>> Femat(
            Fe->data(), dim, dim);  // TODO: Change this to a constant size when possible
        Eigen::Map<Eigen::Matrix<floatType, -1, -1, Eigen::RowMajor>> invFemat(
            invFe.data(), dim, dim);  // TODO: Change this to a constant size when possible
        invFemat = Femat.inverse().eval();

        dJedFe.zero(dim * dim);

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int I = 0; I < dim; ++I) {
                (*dJedFe.value)[dim * i + I] = (*Je) * invFe[dim * I + i];
            }
        }
    }

    /*!
     * Compute the derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setdJedFe() { setdJedFe(false); }

    /*!
     * Compute the previous derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setdPreviousJedPreviousFe() { setdJedFe(true); }

    /*!
     * Compute the second derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true)
     *     value
     */
    void CHIPFoamStrainEnergy::setd2JedFe2(bool isPrevious) {
        constexpr unsigned int dim = 3;

        const floatType *Je;

        const secondOrderTensor *dJedFe;

        SetDataStorageBase<fourthOrderTensor> d2JedFe2;

        if (isPrevious) {
            Je = get_previousJe();

            dJedFe = get_dPreviousJedPreviousFe();

            d2JedFe2 = get_SetDataStorage_d2PreviousJedPreviousFe2();

        } else {
            Je = get_Je();

            dJedFe = get_dJedFe();

            d2JedFe2 = get_SetDataStorage_d2JedFe2();
        }

        d2JedFe2.zero(dim * dim * dim * dim);

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int I = 0; I < dim; ++I) {
                for (unsigned int j = 0; j < dim; ++j) {
                    for (unsigned int J = 0; J < dim; ++J) {
                        (*d2JedFe2.value)[dim * dim * dim * i + dim * dim * I + dim * j + J] =
                            ((*dJedFe)[dim * j + J] * (*dJedFe)[dim * i + I] -
                             (*dJedFe)[dim * j + I] * (*dJedFe)[dim * i + J]) /
                            (*Je);
                    }
                }
            }
        }
    }

    /*!
     * Compute the derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setd2JedFe2() { setd2JedFe2(false); }

    /*!
     * Compute the previous derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setd2PreviousJedPreviousFe2() { setd2JedFe2(true); }

    /*!
     * Compute the first invariant of the isochoric elastic deformation
     *
     * \param isPrevious: A flag for of the current (false) or previous (true) value should be computed
     */
    void CHIPFoamStrainEnergy::setIbar1(bool isPrevious) {
        constexpr unsigned int dim = 3;

        const secondOrderTensor *Fe;

        const floatType *Je;

        SetDataStorageBase<floatType> Ibar1;

        if (isPrevious) {
            Fe = get_previousFe();

            Je = get_previousJe();

            Ibar1 = get_SetDataStorage_previousIbar1();

        } else {
            Fe = get_Fe();

            Je = get_Je();

            Ibar1 = get_SetDataStorage_Ibar1();
        }

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int I = 0; I < dim; ++I) {
                *Ibar1.value += (*Fe)[dim * i + I] * (*Fe)[dim * i + I];
            }
        }

        *Ibar1.value /= std::pow(*Je, 2. / 3);
    }

    /*!
     * Compute the first invariant of the isochoric elastic deformation
     */
    void CHIPFoamStrainEnergy::setIbar1() { setIbar1(false); }

    /*!
     * Compute the previous first invariant of the isochoric elastic deformation
     */
    void CHIPFoamStrainEnergy::setPreviousIbar1() { setIbar1(true); }

    /*!
     * Compute the derivative of the first invariant of the isochoric elastic deformation
     * with respect to the elastic deformation
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true)
     *     value
     */
    void CHIPFoamStrainEnergy::setdIbar1dFe(bool isPrevious) {
        constexpr unsigned int dim = 3;

        const floatType *Ibar1;

        const secondOrderTensor *Fe;

        const floatType *Je;

        const secondOrderTensor *dJedFe;

        SetDataStorageBase<secondOrderTensor> dIbar1dFe;

        if (isPrevious) {
            Ibar1 = get_previousIbar1();

            Fe = get_previousFe();

            Je = get_previousJe();

            dJedFe = get_dPreviousJedPreviousFe();

            dIbar1dFe = get_SetDataStorage_dPreviousIbar1dPreviousFe();

        } else {
            Ibar1 = get_Ibar1();

            Fe = get_Fe();

            Je = get_Je();

            dJedFe = get_dJedFe();

            dIbar1dFe = get_SetDataStorage_dIbar1dFe();
        }

        dIbar1dFe.zero(dim * dim);

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int I = 0; I < dim; ++I) {
                (*dIbar1dFe.value)[dim * i + I] =
                    2 * ((*Fe)[dim * i + I] * std::pow((*Je), 1. / 3) - (*Ibar1) * (*dJedFe)[dim * i + I] / 3.) / (*Je);
            }
        }
    }

    /*!
     * Compute the derivative of the first invariant of the isochoric elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setdIbar1dFe() { setdIbar1dFe(false); }

    /*!
     * Compute the previous derivative of the first invariant of the isochoric elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setdPreviousIbar1dPreviousFe() { setdIbar1dFe(true); }

    /*!
     * Compute the second derivative of the first invariant of the isochoric elastic deformation
     * with respect to the elastic deformation
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true)
     *     value
     */
    void CHIPFoamStrainEnergy::setd2Ibar1dFe2(bool isPrevious) {
        constexpr unsigned int dim = 3;

        const floatType *Ibar1;

        const secondOrderTensor *dIbar1dFe;

        const secondOrderTensor *Fe;

        const floatType *Je;

        const secondOrderTensor *dJedFe;

        const fourthOrderTensor *d2JedFe2;

        SetDataStorageBase<fourthOrderTensor> d2Ibar1dFe2;

        if (isPrevious) {
            Ibar1 = get_previousIbar1();

            dIbar1dFe = get_dPreviousIbar1dPreviousFe();

            Fe = get_previousFe();

            Je = get_previousJe();

            dJedFe = get_dPreviousJedPreviousFe();

            d2JedFe2 = get_d2PreviousJedPreviousFe2();

            d2Ibar1dFe2 = get_SetDataStorage_d2PreviousIbar1dPreviousFe2();

        } else {
            Ibar1 = get_Ibar1();

            dIbar1dFe = get_dIbar1dFe();

            Fe = get_Fe();

            Je = get_Je();

            dJedFe = get_dJedFe();

            d2JedFe2 = get_d2JedFe2();

            d2Ibar1dFe2 = get_SetDataStorage_d2Ibar1dFe2();
        }

        d2Ibar1dFe2.zero(dim * dim * dim * dim);

        auto Je_23 = std::pow((*Je), 2. / 3.);

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int I = 0; I < dim; ++I) {
                (*d2Ibar1dFe2.value)[dim * dim * dim * i + dim * dim * I + dim * i + I] += 2 / Je_23;
                for (unsigned int j = 0; j < dim; ++j) {
                    for (unsigned int J = 0; J < dim; ++J) {
                        (*d2Ibar1dFe2.value)[dim * dim * dim * i + dim * dim * I + dim * j + J] +=
                            2. / 3. *
                                ((*Fe)[dim * i + I] / Je_23 * (*dJedFe)[dim * j + J] -
                                 (*dIbar1dFe)[dim * j + J] * (*dJedFe)[dim * i + I] -
                                 (*Ibar1) * (*d2JedFe2)[dim * dim * dim * i + dim * dim * I + dim * j + J]) /
                                (*Je) -
                            (*dIbar1dFe)[dim * i + I] / (*Je) * (*dJedFe)[dim * j + J];
                    }
                }
            }
        }
    }

    /*!
     * Compute the derivative of the first invariant of the isochoric elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setd2Ibar1dFe2() { setd2Ibar1dFe2(false); }

    /*!
     * Compute the previous derivative of the first invariant of the isochoric elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setd2PreviousIbar1dPreviousFe2() { setd2Ibar1dFe2(true); }

    /*!
     * Set the buckling Neo-Hookean strain energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setWLB(bool isPrevious) {
        const floatType *Je;

        const floatType *Ibar1;

        auto Jb = get_Jb();

        auto Khat = get_Khat();

        auto Ghat = get_Ghat();

        SetDataStorageBase<floatType> WLB;

        if (isPrevious) {
            Je = get_previousJe();

            Ibar1 = get_previousIbar1();

            WLB = get_SetDataStorage_previousWLB();

        } else {
            Je = get_Je();

            Ibar1 = get_Ibar1();

            WLB = get_SetDataStorage_WLB();
        }

        *WLB.value = Ghat * (*Ibar1 - 3) / 2. + Khat * (Jb - 1) * (*Je - Jb / 2. - 0.5);

        if (*Je >= Jb) {
            *WLB.value += Khat * ((*Je - 1) * (*Je - 1) / 2. - (Jb - 1) * (*Je - Jb / 2. - 0.5));
        }
    }

    /*!
     * Set the buckling Neo-Hookean strain energy
     */
    void CHIPFoamStrainEnergy::setWLB() { setWLB(false); }

    /*!
     * Set the previous buckling Neo-Hookean strain energy
     */
    void CHIPFoamStrainEnergy::setPreviousWLB() { setWLB(true); }

    /*!
     * Set the derivatives of the buckling strain-energy with respect to
     * the deformation measures
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setWLBDerivatives(bool isPrevious) {
        const floatType *Je;

        auto Jb = get_Jb();

        auto Khat = get_Khat();

        auto Ghat = get_Ghat();

        SetDataStorageBase<floatVector> dWLBdD;

        if (isPrevious) {
            Je = get_previousJe();

            dWLBdD = get_SetDataStorage_dPreviousWLBdPreviousD();

        } else {
            Je = get_Je();

            dWLBdD = get_SetDataStorage_dWLBdD();
        }

        dWLBdD.zero(2);

        (*dWLBdD.value)[0] = Khat * (Jb - 1);
        (*dWLBdD.value)[1] = 0.5 * Ghat;

        if (*Je >= Jb) {
            (*dWLBdD.value)[0] += Khat * (*Je - Jb);
        }
    }

    /*!
     * Set the derivatives of the buckling strain-energy with respect to
     * the deformation measures
     */
    void CHIPFoamStrainEnergy::setWLBDerivatives() { setWLBDerivatives(false); }

    /*!
     * Set the previous derivatives of the buckling strain-energy with respect to
     * the deformation measures
     */
    void CHIPFoamStrainEnergy::setPreviousWLBDerivatives() { setWLBDerivatives(true); }

    /*!
     * Set the Hessians of the buckling strain-energy with respect to
     * the deformation measures
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setWLBHessians(bool isPrevious) {
        const floatType *Je;

        auto Jb = get_Jb();

        auto Khat = get_Khat();

        SetDataStorageBase<floatVector> d2WLBdD2;

        if (isPrevious) {
            Je = get_previousJe();

            d2WLBdD2 = get_SetDataStorage_d2PreviousWLBdPreviousD2();

        } else {
            Je = get_Je();

            d2WLBdD2 = get_SetDataStorage_d2WLBdD2();
        }

        d2WLBdD2.zero(4);

        if (*Je >= Jb) {
            (*d2WLBdD2.value)[0] += Khat;
        }
    }

    /*!
     * Set the Hessians of the buckling strain-energy with respect to
     * the deformation measures
     */
    void CHIPFoamStrainEnergy::setWLBHessians() { setWLBHessians(false); }

    /*!
     * Set the previous Hessians of the buckling strain-energy with respect to
     * the deformation measures
     */
    void CHIPFoamStrainEnergy::setPreviousWLBHessians() { setWLBHessians(true); }

    /*!
     * Set the value of Jbar
     *
     * \param isPrevious: Whether to set the current (false) or previous (true)
     *     value of Jbar
     */
    void CHIPFoamStrainEnergy::setJbar(bool isPrevious){

        if (isPrevious){

            auto Jbar = get_SetDataStorage_previousJbar();

            *Jbar.value = Jbar_newton(*get_previousJe());

        }else{

            auto Jbar = get_SetDataStorage_Jbar();

            *Jbar.value = Jbar_newton(*get_Je());

        }

    }

    /*!
     * Set the current value of Jbar
     */
    void CHIPFoamStrainEnergy::setJbar(){ setJbar(false); }

    /*!
     * Set the previous value of Jbar
     */
    void CHIPFoamStrainEnergy::setPreviousJbar(){ setJbar(true); }

    /*!
     * Compute the derivative of Jbar with respect to the
     * net elastic relative volume
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true)
     *     value
     */
    void CHIPFoamStrainEnergy::setdJbardJe(bool isPrevious){

        const floatType *Jbar;

        const floatType *Je;

        SetDataStorageBase<floatType> dJbardJe;

        if ( isPrevious ){

            Jbar = get_previousJbar();

            Je = get_previousJe();

            dJbardJe = get_SetDataStorage_dPreviousJbardPreviousJe();

        }else{

            Jbar = get_Jbar();

            Je = get_Je();

            dJbardJe = get_SetDataStorage_dJbardJe();

        }

        auto jacobian = compute_Jbar_dRdJbar(*Jbar, *Je);

        auto dRdX = compute_Jbar_dRdJe(*Jbar, *Je);

        *dJbardJe.value = -dRdX / jacobian;

    }

    /*!
     * Compute the current derivative of Jbar with respect to the
     * net elastic relative volume
     */
    void CHIPFoamStrainEnergy::setdJbardJe(){ setdJbardJe(false); }

    /*!
     * Compute the previous derivative of Jbar with respect to the
     * net elastic relative volume
     */
    void CHIPFoamStrainEnergy::setdPreviousJbardPreviousJe(){ setdJbardJe(true); }

    /*!
     * Compute the second derivative of Jbar with respect to the
     * net elastic relative volume
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true)
     *     value
     */
    void CHIPFoamStrainEnergy::setd2JbardJe2(bool isPrevious){

        const floatType *Jbar;

        const floatType *Je;

        const floatType *dJbardJe;

        SetDataStorageBase<floatType> d2JbardJe2;

        if ( isPrevious ){

            Jbar = get_previousJbar();

            Je = get_previousJe();

            dJbardJe = get_dPreviousJbardPreviousJe();

            d2JbardJe2 = get_SetDataStorage_d2PreviousJbardPreviousJe2();

        }else{

            Jbar = get_Jbar();

            Je = get_Je();

            dJbardJe = get_dJbardJe();

            d2JbardJe2 = get_SetDataStorage_d2JbardJe2();

        }

        auto _J = compute_Jbar_dRdJbar(*Jbar, *Je);

        auto _dJdX = compute_Jbar_d2RdJbar2(*Jbar, *Je);

        auto _dJdTheta = compute_Jbar_d2RdJedJbar(*Jbar, *Je);

        auto _d2RdTheta2 = compute_Jbar_d2RdJe2(*Jbar, *Je);

        auto _d2RdThetadX = compute_Jbar_d2RdJedJbar(*Jbar, *Je);

        *d2JbardJe2.value = -(_dJdTheta * (*dJbardJe) + _dJdX * (*dJbardJe) * (*dJbardJe) + _d2RdTheta2 + _d2RdThetadX * (*dJbardJe)) / _J;

    }

    /*!
     * Compute the current second derivative of Jbar with respect to the
     * net elastic relative volume
     */
    void CHIPFoamStrainEnergy::setd2JbardJe2(){ setd2JbardJe2(false); }

    /*!
     * Compute the previous second derivative of Jbar with respect to the
     * net elastic relative volume
     */
    void CHIPFoamStrainEnergy::setd2PreviousJbardPreviousJe2(){ setd2JbardJe2(true); }

    /*!
     * Set the derivative of Jbar with respect to Je evaluated when Jbar is one
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true)
     *     value
     */
    void CHIPFoamStrainEnergy::setdJbardJe1(bool isPrevious){

        auto K = get_K();

        auto p0 = get_p0();

        auto gamma = get_gamma();

        const floatType *Je;

        SetDataStorageBase<floatType> dJbardJe1;

        if (isPrevious){
            Je = get_previousJe();

            dJbardJe1 = get_SetDataStorage_previousdJbardJe1();
        }else{
            Je = get_Je();

            dJbardJe1 = get_SetDataStorage_dJbardJe1();
        }

        *dJbardJe1.value = std::exp(p0*(-1 + std::pow(*Je,-gamma))/K) - gamma*p0*std::exp(p0*(-1 + std::pow(*Je,-gamma))/K)/(std::pow(*Je,gamma)*K);

    }

    /*!
     * Set the current derivative of Jbar with respect to Je evaluated when Jbar is one
     */
    void CHIPFoamStrainEnergy::setdJbardJe1(){ setdJbardJe1(false); }

    /*!
     * Set the previous derivative of Jbar with respect to Je evaluated when Jbar is one
     */
    void CHIPFoamStrainEnergy::setPreviousdJbardJe1(){ setdJbardJe1(true); }


    /*!
     * Set the derivative of the derivative of Jbar with respect to Je evaluated when Jbar is one
     * with respect to the net elastic relative volume
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true)
     *     value
     */
    void CHIPFoamStrainEnergy::setddJbardJe1dJe(bool isPrevious){

        auto K = get_K();

        auto p0 = get_p0();

        auto gamma = get_gamma();

        const floatType *Je;

        SetDataStorageBase<floatType> ddJbardJe1dJe;

        if (isPrevious){
            Je = get_previousJe();

            ddJbardJe1dJe = get_SetDataStorage_dPreviousdJbardJe1dPreviousJe();
        }else{
            Je = get_Je();

            ddJbardJe1dJe = get_SetDataStorage_ddJbardJe1dJe();
        }

        *ddJbardJe1dJe.value = gamma*gamma*p0*std::exp(p0*(-1 + std::pow(*Je,-gamma))/K)/((*Je)*std::pow(*Je,gamma)*K) - gamma*p0*std::exp(p0*(-1 + std::pow(*Je,-gamma))/K)/((*Je)*std::pow(*Je,gamma)*K) + gamma*gamma*p0*p0*std::exp(p0*(-1 + std::pow(*Je,-gamma))/K)/((*Je)*std::pow(*Je,2*gamma)*K*K);
    }

    /*!
     * Set the current derivative of the derivative of Jbar with respect to Je evaluated when Jbar is one
     * with respect to the net elastic relative volume
     */
    void CHIPFoamStrainEnergy::setddJbardJe1dJe(){ setddJbardJe1dJe(false); }

    /*!
     * Set the previous derivative of the derivative of Jbar with respect to Je evaluated when Jbar is one
     * with respect to the net elastic relative volume
     */
    void CHIPFoamStrainEnergy::setdPreviousdJbardJe1dPreviousJe(){ setddJbardJe1dJe(true); }

    /*!
     * Set the second derivative of the derivative of Jbar with respect to Je evaluated when Jbar is one
     * with respect to the net elastic relative volume
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true)
     *     value
     */
    void CHIPFoamStrainEnergy::setd2dJbardJe1dJe2(bool isPrevious){

        auto K = get_K();

        auto p0 = get_p0();

        auto gamma = get_gamma();

        const floatType *Je;

        SetDataStorageBase<floatType> d2dJbardJe1dJe2;

        if (isPrevious){
            Je = get_previousJe();

            d2dJbardJe1dJe2 = get_SetDataStorage_d2PreviousdJbardJe1dPreviousJe2();
        }else{
            Je = get_Je();

            d2dJbardJe1dJe2 = get_SetDataStorage_d2dJbardJe1dJe2();
        }

        *d2dJbardJe1dJe2.value = -gamma*gamma*gamma*p0*exp(p0*(-1 + std::pow(*Je,-gamma))/K)/((*Je)*(*Je)*std::pow(*Je,gamma)*K) + gamma*p0*exp(p0*(-1 + std::pow(*Je,-gamma))/K)/((*Je)*(*Je)*std::pow(*Je,gamma)*K) - 3*gamma*gamma*gamma*p0*p0*exp(p0*(-1 + std::pow(*Je,-gamma))/K)/((*Je)*(*Je)*std::pow(*Je,2*gamma)*K*K) - gamma*gamma*gamma*p0*p0*p0*exp(p0*(-1 + std::pow(*Je,-gamma))/K)/((*Je)*(*Je)*std::pow(*Je,3*gamma)*K*K*K);

    }

    /*!
     * Set the current second derivative of the derivative of Jbar with respect to Je evaluated when Jbar is one
     * with respect to the net elastic relative volume
     */
    void CHIPFoamStrainEnergy::setd2dJbardJe1dJe2(){ setd2dJbardJe1dJe2(false); }

    /*!
     * Set the previous second derivative of the derivative of Jbar with respect to Je evaluated when Jbar is one
     * with respect to the net elastic relative volume
     */
    void CHIPFoamStrainEnergy::setd2PreviousdJbardJe1dPreviousJe2(){ setd2dJbardJe1dJe2(true); }

    /*!
     * Set the value of the modified Danielsson function
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setWDC(bool isPrevious){

        auto phi0 = get_phi0();

        auto C10 = get_C10();

        const floatType *Je;

        const floatType *Ibar1;

        const floatType *Jbar;

        const floatType *dJbardJe1;

        SetDataStorageBase<floatType> WDC;

        if(isPrevious){

            Je = get_previousJe();

            Ibar1 = get_previousIbar1();

            Jbar = get_previousJbar();

            dJbardJe1 = get_previousdJbardJe1();

            WDC = get_SetDataStorage_previousWDC();

        }else{

            Je = get_Je();

            Ibar1 = get_Ibar1();

            Jbar = get_Jbar();

            dJbardJe1 = get_dJbardJe1();

            WDC = get_SetDataStorage_WDC();

        }

        auto Jm = compute_Jm(*Jbar, *Je);

        auto f = compute_f(*Jbar);

        *WDC.value = C10 * (Jm * ( (*Ibar1) * f - 3 * ( 1 - phi0 ) ) - (*Je) * ( (*Ibar1) - 3 )*(1 - phi0) * (*dJbardJe1));

    }

    /*!
     * Set the current value of the modified Danielsson function
     */
    void CHIPFoamStrainEnergy::setWDC(){ setWDC(false); }

    /*!
     * Set the previous value of the modified Danielsson function
     */
    void CHIPFoamStrainEnergy::setPreviousWDC(){ setWDC(true); }

    /*!
     * Set the derivatives of the modified Danielsson function
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setWDCDerivatives(bool isPrevious){

        auto phi0 = get_phi0();

        auto C10 = get_C10();

        const floatType *Je;

        const floatType *Ibar1;

        const floatType *Jbar;

        const floatType *dJbardJe;

        const floatType *dJbardJe1;

        const floatType *ddJbardJe1dJe;

        SetDataStorageBase<floatVector> dWDCdD;

        if(isPrevious){

            Je = get_previousJe();

            Ibar1 = get_previousIbar1();

            Jbar = get_previousJbar();

            dJbardJe = get_dPreviousJbardPreviousJe();

            dJbardJe1 = get_previousdJbardJe1();

            ddJbardJe1dJe = get_dPreviousdJbardJe1dPreviousJe();

            dWDCdD = get_SetDataStorage_dPreviousWDCdPreviousD();

        }else{

            Je = get_Je();

            Ibar1 = get_Ibar1();

            Jbar = get_Jbar();

            dJbardJe = get_dJbardJe();

            dJbardJe1 = get_dJbardJe1();

            ddJbardJe1dJe = get_ddJbardJe1dJe();

            dWDCdD = get_SetDataStorage_dWDCdD();

        }

        auto Jm = compute_Jm(*Jbar, *Je);

        auto dJmdJe = compute_dJmdJe(*Jbar, *Je) + compute_dJmdJbar(*Jbar, *Je) * (*dJbardJe);

        auto f = compute_f(*Jbar);

        auto dfdJe = compute_dfdJ(*Jbar) * (*dJbardJe);

        dWDCdD.zero(2);

        (*dWDCdD.value)[0] = C10 * (dJmdJe * ( (*Ibar1) * f - 3 * ( 1 - phi0 ) ) + Jm * (*Ibar1) * dfdJe - ( (*Ibar1) - 3 )*(1 - phi0) * (*dJbardJe1) - (*Je) * ( (*Ibar1) - 3 )*(1 - phi0) * (*ddJbardJe1dJe));
        (*dWDCdD.value)[1] = C10 * (Jm * f - (*Je) * (1 - phi0) * (*dJbardJe1));
    }

    /*!
     * Set the current derivatives of the modified Danielsson function
     */
    void CHIPFoamStrainEnergy::setWDCDerivatives(){ setWDCDerivatives(false); }

    /*!
     * Set the previous derivatives of the modified Danielsson function
     */
    void CHIPFoamStrainEnergy::setPreviousWDCDerivatives(){ setWDCDerivatives(true); }

    /*!
     * Set the Hessians of the modified Danielsson function
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setWDCHessians(bool isPrevious){

        auto phi0 = get_phi0();

        auto C10 = get_C10();

        const floatType *Je;

        const floatType *Ibar1;

        const floatType *Jbar;

        const floatType *dJbardJe;

        const floatType *d2JbardJe2;

        const floatType *dJbardJe1;

        const floatType *ddJbardJe1dJe;

        const floatType *d2dJbardJe1dJe2;

        SetDataStorageBase<floatVector> d2WDCdD2;

        if(isPrevious){

            Je = get_previousJe();

            Ibar1 = get_previousIbar1();

            Jbar = get_previousJbar();

            dJbardJe = get_dPreviousJbardPreviousJe();

            d2JbardJe2 = get_d2PreviousJbardPreviousJe2();

            dJbardJe1 = get_previousdJbardJe1();

            ddJbardJe1dJe = get_dPreviousdJbardJe1dPreviousJe();

            d2dJbardJe1dJe2 = get_d2PreviousdJbardJe1dPreviousJe2();

            d2WDCdD2 = get_SetDataStorage_d2PreviousWDCdPreviousD2();

        }else{

            Je = get_Je();

            Ibar1 = get_Ibar1();

            Jbar = get_Jbar();

            dJbardJe = get_dJbardJe();

            d2JbardJe2 = get_d2JbardJe2();

            dJbardJe1 = get_dJbardJe1();

            ddJbardJe1dJe = get_ddJbardJe1dJe();

            d2dJbardJe1dJe2 = get_d2dJbardJe1dJe2();

            d2WDCdD2 = get_SetDataStorage_d2WDCdD2();

        }

        auto Jm = compute_Jm(*Jbar, *Je);

        auto dJmdJe = compute_dJmdJe(*Jbar, *Je) + compute_dJmdJbar(*Jbar, *Je) * (*dJbardJe);

        auto d2JmdJe2 = compute_d2JmdJe2(*Jbar, *Je) + 2 * compute_d2JmdJedJbar(*Jbar, *Je) * (*dJbardJe) + compute_d2JmdJbar2(*Jbar, *Je) * (*dJbardJe) * (*dJbardJe) + compute_dJmdJbar(*Jbar, *Je) * (*d2JbardJe2);

        auto f = compute_f(*Jbar);

        auto dfdJe = compute_dfdJ(*Jbar) * (*dJbardJe);

        auto d2fdJe2 = compute_d2fdJ2(*Jbar) * (*dJbardJe) * (*dJbardJe) + compute_dfdJ(*Jbar) * (*d2JbardJe2);

        d2WDCdD2.zero(4);

        (*d2WDCdD2.value)[0] = C10 * (d2JmdJe2 * ( (*Ibar1) * f - 3 * ( 1 - phi0 ) ) + dJmdJe * (*Ibar1) * dfdJe + dJmdJe * (*Ibar1) * dfdJe + Jm * (*Ibar1) * d2fdJe2 - ( (*Ibar1) - 3 )*(1 - phi0) * (*ddJbardJe1dJe) - ( (*Ibar1) - 3 )*(1 - phi0) * (*ddJbardJe1dJe) - (*Je) * ( (*Ibar1) - 3 )*(1 - phi0) * (*d2dJbardJe1dJe2));
        (*d2WDCdD2.value)[1] = C10 * (dJmdJe * f + Jm * dfdJe - (1 - phi0) * (*dJbardJe1) - (*Je) * (1 - phi0) * (*ddJbardJe1dJe));
        (*d2WDCdD2.value)[2] = (*d2WDCdD2.value)[1];
        (*d2WDCdD2.value)[3] = 0.;

    }

    /*!
     * Set the current Hessians of the modified Danielsson function
     */
    void CHIPFoamStrainEnergy::setWDCHessians(){ setWDCHessians(false); }

    /*!
     * Set the previous Hessians of the modified Danielsson function
     */
    void CHIPFoamStrainEnergy::setPreviousWDCHessians(){ setWDCHessians(true); }

    /*!
     * Set the value of the gas energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setWG(bool isPrevious){

        auto phi0 = get_phi0();

        auto p0   = get_p0();

        auto gamma = get_gamma();

        const floatType *Je;

        const floatType *Jbar;

        SetDataStorageBase<floatType> WG;

        if(isPrevious){

            Je = get_previousJe();

            Jbar = get_previousJbar();

            WG = get_SetDataStorage_previousWG();

        }else{

            Je = get_Je();

            Jbar = get_Jbar();

            WG = get_SetDataStorage_WG();

        }

        auto Jm = (*Je) / (*Jbar);

        auto Jg = Jm * ( (*Jbar) - 1 + phi0 ) / phi0;

        if ( gas_isothermal_compression ){
            *WG.value = Jg - std::log(Jg) - 1;
        }else{
            *WG.value = Jg - 1./(gamma - 1)*(gamma - std::pow(Jg, 1 - gamma));
        }

        *WG.value *= p0 * phi0;

    }

    /*!
     * Set the current value of the gas energy
     */
    void CHIPFoamStrainEnergy::setWG(){ setWG(false); }

    /*!
     * Set the previous value of the gas energy
     */
    void CHIPFoamStrainEnergy::setPreviousWG(){ setWG(true); }

    /*!
     * Set the derivatives of the gas energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setWGDerivatives(bool isPrevious){

        auto phi0 = get_phi0();

        auto p0   = get_p0();

        auto gamma = get_gamma();

        const floatType *Je;

        const floatType *Jbar;

        const floatType *dJbardJe;

        SetDataStorageBase<floatVector> dWGdD;

        if(isPrevious){

            Je = get_previousJe();

            Jbar = get_previousJbar();

            dJbardJe = get_dPreviousJbardPreviousJe();

            dWGdD = get_SetDataStorage_dPreviousWGdPreviousD();

        }else{

            Je = get_Je();

            Jbar = get_Jbar();

            dJbardJe = get_dJbardJe();

            dWGdD = get_SetDataStorage_dWGdD();

        }

        auto Jm = (*Je) / (*Jbar);

        auto dJmdJe = 1. / (*Jbar) - (*Je) / ((*Jbar)*(*Jbar)) * (*dJbardJe);

        auto Jg = Jm * ( (*Jbar) - 1 + phi0 ) / phi0;

        auto dJgdJe = dJmdJe * ( (*Jbar) - 1 + phi0 ) / phi0 + Jm * (*dJbardJe) / phi0;

        dWGdD.zero(2);

        if ( gas_isothermal_compression ){
            (*dWGdD.value)[0] = dJgdJe * ( 1 - 1 / Jg );
        }else{
            (*dWGdD.value)[0] = dJgdJe * (1 - std::pow(Jg, -gamma));
        }

        (*dWGdD.value)[0] *= p0 * phi0;

    }

    /*!
     * Set the current derivatives of the gas energy
     */
    void CHIPFoamStrainEnergy::setWGDerivatives(){ setWGDerivatives(false); }

    /*!
     * Set the previous derivatives of the gas energy
     */
    void CHIPFoamStrainEnergy::setPreviousWGDerivatives(){ setWGDerivatives(true); }

    /*!
     * Set the Hessians of the gas energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setWGHessians(bool isPrevious){

        auto phi0 = get_phi0();

        auto p0   = get_p0();

        auto gamma = get_gamma();

        const floatType *Je;

        const floatType *Jbar;

        const floatType *dJbardJe;

        const floatType *d2JbardJe2;

        SetDataStorageBase<floatVector> d2WGdD2;

        if(isPrevious){

            Je = get_previousJe();

            Jbar = get_previousJbar();

            dJbardJe = get_dPreviousJbardPreviousJe();

            d2JbardJe2 = get_d2PreviousJbardPreviousJe2();

            d2WGdD2 = get_SetDataStorage_d2PreviousWGdPreviousD2();

        }else{

            Je = get_Je();

            Jbar = get_Jbar();

            dJbardJe = get_dJbardJe();

            d2JbardJe2 = get_d2JbardJe2();

            d2WGdD2 = get_SetDataStorage_d2WGdD2();

        }

        auto Jm = (*Je) / (*Jbar);

        auto dJmdJe = 1. / (*Jbar) - (*Je) / ((*Jbar)*(*Jbar)) * (*dJbardJe);

        auto d2JmdJe2 = -(*dJbardJe) / ((*Jbar)*(*Jbar)) - (*dJbardJe) / ((*Jbar)*(*Jbar)) + 2 * (*Je) / ((*Jbar)*(*Jbar)*(*Jbar)) * (*dJbardJe) * (*dJbardJe) - (*Je) / ((*Jbar)*(*Jbar)) * (*d2JbardJe2);

        auto Jg = Jm * ( (*Jbar) - 1 + phi0 ) / phi0;

        auto dJgdJe = dJmdJe * ( (*Jbar) - 1 + phi0 ) / phi0 + Jm * (*dJbardJe) / phi0;

        auto d2JgdJe2 = d2JmdJe2 * ( (*Jbar) - 1 + phi0 ) / phi0 + dJmdJe * (*dJbardJe) / phi0 + dJmdJe * (*dJbardJe) / phi0 + Jm * (*d2JbardJe2) / phi0;

        d2WGdD2.zero(4);

        if ( gas_isothermal_compression ){
            (*d2WGdD2.value)[0] = d2JgdJe2 * ( 1 - 1 / Jg ) + dJgdJe * dJgdJe / (Jg * Jg);
        }else{
            (*d2WGdD2.value)[0] = d2JgdJe2 * (1 - std::pow(Jg, -gamma)) + dJgdJe * dJgdJe * gamma * std::pow(Jg, -(gamma+1));
        }

        (*d2WGdD2.value)[0] *= p0 * phi0;
    }

    /*!
     * Set the current Hessians of the gas energy
     */
    void CHIPFoamStrainEnergy::setWGHessians(){ setWGHessians(false); }

    /*!
     * Set the previous Hessians of the gas energy
     */
    void CHIPFoamStrainEnergy::setPreviousWGHessians(){ setWGHessians(true); }

    /*!
     * Set the value of the parent material energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setWM(bool isPrevious){

        auto phi0 = get_phi0();

        auto K = get_K();

        const floatType *Je;

        const floatType *Jbar;

        SetDataStorageBase<floatType> WM;

        if(isPrevious){

            Je = get_previousJe();

            Jbar = get_previousJbar();

            WM = get_SetDataStorage_previousWM();

        }else{

            Je = get_Je();

            Jbar = get_Jbar();

            WM = get_SetDataStorage_WM();

        }

        auto Jm = (*Je) / (*Jbar);

        *WM.value = (1 - phi0) * K * ( Jm * std::log(Jm) - Jm + 1);

    }

    /*!
     * Set the current value of the parent material energy
     */
    void CHIPFoamStrainEnergy::setWM(){ setWM(false); }

    /*!
     * Set the previous value of the parent material energy
     */
    void CHIPFoamStrainEnergy::setPreviousWM(){ setWM(true); }

    /*!
     * Set the derivatives of the parent energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setWMDerivatives(bool isPrevious){

        auto phi0 = get_phi0();

        auto K = get_K();

        const floatType *Je;

        const floatType *Jbar;

        const floatType *dJbardJe;

        SetDataStorageBase<floatVector> dWMdD;

        if(isPrevious){

            Je = get_previousJe();

            Jbar = get_previousJbar();

            dJbardJe = get_dPreviousJbardPreviousJe();

            dWMdD = get_SetDataStorage_dPreviousWMdPreviousD();

        }else{

            Je = get_Je();

            Jbar = get_Jbar();

            dJbardJe = get_dJbardJe();

            dWMdD = get_SetDataStorage_dWMdD();

        }

        auto Jm = (*Je) / (*Jbar);

        auto dJmdJe = 1. / (*Jbar) - (*Je) / ((*Jbar)*(*Jbar)) * (*dJbardJe);

        dWMdD.zero(2);

        (*dWMdD.value)[0] = (1 - phi0) * K * std::log(Jm) * dJmdJe;

    }

    /*!
     * Set the current derivatives of the parent material energy
     */
    void CHIPFoamStrainEnergy::setWMDerivatives(){ setWMDerivatives(false); }

    /*!
     * Set the previous derivatives of the parent material energy
     */
    void CHIPFoamStrainEnergy::setPreviousWMDerivatives(){ setWMDerivatives(true); }

    /*!
     * Set the Hessians of the parent energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setWMHessians(bool isPrevious){

        auto phi0 = get_phi0();

        auto K = get_K();

        const floatType *Je;

        const floatType *Jbar;

        const floatType *dJbardJe;

        const floatType *d2JbardJe2;

        SetDataStorageBase<floatVector> d2WMdD2;

        if(isPrevious){

            Je = get_previousJe();

            Jbar = get_previousJbar();

            dJbardJe = get_dPreviousJbardPreviousJe();

            d2JbardJe2 = get_d2PreviousJbardPreviousJe2();

            d2WMdD2 = get_SetDataStorage_d2PreviousWMdPreviousD2();

        }else{

            Je = get_Je();

            Jbar = get_Jbar();

            dJbardJe = get_dJbardJe();

            d2JbardJe2 = get_d2JbardJe2();

            d2WMdD2 = get_SetDataStorage_d2WMdD2();

        }

        auto Jm = (*Je) / (*Jbar);

        auto dJmdJe = 1. / (*Jbar) - (*Je) / ((*Jbar)*(*Jbar)) * (*dJbardJe);

        auto d2JmdJe2 = -(*dJbardJe) / ((*Jbar)*(*Jbar)) - (*dJbardJe) / ((*Jbar)*(*Jbar)) + 2 * (*Je) / ((*Jbar)*(*Jbar)*(*Jbar)) * (*dJbardJe) * (*dJbardJe) - (*Je) / ((*Jbar)*(*Jbar)) * (*d2JbardJe2);

        d2WMdD2.zero(4);

        (*d2WMdD2.value)[0] = (1 - phi0) * K * dJmdJe * dJmdJe / Jm + (1 - phi0) * K * std::log(Jm) * d2JmdJe2;
    }

    /*!
     * Set the current Hessians of the parent material energy
     */
    void CHIPFoamStrainEnergy::setWMHessians(){ setWMHessians(false); }

    /*!
     * Set the previous Hessians of the parent material energy
     */
    void CHIPFoamStrainEnergy::setPreviousWMHessians(){ setWMHessians(true); }

    /*!
     * Set the value of the strain energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setStrainEnergy(const bool isPrevious){

        const floatType *WLB;

        const floatType *WDC;

        const floatType *WM;

        const floatType *WG;

        SetDataStorageBase<floatType> strainEnergy;

        if (isPrevious){

            WLB = get_previousWLB();

            WDC = get_previousWDC();

            WM = get_previousWM();

            WG = get_previousWG();

            strainEnergy = get_SetDataStorage_previousStrainEnergy();

        }
        else{

            WLB = get_WLB();

            WDC = get_WDC();

            WM = get_WM();

            WG = get_WG();

            strainEnergy = get_SetDataStorage_strainEnergy();

        }

        *strainEnergy.value = (*WLB) + (*WDC) + (*WM) + (*WG);

    }

    /*!
     * Set the derivatives of the strain energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setStrainEnergyJacobians(const bool isPrevious){

        constexpr unsigned int dim = 3;

        const floatVector *dWLBdD;

        const floatVector *dWDCdD;

        const floatVector *dWMdD;

        const floatVector *dWGdD;

        const floatVector *dJedFe;

        const floatVector *dIbar1dFe;

        SetDataStorageBase<floatVector> dStrainEnergydFe;

        if (isPrevious){

            dJedFe = get_dPreviousJedPreviousFe();

            dIbar1dFe = get_dPreviousIbar1dPreviousFe();

            dWLBdD = get_dPreviousWLBdPreviousD();
                                  
            dWDCdD = get_dPreviousWDCdPreviousD();

            dWMdD = get_dPreviousWMdPreviousD();

            dWGdD = get_dPreviousWGdPreviousD();

            dStrainEnergydFe = get_SetDataStorage_dPreviousStrainEnergydPreviousFe();

        }
        else{

            dJedFe = get_dJedFe();

            dIbar1dFe = get_dIbar1dFe();

            dWLBdD = get_dWLBdD();

            dWDCdD = get_dWDCdD();

            dWMdD = get_dWMdD();

            dWGdD = get_dWGdD();

            dStrainEnergydFe = get_SetDataStorage_dStrainEnergydFe();

        }

        dStrainEnergydFe.zero(dim*dim);

        for ( unsigned int iI = 0; iI < dim * dim; ++iI ){

            (*dStrainEnergydFe.value)[iI] += ((*dWLBdD)[0] + (*dWDCdD)[0] + (*dWMdD)[0] + (*dWGdD)[0]) * (*dJedFe)[iI];
            (*dStrainEnergydFe.value)[iI] += ((*dWLBdD)[1] + (*dWDCdD)[1] + (*dWMdD)[1] + (*dWGdD)[1]) * (*dIbar1dFe)[iI];

        }

    }

    /*!
     * Set the Hessians of the strain energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) value
     */
    void CHIPFoamStrainEnergy::setStrainEnergyHessians(const bool isPrevious){

        constexpr unsigned int dim = 3;

        const floatVector *dWLBdD;

        const floatVector *dWDCdD;

        const floatVector *dWMdD;

        const floatVector *dWGdD;

        const floatVector *dJedFe;

        const floatVector *dIbar1dFe;

        const floatVector *d2WLBdD2;

        const floatVector *d2WDCdD2;

        const floatVector *d2WMdD2;

        const floatVector *d2WGdD2;

        const floatVector *d2JedFe2;

        const floatVector *d2Ibar1dFe2;

        SetDataStorageBase<floatVector> d2StrainEnergydFe2;

        SetDataStorageBase<floatVector> d2StrainEnergydFedT;

        if (isPrevious){

            dJedFe = get_dPreviousJedPreviousFe();

            dIbar1dFe = get_dPreviousIbar1dPreviousFe();

            dWLBdD = get_dPreviousWLBdPreviousD();
                                  
            dWDCdD = get_dPreviousWDCdPreviousD();

            dWMdD = get_dPreviousWMdPreviousD();

            dWGdD = get_dPreviousWGdPreviousD();

            d2JedFe2 = get_d2PreviousJedPreviousFe2();

            d2Ibar1dFe2 = get_d2PreviousIbar1dPreviousFe2();

            d2WLBdD2 = get_d2PreviousWLBdPreviousD2();

            d2WDCdD2 = get_d2PreviousWDCdPreviousD2();

            d2WMdD2 = get_d2PreviousWMdPreviousD2();

            d2WGdD2 = get_d2PreviousWGdPreviousD2();

            d2StrainEnergydFe2 = get_SetDataStorage_d2PreviousStrainEnergydPreviousFe2();

            d2StrainEnergydFedT = get_SetDataStorage_d2PreviousStrainEnergydPreviousFedPreviousT();

        }
        else{

            dJedFe = get_dJedFe();

            dIbar1dFe = get_dIbar1dFe();

            dWLBdD = get_dWLBdD();

            dWDCdD = get_dWDCdD();

            dWMdD = get_dWMdD();

            dWGdD = get_dWGdD();

            d2JedFe2 = get_d2JedFe2();

            d2Ibar1dFe2 = get_d2Ibar1dFe2();

            d2WLBdD2 = get_d2WLBdD2();

            d2WDCdD2 = get_d2WDCdD2();

            d2WMdD2 = get_d2WMdD2();

            d2WGdD2 = get_d2WGdD2();

            d2StrainEnergydFe2 = get_SetDataStorage_d2StrainEnergydFe2();

            d2StrainEnergydFedT = get_SetDataStorage_d2StrainEnergydFedT();

        }

        d2StrainEnergydFe2.zero(dim*dim*dim*dim);
        d2StrainEnergydFedT.zero(dim*dim);

        for ( unsigned int iI = 0; iI < dim * dim; ++iI ){

            for ( unsigned int aA = 0; aA < dim * dim; ++aA ){

                (*d2StrainEnergydFe2.value)[dim * dim * iI + aA] += ((*d2WLBdD2)[0] + (*d2WDCdD2)[0] + (*d2WMdD2)[0] + (*d2WGdD2)[0]) * (*dJedFe)[iI] * (*dJedFe)[aA]
                                                                  + ((*d2WLBdD2)[1] + (*d2WDCdD2)[1] + (*d2WMdD2)[1] + (*d2WGdD2)[1]) * (*dJedFe)[iI] * (*dIbar1dFe)[aA]
                                                                  + ((*d2WLBdD2)[2] + (*d2WDCdD2)[2] + (*d2WMdD2)[2] + (*d2WGdD2)[2]) * (*dIbar1dFe)[iI] * (*dJedFe)[aA]
                                                                  + ((*d2WLBdD2)[3] + (*d2WDCdD2)[3] + (*d2WMdD2)[3] + (*d2WGdD2)[3]) * (*dIbar1dFe)[iI] * (*dIbar1dFe)[aA]
                                                                  + ((*dWLBdD)[0] + (*dWDCdD)[0] + (*dWMdD)[0] + (*dWGdD)[0]) * (*d2JedFe2)[dim * dim * iI + aA]
                                                                  + ((*dWLBdD)[1] + (*dWDCdD)[1] + (*dWMdD)[1] + (*dWGdD)[1]) * (*d2Ibar1dFe2)[dim * dim * iI + aA];

            }

        }

    }

}  // namespace tardigradeHydra
