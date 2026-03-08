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
    void CHIPFoamStrainEnergy::setdJedFe(bool isPrevious){

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

        secondOrderTensor invFe(dim*dim,0);
        Eigen::Map<const Eigen::Matrix<floatType, -1, -1, Eigen::RowMajor>> Femat(
            Fe->data(), dim, dim);  // TODO: Change this to a constant size when possible
        Eigen::Map<Eigen::Matrix<floatType, -1, -1, Eigen::RowMajor>> invFemat(
            invFe.data(), dim, dim);  // TODO: Change this to a constant size when possible
        invFemat = Femat.inverse().eval();

        dJedFe.zero(dim*dim);

        for ( unsigned int i = 0; i < dim; ++i){
            for ( unsigned int I = 0; I < dim; ++I){
                (*dJedFe.value)[dim*i+I] = (*Je)*invFe[dim*I+i];
            }
        }

    }

    /*!
     * Compute the derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setdJedFe(){ setdJedFe(false); }

    /*!
     * Compute the previous derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setdPreviousJedPreviousFe(){ setdJedFe(true); }

    /*!
     * Compute the second derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true)
     *     value
     */
    void CHIPFoamStrainEnergy::setd2JedFe2(bool isPrevious){

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

        d2JedFe2.zero(dim*dim*dim*dim);

        for ( unsigned int i = 0; i < dim; ++i){
            for ( unsigned int I = 0; I < dim; ++I){
                for ( unsigned int j = 0; j < dim; ++j){
                    for ( unsigned int J = 0; J < dim; ++J){
                        (*d2JedFe2.value)[dim*dim*dim*i+dim*dim*I+dim*j+J] = ((*dJedFe)[dim*j+J]*(*dJedFe)[dim*i+I] - (*dJedFe)[dim*j+I]*(*dJedFe)[dim*i+J])/(*Je);
                    }
                }
            }
        }
    } 

    /*!
     * Compute the derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setd2JedFe2(){ setd2JedFe2(false); }

    /*!
     * Compute the previous derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setd2PreviousJedPreviousFe2(){ setd2JedFe2(true); }

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
    void CHIPFoamStrainEnergy::setdIbar1dFe(bool isPrevious){

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

        dIbar1dFe.zero(dim*dim);

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int I = 0; I < dim; ++I) {
                (*dIbar1dFe.value)[dim*i+I] = 2*((*Fe)[dim * i + I] * std::pow((*Je),1./3) - (*Ibar1) * (*dJedFe)[dim * i + I] / 3.) / (*Je);
            }
        }

    }

    /*!
     * Compute the derivative of the first invariant of the isochoric elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setdIbar1dFe(){ setdIbar1dFe(false); }

    /*!
     * Compute the previous derivative of the first invariant of the isochoric elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setdPreviousIbar1dPreviousFe(){ setdIbar1dFe(true); }

    /*!
     * Compute the second derivative of the first invariant of the isochoric elastic deformation
     * with respect to the elastic deformation
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true)
     *     value
     */
    void CHIPFoamStrainEnergy::setd2Ibar1dFe2(bool isPrevious){

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

        d2Ibar1dFe2.zero(dim*dim*dim*dim);

        auto Je_23 = std::pow((*Je),2./3.);

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int I = 0; I < dim; ++I) {
                (*d2Ibar1dFe2.value)[dim*dim*dim*i+dim*dim*I+dim*i+I] += 2 / Je_23;
                for ( unsigned int j = 0; j < dim; ++j){

                    for ( unsigned int J = 0; J < dim; ++J){
                        (*d2Ibar1dFe2.value)[dim*dim*dim*i+dim*dim*I+dim*j+J] += 2./3. * ((*Fe)[dim * i + I] / Je_23 * (*dJedFe)[dim * j + J] - (*dIbar1dFe)[dim*j+J] * (*dJedFe)[dim * i + I] - (*Ibar1) * (*d2JedFe2)[dim*dim*dim*i+dim*dim*I+dim*j+J]) / (*Je) - (*dIbar1dFe)[dim*i+I] / (*Je) * (*dJedFe)[dim*j+J];
                    }
                }
            }
        }

    }

    /*!
     * Compute the derivative of the first invariant of the isochoric elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setd2Ibar1dFe2(){ setd2Ibar1dFe2(false); }

    /*!
     * Compute the previous derivative of the first invariant of the isochoric elastic deformation
     * with respect to the elastic deformation
     */
    void CHIPFoamStrainEnergy::setd2PreviousIbar1dPreviousFe2(){ setd2Ibar1dFe2(true); }

}  // namespace tardigradeHydra
