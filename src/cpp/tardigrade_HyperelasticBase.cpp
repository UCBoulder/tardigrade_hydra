/**
 ******************************************************************************
 * \file tardigrade_HyperelasticBase.cpp
 ******************************************************************************
 * The base class for Hyperelastic materials
 ******************************************************************************
 */

#include "tardigrade_HyperelasticBase.h"

namespace tardigradeHydra {

    /*!
     * Set the value of the elastic deformation gradient
     *
     * \param isPrevious: Flag for whether to set the current (false) or previous (true) elastic deformation
     * gradient
     */
    void HyperelasticBase::setFe(const bool isPrevious) {
        if (isPrevious) {
            auto previousFe = get_SetDataStorage_previousFe();

            *previousFe.value = secondOrderTensor(hydra->deformation->get_previousConfigurations()->begin(),
                                                  hydra->deformation->get_previousConfigurations()->begin() + dimension * dimension);

        } else {
            auto Fe = get_SetDataStorage_Fe();

            *Fe.value = secondOrderTensor(hydra->deformation->get_configurations()->begin(),
                                          hydra->deformation->get_configurations()->begin() + dimension * dimension);
        }
    }

    /*!
     * Set the value of the elastic deformation gradient
     */
    void HyperelasticBase::setFe() { setFe(false); }

    /*!
     * Set the value of the previous elastic deformation gradient
     */
    void HyperelasticBase::setPreviousFe() { setFe(true); }

    /*!
     * Set the value of the derivatives of the elastic strain
     *
     * \param isPrevious: Whether this is the previous or current elastic deforamtion
     */
    void HyperelasticBase::setFeDerivatives(const bool isPrevious) {
        if (isPrevious) {
            auto previousdFedF = get_SetDataStorage_previousdFedF();

            auto previousdFedFn = get_SetDataStorage_previousdFedFn();

            *previousdFedF.value = *hydra->deformation->get_previousdF1dF();

            *previousdFedFn.value = *hydra->deformation->get_previousdF1dFn();

        } else {
            auto dFedF = get_SetDataStorage_dFedF();

            auto dFedFn = get_SetDataStorage_dFedFn();

            *dFedF.value = *hydra->deformation->get_dF1dF();

            *dFedFn.value = *hydra->deformation->get_dF1dFn();
        }
    }

    /*!
     * Set the value of the derivative of the elastic deformation gradient w.r.t. the deformation gradient
     */
    void HyperelasticBase::setdFedF() { setFeDerivatives(false); }

    /*!
     * Set the value of the derivative of the elastic deformation gradient w.r.t. the sub-deformation gradients
     */
    void HyperelasticBase::setdFedFn() { setFeDerivatives(false); }

    /*!
     * Set the value of the previous derivative of the elastic deformation gradient w.r.t. the deformation
     * gradient
     */
    void HyperelasticBase::setPreviousdFedF() { setFeDerivatives(true); }

    /*!
     * Set the value of the previous derivative of the elastic deformation gradient w.r.t. the sub-deformation
     * gradients
     */
    void HyperelasticBase::setPreviousdFedFn() { setFeDerivatives(true); }

    /*!
     * Set the strain energy
     *
     * \param isPrevious: A flag for whether to set the current (false) or previous (true) strain energy
     */
    void HyperelasticBase::setStrainEnergy(const bool isPrevious) {
        throw std::runtime_error("The strain energy calculation must be defined");
    }

    /*!
     * Compute the Jacobian of the elastic deformation
     *
     * \param isPrevious: A flag for of the current (false) or previous (true) value should be computed
     */
    void HyperelasticBase::setJe(bool isPrevious) {
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
    void HyperelasticBase::setJe() { setJe(false); }

    /*!
     * Compute the previous Jacobian of the elastic deformation
     */
    void HyperelasticBase::setPreviousJe() { setJe(true); }

    /*!
     * Compute the derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true)
     *     value
     */
    void HyperelasticBase::setdJedFe(bool isPrevious) {
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
    void HyperelasticBase::setdJedFe() { setdJedFe(false); }

    /*!
     * Compute the previous derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     */
    void HyperelasticBase::setdPreviousJedPreviousFe() { setdJedFe(true); }

    /*!
     * Compute the second derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     *
     * \param isPrevious: Whether to compute the current (false) or previous (true)
     *     value
     */
    void HyperelasticBase::setd2JedFe2(bool isPrevious) {
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
    void HyperelasticBase::setd2JedFe2() { setd2JedFe2(false); }

    /*!
     * Compute the previous derivative of the Jacobian of the elastic deformation
     * with respect to the elastic deformation
     */
    void HyperelasticBase::setd2PreviousJedPreviousFe2() { setd2JedFe2(true); }

    /*!
     * Compute the current strain energy
     */
    void HyperelasticBase::setStrainEnergy() { setStrainEnergy(false); }

    /*!
     * Compute the previous strain energy
     */
    void HyperelasticBase::setPreviousStrainEnergy() { setStrainEnergy(true); }

    /*!
     * Set the Jacobians of the strain energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) Jacobians of the strain energy
     */
    void HyperelasticBase::setStrainEnergyJacobians(const bool isPrevious) {
        throw std::runtime_error("The Jacobians of the strain energy must be defined");
    }

    /*!
     * Set the Jacobians of the strain energy
     */
    void HyperelasticBase::setStrainEnergyJacobians() { setStrainEnergyJacobians(false); }

    /*!
     * Set the previous Jacobians of the strain energy
     */
    void HyperelasticBase::setPreviousStrainEnergyJacobians() { setStrainEnergyJacobians(true); }

    /*!
     * Set the Hessians of the strain energy
     *
     * \param isPrevious: Whether to set the current (false) or previous (true) Hessians of the strain energy
     */
    void HyperelasticBase::setStrainEnergyHessians(const bool isPrevious) {
        throw std::runtime_error("The Hessians of the strain energy must be defined");
    }

    /*!
     * Set the Hessians of the strain energy
     */
    void HyperelasticBase::setStrainEnergyHessians() { setStrainEnergyHessians(false); }

    /*!
     * Set the previous Hessians of the strain energy
     */
    void HyperelasticBase::setPreviousStrainEnergyHessians() { setStrainEnergyHessians(true); }

    /*!
     * Set the Cauchy stress from the strain-energy function
     *
     * \param isPrevious: Flag for whether to set the current (false) or previous (true) strain energy function
     */
    void HyperelasticBase::setCauchyStress(const bool isPrevious) {
        const secondOrderTensor *Fe;

        const floatType *Je;

        const secondOrderTensor *dStrainEnergydFe;

        SetDataStorageBase<secondOrderTensor> cauchyStress;

        if (isPrevious) {
            Fe = get_previousFe();

            Je = get_previousJe();

            dStrainEnergydFe = get_dPreviousStrainEnergydPreviousFe();

            cauchyStress = get_SetDataStorage_previousCauchyStress();

        } else {
            Fe = get_Fe();

            Je = get_Je();

            dStrainEnergydFe = get_dStrainEnergydFe();

            cauchyStress = get_SetDataStorage_cauchyStress();
        }

        TARDIGRADE_ERROR_TOOLS_CHECK((std::end(*Fe) - std::begin(*Fe)) == dimension * dimension,
                                     "The elastic deformation must have a size of " + std::to_string(dimension * dimension));

        TARDIGRADE_ERROR_TOOLS_CHECK(
            (std::end(*dStrainEnergydFe) - std::begin(*dStrainEnergydFe)) == dimension * dimension,
            "The derivative of the strain energy with respect to the elastic deformation must have a size of " +
                std::to_string(dimension * dimension));

        cauchyStress.zero(dimension * dimension);

        for (unsigned int i = 0; i < dimension; ++i) {
            for (unsigned int j = 0; j < dimension; ++j) {
                for (unsigned int I = 0; I < dimension; ++I) {
                    (*cauchyStress.value)[dimension * i + j] += (*dStrainEnergydFe)[dimension * i + I] * (*Fe)[dimension * j + I];
                }

                (*cauchyStress.value)[dimension * i + j] /= *Je;
            }
        }
    }

    /*!
     * Set the current Cauchy stress
     */
    void HyperelasticBase::setCauchyStress() { setCauchyStress(false); }

    /*!
     * Set the previous Cauchy stress
     */
    void HyperelasticBase::setPreviousCauchyStress() { setCauchyStress(true); }

    /*!
     * Set the Cauchy stress jacobians
     *
     * \param isPrevious: Whether we should compute the previous or current Jacobians
     */
    void HyperelasticBase::setCauchyStressJacobians(const bool isPrevious) {
        const secondOrderTensor *Fe;

        const floatType *Je;

        const secondOrderTensor *dStrainEnergydFe;

        const fourthOrderTensor *dFedF;

        const floatVector *dFedFn;

        const fourthOrderTensor *d2StrainEnergydFe2;

        const secondOrderTensor *d2StrainEnergydFedT;

        SetDataStorageBase<fourthOrderTensor> dCauchyStressdF;

        SetDataStorageBase<floatVector> dCauchyStressdFn;

        SetDataStorageBase<secondOrderTensor> dCauchyStressdT;

        if (isPrevious) {
            Fe = get_previousFe();

            Je = get_previousJe();

            dStrainEnergydFe = get_dPreviousStrainEnergydPreviousFe();

            dFedF = get_previousdFedF();

            dFedFn = get_previousdFedFn();

            d2StrainEnergydFe2 = get_d2PreviousStrainEnergydPreviousFe2();

            d2StrainEnergydFedT = get_d2PreviousStrainEnergydPreviousFedPreviousT();

            dCauchyStressdF = get_SetDataStorage_dPreviousCauchyStressdPreviousF();

            dCauchyStressdFn = get_SetDataStorage_dPreviousCauchyStressdPreviousFn();

            dCauchyStressdT = get_SetDataStorage_dPreviousCauchyStressdPreviousT();

        } else {
            Fe = get_Fe();

            Je = get_Je();

            dStrainEnergydFe = get_dStrainEnergydFe();

            dFedF = get_dFedF();

            dFedFn = get_dFedFn();

            d2StrainEnergydFe2 = get_d2StrainEnergydFe2();

            d2StrainEnergydFedT = get_d2StrainEnergydFedT();

            dCauchyStressdF = get_SetDataStorage_dCauchyStressdF();

            dCauchyStressdFn = get_SetDataStorage_dCauchyStressdFn();

            dCauchyStressdT = get_SetDataStorage_dCauchyStressdT();
        }

        TARDIGRADE_ERROR_TOOLS_CHECK((std::end(*Fe) - std::begin(*Fe)) == dimension * dimension,
                                     "The elastic deformation must have a size of " + std::to_string(dimension * dimension));

        TARDIGRADE_ERROR_TOOLS_CHECK(
            (std::end(*dStrainEnergydFe) - std::begin(*dStrainEnergydFe)) == dimension * dimension,
            "The derivative of the strain energy with respect to the elastic deformation must have a size of " +
                std::to_string(dimension * dimension));

        dCauchyStressdF.zero(dimension * dimension * dimension * dimension);

        dCauchyStressdFn.zero(dimension * dimension * dimension * dimension * (hydra->getNumConfigurations() - 1));

        dCauchyStressdT.zero(dimension * dimension);

        secondOrderTensor invFe(dimension * dimension, 0);

        Eigen::Map<const Eigen::Matrix<floatType, dimension, dimension, Eigen::RowMajor>> Femat(
            Fe->data(), dimension, dimension);
        Eigen::Map<Eigen::Matrix<floatType, dimension, dimension, Eigen::RowMajor>> invFemat(
            invFe.data(), dimension, dimension);
        invFemat = Femat.inverse().eval();

        secondOrderTensor JecauchyStress(dimension * dimension, 0);
        secondOrderTensor invFedFedF(dimension * dimension, 0);
        secondOrderTensor invFedFedFn(dimension * dimension * (hydra->getNumConfigurations() - 1), 0);
        fourthOrderTensor d2StrainEnergydFedF(dimension * dimension * dimension * dimension, 0);
        fourthOrderTensor d2StrainEnergydFedFn(dimension * dimension * dimension * dimension * (hydra->getNumConfigurations() - 1), 0);

        for (unsigned int iI = 0; iI < dimension * dimension; ++iI) {
            for (unsigned int aA = 0; aA < dimension * dimension; ++aA) {
                for (unsigned int bB = 0; bB < dimension * dimension; ++bB) {
                    d2StrainEnergydFedF[dimension * dimension * iI + bB] +=
                        (*d2StrainEnergydFe2)[dimension * dimension * iI + aA] * (*dFedF)[dimension * dimension * aA + bB];
                }
                for (unsigned int bB = 0; bB < dimension * dimension * (hydra->getNumConfigurations() - 1); ++bB) {
                    d2StrainEnergydFedFn[dimension * dimension * (hydra->getNumConfigurations() - 1) * iI + bB] +=
                        (*d2StrainEnergydFe2)[dimension * dimension * iI + aA] *
                        (*dFedFn)[dimension * dimension * (hydra->getNumConfigurations() - 1) * aA + bB];
                }
            }
        }

        for (unsigned int i = 0; i < dimension; ++i) {
            for (unsigned int I = 0; I < dimension; ++I) {
                for (unsigned int aA = 0; aA < dimension * dimension; ++aA) {
                    invFedFedF[aA] += invFe[dimension * I + i] * (*dFedF)[dimension * dimension * dimension * i + dimension * dimension * I + aA];
                }

                for (unsigned int aA = 0; aA < dimension * dimension * (hydra->getNumConfigurations() - 1); ++aA) {
                    invFedFedFn[aA] +=
                        invFe[dimension * I + i] * (*dFedFn)[dimension * dimension * dimension * (hydra->getNumConfigurations() - 1) * i +
                                                       dimension * dimension * (hydra->getNumConfigurations() - 1) * I + aA];
                }
            }
        }

        for (unsigned int i = 0; i < dimension; ++i) {
            for (unsigned int j = 0; j < dimension; ++j) {
                for (unsigned int I = 0; I < dimension; ++I) {
                    (*dCauchyStressdT.value)[dimension * i + j] += (*d2StrainEnergydFedT)[dimension * i + I] * (*Fe)[dimension * j + I];
                    JecauchyStress[dimension * i + j] += (*dStrainEnergydFe)[dimension * i + I] * (*Fe)[dimension * j + I];
                }

                (*dCauchyStressdT.value)[dim * i + j] /= *Je;
            }
        }

        for (unsigned int i = 0; i < dimension; ++i) {
            for (unsigned int j = 0; j < dimension; ++j) {
                for (unsigned int I = 0; I < dimension; ++I) {
                    for (unsigned int aA = 0; aA < dimension * dimension; ++aA) {
                        (*dCauchyStressdF.value)[dimension * dimension * dimension * i + dimension * dimension * j + aA] +=
                            d2StrainEnergydFedF[dimension * dimension * dimension * i + dimension * dimension * I + aA] * (*Fe)[dimension * j + I] +
                            (*dStrainEnergydFe)[dimension * i + I] * (*dFedF)[dimension * dimension * dimension * j + dimension * dimension * I + aA];
                    }

                    for (unsigned int aA = 0; aA < dimension * dimension * (hydra->getNumConfigurations() - 1); ++aA) {
                        (*dCauchyStressdFn.value)[dimension * dimension * dimension * (hydra->getNumConfigurations() - 1) * i +
                                                  dimension * dimension * (hydra->getNumConfigurations() - 1) * j + aA] +=
                            d2StrainEnergydFedFn[dimension * dimension * dimension * (hydra->getNumConfigurations() - 1) * i +
                                                 dimension * dimension * (hydra->getNumConfigurations() - 1) * I + aA] *
                                (*Fe)[dimension * j + I] +
                            (*dStrainEnergydFe)[dimension * i + I] *
                                (*dFedFn)[dimension * dimension * dimension * (hydra->getNumConfigurations() - 1) * j +
                                          dimension * dimension * (hydra->getNumConfigurations() - 1) * I + aA];
                    }
                }
            }
        }

        for (unsigned int ij = 0; ij < dimension * dimension; ++ij) {
            for (unsigned int aA = 0; aA < dimension * dimension; ++aA) {
                (*dCauchyStressdF.value)[dimension * dimension * ij + aA] -= JecauchyStress[ij] * invFedFedF[aA];
                (*dCauchyStressdF.value)[dimension * dimension * ij + aA] /= *Je;
            }

            for (unsigned int aA = 0; aA < dimension * dimension * (hydra->getNumConfigurations() - 1); ++aA) {
                (*dCauchyStressdFn.value)[dimension * dimension * (hydra->getNumConfigurations() - 1) * ij + aA] -=
                    JecauchyStress[ij] * invFedFedFn[aA];
                (*dCauchyStressdFn.value)[dimension * dimension * (hydra->getNumConfigurations() - 1) * ij + aA] /= *Je;
            }
        }
    }

    /*!
     * Set the Jacobian of the Cauchy stress w.r.t. the deformation gradient
     */
    void HyperelasticBase::setdCauchyStressdF() { setCauchyStressJacobians(false); }

    /*!
     * Set the Jacobian of the Cauchy stress w.r.t. the temperature
     */
    void HyperelasticBase::setdCauchyStressdT() { setCauchyStressJacobians(false); }

    /*!
     * Set the Jacobian of the Cauchy stress w.r.t. the other deformation gradients
     */
    void HyperelasticBase::setdCauchyStressdFn() { setCauchyStressJacobians(false); }

    /*!
     * Set the previous Jacobian of the Cauchy stress w.r.t. the deformation gradient
     */
    void HyperelasticBase::setdPreviousCauchyStressdPreviousF() { setCauchyStressJacobians(true); }

    /*!
     * Set the previous Jacobian of the Cauchy stress w.r.t. the temperature
     */
    void HyperelasticBase::setdPreviousCauchyStressdPreviousT() { setCauchyStressJacobians(true); }

    /*!
     * Set the previous Jacobian of the Cauchy stress w.r.t. the other deformation gradients
     */
    void HyperelasticBase::setdPreviousCauchyStressdPreviousFn() { setCauchyStressJacobians(true); }

    void HyperelasticBase::setStress() {
        /*!
         * Set the stress
         *
         * Currently uses the Cauchy stress
         */

        setCauchyStress(false);

        auto stress = get_SetDataStorage_stress();

        *stress.value = *get_cauchyStress();
    }

    void HyperelasticBase::setPreviousStress() {
        /*!
         * Set the previous stress
         *
         * Currently uses the Cauchy stress
         */

        setCauchyStress(true);

        auto previousStress = get_SetDataStorage_previousStress();

        *previousStress.value = *get_previousCauchyStress();
    }

    void HyperelasticBase::setResidual() {
        /*!
         * Set the residual value
         */

        auto residual = get_SetDataStorage_residual();

        const secondOrderTensor *stress = getStress();

        TARDIGRADE_ERROR_TOOLS_CATCH(*residual.value = *stress - *hydra->getStress());
    }

    void HyperelasticBase::setJacobian() {
        /*!
         * Set the Jacobian value
         */

        auto num_configs = hydra->getNumConfigurations();

        auto num_unknown_config_vars = (num_configs - 1) * dimension * dimension;

        auto num_unknowns = hydra->getNumUnknowns();

        // Form the Jacobian
        auto jacobian = get_SetDataStorage_jacobian();

        jacobian.zero(dimension * dimension * num_unknowns);

        for (unsigned int i = 0; i < dimension; i++) {
            for (unsigned int j = 0; j < dimension; j++) {
                (*jacobian.value)[num_unknowns * dimension * i + num_unknowns * j + dimension * i + j] = -1;

                for (unsigned int I = 0; I < num_unknown_config_vars; I++) {
                    (*jacobian.value)[num_unknowns * dimension * i + num_unknowns * j + getStress()->size() + I] =
                        (*get_dCauchyStressdFn())[dimension * num_unknown_config_vars * i + num_unknown_config_vars * j + I];
                }
            }
        }
    }

    void HyperelasticBase::setdRdT() {
        /*!
         * Set the derivative of the residual w.r.t. the temperature
         */

        auto dRdT = get_SetDataStorage_dRdT();

        *dRdT.value = *get_dCauchyStressdT();
    }

    void HyperelasticBase::setdRdF() {
        /*!
         * Set the derivative of the residual w.r.t. the deformation gradient
         */

        auto dRdF = get_SetDataStorage_dRdF();

        *dRdF.value = *get_dCauchyStressdF();
    }

}  // namespace tardigradeHydra
