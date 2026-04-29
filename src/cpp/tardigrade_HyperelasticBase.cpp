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
        constexpr unsigned int dim     = 3;
        constexpr unsigned int sot_dim = dim * dim;

        if (isPrevious) {
            auto previousFe = get_SetDataStorage_previousFe();

            *previousFe.value = secondOrderTensor(hydra->deformation->get_previousConfigurations()->begin(),
                                                  hydra->deformation->get_previousConfigurations()->begin() + sot_dim);

        } else {
            auto Fe = get_SetDataStorage_Fe();

            *Fe.value = secondOrderTensor(hydra->deformation->get_configurations()->begin(),
                                          hydra->deformation->get_configurations()->begin() + sot_dim);
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
        auto dim = hydra->deformation->dimension;

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

        TARDIGRADE_ERROR_TOOLS_CHECK((std::end(*Fe) - std::begin(*Fe)) == dim * dim,
                                     "The elastic deformation must have a size of " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK(
            (std::end(*dStrainEnergydFe) - std::begin(*dStrainEnergydFe)) == dim * dim,
            "The derivative of the strain energy with respect to the elastic deformation must have a size of " +
                std::to_string(dim * dim));

        cauchyStress.zero(dim * dim);

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int I = 0; I < dim; ++I) {
                    (*cauchyStress.value)[dim * i + j] += (*dStrainEnergydFe)[dim * i + I] * (*Fe)[dim * j + I];
                }

                (*cauchyStress.value)[dim * i + j] /= *Je;
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
        auto dim = hydra->deformation->dimension;

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

        TARDIGRADE_ERROR_TOOLS_CHECK((std::end(*Fe) - std::begin(*Fe)) == dim * dim,
                                     "The elastic deformation must have a size of " + std::to_string(dim * dim));

        TARDIGRADE_ERROR_TOOLS_CHECK(
            (std::end(*dStrainEnergydFe) - std::begin(*dStrainEnergydFe)) == dim * dim,
            "The derivative of the strain energy with respect to the elastic deformation must have a size of " +
                std::to_string(dim * dim));

        dCauchyStressdF.zero(dim * dim * dim * dim);

        dCauchyStressdFn.zero(dim * dim * dim * dim * (hydra->getNumConfigurations() - 1));

        dCauchyStressdT.zero(dim * dim);

        secondOrderTensor invFe(dim * dim, 0);

        Eigen::Map<const Eigen::Matrix<floatType, -1, -1, Eigen::RowMajor>> Femat(
            Fe->data(), dim, dim);  // TODO: Change this to a constant size when possible
        Eigen::Map<Eigen::Matrix<floatType, -1, -1, Eigen::RowMajor>> invFemat(
            invFe.data(), dim, dim);  // TODO: Change this to a constant size when possible
        invFemat = Femat.inverse().eval();

        secondOrderTensor JecauchyStress(dim * dim, 0);
        secondOrderTensor invFedFedF(dim * dim, 0);
        secondOrderTensor invFedFedFn(dim * dim * (hydra->getNumConfigurations() - 1), 0);
        fourthOrderTensor d2StrainEnergydFedF(dim * dim * dim * dim, 0);
        fourthOrderTensor d2StrainEnergydFedFn(dim * dim * dim * dim * (hydra->getNumConfigurations() - 1), 0);

        for (unsigned int iI = 0; iI < dim * dim; ++iI) {
            for (unsigned int aA = 0; aA < dim * dim; ++aA) {
                for (unsigned int bB = 0; bB < dim * dim; ++bB) {
                    d2StrainEnergydFedF[dim * dim * iI + bB] +=
                        (*d2StrainEnergydFe2)[dim * dim * iI + aA] * (*dFedF)[dim * dim * aA + bB];
                }
                for (unsigned int bB = 0; bB < dim * dim * (hydra->getNumConfigurations() - 1); ++bB) {
                    d2StrainEnergydFedFn[dim * dim * (hydra->getNumConfigurations() - 1) * iI + bB] +=
                        (*d2StrainEnergydFe2)[dim * dim * iI + aA] *
                        (*dFedFn)[dim * dim * (hydra->getNumConfigurations() - 1) * aA + bB];
                }
            }
        }

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int I = 0; I < dim; ++I) {
                for (unsigned int aA = 0; aA < dim * dim; ++aA) {
                    invFedFedF[aA] += invFe[dim * I + i] * (*dFedF)[dim * dim * dim * i + dim * dim * I + aA];
                }

                for (unsigned int aA = 0; aA < dim * dim * (hydra->getNumConfigurations() - 1); ++aA) {
                    invFedFedFn[aA] +=
                        invFe[dim * I + i] * (*dFedFn)[dim * dim * dim * (hydra->getNumConfigurations() - 1) * i +
                                                       dim * dim * (hydra->getNumConfigurations() - 1) * I + aA];
                }
            }
        }

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int I = 0; I < dim; ++I) {
                    (*dCauchyStressdT.value)[dim * i + j] += (*d2StrainEnergydFedT)[dim * i + I] * (*Fe)[dim * j + I];
                    JecauchyStress[dim * i + j] += (*dStrainEnergydFe)[dim * i + I] * (*Fe)[dim * j + I];
                }

                (*dCauchyStressdT.value)[dim * i + j] /= *Je;
            }
        }

        for (unsigned int i = 0; i < dim; ++i) {
            for (unsigned int j = 0; j < dim; ++j) {
                for (unsigned int I = 0; I < dim; ++I) {
                    for (unsigned int aA = 0; aA < dim * dim; ++aA) {
                        (*dCauchyStressdF.value)[dim * dim * dim * i + dim * dim * j + aA] +=
                            d2StrainEnergydFedF[dim * dim * dim * i + dim * dim * I + aA] * (*Fe)[dim * j + I] +
                            (*dStrainEnergydFe)[dim * i + I] * (*dFedF)[dim * dim * dim * j + dim * dim * I + aA];
                    }

                    for (unsigned int aA = 0; aA < dim * dim * (hydra->getNumConfigurations() - 1); ++aA) {
                        (*dCauchyStressdFn.value)[dim * dim * dim * (hydra->getNumConfigurations() - 1) * i +
                                                  dim * dim * (hydra->getNumConfigurations() - 1) * j + aA] +=
                            d2StrainEnergydFedFn[dim * dim * dim * (hydra->getNumConfigurations() - 1) * i +
                                                 dim * dim * (hydra->getNumConfigurations() - 1) * I + aA] *
                                (*Fe)[dim * j + I] +
                            (*dStrainEnergydFe)[dim * i + I] *
                                (*dFedFn)[dim * dim * dim * (hydra->getNumConfigurations() - 1) * j +
                                          dim * dim * (hydra->getNumConfigurations() - 1) * I + aA];
                    }
                }
            }
        }

        for (unsigned int ij = 0; ij < dim * dim; ++ij) {
            for (unsigned int aA = 0; aA < dim * dim; ++aA) {
                (*dCauchyStressdF.value)[dim * dim * ij + aA] -= JecauchyStress[ij] * invFedFedF[aA];
                (*dCauchyStressdF.value)[dim * dim * ij + aA] /= *Je;
            }

            for (unsigned int aA = 0; aA < dim * dim * (hydra->getNumConfigurations() - 1); ++aA) {
                (*dCauchyStressdFn.value)[dim * dim * (hydra->getNumConfigurations() - 1) * ij + aA] -=
                    JecauchyStress[ij] * invFedFedFn[aA];
                (*dCauchyStressdFn.value)[dim * dim * (hydra->getNumConfigurations() - 1) * ij + aA] /= *Je;
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

        auto dim = hydra->getDimension();

        auto sot_dim = hydra->getSOTDimension();

        auto num_configs = hydra->getNumConfigurations();

        auto num_unknown_config_vars = (num_configs - 1) * sot_dim;

        auto num_unknowns = hydra->getNumUnknowns();

        // Form the Jacobian
        auto jacobian = get_SetDataStorage_jacobian();

        jacobian.zero(sot_dim * num_unknowns);

        for (unsigned int i = 0; i < dim; i++) {
            for (unsigned int j = 0; j < dim; j++) {
                (*jacobian.value)[num_unknowns * dim * i + num_unknowns * j + dim * i + j] = -1;

                for (unsigned int I = 0; I < num_unknown_config_vars; I++) {
                    (*jacobian.value)[num_unknowns * dim * i + num_unknowns * j + getStress()->size() + I] =
                        (*get_dCauchyStressdFn())[dim * num_unknown_config_vars * i + num_unknown_config_vars * j + I];
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
