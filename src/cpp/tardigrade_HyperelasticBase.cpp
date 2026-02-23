/**
 ******************************************************************************
 * \file tardigrade_HyperelasticBase.cpp
 ******************************************************************************
 * The base class for Hyperelastic materials
 ******************************************************************************
 */

#include "tardigrade_HyperelasticBase.h"

namespace tardigradeHydra {

        void HyperelasticBase::setFe(const bool isPrevious) {
            /*!
             * Set the value of the elastic deformation gradient
             *
             * \param isPrevious: Flag for whether to set the current (false) or previous (true) elastic deformation
             * gradient
             */

            constexpr unsigned int dim     = 3;
            constexpr unsigned int sot_dim = dim * dim;

            if (isPrevious) {
                auto previousFe = get_SetDataStorage_previousFe();

                *previousFe.value =
                    secondOrderTensor(hydra->deformation->get_previousConfigurations()->begin(),
                                      hydra->deformation->get_previousConfigurations()->begin() + sot_dim);

            } else {
                auto Fe = get_SetDataStorage_Fe();

                *Fe.value = secondOrderTensor(hydra->deformation->get_configurations()->begin(),
                                              hydra->deformation->get_configurations()->begin() + sot_dim);
            }
        }

        void HyperelasticBase::setFe() {
            /*!
             * Set the value of the elastic deformation gradient
             */

            setFe(false);
        }

        void HyperelasticBase::setPreviousFe() {
            /*!
             * Set the value of the previous elastic deformation gradient
             */

            setFe(true);
        }

        void HyperelasticBase::setFeDerivatives(const bool isPrevious) {
            /*!
             * Set the value of the derivatives of the elastic strain
             */

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

        void HyperelasticBase::setdFedF() {
            /*!
             * Set the value of the derivative of the elastic deformation gradient w.r.t. the deformation gradient
             */

            setFeDerivatives(false);
        }

        void HyperelasticBase::setdFedFn() {
            /*!
             * Set the value of the derivative of the elastic deformation gradient w.r.t. the sub-deformation gradients
             */

            setFeDerivatives(false);
        }

        void HyperelasticBase::setPreviousdFedF() {
            /*!
             * Set the value of the previous derivative of the elastic deformation gradient w.r.t. the deformation
             * gradient
             */

            setFeDerivatives(true);
        }

        void HyperelasticBase::setPreviousdFedFn() {
            /*!
             * Set the value of the previous derivative of the elastic deformation gradient w.r.t. the sub-deformation
             * gradients
             */

            setFeDerivatives(true);
        }

}
