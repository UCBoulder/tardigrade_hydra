/**
 ******************************************************************************
 * \file tardigrade_hydraDOFVelocityGradientDeformation.cpp
 ******************************************************************************
 * An implementation of a deformation where the velocity gradient is in the
 * additionalDOF vector.
 ******************************************************************************
 */

#include <tardigrade_constitutive_tools.h>
#include <tardigrade_hydraDOFVelocityGradientDeformation.h>

namespace tardigradeHydra {

    namespace dofVelocityGradientDeformation {

        void residual::decomposeParameters(const floatType *parameters, const unsigned int parameters_size) {
            /*!
             * Decompose the parameter vector
             *
             * \param *parameters: The starting iterator of the parameter vector
             * \param parameters_size: The size of the parameter vector
             */

            TARDIGRADE_ERROR_TOOLS_CHECK(parameters_size == 2, "The parameter vector must have a size of 2")

            *get_SetDataStorage_massChangeRateFactor().value             = *(parameters + 0);
            *get_SetDataStorage_internalHeatGenerationRateFactor().value = *(parameters + 1);
        }

        void residual::decomposeAdditionalDOF() {
            /*!
             * Decompose the additional DOF vectors
             */

            constexpr unsigned int dim = 3;

            TARDIGRADE_ERROR_TOOLS_CHECK(hydra->getAdditionalDOF()->size() >= getDOFVelocityGradientIndex() + dim * dim,
                                         "The additional DOF vector is of size " +
                                             std::to_string(hydra->getAdditionalDOF()->size()) +
                                             " which is less than the required size of " +
                                             std::to_string(getDOFVelocityGradientIndex() + dim * dim));

            TARDIGRADE_ERROR_TOOLS_CHECK(hydra->getPreviousAdditionalDOF()->size() >=
                                             getDOFVelocityGradientIndex() + dim * dim,
                                         "The additional DOF vector is of size " +
                                             std::to_string(hydra->getPreviousAdditionalDOF()->size()) +
                                             " which is less than the required size of " +
                                             std::to_string(getDOFVelocityGradientIndex() + dim * dim));

            auto density = get_SetDataStorage_density();

            auto internalEnergy = get_SetDataStorage_internalEnergy();

            auto dofVelocityGradient = get_SetDataStorage_dofVelocityGradient();

            auto previousDOFVelocityGradient = get_SetDataStorage_previousDOFVelocityGradient();

            *density.value = *(hydra->getAdditionalDOF()->begin() + getDensityIndex());

            *internalEnergy.value = *(hydra->getAdditionalDOF()->begin() + getInternalEnergyIndex());

            *dofVelocityGradient.value =
                secondOrderTensor(hydra->getAdditionalDOF()->begin() + getDOFVelocityGradientIndex(),
                                  hydra->getAdditionalDOF()->begin() + getDOFVelocityGradientIndex() + dim * dim);

            *previousDOFVelocityGradient.value =
                secondOrderTensor(hydra->getPreviousAdditionalDOF()->begin() + getDOFVelocityGradientIndex(),
                                  hydra->getPreviousAdditionalDOF()->begin() + getDOFVelocityGradientIndex() +
                                      dim * dim);
        }

        void residual::setPrecedingDeformationGradient(const bool &isPrevious) {
            /*!
             * Set the preceding deformation gradient
             *
             * \param &isPrevious: Flag for whether to set the current (false) or previous (true) value
             */

            if (isPrevious) {
                auto precedingDeformationGradient = get_SetDataStorage_previousPrecedingDeformationGradient();
                *precedingDeformationGradient.value =
                    hydra->deformation->getPreviousPrecedingConfiguration(getDOFConfigurationIndex());

            } else {
                auto precedingDeformationGradient   = get_SetDataStorage_precedingDeformationGradient();
                *precedingDeformationGradient.value = hydra->deformation->getPrecedingConfiguration(getDOFConfigurationIndex());
            }
        }

        void residual::setPrecedingDeformationGradientDerivatives(const bool &isPrevious) {
            /*!
             * Set the derivatives of the preceding deformation gradient
             *
             * \param &isPrevious Flag for whether to set the current (false) or previous (true) values
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            auto num_configs = hydra->getNumConfigurations();

            const fourthOrderTensor *dF1dF;

            const floatVector *dF1dFn;

            floatVector dpFdFs;

            SetDataStorageBase<secondOrderTensor> precedingDeformationGradient;

            SetDataStorageBase<fourthOrderTensor> dpFdF;

            SetDataStorageBase<floatVector> dpFdFn;

            if (isPrevious) {
                TARDIGRADE_ERROR_TOOLS_CATCH(dF1dF = hydra->get_previousdF1dF())

                TARDIGRADE_ERROR_TOOLS_CATCH(dF1dFn = hydra->get_previousdF1dFn())

                TARDIGRADE_ERROR_TOOLS_CATCH(
                    dpFdFs = hydra->getPreviousPrecedingConfigurationJacobian(getDOFConfigurationIndex()))

                auto precedingDeformationGradient = get_SetDataStorage_previousPrecedingDeformationGradient();
                *precedingDeformationGradient.value =
                    hydra->deformation->getPreviousPrecedingConfiguration(getDOFConfigurationIndex());

                dpFdF = get_SetDataStorage_dPreviousPrecedingDeformationGradientdPreviousDeformationGradient();

                dpFdFn = get_SetDataStorage_dPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients();

            } else {
                TARDIGRADE_ERROR_TOOLS_CATCH(dF1dF = hydra->get_dF1dF())

                TARDIGRADE_ERROR_TOOLS_CATCH(dF1dFn = hydra->get_dF1dFn())

                TARDIGRADE_ERROR_TOOLS_CATCH(dpFdFs =
                                                 hydra->getPrecedingConfigurationJacobian(getDOFConfigurationIndex()))

                auto precedingDeformationGradient   = get_SetDataStorage_precedingDeformationGradient();
                *precedingDeformationGradient.value = hydra->deformation->getPrecedingConfiguration(getDOFConfigurationIndex());

                dpFdF = get_SetDataStorage_dPrecedingDeformationGradientdDeformationGradient();

                dpFdFn = get_SetDataStorage_dPrecedingDeformationGradientdSubDeformationGradients();
            }

            dpFdF.zero(sot_dim * sot_dim);

            dpFdFn.zero(sot_dim * sot_dim * (num_configs - 1));

            for (unsigned int i = 0; i < sot_dim; i++) {
                for (unsigned int j = 0; j < sot_dim; j++) {
                    for (unsigned int k = 0; k < sot_dim; k++) {
                        (*dpFdF.value)[sot_dim * i + k] +=
                            dpFdFs[num_configs * sot_dim * i + j] * (*dF1dF)[sot_dim * j + k];
                    }
                }
            }

            for (unsigned int i = 0; i < sot_dim; i++) {
                for (unsigned int j = 0; j < (num_configs - 1) * sot_dim; j++) {
                    (*dpFdFn.value)[(num_configs - 1) * sot_dim * i + j] +=
                        dpFdFs[num_configs * sot_dim * i + j + sot_dim];

                    for (unsigned int k = 0; k < sot_dim; k++) {
                        (*dpFdFn.value)[(num_configs - 1) * sot_dim * i + j] +=
                            dpFdFs[num_configs * sot_dim * i + k] * (*dF1dFn)[(num_configs - 1) * sot_dim * k + j];
                    }
                }
            }
        }

        void residual::setPrecedingDeformationGradient() {
            /*!
             * Set the value of the preceding deformation gradient
             */

            setPrecedingDeformationGradient(false);
        }

        void residual::setPreviousPrecedingDeformationGradient() {
            /*!
             * Set the value of the previous preceding deformation gradient
             */

            setPrecedingDeformationGradient(true);
        }

        void residual::setdPrecedingDeformationGradientdDeformationGradient() {
            /*!
             * Set the derivative of the preceding deformation gradient w.r.t. the total deformation gradient
             */

            setPrecedingDeformationGradientDerivatives(false);
        }

        void residual::setdPrecedingDeformationGradientdSubDeformationGradients() {
            /*!
             * Set the derivative of the preceding deformation gradient w.r.t. the sub-deformation gradients
             */

            setPrecedingDeformationGradientDerivatives(false);
        }

        void residual::setdPreviousPrecedingDeformationGradientdPreviousDeformationGradient() {
            /*!
             * Set the derivative of the previous preceding deformation gradient w.r.t. the previous total deformation
             * gradient
             */

            setPrecedingDeformationGradientDerivatives(true);
        }

        void residual::setdPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients() {
            /*!
             * Set the derivative of the previous preceding deformation gradient w.r.t. the previous sub-deformation
             * gradients
             */

            setPrecedingDeformationGradientDerivatives(true);
        }

        void residual::setDOFIntermediateVelocityGradient(const bool &isPrevious) {
            /*!
             * Set the velocity gradient in the intermediate configuration
             *
             * \param &isPrevious: Flag for whether this is being computed for the current or previous timestep
             */

            const secondOrderTensor *velocityGradient;

            const secondOrderTensor *precedingDeformationGradient;

            SetDataStorageBase<secondOrderTensor> intermediateVelocityGradient;

            if (isPrevious) {
                TARDIGRADE_ERROR_TOOLS_CATCH(velocityGradient = get_previousDOFVelocityGradient())

                TARDIGRADE_ERROR_TOOLS_CATCH(precedingDeformationGradient = get_previousPrecedingDeformationGradient())

                intermediateVelocityGradient = get_SetDataStorage_previousDOFIntermediateVelocityGradient();

            } else {
                TARDIGRADE_ERROR_TOOLS_CATCH(velocityGradient = get_dofVelocityGradient())

                TARDIGRADE_ERROR_TOOLS_CATCH(precedingDeformationGradient = get_precedingDeformationGradient())

                intermediateVelocityGradient = get_SetDataStorage_dofIntermediateVelocityGradient();
            }

            tardigradeConstitutiveTools::pullBackVelocityGradient(*velocityGradient, *precedingDeformationGradient,
                                                                  *intermediateVelocityGradient.value);
        }

        void residual::setDOFIntermediateVelocityGradientDerivatives(const bool &isPrevious) {
            /*!
             * Set the derivatives of the velocity gradient in the intermediate configuration
             *
             * \param &isPrevious: Flag for whether this is being computed for the current or previous timestep
             */

            constexpr unsigned int dim = 3;

            auto sot_dim = dim * dim;

            auto num_configs = hydra->getNumConfigurations();

            const secondOrderTensor *velocityGradient;

            const secondOrderTensor *precedingDeformationGradient;

            const fourthOrderTensor *dPFdF;

            const floatVector *dPFdFn;

            SetDataStorageBase<secondOrderTensor> intermediateVelocityGradient;

            SetDataStorageBase<secondOrderTensor> dILdL;

            SetDataStorageBase<fourthOrderTensor> dILdF;

            SetDataStorageBase<fourthOrderTensor> dILdFn;

            if (isPrevious) {
                TARDIGRADE_ERROR_TOOLS_CATCH(
                    dPFdF = get_dPreviousPrecedingDeformationGradientdPreviousDeformationGradient())

                TARDIGRADE_ERROR_TOOLS_CATCH(
                    dPFdFn = get_dPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients())

                TARDIGRADE_ERROR_TOOLS_CATCH(velocityGradient = get_previousDOFVelocityGradient())

                TARDIGRADE_ERROR_TOOLS_CATCH(precedingDeformationGradient = get_previousPrecedingDeformationGradient())

                intermediateVelocityGradient = get_SetDataStorage_previousDOFIntermediateVelocityGradient();

                dILdL = get_SetDataStorage_dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient();

                dILdF = get_SetDataStorage_dPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient();

                dILdFn = get_SetDataStorage_dPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients();

            } else {
                TARDIGRADE_ERROR_TOOLS_CATCH(dPFdF = get_dPrecedingDeformationGradientdDeformationGradient())

                TARDIGRADE_ERROR_TOOLS_CATCH(dPFdFn = get_dPrecedingDeformationGradientdSubDeformationGradients())

                TARDIGRADE_ERROR_TOOLS_CATCH(velocityGradient = get_dofVelocityGradient())

                TARDIGRADE_ERROR_TOOLS_CATCH(precedingDeformationGradient = get_precedingDeformationGradient())

                intermediateVelocityGradient = get_SetDataStorage_dofIntermediateVelocityGradient();

                dILdL = get_SetDataStorage_dDOFIntermediateVelocityGradientdDOFVelocityGradient();

                dILdF = get_SetDataStorage_dDOFIntermediateVelocityGradientdDeformationGradient();

                dILdFn = get_SetDataStorage_dDOFIntermediateVelocityGradientdSubDeformationGradients();
            }

            fourthOrderTensor dILdPF;

            tardigradeConstitutiveTools::pullBackVelocityGradient(*velocityGradient, *precedingDeformationGradient,
                                                                  *intermediateVelocityGradient.value, *dILdL.value,
                                                                  dILdPF);

            dILdF.zero(sot_dim * sot_dim);

            dILdFn.zero((num_configs - 1) * sot_dim * sot_dim);

            for (unsigned int i = 0; i < sot_dim; i++) {
                for (unsigned int j = 0; j < sot_dim; j++) {
                    for (unsigned int k = 0; k < sot_dim; k++) {
                        (*dILdF.value)[sot_dim * i + k] += dILdPF[sot_dim * i + j] * (*dPFdF)[sot_dim * j + k];
                    }

                    for (unsigned int k = 0; k < (num_configs - 1) * sot_dim; k++) {
                        (*dILdFn.value)[(num_configs - 1) * sot_dim * i + k] +=
                            dILdPF[sot_dim * i + j] * (*dPFdFn)[(num_configs - 1) * sot_dim * j + k];
                    }
                }
            }
        }

        void residual::setDOFIntermediateVelocityGradient() {
            /*!
             * Set the current intermediate velocity gradient
             */

            setDOFIntermediateVelocityGradient(false);
        }

        void residual::setPreviousDOFIntermediateVelocityGradient() {
            /*!
             * Set the previous intermediate velocity gradient
             */

            setDOFIntermediateVelocityGradient(true);
        }

        void residual::setdDOFIntermediateVelocityGradientdDOFVelocityGradient() {
            /*!
             * Set the derivative of the current intermediate velocity gradient w.r.t. the mass change velocity gradient
             */

            setDOFIntermediateVelocityGradientDerivatives(false);
        }

        void residual::setdDOFIntermediateVelocityGradientdDeformationGradient() {
            /*!
             * Set the derivative of the current intermediate velocity gradient w.r.t. the deformation gradient
             */

            setDOFIntermediateVelocityGradientDerivatives(false);
        }

        void residual::setdDOFIntermediateVelocityGradientdSubDeformationGradients() {
            /*!
             * Set the derivative of the current intermediate velocity gradient w.r.t. the sub-deformation gradients
             */

            setDOFIntermediateVelocityGradientDerivatives(false);
        }

        void residual::setdPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient() {
            /*!
             * Set the derivative of the previous intermediate velocity gradient w.r.t. the previous mass change
             * velocity gradient
             */

            setDOFIntermediateVelocityGradientDerivatives(true);
        }

        void residual::setdPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient() {
            /*!
             * Set the derivative of the previous intermediate velocity gradient w.r.t. the previous deformation
             * gradient
             */

            setDOFIntermediateVelocityGradientDerivatives(true);
        }

        void residual::setdPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients() {
            /*!
             * Set the derivative of the previous intermediate velocity gradient w.r.t. the previous sub-deformation
             * gradients
             */

            setDOFIntermediateVelocityGradientDerivatives(true);
        }

        void residual::setDOFDeformationGradient() {
            /*!
             * Set the mass-change deformation gradient
             */

            const secondOrderTensor *intermediateVelocityGradient = get_dofIntermediateVelocityGradient();

            const secondOrderTensor *previousIntermediateVelocityGradient =
                get_previousDOFIntermediateVelocityGradient();

            const secondOrderTensor previousDOFDeformationGradient =
                hydra->getPreviousConfiguration(getDOFConfigurationIndex());

            auto dofDeformationGradient = get_SetDataStorage_dofDeformationGradient();

            if (getUseTrapezoidalIntegration()) {
                secondOrderTensor temp;

                TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeConstitutiveTools::evolveF(
                    hydra->getDeltaTime(), previousDOFDeformationGradient, *previousIntermediateVelocityGradient,
                    *intermediateVelocityGradient, temp, *dofDeformationGradient.value, 1 - getIntegrationParameter())

                )

            } else {
                TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeConstitutiveTools::evolveFExponentialMap(
                    hydra->getDeltaTime(), previousDOFDeformationGradient, *previousIntermediateVelocityGradient,
                    *intermediateVelocityGradient, *dofDeformationGradient.value, getIntegrationParameter()))
            }
        }

        void residual::setDOFDeformationGradientDerivatives(const bool &computePrevious) {
            /*!
             * Compute the derivatives of the mass-change deformation gradient
             *
             * \param &computePrevious: Compute the gradients w.r.t. previous values
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            auto num_configs = hydra->getNumConfigurations();

            const secondOrderTensor *intermediateVelocityGradient = get_dofIntermediateVelocityGradient();

            const fourthOrderTensor *dILdL = get_dDOFIntermediateVelocityGradientdDOFVelocityGradient();

            const fourthOrderTensor *dILdF = get_dDOFIntermediateVelocityGradientdDeformationGradient();

            const floatVector *dILdFn = get_dDOFIntermediateVelocityGradientdSubDeformationGradients();

            const secondOrderTensor *previousIntermediateVelocityGradient =
                get_previousDOFIntermediateVelocityGradient();

            const secondOrderTensor previousDOFDeformationGradient =
                hydra->getPreviousConfiguration(getDOFConfigurationIndex());

            auto dofDeformationGradient = get_SetDataStorage_dofDeformationGradient();

            fourthOrderTensor dFmdIL;

            if (computePrevious) {
                const fourthOrderTensor *dILpdL =
                    get_dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient();

                const fourthOrderTensor *dILpdF =
                    get_dPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient();

                const floatVector *dILpdFn =
                    get_dPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients();

                fourthOrderTensor dFmdFp;

                fourthOrderTensor dFmdILp;

                if (getUseTrapezoidalIntegration()) {
                    secondOrderTensor temp;
                    fourthOrderTensor temp2;

                    TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeConstitutiveTools::evolveFFlatJ(
                        hydra->getDeltaTime(), previousDOFDeformationGradient, *previousIntermediateVelocityGradient,
                        *intermediateVelocityGradient, temp, *dofDeformationGradient.value, dFmdIL, temp2, dFmdFp,
                        dFmdILp, 1 - getIntegrationParameter()))

                } else {
                    TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeConstitutiveTools::evolveFExponentialMap(
                        hydra->getDeltaTime(), previousDOFDeformationGradient, *previousIntermediateVelocityGradient,
                        *intermediateVelocityGradient, *dofDeformationGradient.value, dFmdIL, dFmdFp, dFmdILp,
                        getIntegrationParameter()))
                }

                auto dFmdPreviousL = get_SetDataStorage_dDOFDeformationGradientdPreviousDOFVelocityGradient();
                dFmdPreviousL.zero(sot_dim * sot_dim);

                auto dFmdPreviousF = get_SetDataStorage_dDOFDeformationGradientdPreviousDeformationGradient();
                dFmdPreviousF.zero(sot_dim * sot_dim);

                auto dFmdPreviousFn = get_SetDataStorage_dDOFDeformationGradientdPreviousSubDeformationGradients();
                dFmdPreviousFn.zero(sot_dim * sot_dim * (num_configs - 1));

                for (unsigned int i = 0; i < sot_dim; i++) {
                    for (unsigned int j = 0; j < sot_dim; j++) {
                        (*dFmdPreviousFn.value)[(num_configs - 1) * sot_dim * i + j +
                                                (getDOFConfigurationIndex() - 1) * sot_dim] += dFmdFp[sot_dim * i + j];

                        for (unsigned int k = 0; k < sot_dim; k++) {
                            (*dFmdPreviousL.value)[sot_dim * i + k] +=
                                dFmdILp[sot_dim * i + j] * (*dILpdL)[sot_dim * j + k];

                            (*dFmdPreviousF.value)[sot_dim * i + k] +=
                                dFmdILp[sot_dim * i + j] * (*dILpdF)[sot_dim * j + k];
                        }

                        for (unsigned int k = 0; k < (num_configs - 1) * sot_dim; k++) {
                            (*dFmdPreviousFn.value)[(num_configs - 1) * sot_dim * i + k] +=
                                dFmdILp[sot_dim * i + j] * (*dILpdFn)[(num_configs - 1) * sot_dim * j + k];
                        }
                    }
                }

            } else {
                if (getUseTrapezoidalIntegration()) {
                    TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeConstitutiveTools::evolveFFlatJ(
                        hydra->getDeltaTime(), previousDOFDeformationGradient, *previousIntermediateVelocityGradient,
                        *intermediateVelocityGradient, *dofDeformationGradient.value, dFmdIL,
                        1 - getIntegrationParameter()))

                } else {
                    TARDIGRADE_ERROR_TOOLS_CATCH(tardigradeConstitutiveTools::evolveFExponentialMap(
                        hydra->getDeltaTime(), previousDOFDeformationGradient, *previousIntermediateVelocityGradient,
                        *intermediateVelocityGradient, *dofDeformationGradient.value, dFmdIL,
                        getIntegrationParameter()))
                }
            }

            auto dFmdL = get_SetDataStorage_dDOFDeformationGradientdDOFVelocityGradient();
            dFmdL.zero(sot_dim * sot_dim);

            auto dFmdF = get_SetDataStorage_dDOFDeformationGradientdDeformationGradient();
            dFmdF.zero(sot_dim * sot_dim);

            auto dFmdFn = get_SetDataStorage_dDOFDeformationGradientdSubDeformationGradients();
            dFmdFn.zero(sot_dim * sot_dim * (num_configs - 1));

            for (unsigned int i = 0; i < sot_dim; i++) {
                for (unsigned int j = 0; j < sot_dim; j++) {
                    for (unsigned int k = 0; k < sot_dim; k++) {
                        (*dFmdL.value)[sot_dim * i + k] += dFmdIL[sot_dim * i + j] * (*dILdL)[sot_dim * j + k];

                        (*dFmdF.value)[sot_dim * i + k] += dFmdIL[sot_dim * i + j] * (*dILdF)[sot_dim * j + k];
                    }

                    for (unsigned int k = 0; k < (num_configs - 1) * sot_dim; k++) {
                        (*dFmdFn.value)[(num_configs - 1) * sot_dim * i + k] +=
                            dFmdIL[sot_dim * i + j] * (*dILdFn)[(num_configs - 1) * sot_dim * j + k];
                    }
                }
            }
        }

        void residual::setdDOFDeformationGradientdDOFVelocityGradient() {
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the mass change velocity gradient
             */

            setDOFDeformationGradientDerivatives(false);
        }

        void residual::setdDOFDeformationGradientdDeformationGradient() {
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the deformation gradient
             */

            setDOFDeformationGradientDerivatives(false);
        }

        void residual::setdDOFDeformationGradientdSubDeformationGradients() {
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the sub deformation gradients
             */

            setDOFDeformationGradientDerivatives(false);
        }

        void residual::setdDOFDeformationGradientdPreviousDOFVelocityGradient() {
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the previous mass change velocity
             * gradient
             */

            setDOFDeformationGradientDerivatives(true);
        }

        void residual::setdDOFDeformationGradientdPreviousDeformationGradient() {
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the previous deformation gradient
             */

            setDOFDeformationGradientDerivatives(true);
        }

        void residual::setdDOFDeformationGradientdPreviousSubDeformationGradients() {
            /*!
             * Compute the derivative of the mass-change deformation gradient w.r.t. the previous sub deformation
             * gradients
             */

            setDOFDeformationGradientDerivatives(true);
        }

        void residual::setResidual() {
            /*!
             * Set the value of the residual
             *
             * Defined as the residual's computed thermal deformation gradient minus the value stored in hydra's
             * configurations.
             */

            auto dim = hydra->getDimension();

            auto sot_dim = dim * dim;

            auto dofConfigurationIndex = getDOFConfigurationIndex();

            auto residual = get_SetDataStorage_residual();

            residual.zero(sot_dim + 2);

            std::transform(std::begin(*get_dofDeformationGradient()), std::end(*get_dofDeformationGradient()),
                           hydra->deformation->get_configurations()->begin() + dofConfigurationIndex * sot_dim, residual.begin(),
                           std::minus<>());

            TARDIGRADE_ERROR_TOOLS_CHECK(getStateVariableIndices()->size() == 2,
                                         "The state variable indices must have a size of 2 instead of " +
                                             std::to_string(getStateVariableIndices()->size()));

            (*residual.value)[sot_dim + 0] =
                *get_massChangeRate() - (*hydra->get_nonLinearSolveStateVariables())[(*getStateVariableIndices())[0]];

            (*residual.value)[sot_dim + 1] =
                *get_internalHeatGenerationRate() -
                (*hydra->get_nonLinearSolveStateVariables())[(*getStateVariableIndices())[1]];
        }

        void residual::setJacobian() {
            /*!
             * Set the values of the jacobian
             */

            auto sot_dim = hydra->getSOTDimension();

            auto num_unknowns = hydra->getNumUnknowns();

            auto num_equations = getNumEquations();

            auto num_configs = hydra->getNumConfigurations();

            auto jacobian = get_SetDataStorage_jacobian();
            jacobian.zero(num_equations * num_unknowns);

            const floatVector *dFmdFn = get_dDOFDeformationGradientdSubDeformationGradients();

            for (unsigned int i = 0; i < sot_dim; i++) {
                (*jacobian.value)[num_unknowns * i + sot_dim * getDOFConfigurationIndex() + i] += -1;

                for (unsigned int j = 0; j < (num_configs - 1) * sot_dim; j++) {
                    (*jacobian.value)[num_unknowns * i + j + sot_dim] += (*dFmdFn)[(num_configs - 1) * sot_dim * i + j];
                }
            }

            unsigned int offset = sot_dim * num_configs;

            (*jacobian.value)[num_unknowns * (sot_dim + 0) + offset + (*getStateVariableIndices())[0]] += -1;

            (*jacobian.value)[num_unknowns * (sot_dim + 1) + offset + (*getStateVariableIndices())[1]] += -1;
        }

        void residual::setdRdT() {
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            const unsigned int sot_dim = hydra->getSOTDimension();

            auto dRdT = get_SetDataStorage_dRdT();

            dRdT.zero(sot_dim + getStateVariableIndices()->size());
        }

        void residual::setdRdF() {
            /*!
             * Set the derivative of the residual w.r.t. the deformation gradient
             */

            auto sot_dim = hydra->getSOTDimension();
            auto dRdF    = get_SetDataStorage_dRdF();
            dRdF.zero(getNumEquations() * sot_dim);
            std::copy(std::begin(*get_dDOFDeformationGradientdDeformationGradient()),
                      std::end(*get_dDOFDeformationGradientdDeformationGradient()), dRdF.begin());
        }

        void residual::setdRdAdditionalDOF() {
            /*!
             * Set the additional derivatives
             */

            auto sot_dim = hydra->getSOTDimension();

            auto num_equations = getNumEquations();

            auto num_additional_dof = hydra->getAdditionalDOF()->size();

            const fourthOrderTensor *dDOFDeformationdDOFVelocityGradient =
                get_dDOFDeformationGradientdDOFVelocityGradient();

            auto dRdAdditionalDOF = get_SetDataStorage_dRdAdditionalDOF();
            dRdAdditionalDOF.zero(num_equations * num_additional_dof);

            auto offset = getDOFVelocityGradientIndex();

            for (unsigned int i = 0; i < sot_dim; i++) {
                for (unsigned int j = 0; j < sot_dim; j++) {
                    (*dRdAdditionalDOF.value)[num_additional_dof * i + j + offset] +=
                        (*dDOFDeformationdDOFVelocityGradient)[sot_dim * i + j];
                }
            }

            (*dRdAdditionalDOF.value)[num_additional_dof * (sot_dim + 0) + getDensityIndex()] +=
                *get_dMassChangeRatedDensity();

            std::transform(std::begin(*get_dMassChangeRatedDOFVelocityGradient()),
                           std::end(*get_dMassChangeRatedDOFVelocityGradient()),
                           dRdAdditionalDOF.begin() + num_additional_dof * (sot_dim + 0) +
                               getDOFVelocityGradientIndex(),
                           dRdAdditionalDOF.begin() + num_additional_dof * (sot_dim + 0) +
                               getDOFVelocityGradientIndex(),
                           std::plus<>());

            (*dRdAdditionalDOF.value)[num_additional_dof * (sot_dim + 1) + getDensityIndex()] +=
                *get_dInternalHeatGenerationRatedDensity();

            std::transform(std::begin(*get_dInternalHeatGenerationRatedDOFVelocityGradient()),
                           std::end(*get_dInternalHeatGenerationRatedDOFVelocityGradient()),
                           dRdAdditionalDOF.begin() + num_additional_dof * (sot_dim + 1) +
                               getDOFVelocityGradientIndex(),
                           dRdAdditionalDOF.begin() + num_additional_dof * (sot_dim + 1) +
                               getDOFVelocityGradientIndex(),
                           std::plus<>());

            (*dRdAdditionalDOF.value)[num_additional_dof * (sot_dim + 1) + getInternalEnergyIndex()] +=
                *get_dInternalHeatGenerationRatedInternalEnergy();
        }

        void residual::suggestInitialIterateValues(std::vector<unsigned int> &indices, std::vector<floatType> &values) {
            /*!
             * Suggest initial iterate values to try and improve convergence
             *
             * \param &indices: The indices of the unknown vector to suggest initial values
             * \param &values: The values to suggest
             */

            auto sot_dim = hydra->getSOTDimension();

            auto configuration = getDOFConfigurationIndex();

            auto num_configurations = hydra->getNumConfigurations();

            const secondOrderTensor *dofDeformationGradient = get_dofDeformationGradient();

            indices = std::vector<unsigned int>(sot_dim + 2, sot_dim * configuration);

            for (unsigned int i = 0; i < sot_dim; i++) {
                indices[i] += i;
            }

            indices[sot_dim + 0] = sot_dim * num_configurations + (*getStateVariableIndices())[0];
            indices[sot_dim + 1] = sot_dim * num_configurations + (*getStateVariableIndices())[1];

            values = std::vector<floatType>(sot_dim + 2, 0);

            std::copy(std::begin(*dofDeformationGradient), std::end(*dofDeformationGradient), std::begin(values));

            values[sot_dim + 0] = *get_massChangeRate();
            values[sot_dim + 1] = *get_internalHeatGenerationRate();
        }

        void residual::setMassChangeRate() {
            /*!
             * Set the mass-change rate from the external velocity gradient
             */

            auto dim = hydra->getDimension();

            auto mass_change_rate_factor = get_massChangeRateFactor();

            auto density = get_density();

            auto velocity_gradient = get_dofVelocityGradient();

            auto mass_change_rate = get_SetDataStorage_massChangeRate();

            (*mass_change_rate.value) = 0.;
            for (unsigned int i = 0; i < dim; ++i) {
                *mass_change_rate.value += (*mass_change_rate_factor) * (*density) * (*velocity_gradient)[dim * i + i];
            }
        }

        void residual::setMassChangeRateGradients() {
            /*!
             * Compute the gradients of the mass change rate
             */

            auto dim = hydra->getDimension();

            auto sot_dim = dim * dim;

            auto mass_change_rate_factor = get_massChangeRateFactor();

            auto density = get_density();

            auto velocity_gradient = get_dofVelocityGradient();

            auto dCdRho = get_SetDataStorage_dMassChangeRatedDensity();

            auto dCdV = get_SetDataStorage_dMassChangeRatedDOFVelocityGradient();

            *dCdRho.value = 0.;
            *dCdV.value   = secondOrderTensor(sot_dim, 0);

            for (unsigned int i = 0; i < dim; ++i) {
                *dCdRho.value += (*mass_change_rate_factor) * (*velocity_gradient)[dim * i + i];
                (*dCdV.value)[dim * i + i] += (*mass_change_rate_factor) * (*density);
            }
        }

        void residual::setInternalHeatGenerationRate() {
            /*!
             * Compute the internal heat generation rate
             */

            auto internal_heat_generation_rate_factor = get_internalHeatGenerationRateFactor();

            auto normalize_by_density = getInternalEnergyScaledByDensity();

            auto dim = hydra->getDimension();

            auto density = get_density();

            auto internal_energy = get_internalEnergy();

            auto velocity_gradient = get_dofVelocityGradient();

            auto internal_heat_generation_rate = get_SetDataStorage_internalHeatGenerationRate();

            *internal_heat_generation_rate.value = 0.;
            if (normalize_by_density) {
                for (unsigned int i = 0; i < dim; ++i) {
                    *internal_heat_generation_rate.value += (*internal_heat_generation_rate_factor) *
                                                            (*internal_energy) / (*density) *
                                                            (*velocity_gradient)[dim * i + i];
                }

            } else {
                for (unsigned int i = 0; i < dim; ++i) {
                    *internal_heat_generation_rate.value += (*internal_heat_generation_rate_factor) *
                                                            (*internal_energy) * (*velocity_gradient)[dim * i + i];
                }
            }
        }

        void residual::setInternalHeatGenerationRateGradients() {
            /*!
             * Set the gradients of the internal heat generation rate
             */

            auto internal_heat_generation_rate_factor = get_internalHeatGenerationRateFactor();

            auto normalize_by_density = getInternalEnergyScaledByDensity();

            auto dim = hydra->getDimension();

            auto density = get_density();

            auto internal_energy = get_internalEnergy();

            auto velocity_gradient = get_dofVelocityGradient();

            auto drdRho = get_SetDataStorage_dInternalHeatGenerationRatedDensity();

            auto drdE = get_SetDataStorage_dInternalHeatGenerationRatedInternalEnergy();

            auto drdV = get_SetDataStorage_dInternalHeatGenerationRatedDOFVelocityGradient();

            *drdRho.value = 0.;
            *drdE.value   = 0.;
            *drdV.value   = secondOrderTensor(dim * dim, 0.);
            if (normalize_by_density) {
                for (unsigned int i = 0; i < dim; ++i) {
                    *drdRho.value -= (*internal_heat_generation_rate_factor) * (*internal_energy) /
                                     ((*density) * (*density)) * (*velocity_gradient)[dim * i + i];
                    *drdE.value +=
                        (*internal_heat_generation_rate_factor) / (*density) * (*velocity_gradient)[dim * i + i];
                    (*drdV.value)[dim * i + i] +=
                        (*internal_heat_generation_rate_factor) * (*internal_energy) / (*density);
                }

            } else {
                for (unsigned int i = 0; i < dim; ++i) {
                    *drdE.value += (*internal_heat_generation_rate_factor) * (*velocity_gradient)[dim * i + i];
                    (*drdV.value)[dim * i + i] += (*internal_heat_generation_rate_factor) * (*internal_energy);
                }
            }
        }

    }  // namespace dofVelocityGradientDeformation

}  // namespace tardigradeHydra
