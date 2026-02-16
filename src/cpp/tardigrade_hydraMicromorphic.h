/**
 ******************************************************************************
 * \file tardigrade_hydraMicromorphic.h
 ******************************************************************************
 * A C++ utility for constructing finite deformation micromorphic constitutive
 * models.
 ******************************************************************************
 */

#ifndef TARDIGRADE_HYDRA_MICROMORPHIC_H
#define TARDIGRADE_HYDRA_MICROMORPHIC_H

#include <tardigrade_hydra.h>

namespace tardigradeHydra {

    // Define tensors of known size
    typedef std::vector<floatType> fifthOrderTensor;  //!< Fifth order tensors
    typedef std::vector<floatType> sixthOrderTensor;  //!< Sixth order tensors

    /*!
     * HydraClassicalConfiguration: A class which defines a classical deformation problem
     */
    class HydraMicromorphicConfiguration : public HydraConfigurationBase {
       public:
        HydraMicromorphicConfiguration() { configuration_unknown_count = 45; }
    };

    /*!
     * A storage class for the micromorphic degrees of freedom
     */
    class MicromorphicDOFStorage : public DOFStorageBase {
       public:
        /*!
         * Constructor which sets the dof information
         *
         * TODO: We're eventually going to store the degrees of freedom into a single storage
         *       array but, for now, I'm leaving them as discrete just to help with the
         *       transition.
         *
         * \param &time: The current time
         * \param &deltaTime: The change in time from the previous time
         * \param &temperature: The current temperature
         * \param &previous_temperature: The previous temperature
         * \param &deformation_gradient: The deformation gradient
         * \param &previous_deformation_gradient: The previous deformation gradient
         * \param &micro_deformation: The micro-deformation
         * \param &previous_micro_deformation: The previous micro-deformation
         * version of DOFStorage
         * \param &gradient_micro_deformation: The spatial gradient of the micro-deformation
         * \param &previous_gradient_micro_deformation: The previous spatial gradient of the micro-deformation
         * \param &additional_dof: The additional degrees of freedom
         * \param &previous_additional_dof: The previous additional degrees of freedom
         */
        MicromorphicDOFStorage(const floatType &time, const floatType &deltaTime, const floatType &temperature,
                               const floatType &previous_temperature, const floatVector &deformation_gradient,
                               const floatVector &previous_deformation_gradient, const floatVector &micro_deformation,
                               const floatVector &previous_micro_deformation,
                               const floatVector &gradient_micro_deformation,
                               const floatVector &previous_gradient_micro_deformation,
                               const floatVector &additional_dof, const floatVector &previous_additional_dof)
            : DOFStorageBase(time, deltaTime, temperature, previous_temperature, deformation_gradient,
                             previous_deformation_gradient, additional_dof, previous_additional_dof),
              _micro_deformation(micro_deformation),
              _previous_micro_deformation(previous_micro_deformation),
              _gradient_micro_deformation(gradient_micro_deformation),
              _previous_gradient_micro_deformation(previous_gradient_micro_deformation) {}

        //! The micro deformation
        const floatVector _micro_deformation;

        //! The previous micro deformation
        const floatVector _previous_micro_deformation;

        //! The gradient of the micro deformation
        const floatVector _gradient_micro_deformation;

        //! The previous spatial gradient of the micro deformation
        const floatVector _previous_gradient_micro_deformation;

       protected:
    };

    //! The base class for hydra framework micromorphic material models
    class hydraBaseMicromorphic : public hydraBase {
       public:
        hydraBaseMicromorphic() {}

        hydraBaseMicromorphic(const MicromorphicDOFStorage &DOFStorage, const floatVector &previousStateVariables,
                              const floatVector &parameters, const unsigned int numConfigurations,
                              const unsigned int     numNonLinearSolveStateVariables,
                              HydraConfigurationBase _hydra_configuration = HydraMicromorphicConfiguration());

        virtual void initialize() override;

        //! Get the current micro-deformation tensor
        const secondOrderTensor *getMicroDeformation() { return getScaledMicroDeformation(); }

        //! Get the previous micro-deformation tensor
        const secondOrderTensor *getPreviousMicroDeformation() {
            auto local_dof = static_cast<const tardigradeHydra::MicromorphicDOFStorage *>(dof);
            return &local_dof->_previous_micro_deformation;
        }

        //! Get the current spatial gradient w.r.t. the reference configuration of the micro-deformation tensor
        const thirdOrderTensor *getGradientMicroDeformation() { return getScaledGradientMicroDeformation(); }

        //! Get the previous spatial gradient w.r.t. the reference configuration of the micro-deformation tensor
        const thirdOrderTensor *getPreviousGradientMicroDeformation() {
            auto local_dof = static_cast<const tardigradeHydra::MicromorphicDOFStorage *>(dof);
            return &local_dof->_previous_gradient_micro_deformation;
        }

        secondOrderTensor getSubMicroConfiguration(const unsigned int &lowerIndex, const unsigned int &upperIndex);

        secondOrderTensor getPrecedingMicroConfiguration(const unsigned int &index);

        secondOrderTensor getFollowingMicroConfiguration(const unsigned int &index);

        secondOrderTensor getMicroConfiguration(const unsigned int &index);

        secondOrderTensor getPreviousSubMicroConfiguration(const unsigned int &lowerIndex,
                                                           const unsigned int &upperIndex);

        secondOrderTensor getPreviousPrecedingMicroConfiguration(const unsigned int &index);

        secondOrderTensor getPreviousFollowingMicroConfiguration(const unsigned int &index);

        secondOrderTensor getPreviousMicroConfiguration(const unsigned int &index);

        floatVector getSubMicroConfigurationJacobian(const unsigned int &lowerIndex, const unsigned int &upperIndex);

        floatVector getPrecedingMicroConfigurationJacobian(const unsigned int &index);

        floatVector getFollowingMicroConfigurationJacobian(const unsigned int &index);

        floatVector getPreviousSubMicroConfigurationJacobian(const unsigned int &lowerIndex,
                                                             const unsigned int &upperIndex);

        floatVector getPreviousPrecedingMicroConfigurationJacobian(const unsigned int &index);

        floatVector getPreviousFollowingMicroConfigurationJacobian(const unsigned int &index);

        //! Get a reference to the scaled value of the micro deformation
        const secondOrderTensor *getScaledMicroDeformation() { return &_scaled_microDeformation; }

        //! Get a reference to the scaled value of the gradient of the micro deformation
        const thirdOrderTensor *getScaledGradientMicroDeformation() { return &_scaled_gradientMicroDeformation; }

        const floatVector *getFlatdXdD() {
            /*!
             * Get the total derivative of the unknown vector w.r.t. the deformation.
             *
             * Pass-through to getFlatdXdF just changing the naming convention
             */

            return hydraBase::getFlatdXdF();
        }

       protected:
        // Utility functions
        virtual void initializeUnknownVector() override;

        virtual void updateConfigurationsFromUnknownVector() override;

        //            virtual void decomposeUnknownVector( ) override;

        virtual void decomposeUnknownVectorMicroConfigurations();

        virtual void decomposeStateVariableVector() override;

        virtual void decomposeStateVariableVectorMicroConfigurations();

        virtual void setScaledQuantities() override;

       private:
        secondOrderTensor _scaled_microDeformation;  //!< The current micro-deformation scaled by the scale factor

        thirdOrderTensor _scaled_gradientMicroDeformation;  //!< The spatial gradient of the micro-deformation w.r.t.
                                                            //!< the reference coordinates scaled by the scale factor

        void setFirstMicroConfigurationJacobians();

        void setPreviousFirstMicroConfigurationJacobians();

        void setFirstGradientMicroConfigurationJacobians();

        void setPreviousFirstGradientMicroConfigurationJacobians();

        void computeGradientMicroConfigurations(const floatVector *data_vector, unsigned int start_index,
                                                const floatVector &configurations,
                                                const floatVector &microConfigurations,
                                                const floatVector &gradientMicroConfiguration,
                                                floatVector       &gradientMicroConfigurations);

        void calculateFirstConfigurationGradChi(const floatVector &configurations,
                                                const floatVector &microConfigurations,
                                                const floatVector &gradientMicroConfiguration,
                                                floatVector       &gradientMicroConfigurations);

        void calculateFirstConfigurationGradChiJacobian(
            const floatVector &configurations, const floatVector &microConfigurations,
            const floatVector &gradientMicroConfiguration, const floatVector &gradientMicroConfigurations,
            const floatVector &dChi1dChi, const floatVector &dChi1dChin, floatVector &dGradChi1dCn,
            fifthOrderTensor &dGradChi1dChi, floatVector &dGradChi1dChin, sixthOrderTensor &dGradChi1dGradChi,
            floatVector &dGradChi1dGradChin);

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, microConfigurations, floatVector, unexpectedError)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, inverseMicroConfigurations, floatVector, unexpectedError)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, gradientMicroConfigurations, floatVector, unexpectedError)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousMicroConfigurations, floatVector, unexpectedError)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousInverseMicroConfigurations, floatVector,
                                                  unexpectedError)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousGradientMicroConfigurations, floatVector,
                                                  unexpectedError)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dChi1dChi, fourthOrderTensor,
                                                   setFirstMicroConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dChi1dChin, floatVector,
                                                   setFirstMicroConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousdChi1dChi, fourthOrderTensor,
                                                  setPreviousFirstMicroConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousdChi1dChin, floatVector,
                                                  setPreviousFirstMicroConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dGradChi1dFn, floatVector,
                                                   setFirstGradientMicroConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dGradChi1dChi, fifthOrderTensor,
                                                   setFirstGradientMicroConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dGradChi1dChin, floatVector,
                                                   setFirstGradientMicroConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dGradChi1dGradChi, sixthOrderTensor,
                                                   setFirstGradientMicroConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, dGradChi1dGradChin, floatVector,
                                                   setFirstGradientMicroConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousdGradChi1dFn, floatVector,
                                                  setPreviousFirstGradientMicroConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousdGradChi1dChi, fifthOrderTensor,
                                                  setPreviousFirstGradientMicroConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousdGradChi1dChin, floatVector,
                                                  setPreviousFirstGradientMicroConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousdGradChi1dGradChi, sixthOrderTensor,
                                                  setPreviousFirstGradientMicroConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousdGradChi1dGradChin, floatVector,
                                                  setPreviousFirstGradientMicroConfigurationJacobians)
    };

    //! The base class for micromorphic residuals
    template <class container = hydraBaseMicromorphic>
    class ResidualBaseMicromorphic : public ResidualBase<hydraBase> {
       public:
        using tardigradeHydra::ResidualBase<hydraBase>::ResidualBase;

        /*!
         * Base class for micromorphic residuals
         *
         * \param *_hydra: A pointer to the containing hydra object
         * \param _numEquations: The number of equations the residual defines
         */
        ResidualBaseMicromorphic(hydraBaseMicromorphic *_hydra, unsigned int _numEquations)
            : ResidualBase(_hydra, _numEquations), hydra(_hydra) {}

        container *hydra;  //!< A pointer to the containing hydra object

        virtual void setdRdD() {
            /*!
             * Set the derivative of the residual w.r.t. the deformation.
             */

            TARDIGRADE_ERROR_TOOLS_CATCH(
                throw std::logic_error("The derivative of the residual w.r.t. the deformation is not implemented"));
        }

        virtual void setdRdF() override {
            /*!
             * Rename setdRdF to setdRdD because we will use it for all of the deformations
             */

            TARDIGRADE_ERROR_TOOLS_CATCH(setdRdD());
        }

        void setdRdD(const floatVector &dRdD) {
            /*!
             * Set the derivative of the residual w.r.t. the deformation.
             *
             * Pass-through to setdRdF just changing the naming convention
             *
             * \param &dRdD: The derivative of the resdual with respect to the deformation (F, chi, gradChi )
             */

            TARDIGRADE_ERROR_TOOLS_CATCH(ResidualBase<hydraBase>::setdRdF(dRdD));
        }

        const floatVector *getdRdD() {
            /*!
             * Get the derivative of the residual w.r.t. the deformation.
             *
             * Pass-through to getdRdF just changing the naming convention
             */

            return ResidualBase<hydraBase>::getdRdF();
        }

        tardigradeHydra::ResidualBase<hydraBase>::SetDataStorageIteration<floatVector> get_SetDataStorage_dRdD() {
            /*!
             * Get the setting term for dRdD
             *
             * Pass-through to get_SetDataStorage_dRdF just changing the naming convention
             */

            return ResidualBase<hydraBase>::get_SetDataStorage_dRdF();
        }
    };

}  // namespace tardigradeHydra

#include "tardigrade_hydraMicromorphic.tpp"

#endif
