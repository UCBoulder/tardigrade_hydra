/**
 ******************************************************************************
 * \file tardigrade_hydraMicromorphic.cpp
 ******************************************************************************
 * A C++ utility for constructing finite deformation micromorphic constitutive
 * models.
 ******************************************************************************
 */

#include <tardigrade_hydraMicromorphic.h>

namespace tardigradeHydra {

    /*!
     * The main constructor for the micromorphic hydra base class. Inputs are all the required values for most
     * solves.
     *
     * \param &DOFStorage: The degree of freedom storage class
     * \param &ModelConfiguration: The model configuration class
      \param &_hydra_configuration: Class which defines the hydra configuration
     */
    hydraBaseMicromorphic::hydraBaseMicromorphic(const MicromorphicDOFStorage &DOFStorage,
                                                 const ModelConfigurationBase &ModelConfiguration,
                                                 HydraConfigurationBase        _hydra_configuration)
        : hydraBase(DOFStorage, ModelConfiguration, _hydra_configuration) {}

    /*!
     * Initialize the hydra object
     */
    void hydraBaseMicromorphic::initialize() {
        tardigradeHydra::hydraBase::hydraBase::initialize();

        setScaledQuantities();

        decomposeStateVariableVectorMicroConfigurations();
    }

    void hydraBaseMicromorphic::setScaledQuantities() {
        /*!
         * Scale the current values by the scale factor
         */

        auto local_dof =
            static_cast<const tardigradeHydra::MicromorphicDOFStorage *>(dof);  // TODO: Avoid this static cast

        hydraBase::hydraBase::setScaledQuantities();

        _scaled_microDeformation =
            getScaleFactor() * (local_dof->_micro_deformation - local_dof->_previous_micro_deformation) +
            local_dof->_previous_micro_deformation;

        _scaled_gradientMicroDeformation = getScaleFactor() * (local_dof->_gradient_micro_deformation -
                                                               local_dof->_previous_gradient_micro_deformation) +
                                           local_dof->_previous_gradient_micro_deformation;
    }

    void hydraBaseMicromorphic::initializeUnknownVector() {
        /*!
         * Initialize the unknown vector for the non-linear solve.
         *
         * \f$X = \left\{ \bf{\sigma}, \bf{F}^2, \bf{F}^3, ..., \bf{F}n, \bf{\chi}^2, \bf{\chi}^3, ...,
         * \frac{\partial}{\partial \bf{X}^2} \bf{\chi}^2, \frac{\partial}{\partial \bf{X}^3} \bf{\chi}^3, ..., \xi^1,
         * \xi^2, ..., \xi^m \right\} \f$
         *
         * It is assumed that the first residual calculation also has a method `void getStress( )`
         * which returns a pointer to the current value of the stress.
         */

        constexpr unsigned int sot_dimension = configuration::dimension * configuration::dimension;
        constexpr unsigned int tot_dimension = configuration::dimension * configuration::dimension * configuration::dimension;

        const floatVector *stress;
        TARDIGRADE_ERROR_TOOLS_CATCH(stress = getStress());

        const floatVector *configurations = deformation->get_configurations();

        const floatVector *microConfigurations = get_microConfigurations();

        const floatVector *gradientMicroConfigurations = get_gradientMicroConfigurations();

        const floatVector *nonLinearSolveStateVariables = get_nonLinearSolveStateVariables();

        floatVector X(getNumUnknowns(), 0);

        unsigned int offset = 0;

        std::copy(std::begin(*stress), std::end(*stress), std::begin(X));

        offset += stress->size();

        // Set the initial values of the macro configurations
        std::copy(std::begin(*configurations) + sot_dimension, std::end(*configurations), std::begin(X) + offset);

        offset += configurations->size() - sot_dimension;

        // Set the values of the micro configurations
        std::copy(std::begin(*microConfigurations) + sot_dimension, std::end(*microConfigurations), std::begin(X) + offset);

        offset += microConfigurations->size() - sot_dimension;

        // Set the values of the micro-gradient configurations
        std::copy(std::begin(*gradientMicroConfigurations) + tot_dimension, std::end(*gradientMicroConfigurations),
                  std::begin(X) + offset);

        offset += gradientMicroConfigurations->size() - tot_dimension;

        std::copy(std::begin(*nonLinearSolveStateVariables), std::end(*nonLinearSolveStateVariables),
                  std::begin(X) + offset);

        bool resetRequired = false;

        setCurrentResidualIndexMeaningful(true);
        for (auto residual_ptr = getResidualClasses()->begin(); residual_ptr != getResidualClasses()->end();
             residual_ptr++) {
            setCurrentResidualIndex(residual_ptr - getResidualClasses()->begin());

            std::vector<unsigned int> indices;

            std::vector<floatType> values;

            (*residual_ptr)->suggestInitialIterateValues(indices, values);

            if (indices.size() > 0) {
                resetRequired = true;
            }

            for (auto i = indices.begin(); i != indices.end(); i++) {
                X[*i] = values[(unsigned int)(i - indices.begin())];
            }
        }
        setCurrentResidualIndexMeaningful(false);

        if (resetRequired) {
            updateUnknownVector(X);

        } else {
            setX(X);
        }
    }

    void hydraBaseMicromorphic::decomposeStateVariableVector() {
        /*!
         * Decompose the incoming state variable vector setting the different configurations along the way
         *
         * The state variable vector is assumed to be of the form:
         *
         * \f$ \text{ISV} = \left\{\bf{F}^2 - \bf{I}, \bf{F}^3 - \bf{I}, \cdots, \bf{F}^n - \bf{I}, \bf{\chi}^2 -
         * \bf{I}, \bf{\chi}^3 - \bf{I}, \cdots, \bf{\chi}^n - \bf{I} \frac{\partial}{\partial \bf{X}} \bf{\chi}^2,
         * \frac{\partial}{\partial \bf{X}} \bf{\chi}^3, \cdots, \frac{\partial}{\partial \bf{X}} \bf{\chi}^n, \xi^1,
         * \xi^2, \cdots, \xi^m, \eta^1, \cdots\right\} \f$
         *
         * where the \f$\bf{F}\f$ are the different deformation gradients (configurations), \f$\bf{\chi}\f$ are the
         * micro-deformations, \f$\xi^y\f$ are the other variables to be solved during the non-linear solve, and
         * \f$\eta^z\f$ are other state variables. Note that we decompose the deformation gradient and micro-deformation
         * as
         *
         * \f$\bf{F} = \bf{F}^1 \bf{F}^2 \cdots \bf{F}^n\f$
         *
         * \f$\bf{\chi} = \bf{\chi}^1 \bf{\chi}^2 \cdots \bf{\chi}^n\f$
         *
         * and so because \f$\bf{F}\f$ and \f$\bf{\chi}\f$ are provided we can solve for \f$\bf{F}^1\f$ and
         * \f$\bf{\chi}\f$. Typically, this configuration would be the elastic configuration (i.e., the configuration
         * that generates the stress) though we do not insist that users follow convention.
         *
         * NOTE: Though we overload the decomposeStateVariableVector in hydraBase this function will not be called in
         * hydraBase's constructor because virtual functions do not come into being during the construction of parent
         * classes constructors. We could work around this but instead we will overload and define a local method to do
         * the decomposition of the micro-deformation tensors.
         */

        // Call the parent class decomposition
        hydraBase::decomposeStateVariableVector();

        // Decompose the micro-deformation
        decomposeStateVariableVectorMicroConfigurations();
    }

    void hydraBaseMicromorphic::updateConfigurationsFromUnknownVector() {
        /*!
         * Decompose the incoming unknown vector setting the different configurations along the way
         */

        hydraBase::updateConfigurationsFromUnknownVector();

        decomposeUnknownVectorMicroConfigurations();
    }

    //    void hydraBaseMicromorphic::decomposeUnknownVector( ){
    //        /*!
    //         * Decompose the incoming unknown vector setting the different configurations along the way
    //         *
    //         * The state variable vector is assumed to be of the form:
    //         *
    //         * \f$ \text{ISV} = \left\{\bf{S}^2, \bf{\Sigma}, \bf{M}, \bf{F}^2, \bf{F}^3, \cdots, \bf{F}^n,
    //         \bf{\chi}^2, \bf{\chi}^3, \cdots, \bf{\chi}^n \frac{\partial}{\partial \bf{X}} \bf{\chi}^2,
    //         \frac{\partial}{\partial \bf{X}} \bf{\chi}^3, \cdots, \frac{\partial}{\partial \bf{X}} \bf{\chi}^n,
    //         \xi^1, \xi^2, \cdots, \xi^m\right\} \f$
    //         *
    //         * where \f$\bf{S}^2\f$ is the second Piola Kirchhoff stress, \f$\bf{\Sigma}\f$ is the reference symmetric
    //         micro stress, and \f$\bf{M}\f$
    //         * is the reference higher order stress the \f$\bf{F}\f$ are the different deformation gradients
    //         (configurations), \f$\bf{\chi}\f$ are the micro-deformations,
    //         * \f$\xi^y\f$ are the other variables to be solved during the non-linear solve, and \f$\eta^z\f$ are
    //         other state variables. Note
    //         * that we decompose the deformation gradient and micro-deformation as
    //         *
    //         * \f$\bf{F} = \bf{F}^1 \bf{F}^2 \cdots \bf{F}^n\f$
    //         *
    //         * \f$\bf{\chi} = \bf{\chi}^1 \bf{\chi}^2 \cdots \bf{\chi}^n\f$
    //         *
    //         * and so because \f$\bf{F}\f$ and \f$\bf{\chi}\f$ are provided we can solve for \f$\bf{F}^1\f$ and
    //         \f$\bf{\chi}\f$. Typically,
    //         * this configuration would be the elastic configuration (i.e., the configuration that generates the
    //         stress) though we do not insist that users follow convention.
    //         *
    //         * NOTE: Though we overload the decomposeStateVariableVector in hydraBase this function will not be called
    //         in hydraBase's constructor because
    //         *       virtual functions do not come into being during the construction of parent classes constructors.
    //         We could work around this but instead
    //         *       we will overload and define a local method to do the decomposition of the micro-deformation
    //         tensors.
    //         */
    //
    //        // Call the parent class decomposition
    //        hydraBase::decomposeUnknownVector( );
    //
    //        // Decompose the micro-deformation
    //        decomposeUnknownVectorMicroConfigurations( );
    //
    //    }

    void hydraBaseMicromorphic::computeGradientMicroConfigurations(const floatVector *data_vector,
                                                                   unsigned int       start_index,
                                                                   const floatVector &configurations,
                                                                   const floatVector &microConfigurations,
                                                                   const floatVector &gradientMicroConfiguration,
                                                                   floatVector       &gradientMicroConfigurations) {
        /*!
         * Compute the gradient of the micro-configurations in their reference configurations
         *
         * \param *data_vector: The vector of data
         * \param &start_index: The index at which to start reading the data from data_vector
         * \param &configurations: The macro-scale configurations
         * \param &microConfigurations: The micro-configurations
         * \param &gradientMicroConfiguration: The gradient of the total micro-configuration w.r.t. the reference
         * configuration \param &gradientMicroConfigurations: The resulting gradients of the micro configurations
         */

        constexpr unsigned int tot_dimension     = configuration::dimension * configuration::dimension * configuration::dimension;
        auto num_configs = getNumConfigurations();

        gradientMicroConfigurations =
            tardigradeVectorTools::appendVectors({thirdOrderTensor(tot_dimension, 0),
                                                  floatVector(data_vector->begin() + start_index,
                                                              data_vector->begin() + start_index +
                                                                  (num_configs - 1) * tot_dimension)});

        calculateFirstConfigurationGradChi(configurations, microConfigurations, gradientMicroConfiguration,
                                           gradientMicroConfigurations);
    }

    void hydraBaseMicromorphic::decomposeUnknownVectorMicroConfigurations() {
        /*!
         * Decompose the micro-deformation parts of the unknown vector
         */

        constexpr unsigned int sot_dimension     = configuration::dimension * configuration::dimension;
        auto num_configs = getNumConfigurations();

        unsigned int start_index = getStressSize() + (num_configs - 1) * sot_dimension;

        auto microConfigurations = get_SetDataStorage_microConfigurations();

        auto inverseMicroConfigurations = get_SetDataStorage_inverseMicroConfigurations();

        auto gradientMicroConfigurations = get_SetDataStorage_gradientMicroConfigurations();

        // Compute the micro-configurations

        computeConfigurations(getUnknownVector(), start_index, *getMicroDeformation(), *microConfigurations.value,
                              *inverseMicroConfigurations.value);

        start_index += (num_configs - 1) * sot_dimension;

        computeGradientMicroConfigurations(getUnknownVector(), start_index, *deformation->get_configurations(),
                                           *microConfigurations.value, *getGradientMicroDeformation(),
                                           *gradientMicroConfigurations.value);
    }

    void hydraBaseMicromorphic::decomposeStateVariableVectorMicroConfigurations() {
        /*!
         * Decompose the micro-deformation parts of the state variable vector
         */

        constexpr unsigned int sot_dimension     = configuration::dimension * configuration::dimension;
        auto num_configs = getNumConfigurations();

        unsigned int start_index = (num_configs - 1) * sot_dimension;

        auto microConfigurations = get_SetDataStorage_microConfigurations();

        auto inverseMicroConfigurations = get_SetDataStorage_inverseMicroConfigurations();

        auto gradientMicroConfigurations = get_SetDataStorage_gradientMicroConfigurations();

        auto previousMicroConfigurations = get_SetDataStorage_previousMicroConfigurations();

        auto previousInverseMicroConfigurations = get_SetDataStorage_previousInverseMicroConfigurations();

        auto previousGradientMicroConfigurations = get_SetDataStorage_previousGradientMicroConfigurations();

        // Compute the micro-configurations

        computeConfigurations(getPreviousStateVariables(), start_index, *getMicroDeformation(),
                              *microConfigurations.value, *inverseMicroConfigurations.value, true);

        computeConfigurations(getPreviousStateVariables(), start_index, *getPreviousMicroDeformation(),
                              *previousMicroConfigurations.value, *previousInverseMicroConfigurations.value, true);

        start_index += (num_configs - 1) * sot_dimension;

        computeGradientMicroConfigurations(getPreviousStateVariables(), start_index, *deformation->get_configurations(),
                                           *microConfigurations.value, *getGradientMicroDeformation(),
                                           *gradientMicroConfigurations.value);

        computeGradientMicroConfigurations(getPreviousStateVariables(), start_index,
                                           *deformation->get_previousConfigurations(),
                                           *previousMicroConfigurations.value, *getPreviousGradientMicroDeformation(),
                                           *previousGradientMicroConfigurations.value);
    }

    secondOrderTensor hydraBaseMicromorphic::getSubMicroConfiguration(const unsigned int &lowerIndex,
                                                                      const unsigned int &upperIndex) {
        /*!
         * Get a sub-micro configuration \f$\bf{\chi}^{sc}\f$ defined as
         *
         * \f$ \chi^{sc}_{iI} = \chi^{\text{lowerIndex}}_{i\hat{I}} \chi^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}}
         * \cdots \chi^{\text{upperIndex-1}}_{\bar{I}I} \f$ \param &lowerIndex: The index of the lower configuration
         * (starts at 0 and goes to numConfigurations - 1) \param &upperIndex: The index of the upper configuration
         * (starts at 0 and goes to numConfigurations) Note, the configuration indicated by the index is NOT included in
         * the sub-configuration
         */

        return deformation->getSubConfiguration<3, 3, 3>(*get_microConfigurations(), lowerIndex, upperIndex);
    }

    secondOrderTensor hydraBaseMicromorphic::getPrecedingMicroConfiguration(const unsigned int &index) {
        /*!
         * Get the sub-micro configuration preceding but not including the index
         *
         * \param &index: The index of the configuration immediately following the sub-micro configuration
         */

        return getSubMicroConfiguration(0, index);
    }

    secondOrderTensor hydraBaseMicromorphic::getFollowingMicroConfiguration(const unsigned int &index) {
        /*!
         * Get the sub-micro configuration following but not including the index
         *
         * \param &index: The index of the current configuration immediately before the sub-micro configuration
         */

        return getSubMicroConfiguration(index + 1, getNumConfigurations());
    }

    secondOrderTensor hydraBaseMicromorphic::getMicroConfiguration(const unsigned int &index) {
        /*!
         * Get the micro configuration indicated by the provided index
         *
         * \param &index: The index of the current configuration to be extracted
         */

        return getSubMicroConfiguration(index, index + 1);
    }

    secondOrderTensor hydraBaseMicromorphic::getPreviousSubMicroConfiguration(const unsigned int &lowerIndex,
                                                                              const unsigned int &upperIndex) {
        /*!
         * Get a previous sub-micro configuration \f$\bf{\chi}^{sc}\f$ defined as
         *
         * \f$ \chi^{sc}_{iI} = \chi^{\text{lowerIndex}}_{i\hat{I}} \chi^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}}
         * \cdots \chi^{\text{upperIndex-1}}_{\bar{I}I} \f$ \param &lowerIndex: The index of the lower configuration
         * (starts at 0 and goes to numConfigurations - 1) \param &upperIndex: The index of the upper configuration
         * (starts at 0 and goes to numConfigurations) Note, the configuration indicated by the index is NOT included in
         * the sub-configuration
         */

        return deformation->getSubConfiguration<3, 3, 3>(*get_previousMicroConfigurations(), lowerIndex, upperIndex);
    }

    secondOrderTensor hydraBaseMicromorphic::getPreviousPrecedingMicroConfiguration(const unsigned int &index) {
        /*!
         * Get the previous sub-micro configuration preceding but not including the index
         *
         * \param &index: The index of the configuration immediately following the sub-micro configuration
         */

        return getPreviousSubMicroConfiguration(0, index);
    }

    secondOrderTensor hydraBaseMicromorphic::getPreviousFollowingMicroConfiguration(const unsigned int &index) {
        /*!
         * Get the previous sub-micro configuration following but not including the index
         *
         * \param &index: The index of the current configuration immediately before the sub-micro configuration
         */

        return getPreviousSubMicroConfiguration(index + 1, getNumConfigurations());
    }

    secondOrderTensor hydraBaseMicromorphic::getPreviousMicroConfiguration(const unsigned int &index) {
        /*!
         * Get the previous micro configuration indicated by the provided index
         *
         * \param &index: The index of the current configuration to be extracted
         */

        return getPreviousSubMicroConfiguration(index, index + 1);
    }

    floatVector hydraBaseMicromorphic::getSubMicroConfigurationJacobian(const unsigned int &lowerIndex,
                                                                        const unsigned int &upperIndex) {
        /*!
         * Get the jacobian of a sub-micro configuration \f$\bf{\chi}^{sc}\f$ defined as
         *
         * \f$ \chi^{sc}_{iI} = \chi^{\text{lowerIndex}}_{i\hat{I}} \chi^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}}
         * \cdots \chi^{\text{upperIndex-1}}_{\bar{I}I} \f$
         *
         * with respect to the current configurations.
         *
         * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
         * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
         *   Note, the configuration indicated by the index is NOT included in the sub-configuration
         */

        return deformation->getSubConfigurationJacobian<3, 3, 3>(*get_microConfigurations(), lowerIndex, upperIndex);
    }

    floatVector hydraBaseMicromorphic::getPrecedingMicroConfigurationJacobian(const unsigned int &index) {
        /*!
         * Get the jacobian of the sub-micro configuration preceding but not including the index with respect to the
         * current configurations.
         *
         * \param &index: The index of the configuration immediately following the sub-configuration
         */

        return getSubMicroConfigurationJacobian(0, index);
    }

    floatVector hydraBaseMicromorphic::getFollowingMicroConfigurationJacobian(const unsigned int &index) {
        /*!
         * Get the jacobian of the sub-micro configuration following but not including the index with respect to the
         * current configurations.
         *
         * \param &index: The index of the current configuration immediately before the sub-configuration
         */

        return getSubMicroConfigurationJacobian(index + 1, getNumConfigurations());
    }

    floatVector hydraBaseMicromorphic::getPreviousSubMicroConfigurationJacobian(const unsigned int &lowerIndex,
                                                                                const unsigned int &upperIndex) {
        /*!
         * Get the jacobian of a previous sub-micro configuration \f$\bf{\chi}^{sc}\f$ defined as
         *
         * \f$ \chi^{sc}_{iI} = \chi^{\text{lowerIndex}}_{i\hat{I}} \chi^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}}
         * \cdots \chi^{\text{upperIndex-1}}_{\bar{I}I} \f$
         *
         * with respect to the current configurations.
         *
         * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
         * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
         *   Note, the configuration indicated by the index is NOT included in the sub-configuration
         */

        return deformation->getSubConfigurationJacobian<3, 3, 3>(*get_previousMicroConfigurations(), lowerIndex,
                                                                 upperIndex);
    }

    floatVector hydraBaseMicromorphic::getPreviousPrecedingMicroConfigurationJacobian(const unsigned int &index) {
        /*!
         * Get the jacobian of the previous sub-micro configuration preceding but not including the index with respect
         * to the current configurations.
         *
         * \param &index: The index of the configuration immediately following the sub-configuration
         */

        return getPreviousSubMicroConfigurationJacobian(0, index);
    }

    floatVector hydraBaseMicromorphic::getPreviousFollowingMicroConfigurationJacobian(const unsigned int &index) {
        /*!
         * Get the jacobian of the previous sub-micro configuration following but not including the index with respect
         * to the current configurations.
         *
         * \param &index: The index of the current configuration immediately before the sub-configuration
         */

        return getPreviousSubMicroConfigurationJacobian(index + 1, getNumConfigurations());
    }

    void hydraBaseMicromorphic::setFirstMicroConfigurationJacobians() {
        /*!
         * Set the Jacobians of the first micro configuration w.r.t. the total micro configuration and the remaining
         * sub-micro configurations
         */

        auto dChi1dChi = get_SetDataStorage_dChi1dChi();

        auto dChi1dChin = get_SetDataStorage_dChi1dChin();

        deformation->calculateFirstConfigurationJacobians(*get_microConfigurations(), *dChi1dChi.value,
                                                          *dChi1dChin.value);
    }

    void hydraBaseMicromorphic::setPreviousFirstMicroConfigurationJacobians() {
        /*!
         * Set the Jacobians of the previous first micro configuration w.r.t. the total micro configuration and the
         * remaining sub-micro configurations
         */

        auto previousdChi1dChi = get_SetDataStorage_previousdChi1dChi();

        auto previousdChi1dChin = get_SetDataStorage_previousdChi1dChin();

        deformation->calculateFirstConfigurationJacobians(*get_previousMicroConfigurations(), *previousdChi1dChi.value,
                                                          *previousdChi1dChin.value);
    }

    void hydraBaseMicromorphic::setFirstGradientMicroConfigurationJacobians() {
        /*!
         * Set the Jacobians of the gradient of the first micro configuration w.r.t. the total micro configuration and
         * the remaining sub-micro configurations
         */

        auto dGradChi1dFn = get_SetDataStorage_dGradChi1dFn();

        auto dGradChi1dChi = get_SetDataStorage_dGradChi1dChi();

        auto dGradChi1dChin = get_SetDataStorage_dGradChi1dChin();

        auto dGradChi1dGradChi = get_SetDataStorage_dGradChi1dGradChi();

        auto dGradChi1dGradChin = get_SetDataStorage_dGradChi1dGradChin();

        calculateFirstConfigurationGradChiJacobian(*deformation->get_configurations(), *get_microConfigurations(),
                                                   *getGradientMicroDeformation(), *get_gradientMicroConfigurations(),
                                                   *get_dChi1dChi(), *get_dChi1dChin(), *dGradChi1dFn.value,
                                                   *dGradChi1dChi.value, *dGradChi1dChin.value,
                                                   *dGradChi1dGradChi.value, *dGradChi1dGradChin.value);
    }

    void hydraBaseMicromorphic::setPreviousFirstGradientMicroConfigurationJacobians() {
        /*!
         * Set the Jacobians of the previous gradient of the first micro configuration w.r.t. the total micro
         * configuration and the remaining sub-micro configurations
         */

        auto previousdGradChi1dFn = get_SetDataStorage_previousdGradChi1dFn();

        auto previousdGradChi1dChi = get_SetDataStorage_previousdGradChi1dChi();

        auto previousdGradChi1dChin = get_SetDataStorage_previousdGradChi1dChin();

        auto previousdGradChi1dGradChi = get_SetDataStorage_previousdGradChi1dGradChi();

        auto previousdGradChi1dGradChin = get_SetDataStorage_previousdGradChi1dGradChin();

        calculateFirstConfigurationGradChiJacobian(*deformation->get_previousConfigurations(),
                                                   *get_previousMicroConfigurations(),
                                                   *getPreviousGradientMicroDeformation(),
                                                   *get_previousGradientMicroConfigurations(), *get_previousdChi1dChi(),
                                                   *get_previousdChi1dChin(), *previousdGradChi1dFn.value,
                                                   *previousdGradChi1dChi.value, *previousdGradChi1dChin.value,
                                                   *previousdGradChi1dGradChi.value, *previousdGradChi1dGradChin.value);
    }

    void hydraBaseMicromorphic::calculateFirstConfigurationGradChi(const floatVector      &configurations,
                                                                   const floatVector      &microConfigurations,
                                                                   const thirdOrderTensor &gradientMicroConfiguration,
                                                                   floatVector &gradientMicroConfigurations) {
        /*!
         * Calculate the value of the gradient of the first micro-configuration given all of the configurations, the
         * micro-configurations, the spatial gradient of the micro deformation in the reference configuration, and the
         * gradients of the micro-configurations other than the first in their own reference configurations.
         *
         * \param &configurations: The configuration matrix
         * \param &microConfigurations: The micro-configuration matrix
         * \param &gradientMicroConfiguration: The gradient of the micro-deformation in the reference configuration
         * \param &gradientMicroConfigurations: The matrix of gradients of the micro-configurations in their reference
         * configurations
         */

        constexpr unsigned int sot_dimension     = configuration::dimension * configuration::dimension;
        constexpr unsigned int tot_dimension     = configuration::dimension * configuration::dimension * configuration::dimension;
        auto num_configs = getNumConfigurations();

        // Compute the gradient in the reference configuration
        thirdOrderTensor gradientChi1Reference(
            tot_dimension,
            0);  // = gradientMicroConfiguration; // Initialize to the total gradient in the reference configuration

        thirdOrderTensor temp_tot1(tot_dimension, 0);

        secondOrderTensor chiPrecede(sot_dimension, 0.);
        for (unsigned int i = 0; i < configuration::dimension; i++) {
            chiPrecede[configuration::dimension * i + i] = 1.;
        }

        Eigen::Map<Eigen::Matrix<floatType, configuration::dimension, configuration::dimension, Eigen::RowMajor> >       chip_map(chiPrecede.data(), configuration::dimension, configuration::dimension);
        Eigen::Map<const Eigen::Matrix<floatType, configuration::dimension, configuration::dimension, Eigen::RowMajor> > map1(NULL, configuration::dimension, configuration::dimension);

        for (unsigned int index = 1; index < num_configs; index++) {
            new (&map1) Eigen::Map<const Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> >(microConfigurations.data() +
                                                                                               (index - 1) * sot_dimension,
                                                                                           3, 3);
            chip_map *= map1;

            // Add the contribution of the term
            for (unsigned int i = 0; i < configuration::dimension; i++) {
                for (unsigned int j = 0; j < configuration::dimension; j++) {
                    for (unsigned int I = 0; I < configuration::dimension; I++) {
                        for (unsigned int J = 0; J < configuration::dimension; J++) {
                            gradientChi1Reference[configuration::dimension * configuration::dimension * i + configuration::dimension * I + J] -=
                                chiPrecede[configuration::dimension * i + j] *
                                gradientMicroConfigurations[tot_dimension * index + configuration::dimension * configuration::dimension * j + configuration::dimension * I + J];
                        }
                    }
                }
            }

            if (index == (num_configs - 1)) {
                break;
            }

            //            std::copy( gradientChi1Reference.begin( ), gradientChi1Reference.end( ), temp_tot1.begin( ) );
            std::fill(temp_tot1.begin(), temp_tot1.end(), 0);

            secondOrderTensor chiFollow =
                deformation->getSubConfiguration<3, 3, 3>(microConfigurations, index + 1, getNumConfigurations());

            for (unsigned int i = 0; i < configuration::dimension; i++) {
                for (unsigned int I = 0; I < configuration::dimension; I++) {
                    for (unsigned int J = 0; J < configuration::dimension; J++) {
                        for (unsigned int k = 0; k < configuration::dimension; k++) {
                            temp_tot1[configuration::dimension * configuration::dimension * i + configuration::dimension * I + J] +=
                                chiFollow[configuration::dimension * k + I] * gradientChi1Reference[configuration::dimension * configuration::dimension * i + configuration::dimension * k + J];
                        }
                    }
                }
            }

            //            std::copy( gradientChi1Reference.begin( ), gradientChi1Reference.end( ), temp_tot1.begin( ) );
            std::fill(gradientChi1Reference.begin(), gradientChi1Reference.end(), 0);

            secondOrderTensor FFollow =
                deformation->getSubConfiguration<3, 3, 3>(configurations, index + 1, getNumConfigurations());

            for (unsigned int i = 0; i < configuration::dimension; i++) {
                for (unsigned int I = 0; I < configuration::dimension; I++) {
                    for (unsigned int l = 0; l < configuration::dimension; l++) {
                        for (unsigned int J = 0; J < configuration::dimension; J++) {
                            gradientChi1Reference[configuration::dimension * configuration::dimension * i + configuration::dimension * I + J] +=
                                FFollow[configuration::dimension * l + J] * temp_tot1[configuration::dimension * configuration::dimension * i + configuration::dimension * I + l];
                        }
                    }
                }
            }
        }

        gradientChi1Reference += gradientMicroConfiguration;

        // Map the gradient of the micro-configuration to the reference of the first configuration
        secondOrderTensor invChiFollow = deformation->getSubConfiguration<3, 3, 3>(microConfigurations, 1, num_configs);
        Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> > mat(invChiFollow.data(), 3, 3);
        mat = mat.inverse().eval();

        secondOrderTensor invFFollow = deformation->getSubConfiguration<3, 3, 3>(configurations, 1, num_configs);
        new (&mat) Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> >(invFFollow.data(), 3, 3);
        mat = mat.inverse().eval();

        std::fill(gradientMicroConfigurations.begin(), gradientMicroConfigurations.begin() + tot_dimension, 0.);

        std::fill(temp_tot1.begin(), temp_tot1.end(), 0);

        for (unsigned int i = 0; i < configuration::dimension; i++) {
            for (unsigned int a = 0; a < configuration::dimension; a++) {
                for (unsigned int I = 0; I < configuration::dimension; I++) {
                    for (unsigned int J = 0; J < configuration::dimension; J++) {
                        temp_tot1[configuration::dimension * configuration::dimension * i + configuration::dimension * I + J] +=
                            gradientChi1Reference[sot_dimension * i + configuration::dimension * a + J] * invChiFollow[configuration::dimension * a + I];
                    }
                }
            }
        }

        for (unsigned int i = 0; i < configuration::dimension; i++) {
            for (unsigned int I = 0; I < configuration::dimension; I++) {
                for (unsigned int b = 0; b < configuration::dimension; b++) {
                    for (unsigned int J = 0; J < configuration::dimension; J++) {
                        gradientMicroConfigurations[configuration::dimension * configuration::dimension * i + configuration::dimension * I + J] +=
                            temp_tot1[sot_dimension * i + configuration::dimension * I + b] * invFFollow[configuration::dimension * b + J];
                    }
                }
            }
        }
    }

    void hydraBaseMicromorphic::calculateFirstConfigurationGradChiJacobian(
        const floatVector &configurations, const floatVector &microConfigurations,
        const thirdOrderTensor &gradientMicroConfiguration, const floatVector &gradientMicroConfigurations,
        const fourthOrderTensor &dChi1dChi, const floatVector &dChi1dChin, floatVector &dGradChi1dCn,
        fifthOrderTensor &dGradChi1dChi, floatVector &dGradChi1dChin, sixthOrderTensor &dGradChi1dGradChi,
        floatVector &dGradChi1dGradChin) {
        /*!
         * Calculate the value of the jacobian of the gradient of the first micro-configuration given all of the
         * configurations, the micro-configurations, the spatial gradient of the micro deformation in the reference
         * configuration, and the gradients of the micro-configurations other than the first in their own reference
         * configurations.
         *
         * \param &configurations: The configuration matrix
         * \param &microConfigurations: The micro-configuration matrix
         * \param &gradientMicroConfiguration: The gradient of the micro-deformation in the reference configuration
         * \param &gradientMicroConfigurations: The gradient of the micro-deformations in their reference configurations
         * \param &dChi1dChi: The gradient of the first micro sub-configuration w.r.t. the total micro deformation
         * \param &dChi1dChin: The gradient of the first micro sub-configuration w.r.t. the remaining sub-micro
         * configurations \param &dGradChi1dCn: The Jacobian of the gradient of the first micro-configuration w.r.t. the
         * remaining configurations \param &dGradChi1dChi: The Jacobian of the gradient of the first micro-configuration
         * w.r.t. the total micro-configuration \param &dGradChi1dChin: The Jacobian of the gradient of the first
         * micro-configuration w.r.t. the remaining micro-configurations \param &dGradChi1dGradChi: The Jacobian of the
         * gradient of the first micro-configuration w.r.t. the gradient of the total micro-configuration \param
         * &dGradChi1dGradChin: The Jacobian of the gradient of the first micro-configuration w.r.t. the gradient of the
         * remaining sub micro-configurations
         */

        constexpr unsigned int sot_dimension     = configuration::dimension * configuration::dimension;
        constexpr unsigned int tot_dimension     = configuration::dimension * configuration::dimension * configuration::dimension;
        auto num_configs = getNumConfigurations();

        // Compute the gradient in the reference configuration
        thirdOrderTensor gradientChi1Reference(tot_dimension, 0);

        floatVector dGradientChi1ReferencedCn(tot_dimension * (num_configs - 1) * sot_dimension, 0);

        fifthOrderTensor dGradientChi1ReferencedChi(tot_dimension * sot_dimension, 0);

        floatVector dGradientChi1ReferencedChin(tot_dimension * (num_configs - 1) * sot_dimension, 0);

        floatVector dGradientChi1ReferencedGradChin(tot_dimension * (num_configs - 1) * tot_dimension, 0);

        thirdOrderTensor temp_tot1(tot_dimension, 0);

        thirdOrderTensor temp_tot2(tot_dimension, 0);

        thirdOrderTensor temp_tot2a(tot_dimension, 0);

        thirdOrderTensor temp_tot3(tot_dimension, 0);

        thirdOrderTensor temp_tot3a(tot_dimension, 0);

        secondOrderTensor chiPrecede(sot_dimension, 0), chiFollow(sot_dimension, 0), FFollow(sot_dimension, 0);

        for (unsigned int index = 1; index < num_configs; index++) {
            chiPrecede = deformation->getSubConfiguration<3, 3, 3>(microConfigurations, 0, index);

            // Set the Jacobians of the mapping terms
            floatVector dFFollowdCs =
                deformation->getSubConfigurationJacobian<3, 3, 3>(configurations, index + 1, num_configs);

            floatVector dChiPrecededChis =
                deformation->getSubConfigurationJacobian<3, 3, 3>(microConfigurations, 0, index);

            fourthOrderTensor dChiPrecededChi(sot_dimension * sot_dimension, 0);

            floatVector dChiPrecededChin(sot_dimension * (num_configs - 1) * sot_dimension, 0);

            for (unsigned int i = 0; i < sot_dimension; i++) {
                for (unsigned int k = 0; k < sot_dimension; k++) {
                    for (unsigned int j = 0; j < sot_dimension; j++) {
                        dChiPrecededChi[sot_dimension * i + j] +=
                            dChiPrecededChis[num_configs * sot_dimension * i + k] * dChi1dChi[sot_dimension * k + j];
                    }
                }
                for (unsigned int j = 0; j < sot_dimension * (num_configs - 1); j++) {
                    dChiPrecededChin[(num_configs - 1) * sot_dimension * i + j] +=
                        dChiPrecededChis[num_configs * sot_dimension * i + j + sot_dimension];

                    for (unsigned int k = 0; k < sot_dimension; k++) {
                        dChiPrecededChin[(num_configs - 1) * sot_dimension * i + j] +=
                            dChiPrecededChis[num_configs * sot_dimension * i + k] *
                            dChi1dChin[(num_configs - 1) * sot_dimension * k + j];
                    }
                }
            }

            floatVector dChiFollowdChis =
                deformation->getSubConfigurationJacobian<3, 3, 3>(microConfigurations, index + 1, num_configs);

            // Add the contribution of the term
            for (unsigned int i = 0; i < configuration::dimension; i++) {
                for (unsigned int j = 0; j < configuration::dimension; j++) {
                    for (unsigned int I = 0; I < configuration::dimension; I++) {
                        for (unsigned int J = 0; J < configuration::dimension; J++) {
                            gradientChi1Reference[configuration::dimension * configuration::dimension * i + configuration::dimension * I + J] -=
                                chiPrecede[configuration::dimension * i + j] *
                                gradientMicroConfigurations[tot_dimension * index + configuration::dimension * configuration::dimension * j + configuration::dimension * I + J];
                        }
                    }
                }
            }

            std::copy(gradientChi1Reference.begin(), gradientChi1Reference.end(), temp_tot3.begin());

            if (index != (num_configs - 1)) {
                std::copy(gradientChi1Reference.begin(), gradientChi1Reference.end(), temp_tot1.begin());
                std::fill(gradientChi1Reference.begin(), gradientChi1Reference.end(), 0);
                std::fill(temp_tot2.begin(), temp_tot2.end(), 0);

                chiFollow =
                    deformation->getSubConfiguration<3, 3, 3>(microConfigurations, index + 1, getNumConfigurations());

                for (unsigned int i = 0; i < configuration::dimension; i++) {
                    for (unsigned int I = 0; I < configuration::dimension; I++) {
                        for (unsigned int J = 0; J < configuration::dimension; J++) {
                            for (unsigned int k = 0; k < configuration::dimension; k++) {
                                gradientChi1Reference[configuration::dimension * configuration::dimension * i + configuration::dimension * I + J] +=
                                    chiFollow[configuration::dimension * k + I] * temp_tot1[configuration::dimension * configuration::dimension * i + configuration::dimension * k + J];
                                temp_tot2[configuration::dimension * configuration::dimension * i + configuration::dimension * I + J] -=
                                    chiFollow[configuration::dimension * k + I] *
                                    gradientMicroConfigurations[tot_dimension * index + configuration::dimension * configuration::dimension * i + configuration::dimension * k + J];
                            }
                        }
                    }
                }

                std::copy(gradientChi1Reference.begin(), gradientChi1Reference.end(), temp_tot1.begin());
                std::fill(gradientChi1Reference.begin(), gradientChi1Reference.end(), 0);
                std::fill(temp_tot2a.begin(), temp_tot2a.end(), 0);
                std::fill(temp_tot3a.begin(), temp_tot3a.end(), 0);

                FFollow = deformation->getSubConfiguration<3, 3, 3>(configurations, index + 1, getNumConfigurations());

                for (unsigned int i = 0; i < configuration::dimension; i++) {
                    for (unsigned int I = 0; I < configuration::dimension; I++) {
                        for (unsigned int l = 0; l < configuration::dimension; l++) {
                            for (unsigned int J = 0; J < configuration::dimension; J++) {
                                gradientChi1Reference[configuration::dimension * configuration::dimension * i + configuration::dimension * I + J] +=
                                    FFollow[configuration::dimension * l + J] * temp_tot1[configuration::dimension * configuration::dimension * i + configuration::dimension * I + l];

                                temp_tot2a[configuration::dimension * configuration::dimension * i + configuration::dimension * I + J] +=
                                    FFollow[configuration::dimension * l + J] * temp_tot2[configuration::dimension * configuration::dimension * i + configuration::dimension * I + l];

                                temp_tot3a[configuration::dimension * configuration::dimension * i + configuration::dimension * I + J] +=
                                    FFollow[configuration::dimension * l + J] * temp_tot3[configuration::dimension * configuration::dimension * i + configuration::dimension * I + l];

                                for (unsigned int A = 0; A < (num_configs - 1) * sot_dimension; A++) {
                                    dGradientChi1ReferencedCn[configuration::dimension * configuration::dimension * (num_configs - 1) * sot_dimension * i +
                                                              configuration::dimension * (num_configs - 1) * sot_dimension * I +
                                                              (num_configs - 1) * sot_dimension * J + A] +=
                                        dFFollowdCs[configuration::dimension * num_configs * sot_dimension * l + num_configs * sot_dimension * J +
                                                    sot_dimension + A] *
                                        temp_tot1[configuration::dimension * configuration::dimension * i + configuration::dimension * I + l];
                                }
                            }
                        }
                    }
                }

                std::copy(temp_tot2a.begin(), temp_tot2a.end(), temp_tot2.begin());
                std::copy(temp_tot3a.begin(), temp_tot3a.end(), temp_tot3.begin());

            } else {
                std::fill(chiFollow.begin(), chiFollow.end(), 0.);

                std::fill(FFollow.begin(), FFollow.end(), 0.);

                std::transform(gradientMicroConfigurations.begin() + tot_dimension * index,
                               gradientMicroConfigurations.begin() + tot_dimension * (index + 1), temp_tot2.begin(),
                               std::negate<floatType>());

                for (unsigned int i = 0; i < configuration::dimension; i++) {
                    chiFollow[configuration::dimension * i + i] = 1.;
                    FFollow[configuration::dimension * i + i]   = 1.;
                }

                for (unsigned int i = 0; i < configuration::dimension; i++) {
                    for (unsigned int I = 0; I < configuration::dimension; I++) {
                        for (unsigned int l = 0; l < configuration::dimension; l++) {
                            for (unsigned int J = 0; J < configuration::dimension; J++) {
                                for (unsigned int A = 0; A < (num_configs - 1) * sot_dimension; A++) {
                                    dGradientChi1ReferencedCn[configuration::dimension * configuration::dimension * (num_configs - 1) * sot_dimension * i +
                                                              configuration::dimension * (num_configs - 1) * sot_dimension * I +
                                                              (num_configs - 1) * sot_dimension * J + A] +=
                                        dFFollowdCs[configuration::dimension * num_configs * sot_dimension * l + num_configs * sot_dimension * J +
                                                    sot_dimension + A] *
                                        gradientChi1Reference[configuration::dimension * configuration::dimension * i + configuration::dimension * I + l];
                                }
                            }
                        }
                    }
                }
            }

            for (unsigned int i = 0; i < configuration::dimension; i++) {
                for (unsigned int I = 0; I < configuration::dimension; I++) {
                    for (unsigned int J = 0; J < configuration::dimension; J++) {
                        for (unsigned int j = 0; j < configuration::dimension; j++) {
                            for (unsigned int A = 0; A < sot_dimension; A++) {
                                dGradientChi1ReferencedChi[configuration::dimension * configuration::dimension * sot_dimension * i + configuration::dimension * sot_dimension * I + sot_dimension * J +
                                                           A] += dChiPrecededChi[configuration::dimension * sot_dimension * i + sot_dimension * j + A] *
                                                                 temp_tot2[configuration::dimension * configuration::dimension * j + configuration::dimension * I + J];
                            }

                            for (unsigned int A = 0; A < (num_configs - 1) * sot_dimension; A++) {
                                dGradientChi1ReferencedChin[configuration::dimension * configuration::dimension * (num_configs - 1) * sot_dimension * i +
                                                            configuration::dimension * (num_configs - 1) * sot_dimension * I +
                                                            (num_configs - 1) * sot_dimension * J + A] +=
                                    dChiPrecededChin[configuration::dimension * (num_configs - 1) * sot_dimension * i +
                                                     (num_configs - 1) * sot_dimension * j + A] *
                                        temp_tot2[configuration::dimension * configuration::dimension * j + configuration::dimension * I + J] +
                                    dChiFollowdChis[configuration::dimension * num_configs * sot_dimension * j + num_configs * sot_dimension * I + A +
                                                    sot_dimension] *
                                        temp_tot3[configuration::dimension * configuration::dimension * i + configuration::dimension * j + J];
                            }

                            for (unsigned int k = 0; k < configuration::dimension; k++) {
                                for (unsigned int l = 0; l < configuration::dimension; l++) {
                                    dGradientChi1ReferencedGradChin[configuration::dimension * configuration::dimension * (num_configs - 1) * tot_dimension * i +
                                                                    configuration::dimension * (num_configs - 1) * tot_dimension * I +
                                                                    (num_configs - 1) * tot_dimension * J +
                                                                    tot_dimension * (index - 1) + configuration::dimension * configuration::dimension * j + configuration::dimension * k +
                                                                    l] -=
                                        chiPrecede[configuration::dimension * i + j] * chiFollow[configuration::dimension * k + I] * FFollow[configuration::dimension * l + J];
                                }
                            }
                        }
                    }
                }
            }
        }

        gradientChi1Reference += gradientMicroConfiguration;

        // Map the gradient of the micro-configuration to the reference of the first configuration
        chiFollow = deformation->getSubConfiguration<3, 3, 3>(microConfigurations, 1, num_configs);

        FFollow = deformation->getSubConfiguration<3, 3, 3>(configurations, 1, num_configs);

        secondOrderTensor                                            invChiFollow = chiFollow;
        Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> > mat(invChiFollow.data(), 3, 3);
        mat                              = mat.inverse().eval();
        secondOrderTensor invChiFollow_T = invChiFollow;
        new (&mat) Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> >(invChiFollow_T.data(), 3, 3);
        mat = mat.transpose().eval();

        secondOrderTensor invFFollow = FFollow;
        new (&mat) Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> >(invFFollow.data(), 3, 3);
        mat                            = mat.inverse().eval();
        secondOrderTensor invFFollow_T = invFFollow;
        new (&mat) Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> >(invFFollow_T.data(), 3, 3);
        mat = mat.transpose().eval();

        floatVector dChiFollowdChis =
            deformation->getSubConfigurationJacobian<3, 3, 3>(microConfigurations, 1, num_configs);

        floatVector dFFollowdFs = deformation->getSubConfigurationJacobian<3, 3, 3>(configurations, 1, num_configs);

        fourthOrderTensor dInvChiFollowdChiFollow = tardigradeVectorTools::computeFlatDInvADA(invChiFollow, configuration::dimension, configuration::dimension);

        fourthOrderTensor dInvFFollowdFFollow = tardigradeVectorTools::computeFlatDInvADA(invFFollow, configuration::dimension, configuration::dimension);

        floatVector dInvChiFollowdChin(sot_dimension * sot_dimension * (num_configs - 1), 0);

        floatVector dInvFFollowdFn(sot_dimension * sot_dimension * (num_configs - 1), 0);

        for (unsigned int i = 0; i < sot_dimension; i++) {
            for (unsigned int k = 0; k < sot_dimension; k++) {
                for (unsigned int j = 0; j < (num_configs - 1) * sot_dimension; j++) {
                    dInvChiFollowdChin[(num_configs - 1) * sot_dimension * i + j] +=
                        dInvChiFollowdChiFollow[sot_dimension * i + k] *
                        dChiFollowdChis[num_configs * sot_dimension * k + j + sot_dimension];

                    dInvFFollowdFn[(num_configs - 1) * sot_dimension * i + j] +=
                        dInvFFollowdFFollow[sot_dimension * i + k] * dFFollowdFs[num_configs * sot_dimension * k + j + sot_dimension];
                }
            }
        }

        dGradChi1dChi = fifthOrderTensor(tot_dimension * sot_dimension, 0);

        dGradChi1dChin = floatVector(tot_dimension * (num_configs - 1) * sot_dimension, 0);

        dGradChi1dGradChi = sixthOrderTensor(tot_dimension * tot_dimension, 0);

        dGradChi1dCn = floatVector(tot_dimension * (num_configs - 1) * sot_dimension, 0);

        dGradChi1dGradChin = floatVector(tot_dimension * (num_configs - 1) * tot_dimension, 0);

        std::fill(temp_tot1.begin(), temp_tot1.end(), 0.);
        std::fill(temp_tot2.begin(), temp_tot2.end(), 0.);

        fifthOrderTensor temp_fiot(tot_dimension * sot_dimension, 0);

        floatVector temp_siot1(tot_dimension * (num_configs - 1) * sot_dimension, 0);
        floatVector temp_siot2(tot_dimension * (num_configs - 1) * sot_dimension, 0);

        floatVector temp_seot(tot_dimension * (num_configs - 1) * tot_dimension, 0);

        for (unsigned int i = 0; i < configuration::dimension; i++) {
            for (unsigned int I = 0; I < configuration::dimension; I++) {
                for (unsigned int J = 0; J < configuration::dimension; J++) {
                    for (unsigned int a = 0; a < configuration::dimension; a++) {
                        temp_tot1[configuration::dimension * configuration::dimension * i + configuration::dimension * I + J] +=
                            gradientChi1Reference[configuration::dimension * configuration::dimension * i + configuration::dimension * a + J] * invChiFollow_T[configuration::dimension * I + a];

                        temp_tot2[configuration::dimension * configuration::dimension * i + configuration::dimension * I + J] +=
                            gradientChi1Reference[configuration::dimension * configuration::dimension * i + configuration::dimension * I + a] * invFFollow_T[configuration::dimension * J + a];

                        for (unsigned int b = 0; b < configuration::dimension; b++) {
                            dGradChi1dGradChi[configuration::dimension * configuration::dimension * tot_dimension * i + configuration::dimension * tot_dimension * I + tot_dimension * J +
                                              configuration::dimension * configuration::dimension * i + configuration::dimension * a + b] +=
                                invChiFollow_T[configuration::dimension * I + a] * invFFollow_T[configuration::dimension * J + b];

                            for (unsigned int k = 0; k < configuration::dimension; k++) {
                                temp_fiot[configuration::dimension * configuration::dimension * sot_dimension * i + configuration::dimension * sot_dimension * I + sot_dimension * J + configuration::dimension * a + b] +=
                                    dGradientChi1ReferencedChi[configuration::dimension * configuration::dimension * sot_dimension * i + configuration::dimension * sot_dimension * k +
                                                               sot_dimension * J + configuration::dimension * a + b] *
                                    invChiFollow_T[configuration::dimension * I + k];
                            }
                        }
                    }
                }
            }
        }
        for (unsigned int i = 0; i < configuration::dimension; i++) {
            for (unsigned int I = 0; I < configuration::dimension; I++) {
                for (unsigned int J = 0; J < configuration::dimension; J++) {
                    for (unsigned int k = 0; k < configuration::dimension; k++) {
                        for (unsigned int a = 0; a < configuration::dimension; a++) {
                            for (unsigned int b = 0; b < configuration::dimension; b++) {
                                dGradChi1dChi[configuration::dimension * configuration::dimension * sot_dimension * i + configuration::dimension * sot_dimension * I + sot_dimension * J + configuration::dimension * a +
                                              b] +=
                                    temp_fiot[configuration::dimension * configuration::dimension * sot_dimension * i + configuration::dimension * sot_dimension * I + sot_dimension * k + configuration::dimension * a + b] *
                                    invFFollow_T[configuration::dimension * J + k];
                            }
                        }
                    }
                }
            }
        }

        for (unsigned int i = 0; i < configuration::dimension; i++) {
            for (unsigned int I = 0; I < configuration::dimension; I++) {
                for (unsigned int J = 0; J < configuration::dimension; J++) {
                    for (unsigned int index = 1; index < num_configs; index++) {
                        for (unsigned int a = 0; a < configuration::dimension; a++) {
                            for (unsigned int b = 0; b < configuration::dimension; b++) {
                                for (unsigned int k = 0; k < configuration::dimension; k++) {
                                    temp_siot1[configuration::dimension * configuration::dimension * (num_configs - 1) * sot_dimension * i +
                                               configuration::dimension * (num_configs - 1) * sot_dimension * I + (num_configs - 1) * sot_dimension * J +
                                               configuration::dimension * configuration::dimension * (index - 1) + configuration::dimension * a + b] +=
                                        dGradientChi1ReferencedChin[configuration::dimension * configuration::dimension * (num_configs - 1) * sot_dimension * i +
                                                                    configuration::dimension * (num_configs - 1) * sot_dimension * k +
                                                                    (num_configs - 1) * sot_dimension * J +
                                                                    configuration::dimension * configuration::dimension * (index - 1) + configuration::dimension * a + b] *
                                        invChiFollow_T[configuration::dimension * I + k];

                                    temp_siot2[configuration::dimension * configuration::dimension * (num_configs - 1) * sot_dimension * i +
                                               configuration::dimension * (num_configs - 1) * sot_dimension * I + (num_configs - 1) * sot_dimension * J +
                                               configuration::dimension * configuration::dimension * (index - 1) + configuration::dimension * a + b] +=
                                        dGradientChi1ReferencedCn[configuration::dimension * configuration::dimension * (num_configs - 1) * sot_dimension * i +
                                                                  configuration::dimension * (num_configs - 1) * sot_dimension * k +
                                                                  (num_configs - 1) * sot_dimension * J +
                                                                  configuration::dimension * configuration::dimension * (index - 1) + configuration::dimension * a + b] *
                                        invChiFollow_T[configuration::dimension * I + k];
                                }
                            }
                        }
                    }
                }
            }
        }

        for (unsigned int i = 0; i < configuration::dimension; i++) {
            for (unsigned int I = 0; I < configuration::dimension; I++) {
                for (unsigned int c = 0; c < configuration::dimension; c++) {
                    for (unsigned int J = 0; J < configuration::dimension; J++) {
                        for (unsigned int index = 1; index < num_configs; index++) {
                            for (unsigned int a = 0; a < configuration::dimension; a++) {
                                for (unsigned int b = 0; b < configuration::dimension; b++) {
                                    dGradChi1dCn[configuration::dimension * configuration::dimension * (num_configs - 1) * sot_dimension * i +
                                                 configuration::dimension * (num_configs - 1) * sot_dimension * I +
                                                 (num_configs - 1) * sot_dimension * J + configuration::dimension * configuration::dimension * (index - 1) + configuration::dimension * a +
                                                 b] += temp_tot1[configuration::dimension * configuration::dimension * i + configuration::dimension * I + c] *
                                                       dInvFFollowdFn[configuration::dimension * (num_configs - 1) * sot_dimension * c +
                                                                      (num_configs - 1) * sot_dimension * J +
                                                                      configuration::dimension * configuration::dimension * (index - 1) + configuration::dimension * a + b];

                                    dGradChi1dChin[configuration::dimension * configuration::dimension * (num_configs - 1) * sot_dimension * i +
                                                   configuration::dimension * (num_configs - 1) * sot_dimension * I +
                                                   (num_configs - 1) * sot_dimension * J + configuration::dimension * configuration::dimension * (index - 1) + configuration::dimension * a +
                                                   b] += temp_tot2[configuration::dimension * configuration::dimension * i + configuration::dimension * c + J] *
                                                         dInvChiFollowdChin[configuration::dimension * (num_configs - 1) * sot_dimension * c +
                                                                            (num_configs - 1) * sot_dimension * I +
                                                                            configuration::dimension * configuration::dimension * (index - 1) + configuration::dimension * a + b];
                                    dGradChi1dChin[configuration::dimension * configuration::dimension * (num_configs - 1) * sot_dimension * i +
                                                   configuration::dimension * (num_configs - 1) * sot_dimension * I +
                                                   (num_configs - 1) * sot_dimension * J + configuration::dimension * configuration::dimension * (index - 1) + configuration::dimension * a +
                                                   b] += temp_siot1[configuration::dimension * configuration::dimension * (num_configs - 1) * sot_dimension * i +
                                                                    configuration::dimension * (num_configs - 1) * sot_dimension * I +
                                                                    (num_configs - 1) * sot_dimension * c +
                                                                    configuration::dimension * configuration::dimension * (index - 1) + configuration::dimension * a + b] *
                                                         invFFollow_T[configuration::dimension * J + c];

                                    dGradChi1dCn[configuration::dimension * configuration::dimension * (num_configs - 1) * sot_dimension * i +
                                                 configuration::dimension * (num_configs - 1) * sot_dimension * I +
                                                 (num_configs - 1) * sot_dimension * J + configuration::dimension * configuration::dimension * (index - 1) + configuration::dimension * a +
                                                 b] += temp_siot2[configuration::dimension * configuration::dimension * (num_configs - 1) * sot_dimension * i +
                                                                  configuration::dimension * (num_configs - 1) * sot_dimension * I +
                                                                  (num_configs - 1) * sot_dimension * c +
                                                                  configuration::dimension * configuration::dimension * (index - 1) + configuration::dimension * a + b] *
                                                       invFFollow[configuration::dimension * c + J];

                                    for (unsigned int k = 0; k < configuration::dimension; k++) {
                                        temp_seot[configuration::dimension * configuration::dimension * (num_configs - 1) * tot_dimension * i +
                                                  configuration::dimension * (num_configs - 1) * tot_dimension * I +
                                                  (num_configs - 1) * tot_dimension * J + tot_dimension * (index - 1) +
                                                  configuration::dimension * configuration::dimension * a + configuration::dimension * b + k] +=
                                            dGradientChi1ReferencedGradChin[configuration::dimension * configuration::dimension * (num_configs - 1) * tot_dimension *
                                                                                i +
                                                                            configuration::dimension * (num_configs - 1) * tot_dimension * c +
                                                                            (num_configs - 1) * tot_dimension * J +
                                                                            tot_dimension * (index - 1) + configuration::dimension * configuration::dimension * a +
                                                                            configuration::dimension * b + k] *
                                            invChiFollow_T[configuration::dimension * I + c];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (unsigned int i = 0; i < configuration::dimension; i++) {
            for (unsigned int I = 0; I < configuration::dimension; I++) {
                for (unsigned int J = 0; J < configuration::dimension; J++) {
                    for (unsigned int c = 0; c < configuration::dimension; c++) {
                        for (unsigned int index = 1; index < num_configs; index++) {
                            for (unsigned int a = 0; a < configuration::dimension; a++) {
                                for (unsigned int b = 0; b < configuration::dimension; b++) {
                                    for (unsigned int k = 0; k < configuration::dimension; k++) {
                                        dGradChi1dGradChin[configuration::dimension * configuration::dimension * (num_configs - 1) * tot_dimension * i +
                                                           configuration::dimension * (num_configs - 1) * tot_dimension * I +
                                                           (num_configs - 1) * tot_dimension * J + tot_dimension * (index - 1) +
                                                           configuration::dimension * configuration::dimension * a + configuration::dimension * b + k] +=
                                            temp_seot[configuration::dimension * configuration::dimension * (num_configs - 1) * tot_dimension * i +
                                                      configuration::dimension * (num_configs - 1) * tot_dimension * I +
                                                      (num_configs - 1) * tot_dimension * c + tot_dimension * (index - 1) +
                                                      configuration::dimension * configuration::dimension * a + configuration::dimension * b + k] *
                                            invFFollow_T[configuration::dimension * J + c];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

}  // namespace tardigradeHydra
