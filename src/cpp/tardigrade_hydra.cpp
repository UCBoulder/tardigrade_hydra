/**
 ******************************************************************************
 * \file tardigrade_hydra.cpp
 ******************************************************************************
 * A C++ library for defining frameworks to solve finite deformation material
 * models.
 ******************************************************************************
 */

#include "tardigrade_hydra.h"

#include "tardigrade_DeformationDecompositionBase.h"  //TEMP

namespace tardigradeHydra {

    /*!
     * The main constructor for the hydra base class. Inputs are all the required values for most solves.
     *
     * \param &DOFStorage: The degrees of freedom storage object
     * \param &ModelConfiguration: The model configuration object
      \param &_hydra_configuration: Class which defines the hydra configuration
     */
    hydraBase::hydraBase(const DOFStorageBase &DOFStorage, const ModelConfigurationBase &ModelConfiguration,
                         HydraConfigurationBase _hydra_configuration)
        : hydra_configuration(_hydra_configuration),
          dof(&DOFStorage),
          model_configuration(&ModelConfiguration),
          _stress_size(_hydra_configuration.configuration_unknown_count) {
        // TEMP
        _solver.hydra                  = this;
        _solver.internal_solver->hydra = this;
        auto local_internal_solver     = dynamic_cast<tardigradeHydra::RelaxedSolverBase *>(_solver.internal_solver);
        local_internal_solver->internal_solver->hydra = this;
        // END TEMP
    }

    /*!
     * Set the value of the stress
     *
     * \param &stress: The stress in row-major form
     */
    void hydraBase::setStress(const floatVector &stress) { setIterationData(stress, _stress); }

    /*!
     * Get a SetDataStorage object for the stress
     */
    hydraBase::SetDataStorageIteration<secondOrderTensor> hydraBase::get_SetDataStorage_stress() {
        return hydraBase::SetDataStorageIteration<secondOrderTensor>(&_stress, this);
    }

    //! Get the number of terms in the unknown vector
    const unsigned int hydraBase::getNumUnknowns() {
        return getNumConfigurations() * getConfigurationUnknownCount() + getNumNonLinearSolveStateVariables();
    }

    //! Get the number of additional degrees of freedom
    const unsigned int hydraBase::getNumAdditionalDOF() { return getAdditionalDOF()->size(); }

    //! Get the current residual index
    const unsigned int hydraBase::getCurrentResidualIndex() {
        TARDIGRADE_ERROR_TOOLS_CHECK(currentResidualIndexMeaningful(), "The current residual index isn't meaningful");
        return _current_residual_index;
    }

    /*!
     * Set the verbosity level for failures
     *
     * \param &value: The verbosity level of the failure (defaults to zero)
     */
    void hydraBase::setFailureVerbosityLevel(const unsigned int &value) { _failure_verbosity_level = value; }

    /*! Get a reference to the full residual that is mutable. Returns NULL if it's not allowed.
     *
     * This should only be called in residual classes that need to modify the full residual in their
     * modifyGlobalResidual methods.
     *
     * Be careful!
     */
    floatVector *hydraBase::getMutableResidual() {
        if (_allow_modify_global_residual) {
            return &_residual.second;
        }

        return NULL;
    }

    /*! Get a reference to the full jacobian that is mutable. Returns NULL if it's not allowed.
     *
     * This should only be called in residual classes that need to modify the full residual in their
     * modifyGlobalJacobian methods.
     *
     * Be careful!
     */
    floatVector *hydraBase::getMutableJacobian() {
        if (_allow_modify_global_jacobian) {
            return &_jacobian.second;
        }

        return NULL;
    }

    /*! Get a reference to the full dRdT that is mutable. Returns NULL if it's not allowed.
     *
     * This should only be called in residual classes that need to modify the full residual in their
     * modifyGlobaldRdT methods.
     *
     * Be careful!
     */
    floatVector *hydraBase::getMutabledRdT() {
        if (_allow_modify_global_dRdT) {
            return &_dRdT.second;
        }

        return NULL;
    }

    /*! Get a reference to the full dRdF that is mutable. Returns NULL if it's not allowed.
     *
     * This should only be called in residual classes that need to modify the full residual in their
     * modifyGlobaldRdF methods.
     *
     * Be careful!
     */
    floatVector *hydraBase::getMutabledRdF() {
        if (_allow_modify_global_dRdF) {
            return &_dRdF.second;
        }

        return NULL;
    }

    /*! Get a reference to the full dRdAdditionalDOF that is mutable. Returns NULL if it's not allowed.
     *
     * This should only be called in residual classes that need to modify the full residual in their
     * modifyGlobaldRdAdditionalDOF methods.
     *
     * Be careful!
     */
    floatVector *hydraBase::getMutabledRdAdditionalDOF() {
        if (_allow_modify_global_dRdAdditionalDOF) {
            return &_dRdAdditionalDOF.second;
        }

        return NULL;
    }

    /*!
     * Return if the current residual index is meaningful or not
     */
    const bool hydraBase::currentResidualIndexMeaningful() { return _current_residual_index_set; }

    /*!
     * Loosen the convergence tolerance for the next iteration
     * Useful if a Residual's form is changing in a non-smooth way
     *
     * \param factor: The scale factor to be applied to the current tolerance
     */
    void hydraBase::setToleranceScaleFactor(floatType factor) {
        if (factor > _residual_scale_factor) {
            _residual_scale_factor = factor;
        }
    }

    /*!
     * Update the additional state variable vector
     */
    void hydraBase::updateAdditionalStateVariables() {
        for (auto v = std::begin(*getResidualClasses()); v != std::end(*getResidualClasses()); ++v) {
            (*v)->updateAdditionalStateVariables(_additionalStateVariables.second);
        }
    }

    /*!
     * Extract the stresses out of the unknown vector
     */
    void hydraBase::extractStress() {
        const floatVector *unknownVector = getUnknownVector();

        auto stress = get_SetDataStorage_stress();

        stress.zero(getConfigurationUnknownCount());

        std::copy(std::begin(*unknownVector), std::begin(*unknownVector) + getConfigurationUnknownCount(),
                  std::begin(*stress.value));
    }

    /*!
     * Compute the configurations from the provided vector. Each configuration is assumed to have a dimension
     * of dimension x dimension
     *
     * \param *data_vector: A pointer to the vector of data which contains the configurations and other information
     * \param &start_index: The starting index for the vector
     * \param &total_transformation: The total transformation from the reference to the current configuration
     * \param &configurations: The resulting collection of configurations
     * \param &inverseConfigurations: The resulting inverse configurations
     * \param add_eye: A flag for whether to add the identity matrix to each of the configurations except for the
     *     total transformation. Defaults to false.
     */
    void hydraBase::computeConfigurations(const floatVector *data_vector, const unsigned int start_index,
                                          const floatVector &total_transformation, floatVector &configurations,
                                          floatVector &inverseConfigurations, const bool add_eye) {
        constexpr unsigned int sot_dimension = configuration::dimension * configuration::dimension;

        auto num_configs = getNumConfigurations();

        // Set the configurations
        configurations = floatVector(num_configs * sot_dimension, 0);

        inverseConfigurations = floatVector(num_configs * sot_dimension, 0);

        auto mat = tardigradeHydra::getFixedSizeMatrixMap<floatType, 3, 3>(inverseConfigurations.data());
#ifdef TARDIGRADE_HYDRA_USE_LLXSMM
        kernel_type kernel(LIBXSMM_GEMM_FLAG_NONE, configuration::dimension, configuration::dimension, configuration::dimension, 1, 0);

        // Initialize the first configuration with the total deformation gradient
        secondOrderTensor temp(sot_dimension, 0);
#else
        auto mat2 = tardigradeHydra::getFixedSizeMatrixMap<floatType, 3, 3>(configurations.data());
#endif

        std::copy(total_transformation.begin(), total_transformation.end(), configurations.begin());

        for (int i = num_configs - 2; i >= 0; i--) {
            // Set the current configuration as being equal to the previous
            std::copy(data_vector->begin() + i * sot_dimension + start_index,
                      data_vector->begin() + (i + 1) * sot_dimension + start_index,
                      configurations.begin() + sot_dimension * (i + 1));

            if (add_eye) {
                for (unsigned int j = 0; j < configuration::dimension; j++) {
                    configurations[sot_dimension * (i + 1) + configuration::dimension * j + j] += 1;
                }
            }

            // Compute the inverse of the current configuration and store it
            std::copy(configurations.begin() + sot_dimension * (i + 1), configurations.begin() + sot_dimension * (i + 2),
                      inverseConfigurations.begin() + sot_dimension * (i + 1));
            new (&mat) Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> >(inverseConfigurations.data() +
                                                                                        sot_dimension * (i + 1),
                                                                                    3, 3);
            mat = mat.inverse().eval();

#ifdef TARDIGRADE_HYDRA_USE_LLXSMM
            std::copy(configurations.begin(), configurations.begin() + sot_dimension, temp.begin());

            kernel(&inverseConfigurations[sot_dimension * (i + 1)], &temp[0], &configurations[0]);
#else
            // Add contribution of deformation gradient to the first configuration

            new (&mat2) Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> >(configurations.data(), 3, 3);

            mat2 *= mat;
#endif
        }

        std::copy(configurations.begin(), configurations.begin() + sot_dimension, inverseConfigurations.begin());

        new (&mat) Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor> >(inverseConfigurations.data(), 3, 3);
        mat = mat.inverse().eval();

        return;
    }

    /*!
     * Update the configurations from the unknown vector
     */
    void hydraBase::updateConfigurationsFromUnknownVector() {
        const floatVector *unknownVector = getUnknownVector();

        // Set the configurations
        auto configurations = deformation->get_SetDataStorage_configurations();

        auto inverseConfigurations = deformation->get_SetDataStorage_inverseConfigurations();

        computeConfigurations(unknownVector, getStressSize(), *getDeformationGradient(), *configurations.value,
                              *inverseConfigurations.value);

        // Extract the remaining state variables required for the non-linear solve
        auto nonLinearSolveStateVariables = get_SetDataStorage_nonLinearSolveStateVariables();

        auto nNLISV = getNumNonLinearSolveStateVariables();

        nonLinearSolveStateVariables.zero(nNLISV);

        std::copy(std::begin(*unknownVector) + getNumConfigurations() * getConfigurationUnknownCount(),
                  std::end(*unknownVector), std::begin(*nonLinearSolveStateVariables.value));
    }

    /*!
     * Decompose the unknown vector into the cauchy stress, configurations, and state variables used for the
     * non-linear solve
     */
    void hydraBase::decomposeUnknownVector() {
        // Set the stress
        extractStress();

        updateConfigurationsFromUnknownVector();
    }

    /*!
     * Decompose the incoming state variable vector setting the different configurations along the way
     *
     * The state variable vector is assumed to be of the form:
     *
     * \f$ \text{ISV} = \left\{\bf{F}^2 - \bf{I}, \bf{F}^3 - \bf{I}, \cdots, \bf{F}^n - \bf{I}, \xi^1, \xi^2,
     * \cdots, \xi^m, \eta^1, \cdots\right\} \f$
     *
     * where the \f$\bf{F}^x\f$ are the different configurations, \f$\xi^y\f$ are the other variables to be solved
     * during the non-linear solve and \f$\eta^z\f$ are other state variables. Note that we decompose the
     * deformation gradient as
     *
     * \f$\bf{F} = \bf{F}^1 \bf{F}^2 \cdots \bf{F}^n\f$
     *
     * and so because \f$\bf{F}\f$ is provided we can solve for \f$\bf{F}^1\f$. Typically, this configuration would
     * be the elastic configuration (i.e., the configuration that generates the stress) though we do not insist that
     * users follow convention.
     */
    void hydraBase::decomposeStateVariableVector() {
        auto nConfig = getNumConfigurations();

        auto nNLISV = getNumNonLinearSolveStateVariables();

        // Extract the previous configurations
        if (getPreviousStateVariables()->size() < ((nConfig - 1) * getConfigurationUnknownCount() + nNLISV)) {
            std::string message = "The number of state variables is less than required for the configurations and ";
            message += "non-linear state variables\n";
            message += "  # previousStateVariables                               : " +
                       std::to_string(getPreviousStateVariables()->size()) + "\n";
            message += "  # ( configurations - 1 ) * configuration_unknown_count : " +
                       std::to_string((nConfig - 1) * getConfigurationUnknownCount()) + "\n";
            message += "  # non-linear solve ISVs                                : " + std::to_string(nNLISV) + "\n";
            message += "  # minimum required ISVs                                : " +
                       std::to_string((nConfig - 1) * getConfigurationUnknownCount() + nNLISV);

            TARDIGRADE_ERROR_TOOLS_CATCH(throw std::runtime_error(message));
        }

        auto configurations = deformation->get_SetDataStorage_configurations();

        auto previousConfigurations = deformation->get_SetDataStorage_previousConfigurations();

        auto inverseConfigurations = deformation->get_SetDataStorage_inverseConfigurations();

        auto previousInverseConfigurations = deformation->get_SetDataStorage_previousInverseConfigurations();

        // Compute the configurations
        computeConfigurations(getPreviousStateVariables(), 0, *getDeformationGradient(), *configurations.value,
                              *inverseConfigurations.value, true);

        computeConfigurations(getPreviousStateVariables(), 0, *getPreviousDeformationGradient(),
                              *previousConfigurations.value, *previousInverseConfigurations.value, true);

        // Extract the remaining state variables required for the non-linear solve
        auto nonLinearSolveStateVariables         = get_SetDataStorage_nonLinearSolveStateVariables();
        auto previousNonLinearSolveStateVariables = get_SetDataStorage_previousNonLinearSolveStateVariables();

        previousNonLinearSolveStateVariables.zero(nNLISV);

        std::copy(std::begin(*getPreviousStateVariables()) + (nConfig - 1) * getConfigurationUnknownCount(),
                  std::begin(*getPreviousStateVariables()) + (nConfig - 1) * getConfigurationUnknownCount() + nNLISV,
                  std::begin(*previousNonLinearSolveStateVariables.value));

        *nonLinearSolveStateVariables.value = *previousNonLinearSolveStateVariables.value;

        // Extract the additional state variables
        auto additionalStateVariables         = get_SetDataStorage_additionalStateVariables();
        auto previousAdditionalStateVariables = get_SetDataStorage_previousAdditionalStateVariables();

        unsigned int nAISV = (unsigned int)(std::end(*getPreviousStateVariables()) -
                                            (std::begin(*getPreviousStateVariables()) +
                                             (nConfig - 1) * getConfigurationUnknownCount() + nNLISV));

        previousAdditionalStateVariables.zero(nAISV);

        std::copy(std::begin(*getPreviousStateVariables()) + (nConfig - 1) * getConfigurationUnknownCount() + nNLISV,
                  std::end(*getPreviousStateVariables()), std::begin(*previousAdditionalStateVariables.value));

        *additionalStateVariables.value = *previousAdditionalStateVariables.value;
    }

    /*!
     * Build an error message for when the upper index is larger than the number of configurations
     *
     * \param upperIndex: The upper index
     * \param num_configurations: The number of configurations
     */
    std::string hydraBase::build_upper_index_out_of_range_error_string(const unsigned int upperIndex,
                                                                       const unsigned int num_configurations) {
        std::string message = "The upper index must be less than or equal to the total number of configurations\n";
        message += "  upperIndex      : " + std::to_string(upperIndex) + "\n";
        message += "  # configurations: " + std::to_string(num_configurations);

        return message;
    }

    /*!
     * Build an error message for when the lower index is larger than the upper index
     *
     * \param lowerIndex: The lower configuration index
     * \param upperIndex: The upper configuration index
     */
    std::string hydraBase::build_lower_index_out_of_range_error_string(const unsigned int lowerIndex,
                                                                       const unsigned int upperIndex) {
        std::string message = "The upper index must be greater than or equal to the lower index\n";
        message += "  lowerIndex: " + std::to_string(lowerIndex) + "\n";
        message += "  upperIndex: " + std::to_string(upperIndex) + "\n";

        return message;
    }

    /*!
     * Reset the iteration data to the new base state
     */
    void hydraBase::resetIterationData() {
        for (auto d = _iterationData.begin(); d != _iterationData.end(); d++) {
            (*d)->clear();
        }

        _iterationData.clear();
    }

    /*!
     * Reset the nonlinear step data to the new base state
     */
    void hydraBase::resetNLStepData() {
        for (auto d = _nlStepData.begin(); d != _nlStepData.end(); d++) {
            (*d)->clear();
        }

        _nlStepData.clear();
    }

    /*!
     * Set the vectors for the residuals.
     *
     * The expected form of the residual vector is cauchy stress, configurations,
     * state variables though only the residual on the cauchy stress must come
     * in this specific order if the deconstructSolutionVector function is redefined.
     *
     * The user should define a vector of ResidualBase objects and use the
     * setResidualClasses( std::vector< ResidualBase > & ) function here.
     *
     * The resulting residual should have the form
     *
     * residual = { cauchyResidual, F2residual, ... Fnresidual, xiresidual1, xiresidual2, ... }
     *
     * and can be formed by any number of residual classes. The first residual class must also
     * have the method `void getStress( )` defined which will return the current value
     * of the stress.
     */
    void hydraBase::setResidualClasses() {}

    /*!
     * Set the residual classes
     *
     * \param &residualClasses: A vector of residual classes which will be used to
     *     populate the residual and jacobian matrices for the non-linear solve
     */
    void hydraBase::setResidualClasses(std::vector<ResidualBase<hydraBase> *> &residualClasses) {
        unsigned int numEquations = 0;

        _residualClasses.second = std::vector<ResidualBase<hydraBase> *>(residualClasses.size());

        for (auto c = residualClasses.begin(); c != residualClasses.end(); c++) {
            numEquations += (*c)->getNumEquations();

            _residualClasses.second[c - residualClasses.begin()] = *c;
        }

        if (numEquations !=
            (getNumConfigurations() * getConfigurationUnknownCount() + getNumNonLinearSolveStateVariables())) {
            std::string message =
                "The number of equations for the non-linear solve is not equal to the number of equations defined\n";
            message += "  expected number of equations: " +
                       std::to_string(getNumConfigurations() * getConfigurationUnknownCount() +
                                      getNumNonLinearSolveStateVariables()) +
                       "\n";
            message += "  number of defined equations:  " + std::to_string(numEquations) + "\n";

            TARDIGRADE_ERROR_TOOLS_CATCH(throw std::runtime_error(message));
        }

        _residualClasses.first = true;
    }

    /*!
     * Get a pointer to the vector of residual class pointers
     */
    std::vector<ResidualBase<hydraBase> *> *hydraBase::getResidualClasses() {
        if (!_residualClasses.first) {
            TARDIGRADE_ERROR_TOOLS_CATCH(setResidualClasses());
        }

        return &_residualClasses.second;
    }

    /*!
     * Form the residual
     */
    void hydraBase::formNonLinearResidual() {
        auto configurationUnknownCount = getConfigurationUnknownCount();

        auto residualSize = getNumConfigurations() * configurationUnknownCount + getNumNonLinearSolveStateVariables();

        auto numAdditionalDOF = getAdditionalDOF()->size();

        _residual.second = floatVector(residualSize, 0);

        _jacobian.second = floatVector(residualSize * residualSize, 0);

        _dRdF.second = floatVector(residualSize * configurationUnknownCount, 0);

        _dRdT.second = floatVector(residualSize, 0);

        _dRdAdditionalDOF.second = floatVector(residualSize * numAdditionalDOF, 0);

        _additionalDerivatives.second.clear();

        unsigned int offset = 0;

        setCurrentResidualIndexMeaningful(true);
        for (auto residual_ptr = getResidualClasses()->begin(); residual_ptr != getResidualClasses()->end();
             ++residual_ptr) {
            setCurrentResidualIndex(residual_ptr - getResidualClasses()->begin());

            // Extract the terms

            const floatVector *localResidual;
            TARDIGRADE_ERROR_TOOLS_CATCH(localResidual = (*residual_ptr)->getResidual());

            // Check the contributions to make sure they are consistent sizes

            TARDIGRADE_ERROR_TOOLS_CHECK(localResidual->size() == (*residual_ptr)->getNumEquations(),
                                         "The residual for residual " +
                                             std::to_string(residual_ptr - getResidualClasses()->begin()) +
                                             " is not the expected length\n" +
                                             "  expected: " + std::to_string((*residual_ptr)->getNumEquations()) +
                                             "\n" + "  actual:   " + std::to_string(localResidual->size()) + "\n")

            // Store the values in the global quantities

            // Copy over the values of the local vector to the global structures
            std::copy(localResidual->begin(), localResidual->end(), _residual.second.begin() + offset);

            offset += (*residual_ptr)->getNumEquations();
        }
        setCurrentResidualIndexMeaningful(false);

        // Allow the residuals to modify the global residual if needed
        setAllowModifyGlobalResidual(true);
        setCurrentResidualIndexMeaningful(true);
        for (auto residual_ptr = getResidualClasses()->begin(); residual_ptr != getResidualClasses()->end();
             ++residual_ptr) {
            setCurrentResidualIndex(residual_ptr - getResidualClasses()->begin());

            (*residual_ptr)->modifyGlobalResidual();
        }
        setCurrentResidualIndexMeaningful(false);
        setAllowModifyGlobalResidual(false);

        _residual.first = true;

        addIterationData(&_residual);
    }

    /*!
     * Form the jacobian and gradient matrices
     */
    void hydraBase::formNonLinearDerivatives() {
        auto configurationUnknownCount = getConfigurationUnknownCount();

        auto residualSize = getNumConfigurations() * configurationUnknownCount + getNumNonLinearSolveStateVariables();

        auto numAdditionalDOF = getAdditionalDOF()->size();

        _jacobian.second = floatVector(residualSize * residualSize, 0);

        _dRdF.second = floatVector(residualSize * configurationUnknownCount, 0);

        _dRdT.second = floatVector(residualSize, 0);

        _dRdAdditionalDOF.second = floatVector(residualSize * numAdditionalDOF, 0);

        _additionalDerivatives.second.clear();

        unsigned int offset = 0;

        unsigned int numAdditionalDerivatives = 0;

        setCurrentResidualIndexMeaningful(true);
        for (auto residual_ptr = getResidualClasses()->begin(); residual_ptr != getResidualClasses()->end();
             ++residual_ptr) {
            setCurrentResidualIndex(residual_ptr - getResidualClasses()->begin());

            // Extract the terms

            const floatVector *localJacobian;
            TARDIGRADE_ERROR_TOOLS_CATCH(localJacobian = (*residual_ptr)->getJacobian());

            const floatVector *localdRdF;
            TARDIGRADE_ERROR_TOOLS_CATCH(localdRdF = (*residual_ptr)->getdRdF());

            const floatVector *localdRdT;
            TARDIGRADE_ERROR_TOOLS_CATCH(localdRdT = (*residual_ptr)->getdRdT());

            const floatVector *localdRdAdditionalDOF;
            TARDIGRADE_ERROR_TOOLS_CATCH(localdRdAdditionalDOF = (*residual_ptr)->getdRdAdditionalDOF());

            const floatVector *localAdditionalDerivatives;
            TARDIGRADE_ERROR_TOOLS_CATCH(localAdditionalDerivatives = (*residual_ptr)->getAdditionalDerivatives());

            // Check the contributions to make sure they are consistent sizes

            TARDIGRADE_ERROR_TOOLS_CHECK(localJacobian->size() == (*residual_ptr)->getNumEquations() * residualSize,
                                         "The jacobian for residual " +
                                             std::to_string(residual_ptr - getResidualClasses()->begin()) +
                                             " is not the expected length\n" + "  expected: " +
                                             std::to_string((*residual_ptr)->getNumEquations() * residualSize) + "\n" +
                                             "  actual:   " + std::to_string(localJacobian->size()) + "\n")

            TARDIGRADE_ERROR_TOOLS_CHECK(
                localdRdF->size() == (*residual_ptr)->getNumEquations() * configurationUnknownCount,
                "dRdF for residual " + std::to_string(residual_ptr - getResidualClasses()->begin()) +
                    " is not the expected length\n" +
                    "  expected: " + std::to_string((*residual_ptr)->getNumEquations() * configurationUnknownCount) +
                    "\n" + "  actual:   " + std::to_string(localdRdF->size()) + "\n")

            TARDIGRADE_ERROR_TOOLS_CHECK(localdRdT->size() == (*residual_ptr)->getNumEquations(),
                                         "dRdT for residual " +
                                             std::to_string(residual_ptr - getResidualClasses()->begin()) +
                                             " is not the expected length\n" +
                                             "  expected: " + std::to_string((*residual_ptr)->getNumEquations()) +
                                             "\n" + "  actual:   " + std::to_string(localdRdT->size()) + "\n")

            if (localdRdAdditionalDOF->size() != 0) {
                TARDIGRADE_ERROR_TOOLS_CHECK(
                    localdRdAdditionalDOF->size() == ((*residual_ptr)->getNumEquations() * numAdditionalDOF),
                    "dRdAdditionalDOF for residual " + std::to_string(residual_ptr - getResidualClasses()->begin()) +
                        " is not the expected length\n" +
                        "  expected: " + std::to_string((*residual_ptr)->getNumEquations() * numAdditionalDOF) + "\n" +
                        "  actual  : " + std::to_string(localdRdAdditionalDOF->size()) + "\n")

                std::copy(localdRdAdditionalDOF->begin(), localdRdAdditionalDOF->end(),
                          _dRdAdditionalDOF.second.begin() + numAdditionalDOF * offset);
            }

            if (localAdditionalDerivatives->size() != 0) {
                if ((*localAdditionalDerivatives).size() !=
                    (*residual_ptr)->getNumEquations() * numAdditionalDerivatives) {
                    if (numAdditionalDerivatives == 0) {
                        numAdditionalDerivatives =
                            (*localAdditionalDerivatives).size() / (*residual_ptr)->getNumEquations();

                        _additionalDerivatives.second = floatVector(residualSize * numAdditionalDerivatives, 0);

                    } else {
                        std::string message = "The additional derivatives for residual " +
                                              std::to_string(residual_ptr - getResidualClasses()->begin()) +
                                              " are not the expected length as determined from the first residual\n";
                        message += "  expected: " + std::to_string(numAdditionalDerivatives) + "\n";
                        message += "  actual:   " + std::to_string((*localAdditionalDerivatives).size()) + "\n";

                        TARDIGRADE_ERROR_TOOLS_CATCH(throw std::runtime_error(message));
                    }
                }
            }

            // Store the values in the global quantities

            // Copy over the values of the local vector to the global structures
            std::copy(localJacobian->begin(), localJacobian->end(), _jacobian.second.begin() + residualSize * offset);

            std::copy(localdRdF->begin(), localdRdF->end(), _dRdF.second.begin() + configurationUnknownCount * offset);

            std::copy(localdRdT->begin(), localdRdT->end(), _dRdT.second.begin() + offset);

            std::copy(localAdditionalDerivatives->begin(), localAdditionalDerivatives->end(),
                      _additionalDerivatives.second.begin() + numAdditionalDerivatives * offset);

            offset += (*residual_ptr)->getNumEquations();
        }
        setCurrentResidualIndexMeaningful(false);

        setAllowModifyGlobalJacobian(true);
        setAllowModifyGlobaldRdT(true);
        setAllowModifyGlobaldRdF(true);
        setAllowModifyGlobaldRdAdditionalDOF(true);
        setCurrentResidualIndexMeaningful(true);
        for (auto residual_ptr = getResidualClasses()->begin(); residual_ptr != getResidualClasses()->end();
             residual_ptr++) {
            setCurrentResidualIndexMeaningful(residual_ptr - getResidualClasses()->begin());

            (*residual_ptr)->modifyGlobalJacobian();
            (*residual_ptr)->modifyGlobaldRdT();
            (*residual_ptr)->modifyGlobaldRdF();
            (*residual_ptr)->modifyGlobaldRdAdditionalDOF();
        }
        setCurrentResidualIndexMeaningful(false);
        setAllowModifyGlobalJacobian(false);
        setAllowModifyGlobaldRdT(false);
        setAllowModifyGlobaldRdF(false);
        setAllowModifyGlobaldRdAdditionalDOF(false);

        _jacobian.first = true;

        _dRdF.first = true;

        _dRdT.first = true;

        _dRdAdditionalDOF.first = true;

        _additionalDerivatives.first = true;

        addIterationData(&_jacobian);

        addIterationData(&_dRdF);

        addIterationData(&_dRdT);

        addIterationData(&_dRdAdditionalDOF);

        addIterationData(&_additionalDerivatives);
    }

    /*!
     * Get the residual vector for the non-linear problem
     */
    const floatVector *hydraBase::getResidual() {
        if (!_residual.first) {
            TARDIGRADE_ERROR_TOOLS_CATCH(formNonLinearResidual());
        }

        return &_residual.second;
    }

    /*!
     * Get the flattened row-major jacobian for the non-linear problem
     */
    const floatVector *hydraBase::getFlatJacobian() {
        if (!_jacobian.first) {
            TARDIGRADE_ERROR_TOOLS_CATCH(formNonLinearDerivatives());
        }

        return &_jacobian.second;
    }

    /*!
     * Get the jacobian for the non-linear problem
     */
    floatMatrix hydraBase::getJacobian() {
        return tardigradeVectorTools::inflate(*getFlatJacobian(), getResidual()->size(), getResidual()->size());
    }

    /*!
     * Get the flattened row-major dRdF for the non-linear problem
     */
    const floatVector *hydraBase::getFlatdRdF() {
        if (!_dRdF.first) {
            TARDIGRADE_ERROR_TOOLS_CATCH(formNonLinearDerivatives());
        }

        return &_dRdF.second;
    }

    /*!
     * Get dRdF for the non-linear problem
     */
    floatMatrix hydraBase::getdRdF() {
        return tardigradeVectorTools::inflate(*getFlatdRdF(), getResidual()->size(), getSOTDimension());
    }

    /*!
     * Get the flattened row-major dRdAdditional for the non-linear problem
     */
    const floatVector *hydraBase::getFlatdRdAdditionalDOF() {
        if (!_dRdAdditionalDOF.first) {
            TARDIGRADE_ERROR_TOOLS_CATCH(formNonLinearDerivatives());
        }

        return &_dRdAdditionalDOF.second;
    }

    /*!
     * Get dRdAdditionalDOF for the non-linear problem
     */
    floatMatrix hydraBase::getdRdAdditionalDOF() {
        return tardigradeVectorTools::inflate(*getFlatdRdAdditionalDOF(), getResidual()->size(),
                                              getAdditionalDOF()->size());
    }

    /*!
     * Get dRdT for the non-linear problem
     */
    const floatVector *hydraBase::getdRdT() {
        if (!_dRdT.first) {
            TARDIGRADE_ERROR_TOOLS_CATCH(formNonLinearDerivatives());
        }

        return &_dRdT.second;
    }

    /*!
     * Get the flattened row-major additional derivatives for the non-linear problem
     */
    const floatVector *hydraBase::getFlatAdditionalDerivatives() {
        if (!_additionalDerivatives.first) {
            TARDIGRADE_ERROR_TOOLS_CATCH(formNonLinearDerivatives());
        }

        return &_additionalDerivatives.second;
    }

    /*!
     * Get the additional derivatives for the non-linear problem
     */
    floatMatrix hydraBase::getAdditionalDerivatives() {
        if (getFlatAdditionalDerivatives()->size() > 0) {
            return tardigradeVectorTools::inflate(*getFlatAdditionalDerivatives(), getResidual()->size(),
                                                  getFlatAdditionalDerivatives()->size() / getResidual()->size());
        }

        return floatMatrix(0, floatVector(0, 0));
    }

    /*!
     * Get the stress
     */
    const floatVector *hydraBase::getStress() {
        if (!_stress.first) {
            if (getResidualClasses()->size() == 0) {
                TARDIGRADE_ERROR_TOOLS_CATCH(throw std::runtime_error("No residual classes are defined."));
            }

            TARDIGRADE_ERROR_TOOLS_CATCH(_stress.second = *(*getResidualClasses())[0]->getStress());

            if (getViscoplasticDampingSet()) {
                auto previouslyConvergedStress = getPreviouslyConvergedStress();

                for (auto v = std::begin(*previouslyConvergedStress); v != std::end(*previouslyConvergedStress); ++v) {
                    _stress.second[v - std::begin(*previouslyConvergedStress)] -= getViscoplasticDamping() * (*v);
                }
            }

            _stress.first = true;

            addIterationData(&_stress);
        }

        return &_stress.second;
    }

    /*!
     * Get the previous value of the stress
     */
    const floatVector *hydraBase::getPreviousStress() {
        if (!_previousStress.first) {
            if (getResidualClasses()->size() == 0) {
                TARDIGRADE_ERROR_TOOLS_CATCH(throw std::runtime_error("No residual classes are defined."));
            }

            TARDIGRADE_ERROR_TOOLS_CATCH(_previousStress.second = *(*getResidualClasses())[0]->getPreviousStress());

            _previousStress.first = true;
        }

        return &_previousStress.second;
    }

    /*!
     * Get the previously converged stress value
     */
    const floatVector *hydraBase::getPreviouslyConvergedStress() {
        if (!_previouslyConvergedStress.first) {
            setPreviouslyConvergedStress(*getPreviousStress());
        }

        return &_previouslyConvergedStress.second;
    }

    /*!
     * Use viscoplastic damping to reduce the current stress levels
     *
     * This should be used with care and likely only in the context of a
     * relaxed solve where it will be removed as the solution is obtained.
     *
     * \param &factor: The fraction of the difference between the last converged
     *     stress (which may be the previous stress) and the trial stress which
     *     will be suppressed
     */
    void hydraBase::setViscoplasticDamping(const floatType &factor) {
        _viscoplastic_damping_factor = factor;
        _viscoplastic_damping_set    = true;
    }

    /*!
     * Clear the viscoplastic damping
     */
    void hydraBase::clearViscoplasticDamping() {
        _viscoplastic_damping_factor = 0.;
        _viscoplastic_damping_set    = false;
    }

    /*!
     * Get whether the viscoplastic damping has been set or not
     */
    const bool hydraBase::getViscoplasticDampingSet() { return _viscoplastic_damping_set; }

    /*!
     * Reset the problem to the initial state
     */
    void hydraBase::resetProblem() {
        decomposeStateVariableVector();

        initializeUnknownVector();
    }

    /*!
     * Initialize the unknown vector for the non-linear solve.
     *
     * \f$X = \left\{ \bf{\sigma}, \bf{F}^2, \bf{F}^3, ..., \bf{F}n, \xi^1, \xi^2, ..., \xi^m \right\} \f$
     *
     * It is assumed that the first residual calculation also has a method `void getStress( )`
     * which returns a pointer to the current value of the stress.
     */
    void hydraBase::initializeUnknownVector() {
        constexpr unsigned int sot_dimension = configuration::dimension * configuration::dimension;

        const floatVector *cauchyStress;
        TARDIGRADE_ERROR_TOOLS_CATCH(cauchyStress = getStress());

        const floatVector *configurations = deformation->get_configurations();

        const unsigned int num_local_configs = configurations->size() / sot_dimension;

        TARDIGRADE_ERROR_TOOLS_CHECK(
            configurations->size() % num_local_configs == 0,
            "The size of the configurations vector must be a scalar multiple of the second order tensor size")

        const floatVector *nonLinearSolveStateVariables = get_nonLinearSolveStateVariables();

        floatVector X(getNumUnknowns(), 0);

        std::copy(std::begin(*cauchyStress), std::end(*cauchyStress), std::begin(X));

        std::copy(std::begin(*configurations) + sot_dimension, std::end(*configurations), std::begin(X) + sot_dimension);

        std::copy(std::begin(*nonLinearSolveStateVariables), std::end(*nonLinearSolveStateVariables),
                  std::begin(X) + num_local_configs * sot_dimension);

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

    /*!
     * Get the unknown vector
     */
    const floatVector *hydraBase::getUnknownVector() {
        if (!_X.first) {
            TARDIGRADE_ERROR_TOOLS_CATCH(initializeUnknownVector());
        }

        return &_X.second;
    }

    /*!
     * Update the unknown vector
     *
     * \param &newUnknownVector: The new unknown vector
     */
    void hydraBase::updateUnknownVector(const floatVector &newUnknownVector) {
        // Project the trial unknown vector to the allowable space
        floatVector trialX = newUnknownVector;
        floatVector Xp;

        if (!_X.first) {
            Xp = trialX;

        } else {
            Xp = *getUnknownVector();
        }

        setCurrentResidualIndexMeaningful(true);
        for (auto residual_ptr = getResidualClasses()->begin(); residual_ptr != getResidualClasses()->end();
             residual_ptr++) {
            setCurrentResidualIndex(residual_ptr - getResidualClasses()->begin());
            if ((*residual_ptr)->getUseProjection()) {
                (*residual_ptr)->projectSuggestedX(trialX, Xp);
            }
        }
        setCurrentResidualIndexMeaningful(false);

        // Reset all of the iteration data
        resetIterationData();

        // Set the unknown vector
        setX(trialX);

        // Decompose the unknown vector and update the state
        TARDIGRADE_ERROR_TOOLS_CATCH(decomposeUnknownVector());
    }

    /*!
     * Set if the current residual index is meaningful
     *
     * \param &value: Set if the current residual index is meaningful or not
     */
    const void hydraBase::setCurrentResidualIndexMeaningful(const bool &value) { _current_residual_index_set = value; }

    /*!
     * Set if the current residual index is meaningful
     *
     * \param value: Set the value of the current residual index
     */
    const void hydraBase::setCurrentResidualIndex(const unsigned int value) { _current_residual_index = value; }

    /*!
     * Set the scaled quantities
     */
    void hydraBase::setScaledQuantities() {
        _scaled_time = (_scale_factor - 1) * dof->_deltaTime + dof->_time;

        _scaled_deltaTime = _scale_factor * dof->_deltaTime;

        _scaled_temperature =
            _scale_factor * (dof->_temperature - dof->_previous_temperature) + dof->_previous_temperature;

        _scaled_deformationGradient =
            _scale_factor * (dof->_deformation_gradient - dof->_previous_deformation_gradient) +
            dof->_previous_deformation_gradient;

        _scaled_additionalDOF =
            _scale_factor * (dof->_additional_dof - dof->_previous_additional_dof) + dof->_previous_additional_dof;
    }

    /*!
     * Set the value of the scale factor. Will automatically re-calculate the deformation and trial stresses
     *
     * \param &value: The value of the scale factor
     */
    void hydraBase::setScaleFactor(const floatType &value) {
        // Update the scale factor
        _scale_factor = value;

        // Update the scaled quantities
        setScaledQuantities();

        // Copy the current unknown vector
        floatVector unknownVector = *getUnknownVector();

        // Reset the iteration data
        resetIterationData();

        // Set the unknown vector
        setX(unknownVector);

        // Update the deformation quantities
        updateConfigurationsFromUnknownVector();

        // Compute the new trial stress
        std::copy(getStress()->begin(), getStress()->end(), unknownVector.begin());

        // Re-set the unknown vector
        setX(unknownVector);

        // Extract the stress
        extractStress();
    }

    /*!
     * Initialize the hydra object
     */
    void hydraBase::initialize() {
        // Initialize the scaled-quantities
        setScaledQuantities();

        // Decompose the state variable vector initializing all of the configurations
        decomposeStateVariableVector();

        // Set the residual classes
        setResidualClasses();
    }

    /*!
     * Solve the non-linear problem and update the variables
     */
    void hydraBase::evaluate() {
        initialize();

        solver->solve();
    }

    /*!
     * Compute the values of the consistent tangents
     */
    void hydraBase::computeTangents() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        TARDIGRADE_ERROR_TOOLS_CHECK(solver->step != nullptr, "The step has not been defined");
        TARDIGRADE_ERROR_TOOLS_CHECK(solver->step->trial_step != nullptr, "The trial step has not been defined");
        auto local_trial_step = dynamic_cast<tardigradeHydra::NonlinearStepBase *>(solver->step->trial_step);
        TARDIGRADE_ERROR_TOOLS_CHECK(local_trial_step != nullptr, "The trial_step is not a NonlinearStepBase object");
        TARDIGRADE_ERROR_TOOLS_CHECK(local_trial_step->preconditioner != nullptr,
                                     "The preconditioner has not been defined");

        // Form the solver based on the current value of the jacobian
        floatVector P_dRdF;
        floatVector P_dRdT;
        floatVector P_A;

        local_trial_step->preconditioner->preconditionVector(*getdRdT(), P_dRdT);
        local_trial_step->preconditioner->preconditionMatrix(*getFlatdRdF(), P_dRdF);
        local_trial_step->preconditioner->preconditionMatrix(*getFlatJacobian(), P_A);

        auto P_dRdFmat = tardigradeHydra::getDynamicSizeMatrixMap(P_dRdF.data(), getResidual()->size(),
                                                                  getConfigurationUnknownCount());
        auto P_dRdTmat = tardigradeHydra::getDynamicSizeVectorMap(P_dRdT.data(), getResidual()->size());
        auto P_Amat =
            tardigradeHydra::getDynamicSizeMatrixMap(P_A.data(), getResidual()->size(), getResidual()->size());

        _flatdXdF.second = floatVector(getNumUnknowns() * getConfigurationUnknownCount());
        auto dXdFmat     = tardigradeHydra::getDynamicSizeMatrixMap(_flatdXdF.second.data(), getNumUnknowns(),
                                                                    getConfigurationUnknownCount());

        _flatdXdT.second = floatVector(getNumUnknowns());
        auto dXdTmat     = tardigradeHydra::getDynamicSizeVectorMap(_flatdXdT.second.data(), getNumUnknowns());

        // Solve
        tardigradeVectorTools::solverType<floatType> linear_solver(P_Amat);

        dXdFmat = -linear_solver.solve(P_dRdFmat);
        dXdTmat = -linear_solver.solve(P_dRdTmat);

        unsigned int rank = linear_solver.rank();

        TARDIGRADE_ERROR_TOOLS_CATCH(

            if (solver->step->getRankDeficientError() && (rank != getResidual()->size())) {
                throw convergence_error("The Jacobian is not full rank");
            }

        )

        _flatdXdF.first = true;

        _flatdXdT.first = true;
    }

    /*!
     * Compute the consistent tangent w.r.t. the additional dof
     */
    void hydraBase::computedXdAdditionalDOF() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        TARDIGRADE_ERROR_TOOLS_CHECK(solver->step != nullptr, "The step has not been defined");
        TARDIGRADE_ERROR_TOOLS_CHECK(solver->step->trial_step != nullptr, "The trial step has not been defined");
        auto local_trial_step = dynamic_cast<tardigradeHydra::NonlinearStepBase *>(solver->step->trial_step);
        TARDIGRADE_ERROR_TOOLS_CHECK(local_trial_step->preconditioner != nullptr,
                                     "The preconditioner has not been defined");

        // Form the solver based on the current value of the jacobian
        floatVector P_A;
        floatVector P_dRdAdditionalDOF;

        local_trial_step->preconditioner->preconditionMatrix(*getFlatJacobian(), P_A);
        local_trial_step->preconditioner->preconditionMatrix(*getFlatdRdAdditionalDOF(), P_dRdAdditionalDOF);

        auto P_Amat =
            tardigradeHydra::getDynamicSizeMatrixMap(P_A.data(), getResidual()->size(), getResidual()->size());

        auto P_dRdAdditionalDOFmat =
            tardigradeHydra::getDynamicSizeMatrixMap(P_dRdAdditionalDOF.data(), getResidual()->size(),
                                                     getAdditionalDOF()->size());

        // Form the map for dXdF
        _flatdXdAdditionalDOF.second = floatVector(getNumUnknowns() * getAdditionalDOF()->size(), 0);
        auto dXdAdditionalDOF = tardigradeHydra::getDynamicSizeMatrixMap(_flatdXdAdditionalDOF.second.data(),
                                                                         getNumUnknowns(), getAdditionalDOF()->size());

        // Solve
        tardigradeVectorTools::solverType<floatType> linear_solver(P_Amat);

        dXdAdditionalDOF = -linear_solver.solve(P_dRdAdditionalDOFmat);

        _flatdXdAdditionalDOF.first = true;
    }

    /*!
     * Get the total derivative of X w.r.t. the deformation in row-major format
     */
    const floatVector *hydraBase::getFlatdXdF() {
        if (!_flatdXdF.first) {
            TARDIGRADE_ERROR_TOOLS_CATCH(computeTangents())
        }

        return &_flatdXdF.second;
    }

    /*!
     * Get the total derivative of X w.r.t. the temperature in row-major format
     */
    const floatVector *hydraBase::getFlatdXdT() {
        if (!_flatdXdT.first) {
            TARDIGRADE_ERROR_TOOLS_CATCH(computeTangents())
        }

        return &_flatdXdT.second;
    }

    /*!
     * Get the total derivative of X w.r.t. the additional degrees of freedom
     */
    const floatVector *hydraBase::getFlatdXdAdditionalDOF() {
        if (!_flatdXdAdditionalDOF.first) {
            TARDIGRADE_ERROR_TOOLS_CATCH(computedXdAdditionalDOF());
        }

        return &_flatdXdAdditionalDOF.second;
    }

    /*!
     * Set the constraint values
     */
    void hydraBase::setConstraints() {
        const unsigned int numConstraints = getNumConstraints();

        auto constraints = get_SetDataStorage_constraints();

        constraints.zero(numConstraints);

        unsigned int offset = 0;

        setCurrentResidualIndexMeaningful(true);
        for (auto v = getResidualClasses()->begin(); v != getResidualClasses()->end(); v++) {
            setCurrentResidualIndexMeaningful(v - getResidualClasses()->begin());

            if ((*v)->getNumConstraints() > 0) {
                std::copy((*v)->getConstraints()->begin(), (*v)->getConstraints()->end(),
                          constraints.value->begin() + offset);

                offset += (*v)->getConstraints()->size();
            }
        }
        setCurrentResidualIndexMeaningful(false);
    }

    /*!
     * Set the constraint Jacobians values
     */
    void hydraBase::setConstraintJacobians() {
        const unsigned int numUnknowns = getNumUnknowns();

        const unsigned int numConstraints = getNumConstraints();

        auto constraintJacobians = get_SetDataStorage_constraintJacobians();

        constraintJacobians.zero(numConstraints * numUnknowns);

        unsigned int offset = 0;

        setCurrentResidualIndexMeaningful(true);
        for (auto v = getResidualClasses()->begin(); v != getResidualClasses()->end(); v++) {
            setCurrentResidualIndex(v - getResidualClasses()->begin());

            if ((*v)->getNumConstraints() > 0) {
                std::copy((*v)->getConstraintJacobians()->begin(), (*v)->getConstraintJacobians()->end(),
                          constraintJacobians.value->begin() + offset);

                offset += (*v)->getConstraintJacobians()->size();
            }
        }
        setCurrentResidualIndexMeaningful(false);
    }

    /*!
     * Get the offset of the current residual
     */
    const unsigned int hydraBase::getCurrentResidualOffset() {
        unsigned int offset = 0;
        for (auto v = getResidualClasses()->begin(); v != getResidualClasses()->begin() + getCurrentResidualIndex();
             ++v) {
            offset += (*v)->getNumEquations();
        }
        return offset;
    }

    /*!
     * Get the value of the number of constraint equations
     */
    const unsigned int hydraBase::getNumConstraints() {
        unsigned int value = 0;

        for (auto v = getResidualClasses()->begin(); v != getResidualClasses()->end(); v++) {
            value += (*v)->getNumConstraints();
        }

        return value;
    }

    /*!
     * Get the parameterization information of the residual classes
     */
    std::string hydraBase::getResidualParameterizationInfo() {
        std::string parameterization_info =
            "########################################\n# RESIDUAL PARAMETERIZATION "
            "INFORMATION#\n########################################\n\n";

        for (auto v = std::begin(*getResidualClasses()); v != std::end(*getResidualClasses()); ++v) {
            parameterization_info += "RESIDUAL CLASS:";
            parameterization_info +=
                " " + std::to_string((unsigned int)(v - std::begin(*getResidualClasses()))) + "\n\n";

            (*v)->addParameterizationInfo(parameterization_info);

            parameterization_info += "\n";
        }

        return parameterization_info;
    }

    /*!
     * Function to throw for an unexpected error. A user should never get here!
     */
    void hydraBase::unexpectedError() {
        TARDIGRADE_ERROR_TOOLS_CATCH(
            throw std::runtime_error("You shouldn't have gotten here. If you aren't developing the code then "
                                     "contact a developer with the stack trace."))
    }

    /*!
     * Set the value of the unknown vector
     *
     * \param &X: The unknown vector
     */
    void hydraBase::setX(const floatVector &X) {
        _X.second = X;

        _X.first = true;
    }

    /*!
     * Set a flag for if the global residual can be modified
     *
     * \param value: The updated value
     */
    void hydraBase::setAllowModifyGlobalResidual(const bool value) { _allow_modify_global_residual = value; }

    /*!
     * Set a flag for if the global jacobian can be modified
     *
     * \param value: The updated value
     */
    void hydraBase::setAllowModifyGlobalJacobian(const bool value) { _allow_modify_global_jacobian = value; }

    /*!
     * Set a flag for if the global dRdT can be modified
     *
     * \param value: The updated value
     */
    void hydraBase::setAllowModifyGlobaldRdT(const bool value) { _allow_modify_global_dRdT = value; }

    /*!
     * Set a flag for if the global dRdF can be modified
     *
     * \param value: The updated value
     */
    void hydraBase::setAllowModifyGlobaldRdF(const bool value) { _allow_modify_global_dRdF = value; }

    /*!
     * Set a flag for if the global dRdAdditionalDOF can be modified
     *
     * \param value: The updated value
     */
    void hydraBase::setAllowModifyGlobaldRdAdditionalDOF(const bool value) {
        _allow_modify_global_dRdAdditionalDOF = value;
    }

    /*!
     * Set the value of the previously converged stress.
     *
     * \param &value: The incoming value
     */
    void hydraBase::setPreviouslyConvergedStress(const floatVector &value) {
        _previouslyConvergedStress.second = value;
        _previouslyConvergedStress.first  = true;
    }
}  // namespace tardigradeHydra
