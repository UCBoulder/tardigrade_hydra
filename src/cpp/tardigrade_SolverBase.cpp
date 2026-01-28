/**
 ******************************************************************************
 * \file tardigrade_SolverBase.cpp
 ******************************************************************************
 * A C++ library for the base classes for solvers
 ******************************************************************************
 */

#include "tardigrade_SolverBase.h"

#include "tardigrade_hydra.h"

namespace tardigradeHydra {

    /*!
     * Default constructor of SolverBase
     */
    SolverBase::SolverBase() : hydra(NULL), step(&_step) { step->setSolver(this); }

    /*!
     * Constructor for SolverBase
     *
     * \param *_hydra: The containing hydraBase object
     */
    SolverBase::SolverBase(hydraBase *_hydra) : hydra(_hydra), step(&_step) { step->setSolver(this); }

    /*!
     * Constructor for SolverBase
     *
     * \param *_hydra: The containing hydraBase object
     * \param *_step_ptr: The SolverStepBase object that the solver will use to try and solve
     *     the problem
     */
    SolverBase::SolverBase(hydraBase *_hydra, SolverStepBase *_step_ptr) : hydra(_hydra), step(_step_ptr) {
        step->setSolver(this);
    }

    /*!
     * Add data to the vector of values which will be cleared after each iteration
     *
     * \param *data: The dataBase object to be cleared
     */
    void SolverBase::addIterationData(dataBase *data) {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        hydra->addIterationData(data);
    }

    /*!
     * Add data to the vector of values which will be cleared after each nonlinear step
     *
     * \param *data: The dataBase object to be cleared
     */
    void SolverBase::addNLStepData(dataBase *data) {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        hydra->addNLStepData(data);
    }

    /*!
     * Get the relative tolerance from hydra
     */
    const floatType SolverBase::getRelativeTolerance() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        return hydra->getRelativeTolerance();
    }

    /*!
     * Get the absolute tolerance from hydra
     */
    const floatType SolverBase::getAbsoluteTolerance() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        return hydra->getAbsoluteTolerance();
    }

    /*!
     * Get the number of unknowns from hydra
     */
    const unsigned int SolverBase::getNumUnknowns() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        return hydra->getNumUnknowns();
    }

    /*!
     * Get the unknown vector from hydra
     */
    const floatVector *SolverBase::getUnknownVector() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        return hydra->getUnknownVector();
    }

    /*!
     * Initialize the unknown vector
     */
    void SolverBase::initializeUnknownVector() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        hydra->initializeUnknownVector();
    }

    /*!
     * Update the unknown vector
     *
     * \param &value: The new value of the unknown vector
     */
    void SolverBase::updateUnknownVector(const floatVector &value) {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        hydra->updateUnknownVector(value);
    }

    /*!
     * Get the residual from hydra
     */
    const floatVector *SolverBase::getResidual() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        return hydra->getResidual();
    }

    /*!
     * Get the flat Jacobian from hydra
     */
    const floatVector *SolverBase::getFlatJacobian() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        return hydra->getFlatJacobian();
    }

    /*!
     * Get the number of constraints
     */
    const unsigned int SolverBase::getNumConstraints() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        return hydra->getNumConstraints();
    }

    /*!
     * Get the constraints
     */
    const floatVector *SolverBase::getConstraints() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        return hydra->getConstraints();
    }

    /*!
     * Get the constraint Jacobians
     */
    const floatVector *SolverBase::getConstraintJacobians() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        return hydra->getConstraintJacobians();
    }

    /*!
     * Get the verbosity level for failure messages
     */
    const unsigned int SolverBase::getFailureVerbosityLevel() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        return hydra->getFailureVerbosityLevel();
    }

    /*!
     * Add the string to the failure output message
     *
     * \param &string: The string to add to the failure output message
     */
    void SolverBase::addToFailureOutput(const std::string &string) {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        hydra->addToFailureOutput(string);
    }

    /*!
     * Add a floatVector to the failure output message
     *
     * \param &value: The floatVector to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void SolverBase::addToFailureOutput(const floatVector &value, bool add_endline) {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        hydra->addToFailureOutput(value, add_endline);
    }

    /*!
     * Add a vector of booleans to the failure output message
     *
     * \param &value: The vector of booleans to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void SolverBase::addToFailureOutput(const std::vector<bool> &value, bool add_endline) {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        hydra->addToFailureOutput(value, add_endline);
    }

    /*!
     * Add a floatType to the failure output message
     *
     * \param &value: The floatType to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void SolverBase::addToFailureOutput(const floatType &value, bool add_endline) {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        hydra->addToFailureOutput(value, add_endline);
    }

    /*!
     * Get the scale factor for the tolerance
     */
    const floatType SolverBase::getToleranceScaleFactor() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        return hydra->getToleranceScaleFactor();
    }

    /*!
     * Reset the tolerance scale factor
     */
    void SolverBase::resetToleranceScaleFactor() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        hydra->resetToleranceScaleFactor();
    }

    /*!
     * Set if the current residual index is meaningful or not
     *
     * \param &value: The boolean indicating if the residual index is or isn't meaningful
     */
    void SolverBase::setCurrentResidualIndexMeaningful(const bool &value) {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        hydra->setCurrentResidualIndexMeaningful(value);
    }

    /*!
     * Set the current residual index
     *
     * \param &value: The value of the current residual's index
     */
    void SolverBase::setCurrentResidualIndex(const unsigned int &value) {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        hydra->setCurrentResidualIndex(value);
    }

    /*!
     * Get the residual classes
     */
    const std::vector<tardigradeHydra::ResidualBase<> *> *SolverBase::getResidualClasses() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        return hydra->getResidualClasses();
    }

    /*!
     * Set that the global residual can be modified
     */
    void SolverBase::setAllowModifyGlobalResidual(const bool &value) {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        return hydra->setAllowModifyGlobalResidual(value);
    }

    /*!
     * Reset the solver counts
     */
    void SolverBase::resetCounts() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        TARDIGRADE_ERROR_TOOLS_CATCH(step->resetCounts();)
    }

    /*!
     * Reset the solver
     */
    void SolverBase::reset() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        TARDIGRADE_ERROR_TOOLS_CATCH(step->reset();)
    }

    /*!
     * The function that is called prior to the first attempt
     * of solving the problem
     */
    void SolverBase::initializeSolve() {}

    /*!
     * The function that is called when first attempting to
     * solve the problem
     */
    void SolverBase::initialSolveAttempt() {
        TARDIGRADE_ERROR_TOOLS_CATCH(
            throw std::runtime_error("initialSolveAttempt must be defined for all classes inheriting from SolverBase"));
    }

    /*!
     * The function that is called if there is a convergence
     * error thrown in the initial solve attempt
     */
    void SolverBase::convergenceErrorFunction() { throw; }

    /*!
     * The function that is called if there is an unexpected
     * error thrown in the initial solve attempt
     */
    void SolverBase::unexpectedErrorFunction() { TARDIGRADE_ERROR_TOOLS_CATCH(throw;) }

    /*!
     * Solve the problem
     */
    void SolverBase::solve() {
        initializeSolve();
        try {
            initialSolveAttempt();

        } catch (convergence_error &e) {
            convergenceErrorFunction();

        } catch (std::exception &e) {
            unexpectedErrorFunction();
        }
    }

    // NONLINEAR FUNCTIONS

    /*!
     * Reset all nonlinear step data
     */
    void SolverBase::resetNLStepData() {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        hydra->resetNLStepData();
    }

    // END NONLINEAR FUNCTIONS

}  // namespace tardigradeHydra
