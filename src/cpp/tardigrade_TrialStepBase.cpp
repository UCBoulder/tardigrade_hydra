/**
 ******************************************************************************
 * \file tardigrade_TrialStepBase.cpp
 ******************************************************************************
 * The base class for step damping operations
 ******************************************************************************
 */

#include "tardigrade_TrialStepBase.h"
#define USE_EIGEN
#include "tardigrade_CustomErrors.h"
#include "tardigrade_SolverStepBase.h"
#include "tardigrade_vector_tools.h"

namespace tardigradeHydra {

    /*!
     * The constructor for TrialStepBase
     */
    TrialStepBase::TrialStepBase() : step(NULL) { }

    /*!
     * The constructor for TrialStepBase
     *
     * \param *_step: The containing step object
     */
    TrialStepBase::TrialStepBase(SolverStepBase *_step) : step(_step) {
        step->trial_step           = this;
    }

    /*!
     * Reset the counters
     */
    void TrialStepBase::resetCounts() {}

    /*!
     * Reset the trial step class
     */
    void TrialStepBase::reset() {
        resetCounts();
    }

    /*!
     * Compute the trial step
     *
     * Must set the containing step's deltaX variable
     */
    void TrialStepBase::computeTrial() {
        TARDIGRADE_ERROR_TOOLS_CATCH(throw std::runtime_error("computeTrial must be defined for any class inheriting from TrialStepBase"));
    }

    /*!
     * Add data to the vector of values which will be cleared after each iteration
     *
     * \param *data: The dataBase object to be cleared
     */
    void TrialStepBase::addIterationData(dataBase *data) {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        step->addIterationData(data);
    }

    /*!
     * Add data to the vector of values which will be cleared after each nonlinear step
     *
     * \param *data: The dataBase object to be cleared
     */
    void TrialStepBase::addNLStepData(dataBase *data) {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        step->addNLStepData(data);
    }

    /*!
     * Get the relative tolerance value
     */
    const floatType TrialStepBase::getRelativeTolerance() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        return step->getRelativeTolerance();
    }

    /*!
     * Get the absolute tolerance value
     */
    const floatType TrialStepBase::getAbsoluteTolerance() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        return step->getAbsoluteTolerance();
    }
    /*!
     * Get the residual vector
     */
    const floatVector *TrialStepBase::getResidual() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        return step->getResidual();
    }

    /*!
     * Get the number of unknowns
     */
    const unsigned int TrialStepBase::getNumUnknowns() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        return step->getNumUnknowns();
    }

    /*!
     * Get the Jacobian in row-major format
     */
    const floatVector *TrialStepBase::getFlatJacobian() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        return step->getFlatJacobian();
    }
    /*!
     * Get the number of constraint equations
     */
    const unsigned int TrialStepBase::getNumConstraints() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        return step->getNumConstraints();
    }

    /*!
     * Get the current constraint values
     */
    const floatVector *TrialStepBase::getConstraints() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        return step->getConstraints();
    }

    /*!
     * Get the constraint Jacobians
     */
    const floatVector *TrialStepBase::getConstraintJacobians() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        return step->getConstraintJacobians();
    }

    /*!
     * Get whether a rank-deficient matrix will throw an error
     */
    bool TrialStepBase::getRankDeficientError() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        return step->getRankDeficientError();
    }
    /*!
     * Get the failure verbosity level
     */
    const unsigned int TrialStepBase::getFailureVerbosityLevel() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        return step->getFailureVerbosityLevel();
    }

    /*!
     * Add the string to the failure output message
     *
     * \param &string: The string to add to the failure output message
     */
    void TrialStepBase::addToFailureOutput(const std::string &string) {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        step->addToFailureOutput(string);
    }

    /*!
     * Add a floatVector to the failure output message
     *
     * \param &value: The floatVector to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void TrialStepBase::addToFailureOutput(const floatVector &value, bool add_endline) {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        step->addToFailureOutput(value, add_endline);
    }

    /*!
     * Add a vector of booleans to the failure output message
     *
     * \param &value: The vector of booleans to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void TrialStepBase::addToFailureOutput(const std::vector<bool> &value, bool add_endline) {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        step->addToFailureOutput(value, add_endline);
    }

    /*!
     * Add a floatType to the failure output message
     *
     * \param &value: The floatType to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void TrialStepBase::addToFailureOutput(const floatType &value, bool add_endline) {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        step->addToFailureOutput(value, add_endline);
    }

    /*!
     * Get a pointer to the damping object
     */
    StepDampingBase *TrialStepBase::getDamping( ){
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        return step->damping;
    }
}  // namespace tardigradeHydra
