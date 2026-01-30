/**
 ******************************************************************************
 * \file tardigrade_SolverStepBase.cpp
 ******************************************************************************
 * The base class for solver steps
 ******************************************************************************
 */

#include "tardigrade_SolverStepBase.h"
#define USE_EIGEN
#include "tardigrade_CustomErrors.h"
#include "tardigrade_SolverBase.h"
#include "tardigrade_StepDampingBase.h"
#include "tardigrade_hydra.h"
#include "tardigrade_vector_tools.h"

namespace tardigradeHydra {

    /*!
     * Constructor for NonlinearStepBase
     */
    SolverStepBase::SolverStepBase() : solver(NULL) { }

    /*!
     * Constructor for NonlinearStepBase
     *
     * \param *_solver: The containing solver object
     */
    SolverStepBase::SolverStepBase(SolverBase *_solver) : solver(_solver) { }

    /*!
     * Reset the counts
     */
    void SolverStepBase::resetCounts() {
        TARDIGRADE_ERROR_TOOLS_CATCH(resetNumUndamped();)
        TARDIGRADE_ERROR_TOOLS_CATCH(trial_step->resetCounts();)
        TARDIGRADE_ERROR_TOOLS_CATCH(damping->resetCounts();)
    }

    /*!
     * Reset the step back to an initial state
     */
    void SolverStepBase::reset() {
        TARDIGRADE_ERROR_TOOLS_CATCH(resetCounts();)
        TARDIGRADE_ERROR_TOOLS_CATCH(trial_step->reset();)
        TARDIGRADE_ERROR_TOOLS_CATCH(damping->reset();)
    }

    /*!
     * Get the relative tolerance value
     */
    const floatType SolverStepBase::getRelativeTolerance() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        return solver->getRelativeTolerance();
    }

    /*!
     * Get the absolute tolerance value
     */
    const floatType SolverStepBase::getAbsoluteTolerance() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        return solver->getAbsoluteTolerance();
    }

    /*!
     * Add data to the vector of values which will be cleared after each iteration
     *
     * \param *data: The dataBase object to be cleared
     */
    void SolverStepBase::addIterationData(dataBase *data) {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        solver->addIterationData(data);
    }

    /*!
     * Add data to the vector of values which will be cleared after each nonlinear step
     *
     * \param *data: The dataBase object to be cleared
     */
    void SolverStepBase::addNLStepData(dataBase *data) {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        solver->addNLStepData(data);
    }

    /*!
     * Get the residual vector
     */
    const floatVector *SolverStepBase::getResidual() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        return solver->getResidual();
    }

    /*!
     * Update the unknown vector
     *
     * \param &value: The new value of the unknown vector
     */
    void SolverStepBase::updateUnknownVector(const floatVector &value) {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        solver->updateUnknownVector(value);
    }

    /*!
     * Get the number of unknowns
     */
    const unsigned int SolverStepBase::getNumUnknowns() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        return solver->getNumUnknowns();
    }

    /*!
     * Get the unknown vector
     */
    const floatVector *SolverStepBase::getUnknownVector() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        return solver->getUnknownVector();
    }

    /*!
     * Get the Jacobian in row-major format
     */
    const floatVector *SolverStepBase::getFlatJacobian() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        return solver->getFlatJacobian();
    }

    /*!
     * Get the number of constraint equations
     */
    const unsigned int SolverStepBase::getNumConstraints() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        return solver->getNumConstraints();
    }

    /*!
     * Get the current constraint values
     */
    const floatVector *SolverStepBase::getConstraints() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        return solver->getConstraints();
    }

    /*!
     * Get the constraint Jacobians
     */
    const floatVector *SolverStepBase::getConstraintJacobians() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        return solver->getConstraintJacobians();
    }

    /*!
     * Get the scale factor for the tolerance
     */
    const floatType SolverStepBase::getToleranceScaleFactor() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        return solver->getToleranceScaleFactor();
    }

    /*!
     * Reset the tolerance scale factor
     */
    void SolverStepBase::resetToleranceScaleFactor() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        solver->resetToleranceScaleFactor();
    }

    /*!
     * Get the failure verbosity level
     */
    const unsigned int SolverStepBase::getFailureVerbosityLevel() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        return solver->getFailureVerbosityLevel();
    }

    /*!
     * Add the string to the failure output message
     *
     * \param &string: The string to add to the failure output message
     */
    void SolverStepBase::addToFailureOutput(const std::string &string) {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        solver->addToFailureOutput(string);
    }

    /*!
     * Add a floatVector to the failure output message
     *
     * \param &value: The floatVector to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void SolverStepBase::addToFailureOutput(const floatVector &value, bool add_endline) {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        solver->addToFailureOutput(value, add_endline);
    }

    /*!
     * Add a vector of booleans to the failure output message
     *
     * \param &value: The vector of booleans to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void SolverStepBase::addToFailureOutput(const std::vector<bool> &value, bool add_endline) {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        solver->addToFailureOutput(value, add_endline);
    }

    /*!
     * Add a floatType to the failure output message
     *
     * \param &value: The floatType to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void SolverStepBase::addToFailureOutput(const floatType &value, bool add_endline) {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        solver->addToFailureOutput(value, add_endline);
    }

    /*!
     * Set if the current residual index is meaningful or not
     *
     * \param &value: The boolean indicating if the residual index is or isn't meaningful
     */
    void SolverStepBase::setCurrentResidualIndexMeaningful(const bool &value) {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        solver->setCurrentResidualIndexMeaningful(value);
    }

    /*!
     * Set the current residual index
     *
     * \param &value: The value of the current residual's index
     */
    void SolverStepBase::setCurrentResidualIndex(const unsigned int &value) {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        solver->setCurrentResidualIndex(value);
    }

    /*!
     * Get the residual classes
     */
    const std::vector<tardigradeHydra::ResidualBase<> *> *SolverStepBase::getResidualClasses() {
        TARDIGRADE_ERROR_TOOLS_CHECK(solver != nullptr, "The solver has not been defined");
        return solver->getResidualClasses();
    }

    /*!
     * Increment the solution of the problem
     */
    void SolverStepBase::incrementSolution() {
        TARDIGRADE_ERROR_TOOLS_CHECK(trial_step != nullptr, "The trial step has not been defined");
        TARDIGRADE_ERROR_TOOLS_CHECK(damping != nullptr, "The damping has not been defined");

        X0     = *getUnknownVector();
        deltaX = floatVector(getNumUnknowns(), 0);

        if (getFailureVerbosityLevel() > 0) {
            addToFailureOutput("  X0:\n");
            addToFailureOutput("  ");
            addToFailureOutput(*getUnknownVector());
        }

        damping->setBaseQuantities();

        trial_step->computeTrial();

        if (!damping->applyDamping()) {
            incrementNumUndamped();
        }
    }

    /*!
     * Get whether the Jacobian being rank-deficient will throw an error
     */
    const bool SolverStepBase::getRankDeficientError() { return _rank_deficient_error; }

    /*!
     * Set whether the Jacobian being rank-deficient will throw an error
     *
     * \param &value: The incoming value
     */
    void SolverStepBase::setRankDeficientError(const bool &value) { _rank_deficient_error = value; }

}  // namespace tardigradeHydra
