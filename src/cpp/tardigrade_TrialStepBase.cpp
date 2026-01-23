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
    TrialStepBase::TrialStepBase() : step(NULL) { preconditioner->trial_step = this; }

    /*!
     * The constructor for TrialStepBase
     *
     * \param *_step: The containing step object
     */
    TrialStepBase::TrialStepBase(SolverStepBase *_step) : step(_step) {
        step->trial_step           = this;
        preconditioner->trial_step = this;
    }

    /*!
     * The constructor for TrialStepBase
     *
     * \param *_step: The containing step object
     * \param *_preconditioner_ptr: The preconditioner object used by the trial step
     */
    TrialStepBase::TrialStepBase(SolverStepBase *_step, PreconditionerBase *_preconditioner_ptr)
        : step(_step), preconditioner(_preconditioner_ptr) {
        step->trial_step           = this;
        preconditioner->trial_step = this;
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
        TARDIGRADE_ERROR_TOOLS_CHECK(preconditioner != nullptr, "The preconditioner has not been defined");
        preconditioner->reset();
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

    // BEGIN NONLINEAR SOLVER FUNCTIONS

    /*!
     * Get the RHS vector for the non-linear problem
     */
    const floatVector *TrialStepBase::getNonlinearRHS() { return getResidual(); }

    /*!
     * Get the flat LHS matrix for the non-linear problem
     */
    const floatVector *TrialStepBase::getFlatNonlinearLHS() { return getFlatJacobian(); }

    // END NONLINEAR SOLVER FUNCTIONS

    // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

    /*!
     * Initialize the active constraint vector
     *
     * \param &active_constraints: The current constraints that are active
     */
    void TrialStepBase::initializeActiveConstraints(std::vector<bool> &active_constraints) {
        active_constraints = std::vector<bool>(getNumConstraints(), false);

        for (auto c = getConstraints()->begin(); c != getConstraints()->end(); c++) {
            unsigned int index = (unsigned int)(c - getConstraints()->begin());

            active_constraints[index] = ((*c) < 0.);
        }
    }

    // END SQP SOLVER FUNCTIONS

    /*!
     * Perform a pre-conditioned solve
     *
     * \param &deltaX_tr: The trial chcange in the unknown vector
     */
    void TrialStepBase::performPreconditionedSolve(floatVector &deltaX_tr) {
        TARDIGRADE_ERROR_TOOLS_CHECK(preconditioner != nullptr,
                                     "The preconditioner has not been defined");  // TODO: Move to the trial_step class
        tardigradeVectorTools::solverType<floatType> linearSolver;

        auto dx_map = tardigradeHydra::getDynamicSizeVectorMap(deltaX_tr.data(), getNumUnknowns());

        auto J_map =
            tardigradeHydra::getDynamicSizeMatrixMap(getFlatNonlinearLHS()->data(), getNumUnknowns(), getNumUnknowns());

        auto R_map = tardigradeHydra::getDynamicSizeVectorMap(getNonlinearRHS()->data(), getNumUnknowns());

        if (preconditioner->getPreconditionerIsDiagonal()) {
            auto p_map = tardigradeHydra::getDynamicSizeVectorMap(preconditioner->getFlatPreconditioner()->data(),
                                                                  getNumUnknowns());

            linearSolver = tardigradeVectorTools::solverType<floatType>(p_map.asDiagonal() * J_map);

            dx_map = -linearSolver.solve(p_map.asDiagonal() * R_map);

        } else {
            auto p_map = tardigradeHydra::getDynamicSizeMatrixMap(preconditioner->getFlatPreconditioner()->data(),
                                                                  getNumUnknowns(), getNumUnknowns());

            linearSolver = tardigradeVectorTools::solverType<floatType>(p_map * J_map);

            dx_map = -linearSolver.solve(p_map * R_map);
        }

        unsigned int rank = linearSolver.rank();

        if (getRankDeficientError() && (rank != getResidual()->size())) {
            TARDIGRADE_ERROR_TOOLS_CATCH(throw convergence_error("The Jacobian is not full rank"));
        }
    }

}  // namespace tardigradeHydra
