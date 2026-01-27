/**
  ******************************************************************************
  * \file tardigrade_NonlinearStepBase.cpp
  ******************************************************************************
  * The source file for the base class to determine the nonlinear step
  ******************************************************************************
  */

#include"tardigrade_NonlinearStepBase.h"
#define USE_EIGEN
#include "tardigrade_CustomErrors.h"
#include "tardigrade_SolverStepBase.h"
#include "tardigrade_vector_tools.h"

namespace tardigradeHydra{

    /*!
     * The constructor for NonlinearStepBase
     */
    NonlinearStepBase::NonlinearStepBase() : TrialStepBase() { preconditioner->trial_step = this; }

    /*!
     * The constructor for NonlinearStepBase
     *
     * \param *_step: The containing step object
     */
    NonlinearStepBase::NonlinearStepBase(SolverStepBase *_step) : TrialStepBase(_step) {
        preconditioner->trial_step = this;
    }

    /*!
     * The constructor for NonlinearStepBase
     *
     * \param *_step: The containing step object
     * \param *_preconditioner_ptr: The preconditioner object used by the trial step
     */
    NonlinearStepBase::NonlinearStepBase(SolverStepBase *_step, PreconditionerBase *_preconditioner_ptr)
        : TrialStepBase(_step), preconditioner(_preconditioner_ptr) {
        step->trial_step           = this;
        preconditioner->trial_step = this;
    }

    /*!
     * Compute the trial step
     *
     * Must set the containing step's deltaX variable
     */
    void NonlinearStepBase::computeTrial() {
        TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error("computeTrial must be defined by inheriting classes") );
    }

    /*!
     * Add output to the failure message
     */
    void NonlinearStepBase::addTrialStepOutput() {

        if (getFailureVerbosityLevel() > 0) {
            TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
            addToFailureOutput("  trial deltaX:\n");
            addToFailureOutput("  ");
            addToFailureOutput(step->deltaX);
        }
    }

    /*!
     * Reset the trial step class
     */
    void NonlinearStepBase::reset() {
        TrialStepBase::reset();
        TARDIGRADE_ERROR_TOOLS_CHECK(preconditioner != nullptr, "The preconditioner has not been defined");
        preconditioner->reset();
    }

    /*!
     * Get the RHS vector for the non-linear problem
     */
    const floatVector *NonlinearStepBase::getNonlinearRHS() { return getResidual(); }

    /*!
     * Get the flat LHS matrix for the non-linear problem
     */
    const floatVector *NonlinearStepBase::getFlatNonlinearLHS() { return getFlatJacobian(); }

}
