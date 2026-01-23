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
     * Compute the trial step
     *
     * Must set the containing step's deltaX variable
     */
    void NonlinearStepBase::computeTrial() {
        if (getUseSQPSolver()) {
            solveConstrainedQP(step->deltaX);

        } else {
            solveNewtonUpdate(step->deltaX);
        }

        if (getFailureVerbosityLevel() > 0) {
            addToFailureOutput("  trial deltaX:\n");
            addToFailureOutput("  ");
            addToFailureOutput(step->deltaX);
        }
    }

    // BEGIN NEWTON SOLVER FUNCTIONS

    /*!
     * Solve the Newton update returning the trial value of the unknown vector
     *
     * \param &deltaX_tr: The trial change in the unknown vector
     */
    void NonlinearStepBase::solveNewtonUpdate(floatVector &deltaX_tr) {
        TARDIGRADE_ERROR_TOOLS_CHECK(preconditioner != nullptr,
                                     "The preconditioner has not been defined");  // TODO: Move to the trial_step class
        if (preconditioner->getUsePreconditioner()) {
            performPreconditionedSolve(deltaX_tr);

        } else {
            auto dx_map = tardigradeHydra::getDynamicSizeVectorMap(deltaX_tr.data(), getNumUnknowns());

            auto J_map = tardigradeHydra::getDynamicSizeMatrixMap(getFlatNonlinearLHS()->data(), getNumUnknowns(),
                                                                  getNumUnknowns());

            auto R_map = tardigradeHydra::getDynamicSizeVectorMap(getNonlinearRHS()->data(), getNumUnknowns());

            tardigradeVectorTools::solverType<floatType> linearSolver(J_map);
            dx_map = -linearSolver.solve(R_map);

            unsigned int rank = linearSolver.rank();

            if (getRankDeficientError() && (rank != getResidual()->size())) {
                TARDIGRADE_ERROR_TOOLS_CATCH(throw convergence_error("The Jacobian is not full rank"));
            }
        }
    }

    // END NEWTON SOLVER FUNCTIONS

}
