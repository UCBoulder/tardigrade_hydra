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
}
