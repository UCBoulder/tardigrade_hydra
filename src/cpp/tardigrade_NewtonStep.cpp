/**
 ******************************************************************************
 * \file tardigrade_NewtonStep.cpp
 ******************************************************************************
 * A class which defines a Newton-Raphson step
 ******************************************************************************
 */

#include "tardigrade_NewtonStep.h"

#include "tardigrade_CustomErrors.h"
#include "tardigrade_SolverStepBase.h"
#define USE_EIGEN
#include "tardigrade_vector_tools.h"

namespace tardigradeHydra {

    /*!
     * Compute the trial step
     */
    void NewtonStep::computeTrial() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");

        TARDIGRADE_ERROR_TOOLS_CHECK(preconditioner != nullptr, "The preconditioner has not been defined");

        auto dx_map = tardigradeHydra::getDynamicSizeVectorMap(step->deltaX.data(), getNumUnknowns());

        floatVector P_R;
        floatVector P_J;

        preconditioner->preconditionVector(*getNonlinearRHS(), P_R);
        preconditioner->preconditionMatrix(*getFlatNonlinearLHS(), P_J);

        auto P_R_map = tardigradeHydra::getDynamicSizeVectorMap(P_R.data(), getNumUnknowns());

        auto P_J_map = tardigradeHydra::getDynamicSizeMatrixMap(P_J.data(), getNumUnknowns(), getNumUnknowns());

        tardigradeVectorTools::solverType<floatType> linearSolver(P_J_map);

        dx_map = -linearSolver.solve(P_R_map);

        unsigned int rank = linearSolver.rank();

        if (getRankDeficientError() && (rank != getResidual()->size())) {
            TARDIGRADE_ERROR_TOOLS_CATCH(throw convergence_error("The Jacobian is not full rank"));
        }

        addTrialStepOutput();
    }

}  // namespace tardigradeHydra
