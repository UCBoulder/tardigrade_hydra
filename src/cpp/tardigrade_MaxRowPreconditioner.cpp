/**
 ******************************************************************************
 * \file tardigrade_MaxRowPreconditioner.cpp
 ******************************************************************************
 * Form a preconditioner based on the maximum value in the LHS matrix's rows
 ******************************************************************************
 */

#include "tardigrade_MaxRowPreconditioner.h"

namespace tardigradeHydra {

    /*!
     * Form a left preconditioner comprised of the inverse of the maximum value of each row
     */
    void MaxRowPreconditioner::formMaxRowPreconditioner() {
        TARDIGRADE_ERROR_TOOLS_CHECK(trial_step != nullptr, "The trial step has not been defined");
        const unsigned int problem_size = getNumUnknowns();

        _preconditioner.second = floatVector(problem_size, 0);

        // Find the absolute maximum value in each row
        for (unsigned int i = 0; i < problem_size; i++) {
            _preconditioner.second[i] =
                1 / std::max(std::fabs(*std::max_element(
                                 getFlatNonlinearLHS()->begin() + problem_size * i,
                                 getFlatNonlinearLHS()->begin() + problem_size * (i + 1),
                                 [](const floatType &a, const floatType &b) { return std::fabs(a) < std::fabs(b); })),
                             1e-15);
        }
    }

    /*!
     * Form the preconditioner matrix
     */
    void MaxRowPreconditioner::formPreconditioner() {
        formMaxRowPreconditioner();

        _preconditioner.first = true;

        addIterationData(&_preconditioner);
    }

}
