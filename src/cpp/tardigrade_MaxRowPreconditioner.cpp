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

    /*!
     * Precondition the incoming vector \f$X\f$ via \f$Y_I = P_{IJ} X_J\f$
     *
     * \param &X: The incoming vector to be preconditioned
     * \param &Y: The preconditioned vector
     */
    void MaxRowPreconditioner::preconditionVector( const floatVector &X, floatVector &Y ){

        Y = floatVector(X.size(),0);

        auto X_map =
            tardigradeHydra::getDynamicSizeVectorMap(X.data(), X.size());

        auto Y_map =
            tardigradeHydra::getDynamicSizeVectorMap(Y.data(), Y.size());

        auto p_map = tardigradeHydra::getDynamicSizeVectorMap(getFlatPreconditioner()->data(),
                                                              getFlatPreconditioner()->size()); //Current preconditioner is flat

        Y_map = p_map.asDiagonal( ) * X_map;

    }

    /*!
     * Precondition the incoming matrix \f$A\f$ via \f$B_{IJ} = P_{IK} A_{KJ}\f$
     *
     * \param &A: The incoming matrix to be preconditioned
     * \param &B: The preconditined matrix
     */
    void MaxRowPreconditioner::preconditionMatrix( const floatVector &A, floatVector &B ){

        auto cols = A.size( ) / getNumUnknowns( );

        B = floatVector(A.size(),0);

        auto A_map =
            tardigradeHydra::getDynamicSizeMatrixMap(A.data(), getNumUnknowns(), cols);

        auto B_map =
            tardigradeHydra::getDynamicSizeMatrixMap(B.data(), getNumUnknowns(), cols);

        auto p_map = tardigradeHydra::getDynamicSizeVectorMap(getFlatPreconditioner()->data(),
                                                              getFlatPreconditioner()->size()); //Current preconditioner is flat

        B_map = p_map.asDiagonal( ) * A_map;

    }
}
