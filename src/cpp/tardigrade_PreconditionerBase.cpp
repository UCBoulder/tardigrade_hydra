
/**
 ******************************************************************************
 * \file tardigrade_PreconditionerBase.cpp
 ******************************************************************************
 * A C++ library for the base classes for preconditioners
 ******************************************************************************
 */

#include "tardigrade_PreconditionerBase.h"
#include "tardigrade_NonlinearStepBase.h"

namespace tardigradeHydra {

    /*!
     * Reset the preconditioner back to an initial state
     */
    void PreconditionerBase::reset() {}

    /*!
     * Add data to the vector of values which will be cleared after each iteration
     *
     * \param *data: The dataBase object to be cleared
     */
    void PreconditionerBase::addIterationData(dataBase *data) {
        TARDIGRADE_ERROR_TOOLS_CHECK(trial_step != nullptr, "The trial step has not been defined");
        trial_step->addIterationData(data);
    }

    /*!
     * Add data to the vector of values which will be cleared after each nonlinear step
     *
     * \param *data: The dataBase object to be cleared
     */
    void PreconditionerBase::addNLStepData(dataBase *data) {
        TARDIGRADE_ERROR_TOOLS_CHECK(trial_step != nullptr, "The trial step has not been defined");
        trial_step->addNLStepData(data);
    }

    /*!
     * Get the number of unknowns
     */
    const unsigned int PreconditionerBase::getNumUnknowns(){
        TARDIGRADE_ERROR_TOOLS_CHECK(trial_step != nullptr, "The trial step has not been defined");
        return trial_step->getNumUnknowns();
    }

    /*!
     * Get the row-major nonlinear left hand side
     */
    const floatVector* PreconditionerBase::getFlatNonlinearLHS(){
        TARDIGRADE_ERROR_TOOLS_CHECK(trial_step != nullptr, "The trial step has not been defined");
        return trial_step->getFlatNonlinearLHS();
    }

    /*!
     * Get the nonlinear right hand side
     */
    const floatVector* PreconditionerBase::getNonlinearRHS(){
        TARDIGRADE_ERROR_TOOLS_CHECK(trial_step != nullptr, "The trial step has not been defined");
        return trial_step->getNonlinearRHS();
    }

    /*!
     * Get the flattened row-major preconditioner for the non-linear problem
     */
    const floatVector *PreconditionerBase::getFlatPreconditioner() {
        if (!_preconditioner.first) {
            TARDIGRADE_ERROR_TOOLS_CATCH(formPreconditioner());
        }

        return &_preconditioner.second;
    }

    /*!
     * Form the preconditioner matrix
     */
    void PreconditionerBase::formPreconditioner() {
        _preconditioner.second = floatVector(trial_step->getNumUnknowns(),1.);

        _preconditioner.first = true;

        addIterationData(&_preconditioner);
    }

    /*!
     * Form a left preconditioner comprised of the inverse of the maximum value of each row
     */
    void PreconditionerBase::formMaxRowPreconditioner() {
        const unsigned int problem_size = trial_step->getNumUnknowns();

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
     * Precondition the incoming vector \f$X\f$ via \f$Y_I = P_{IJ} X_J\f$
     *
     * \param &X: The incoming vector to be preconditioned
     * \param &Y: The preconditioned vector
     */
    void PreconditionerBase::preconditionVector( const floatVector &X, floatVector &Y ){

        Y = floatVector(X.size(),0);
        std::copy( std::begin( X ), std::end( X ), std::begin( Y ) );

    }

    /*!
     * Precondition the incoming matrix \f$A\f$ via \f$B_{IJ} = P_{IK} A_{KJ}\f$
     *
     * \param &A: The incoming matrix to be preconditioned
     * \param &B: The preconditined matrix
     */
    void PreconditionerBase::preconditionMatrix( const floatVector &A, floatVector &B ){

        B = floatVector(A.size(),0);
        std::copy( std::begin( A ), std::end( A ), std::begin( B ) );

    }
}  // namespace tardigradeHydra
