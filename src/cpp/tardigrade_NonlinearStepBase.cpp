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

    // BEGIN SQP SOLVER FUNCTIONS

    /*!
     * Assemble the right hand side vector for the KKT matrix
     *
     * \param &dx: The delta vector being solved for
     * \param &KKTRHSVector: The right hand size vector for the KKT matrix
     * \param &active_constraints: The active constraint vector
     */
    void NonlinearStepBase::assembleKKTRHSVector(const floatVector &dx, floatVector &KKTRHSVector,
                                             const std::vector<bool> &active_constraints) {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        const unsigned int numUnknowns = getNumUnknowns();

        const unsigned int numConstraints = getNumConstraints();

        KKTRHSVector = floatVector(numUnknowns + numConstraints, 0);

        Eigen::Map<const Eigen::Vector<floatType, -1> > _dx(dx.data(), numUnknowns);

        Eigen::Map<Eigen::Vector<floatType, -1> > RHS(KKTRHSVector.data(), (numUnknowns + numConstraints),
                                                      (numUnknowns + numConstraints));

        Eigen::Map<const Eigen::Vector<floatType, -1> > R(getResidual()->data(), numUnknowns);

        Eigen::Map<const Eigen::Matrix<floatType, -1, -1, Eigen::RowMajor> > J(getFlatJacobian()->data(), numUnknowns,
                                                                               numUnknowns);

        RHS.head(numUnknowns) = (J.transpose() * (R + J * _dx) + step->damping->getMuk() * _dx).eval();

        for (unsigned int i = 0; i < numConstraints; i++) {
            if (active_constraints[i]) {
                KKTRHSVector[numUnknowns + i] = (*(getConstraints()))[i];

                for (unsigned int I = 0; I < numUnknowns; I++) {
                    KKTRHSVector[numUnknowns + i] += (*(getConstraintJacobians()))[numUnknowns * i + I] * dx[I];
                }
            }
        }
    }

    /*!
     * Assemble the Karush-Kuhn-Tucker matrix for an inequality constrained Newton-Raphson solve
     *
     * \param &KKTMatrix: The Karush-Kuhn-Tucker matrix
     * \param &active_constraints: The vector of currently active constraints.
     */
    void NonlinearStepBase::assembleKKTMatrix(floatVector &KKTMatrix, const std::vector<bool> &active_constraints) {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        const unsigned int numUnknowns = getNumUnknowns();

        const unsigned int numConstraints = getNumConstraints();

        KKTMatrix = floatVector((numUnknowns + numConstraints) * (numUnknowns + numConstraints), 0);

        Eigen::Map<Eigen::Matrix<floatType, -1, -1, Eigen::RowMajor> > K(KKTMatrix.data(),
                                                                         (numUnknowns + numConstraints),
                                                                         (numUnknowns + numConstraints));

        Eigen::Map<const Eigen::Matrix<floatType, -1, -1, Eigen::RowMajor> > J(getFlatJacobian()->data(), numUnknowns,
                                                                               numUnknowns);

        K.block(0, 0, numUnknowns, numUnknowns) = (J.transpose() * J).eval();

        for (unsigned int I = 0; I < numUnknowns; I++) {
            KKTMatrix[(numUnknowns + numConstraints) * I + I] += step->damping->getMuk();
        }

        for (unsigned int i = 0; i < numConstraints; i++) {
            if (active_constraints[i]) {
                for (unsigned int I = 0; I < numUnknowns; I++) {
                    KKTMatrix[(numUnknowns + numConstraints) * (I) + numUnknowns + i] =
                        (*getConstraintJacobians())[numUnknowns * i + I];
                    KKTMatrix[(numUnknowns + numConstraints) * (numUnknowns + i) + I] =
                        (*getConstraintJacobians())[numUnknowns * i + I];
                }

            } else {
                KKTMatrix[(numUnknowns + numConstraints) * (numUnknowns + i) + numUnknowns + i] = 1;
            }
        }
    }

    /*!
     * Update the KKTMatrix if the active constraints have changed
     *
     * \param &KKTMatrix: The Karush-Kuhn-Tucker matrix
     * \param &active_constraints: The vector of currently active constraints.
     */
    void NonlinearStepBase::updateKKTMatrix(floatVector &KKTMatrix, const std::vector<bool> &active_constraints) {
        const unsigned int numUnknowns = getNumUnknowns();

        const unsigned int numConstraints = getNumConstraints();

        Eigen::Map<Eigen::Matrix<floatType, -1, -1, Eigen::RowMajor> > K(KKTMatrix.data(),
                                                                         (numUnknowns + numConstraints),
                                                                         (numUnknowns + numConstraints));

        K.block(0, numUnknowns, numUnknowns, numConstraints).setZero();

        K.block(numUnknowns, 0, numConstraints, numUnknowns).setZero();

        K.block(numUnknowns, numUnknowns, numConstraints, numConstraints).setZero();

        for (unsigned int i = 0; i < numConstraints; i++) {
            if (active_constraints[i]) {
                for (unsigned int I = 0; I < numUnknowns; I++) {
                    KKTMatrix[(numUnknowns + numConstraints) * (I) + numUnknowns + i] =
                        (*getConstraintJacobians())[numUnknowns * i + I];
                    KKTMatrix[(numUnknowns + numConstraints) * (numUnknowns + i) + I] =
                        (*getConstraintJacobians())[numUnknowns * i + I];
                }

            } else {
                KKTMatrix[(numUnknowns + numConstraints) * (numUnknowns + i) + numUnknowns + i] = 1;
            }
        }
    }

    /*!
     * Initialize the active constraint vector
     *
     * \param &active_constraints: The current constraints that are active
     */
    void NonlinearStepBase::initializeActiveConstraints(std::vector<bool> &active_constraints) {
        active_constraints = std::vector<bool>(getNumConstraints(), false);

        for (auto c = getConstraints()->begin(); c != getConstraints()->end(); c++) {
            unsigned int index = (unsigned int)(c - getConstraints()->begin());

            active_constraints[index] = ((*c) < 0.);
        }
    }

    // END SQP SOLVER FUNCTIONS

}
