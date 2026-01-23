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
     * Solve the constrained QP problem to estimate the desired step size
     *
     * \param &dx: The change in the unknown vector
     * \param kmax: The maximum number of iterations (defaults to 100)
     */
    void NonlinearStepBase::solveConstrainedQP(floatVector &dx, const unsigned int kmax) {
        const unsigned int numUnknowns = getNumUnknowns();

        const unsigned int numConstraints = getNumConstraints();

        floatVector K;

        floatVector RHS;

        std::vector<bool> active_constraints;
        initializeActiveConstraints(active_constraints);

        assembleKKTRHSVector(dx, RHS, active_constraints);

        assembleKKTMatrix(K, active_constraints);

        floatType tol = getRelativeTolerance() * (tardigradeVectorTools::l2norm(RHS)) + getAbsoluteTolerance();

        unsigned int k = 0;

        floatVector y(numUnknowns + numConstraints, 0);

        floatVector ck = *getConstraints();

        floatVector ctilde(numConstraints, 0);

        floatVector negp(numUnknowns, 0);

        floatVector lambda(numConstraints, 0);

        floatVector P(numUnknowns + numConstraints, 0);

        for (unsigned int i = 0; i < (numUnknowns + numConstraints); i++) {
            P[i] = 1 / std::max(std::fabs(*std::max_element(K.begin() + (numUnknowns + numConstraints) * i,
                                                            K.begin() + (numUnknowns + numConstraints) * (i + 1),
                                                            [](const floatType &a, const floatType &b) {
                                                                return std::fabs(a) < std::fabs(b);
                                                            })),
                                1e-15);
        }

        Eigen::Map<const Eigen::Vector<floatType, -1> > _P(P.data(), numUnknowns + numConstraints);

        tardigradeVectorTools::solverType<floatType> linearSolver;

        while (k < kmax) {
            Eigen::Map<const Eigen::Matrix<floatType, -1, -1, Eigen::RowMajor> > _K(K.data(),
                                                                                    numConstraints + numUnknowns,
                                                                                    numConstraints + numUnknowns);
            Eigen::Map<const Eigen::Vector<floatType, -1> > _RHS(RHS.data(), numConstraints + numUnknowns);
            Eigen::Map<Eigen::Vector<floatType, -1> >       _y(y.data(), numConstraints + numUnknowns);

            linearSolver = tardigradeVectorTools::solverType<floatType>(_P.asDiagonal() * _K);

            _y = linearSolver.solve(_P.asDiagonal() * _RHS);

            std::copy(y.begin(), y.begin() + numUnknowns, negp.begin());

            std::copy(y.begin() + numUnknowns, y.end(), lambda.begin());

            if (tardigradeVectorTools::l2norm(negp) <= tol) {
                bool negLambda = false;

                floatType minLambda = 1;

                unsigned int imin = 0;

                for (auto v = std::begin(lambda); v != std::end(lambda); v++) {
                    if (*v < 0) {
                        negLambda = true;

                        if ((*v) < minLambda) {
                            imin = (unsigned int)(v - std::begin(lambda));

                            minLambda = *v;
                        }
                    }
                }

                if (negLambda) {
                    active_constraints[imin] = false;

                } else {
                    return;
                }

            } else {
                ck     = *getConstraints();
                ctilde = *getConstraints();
                for (unsigned int i = 0; i < numConstraints; i++) {
                    for (unsigned int j = 0; j < numUnknowns; j++) {
                        ck[i] += (*getConstraintJacobians())[numUnknowns * i + j] * dx[j];
                        ctilde[i] += (*getConstraintJacobians())[numUnknowns * i + j] * (dx[j] - negp[j]);
                    }
                }

                floatType alpha = 1.0;

                unsigned int iblock = 0;

                bool newBlock = false;

                for (unsigned int i = 0; i < numConstraints; i++) {
                    if (!active_constraints[i]) {
                        if (ctilde[i] < -tol) {
                            floatType alpha_trial = -ck[i] / (ctilde[i] - ck[i]);

                            if (alpha_trial <= alpha) {
                                iblock = i;

                                alpha = alpha_trial;

                                newBlock = true;
                            }
                        }
                    }
                }

                if (newBlock) {
                    active_constraints[iblock] = true;
                }

                dx -= alpha * negp;
            }

            updateKKTMatrix(K, active_constraints);

            assembleKKTRHSVector(dx, RHS, active_constraints);

            for (unsigned int i = 0; i < (numUnknowns + numConstraints); i++) {
                P[i] = 1 / std::max(std::fabs(*std::max_element(K.begin() + (numUnknowns + numConstraints) * i,
                                                                K.begin() + (numUnknowns + numConstraints) * (i + 1),
                                                                [](const floatType &a, const floatType &b) {
                                                                    return std::fabs(a) < std::fabs(b);
                                                                })),
                                    1e-15);
            }

            k++;
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

    // END SQP SOLVER FUNCTIONS

}
