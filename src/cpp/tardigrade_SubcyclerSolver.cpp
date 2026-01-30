/**
 ******************************************************************************
 * \file tardigrade_SubcyclerSolver.cpp
 ******************************************************************************
 * A C++ library for the nonlinear solvers which attempt to relax the problem
 * during its solution
 ******************************************************************************
 */

#include "tardigrade_SubcyclerSolver.h"

#include "tardigrade_hydra.h"

namespace tardigradeHydra {

    /*!
     * Default constructor for SubcyclerSolver
     */
    SubcyclerSolver::SubcyclerSolver() : RelaxedSolver() {//IterativeSolverBase() {
        internal_solver->hydra = NULL;
        step                   = internal_solver->step;
    }

    /*!
     * Constructor for SubcyclerSolver
     *
     * \param *_hydra: The containing hydra object
     */
    SubcyclerSolver::SubcyclerSolver(hydraBase *_hydra) : RelaxedSolver(_hydra) {//IterativeSolverBase(_hydra) {
        internal_solver->hydra = _hydra;
        step                   = internal_solver->step;
    }

    /*!
     * Constructor for SubcyclerSolver
     *
     * \param *_hydra: The containing hydra object
     * \param *_internal_solver_ptr: The pointer for the internal solver
     */
    SubcyclerSolver::SubcyclerSolver(hydraBase *_hydra, SolverBase *_internal_solver_ptr) : RelaxedSolver(_hydra) {//IterativeSolverBase(_hydra) {
        internal_solver        = _internal_solver_ptr;
        internal_solver->hydra = hydra;
        step                   = internal_solver->step;
    }

//    /*!
//     * The initial attempt at solving the problem
//     *
//     * Attempts to use the internal solver directly without subcycling
//     */
//    void SubcyclerSolver::initialSolveAttempt() {
//        TARDIGRADE_ERROR_TOOLS_CHECK( internal_solver != nullptr, "The solver hasn't been defined" );
//        internal_solver->solve();
//    }
//
//    /*!
//     * The function that is called if there is a convergence
//     * error thrown in the initial solve attempt
//     *
//     * Attempts a subcycler step
//     */
//    void SubcyclerSolver::convergenceErrorFunction() {
//        performSubcyclerSolve();
//    }
//
//    /*!
//     * The function that is called if there is a unexpected
//     * error thrown in the initial solve attempt
//     *
//     * Attempts a subcycler step
//     */
//    void SubcyclerSolver::unexpectedErrorFunction() {
//        performSubcyclerSolve();
//    }
//
//    /*!
//     * Reset the solver
//     */
//
//    void RelaxedSolver::reset(){
//        internal_solver->reset();
//        tardigradeHydra::IterativeSolverBase::reset();
//    }
//
//    /*!
//     * Initialize the subcycler
//     */
//    void SubcyclerSolver::initializeSubcycler() {
//        sp = 0.0;
//
//        ds = getCutbackFactor();
//
//        num_good = 0;
//
//        callResidualPreSubcycler();
//
//        resetProblem();
//    }
//
//    /*!
//     * Update the pseudo-timestep
//     */
//    void SubcyclerSolver::updatePseudoTimestep() {
//        sp += ds;  // Update the pseudo-time
//
//        // Grow the step if possible
//        if (allowStepGrowth(num_good)) {
//            ds *= getGrowthFactor();
//        }
//
//        // Make sure s will be less than or equal to 1
//        if (sp + ds > 1.0) {
//            ds = 1.0 - sp;
//        }
//    }
//
//    /*!
//     * Perform a subcycler step
//     */
//    void SubcyclerSolver::performSubcyclerStep() {
//        TARDIGRADE_ERROR_TOOLS_CHECK(internal_solver != nullptr, "The internal solver has not been defined");
//        addSubcyclerStepHeader();
//
//        setScaleFactor(sp + ds);  // Update the scaling factor
//
//        resetIterations();  // Reset the non-linear iteration count
//
//        internal_solver->solve();  // Try to solve the non-linear problem
//    }
//
//    /*!
//     * Post-successful subcycler increment updates
//     */
//    void SubcyclerSolver::subcyclerStepSuccess() {
//        setPreviouslyConvergedStress(*getStress());  // Set the previously converged stress
//
//        callResidualPostSubcyclerSuccess();  // Let the residuals know the subcycle step was successful
//
//        num_good++;  // Update the number of good iterations
//
//        updatePseudoTimestep();
//    }
//
//    /*!
//     * Called when there is a failure in a subcycler step
//     */
//    void SubcyclerSolver::subcyclerStepFailure() {
//        callResidualPostSubcyclerFailure();
//
//        // Reduce the time-step and try again
//        num_good = 0;
//
//        ds *= getCutbackFactor();
//
//        setX(initial_unknown);  // Reset X to the last good point
//
//        if (ds < getMinDS()) {
//            throw;
//        }
//    }
//
//    /*!
//     * Solve the problem using the subcycler
//     */
//    void SubcyclerSolver::performSubcyclerSolve() {
//        initializeSubcycler();
//
//        while (sp < 1.0) {
//            try {
//                performSubcyclerStep();
//
//                subcyclerStepSuccess();
//
//            } catch (std::exception &e) {
//                subcyclerStepFailure();
//            }
//        }
//    }
//
//    /*!
//     * Add the header for the subcycler to the output failure string
//     */
//    void SubcyclerSolver::addSubcyclerHeader() {
//        if (getFailureVerbosityLevel() > 0) {
//            addToFailureOutput("\n\n");
//            addToFailureOutput("#########################################\n");
//            addToFailureOutput("###        ENTERING SUB-CYCLER        ###\n");
//            addToFailureOutput("#########################################\n");
//            addToFailureOutput("\n\n");
//        }
//    }
//
//    /*!
//     * Add the subcycler step header to the output failure string
//     */
//    void SubcyclerSolver::addSubcyclerStepHeader() {
//        if (getFailureVerbosityLevel() > 0) {
//            addToFailureOutput("\n\n");
//            addToFailureOutput("######### PSEUDO-TIME INCREMENT #########\n");
//            addToFailureOutput("\n\n    sp, ds: " + std::to_string(sp) + ", " + std::to_string(ds));
//            addToFailureOutput("\n");
//        }
//    }
//
//    /*!
//     * Function to determine if we can increase the step-size for the sub-cycler
//     *
//     * \param &num_good: The number of good increments since the last failure
//     */
//    const bool SubcyclerSolver::allowStepGrowth(const unsigned int &num_good) {
//        if (num_good >= getNumGoodControl()) {
//            return true;
//        }
//
//        return false;
//    }
//
//    /*!
//     * Signal to the residuals that we are entering the subcycler
//     */
//    void SubcyclerSolver::callResidualPreSubcycler() {
//        setCurrentResidualIndexMeaningful(true);
//
//        for (auto residual_ptr = getResidualClasses()->begin(); residual_ptr != getResidualClasses()->end();
//             ++residual_ptr) {
//            setCurrentResidualIndex(residual_ptr - getResidualClasses()->begin());
//
//            try {
//                (*residual_ptr)->preSubcycler();
//
//            } catch (std::exception &e) {
//                if (getFailureVerbosityLevel() > 0) {
//                    addToFailureOutput("Failure in residual " +
//                                       std::to_string(residual_ptr - getResidualClasses()->begin()) + "\n");
//                    std::string message;
//                    tardigradeErrorTools::captureNestedExceptions(e, message);
//                    addToFailureOutput(message);
//                }
//
//                throw;
//            }
//        }
//
//        setCurrentResidualIndexMeaningful(false);
//    }

    /*!
     * Signal to the residuals that we have a successful subcycle increment
     */
    void SubcyclerSolver::callResidualPostSubcyclerSuccess() {
        setCurrentResidualIndexMeaningful(true);

        for (auto residual_ptr = getResidualClasses()->begin(); residual_ptr != getResidualClasses()->end();
             ++residual_ptr) {
            setCurrentResidualIndex(residual_ptr - getResidualClasses()->begin());

            try {
                (*residual_ptr)->postSubcyclerSuccess();

            } catch (std::exception &e) {
                if (getFailureVerbosityLevel() > 0) {
                    addToFailureOutput("Failure in residual " +
                                       std::to_string(residual_ptr - getResidualClasses()->begin()) + "\n");
                    std::string message;
                    tardigradeErrorTools::captureNestedExceptions(e, message);
                    addToFailureOutput(message);
                }

                throw;
            }
        }

        setCurrentResidualIndexMeaningful(false);
    }

    /*!
     * Signal to the residuals that we have a failed subcycle increment
     */
    void SubcyclerSolver::callResidualPostSubcyclerFailure() {
        setCurrentResidualIndexMeaningful(true);

        for (auto residual_ptr = getResidualClasses()->begin(); residual_ptr != getResidualClasses()->end();
             ++residual_ptr) {
            setCurrentResidualIndex(residual_ptr - getResidualClasses()->begin());

            try {
                (*residual_ptr)->postSubcyclerFailure();

            } catch (std::exception &e) {
                if (getFailureVerbosityLevel() > 0) {
                    addToFailureOutput("Failure in residual " +
                                       std::to_string(residual_ptr - getResidualClasses()->begin()) + "\n");
                    std::string message;
                    tardigradeErrorTools::captureNestedExceptions(e, message);
                    addToFailureOutput(message);
                }

                throw;
            }
        }

        setCurrentResidualIndexMeaningful(false);
    }

}
