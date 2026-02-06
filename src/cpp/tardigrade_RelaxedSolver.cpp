/**
 ******************************************************************************
 * \file tardigrade_RelaxedSolver.cpp
 ******************************************************************************
 * A C++ library for the nonlinear solvers which attempt to relax the problem
 * during its solution
 ******************************************************************************
 */

#include "tardigrade_RelaxedSolver.h"

#include "tardigrade_hydra.h"

namespace tardigradeHydra {

    /*!
     * Default constructor for RelaxedSolver
     */
    RelaxedSolver::RelaxedSolver() : IterativeSolverBase() {
        internal_solver->hydra = NULL;
        step                   = internal_solver->step;
    }

    /*!
     * Constructor for RelaxedSolver
     *
     * \param *_hydra: The containing hydra object
     */
    RelaxedSolver::RelaxedSolver(hydraBase *_hydra) : IterativeSolverBase(_hydra) {
        internal_solver->hydra = _hydra;
        step                   = internal_solver->step;
    }

    /*!
     * Constructor for RelaxedSolver
     *
     * \param *_hydra: The containing hydra object
     * \param *_internal_solver_ptr: The pointer for the internal solver
     */
    RelaxedSolver::RelaxedSolver(hydraBase *_hydra, SolverBase *_internal_solver_ptr) : IterativeSolverBase(_hydra) {
        internal_solver        = _internal_solver_ptr;
        internal_solver->hydra = hydra;
        step                   = internal_solver->step;
    }

    /*!
     * Reset the solver
     */
    void RelaxedSolver::reset() {
        resetRelaxedIteration();
        internal_solver->reset();
        tardigradeHydra::IterativeSolverBase::reset();
    }

    /*!
     * Get the current relaxed iteration
     */
    const unsigned int RelaxedSolver::getRelaxedIteration() { return _relaxedIteration; }

    /*!
     * Get the maximum number of relaxed iterations
     */
    const unsigned int RelaxedSolver::getMaxRelaxedIterations() { return _maxRelaxedIterations; }

    /*!
     * Set the maximum allowable number of relaxed iterations
     *
     * \param &value: The number of relaxed iterations
     */
    const void RelaxedSolver::setMaxRelaxedIterations(const unsigned int &value) { _maxRelaxedIterations = value; }

    /*!
     * Set the relaxed iteration number
     *
     * \param &value: The incoming value
     */
    void RelaxedSolver::setRelaxedIteration(const unsigned int &value) { _relaxedIteration = value; }

    /*!
     * Reset the relaxed iteration number
     */
    void RelaxedSolver::resetRelaxedIteration() { setRelaxedIteration(0); }

    /*!
     * Increment the relaxed iteration number
     */
    void RelaxedSolver::incrementRelaxedIteration() { _relaxedIteration++; }

    /*!
     * Initialize the residuals for a relaxed solve
     */
    void RelaxedSolver::initializeResiduals() {
        setCurrentResidualIndexMeaningful(true);

        for (auto residual = std::begin(*getResidualClasses()); residual != std::end(*getResidualClasses());
             ++residual) {
            setCurrentResidualIndex(residual - std::begin(*getResidualClasses()));

            // Prepare the residuals to take a relaxed step
            (*residual)->setupRelaxedStep(getRelaxedIteration());
        }

        setCurrentResidualIndexMeaningful(false);
    }

    /*!
     * Check if the relaxation iterations have converged
     */
    bool RelaxedSolver::checkRelaxedConvergence() {
        bool relaxedConverged = true;

        setCurrentResidualIndexMeaningful(true);

        for (auto residual = std::begin(*getResidualClasses()); residual != std::end(*getResidualClasses());
             ++residual) {
            setCurrentResidualIndex(residual - std::begin(*getResidualClasses()));

            if (!(*residual)->checkRelaxedConvergence()) {
                relaxedConverged = false;
                break;
            }
        }

        setCurrentResidualIndexMeaningful(false);

        return relaxedConverged;
    }

    /*!
     * Signal to the residuals that we have a failed relaxed solve step and
     * determine if a new relaxed step should be taken
     */
    bool RelaxedSolver::callResidualRelaxedStepFailure() {
        bool attempt_relaxed_step = false;

        setCurrentResidualIndexMeaningful(true);

        for (auto residual_ptr = std::begin(*getResidualClasses()); residual_ptr != std::end(*getResidualClasses());
             ++residual_ptr) {
            setCurrentResidualIndex(residual_ptr - std::begin(*getResidualClasses()));

            try {
                auto val = (*residual_ptr)->relaxedStepFailure();

                attempt_relaxed_step = attempt_relaxed_step || val;

            } catch (std::exception &e) {
                if (getFailureVerbosityLevel() > 0) {
                    addToFailureOutput("Failure in residual " +
                                       std::to_string(residual_ptr - std::begin(*getResidualClasses())) + "\n");
                    std::string message;
                    tardigradeErrorTools::captureNestedExceptions(e, message);
                    addToFailureOutput(message);
                }

                throw;
            }
        }

        setCurrentResidualIndexMeaningful(false);

        return attempt_relaxed_step;
    }

    /*!
     * Attempt to perform a solve of the non-linear problem
     */
    bool RelaxedSolver::attemptInternalSolve() {
        TARDIGRADE_ERROR_TOOLS_CHECK(
            internal_solver != nullptr,
            "The solver which is to be relaxed (i.e., the internal solver) has not been defined");

        try {
            internal_solver->solve();

            // Exit if the relaxed solver has converged
            return checkRelaxedConvergence();

        } catch (convergence_error &e) {
            if (!callResidualRelaxedStepFailure()) {
                throw;
            }

            return false;

        } catch (std::exception &e) {
            throw;
        }
    }

    /*!
     * Setup the next relaxed step
     */
    void RelaxedSolver::setupNextRelaxedStep() {
        // Use the current unknown vector as the initial estimate
        setInitializeUnknownVector(false);

        incrementRelaxedIteration();

        // Re-initialize the residuals
        initializeResiduals();

        // Reset hydra
        updateUnknownVector(
            *getUnknownVector());  // This allows for the relaxed to change the projection and adjust the decomposition
        resetIterations();
    }

    /*!
     * Write the relaxed iteration to the failure string
     */
    void RelaxedSolver::logRelaxedIterationHeader() {
        if (getFailureVerbosityLevel() > 0) {
            addToFailureOutput("\n\n###  relaxed iteration: ");
            addToFailureOutput(getRelaxedIteration());
            addToFailureOutput("\n\n");
        }
    }

    /*!
     * The function that is called for the first solve attempt
     */
    void RelaxedSolver::initialSolveAttempt() {
        TARDIGRADE_ERROR_TOOLS_CHECK(internal_solver != nullptr, "The internal solver points to a null pointer")

        internal_solver->solve();
    }

    /*!
     * The function that is called when there is a convergence error
     * thrown by the initial solve
     */
    void RelaxedSolver::convergenceErrorFunction() {
        try {
            performRelaxedSolve();

        } catch (const convergence_error &e) {
            throw;

        } catch (std::exception *e) {
            TARDIGRADE_ERROR_TOOLS_CATCH(throw;)
        }
    }

    /*!
     * Solve the non-linear problem by relaxing difficult sub-problems
     * to achieve a series of solutions.
     */
    void RelaxedSolver::performRelaxedSolve() {
        TARDIGRADE_ERROR_TOOLS_CHECK(
            internal_solver != nullptr,
            "The solver which is to be relaxed (i.e., the internal solver) has not been defined");

        if (getFailureVerbosityLevel() > 0) {
            addToFailureOutput("Failure in conventional solve. Starting relaxed solve.\n");
        }

        initial_unknown = internal_solver->initial_unknown;
        TARDIGRADE_ERROR_TOOLS_CATCH(
            internal_solver->resetCounts();)  // TODO: Maybe this should be a full reset? It causes errors in the tests
                                              // if it is but I'm of two minds
        TARDIGRADE_ERROR_TOOLS_CATCH(updateUnknownVector(initial_unknown);)

        resetRelaxedIteration();

        // Initialize the residuals
        initializeResiduals();

        while (getRelaxedIteration() < getMaxRelaxedIterations()) {
            logRelaxedIterationHeader();

            // Solve the non-linear problem
            if (attemptInternalSolve()) {
                return;
            }

            setupNextRelaxedStep();
        }

        if (getRelaxedIteration() >= getMaxRelaxedIterations()) {
            throw convergence_error("Failure in relaxed solve\n");
        }
    }

}  // namespace tardigradeHydra
