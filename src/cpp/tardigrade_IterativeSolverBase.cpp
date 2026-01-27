/**
 ******************************************************************************
 * \file tardigrade_IterativeSolverBase.cpp
 ******************************************************************************
 * The base class for iterative solver classes
 ******************************************************************************
 */

#include "tardigrade_IterativeSolverBase.h"

#include "tardigrade_hydra.h"

namespace tardigradeHydra {

    /*!
     * reset the solver
     */
    void IterativeSolverBase::reset() {
        resetIterations();
        tardigradeHydra::SolverBase::reset();
    }

    /*!
     * The function that is called when first attempting to
     * solve the problem
     */
    void IterativeSolverBase::initialSolveAttempt() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");

        // Reset the internal steps
        step->reset();

        setRankDeficientError(false);

        // Form the initial unknown vector
        if (getInitializeUnknownVector()) {
            TARDIGRADE_ERROR_TOOLS_CATCH(initializeUnknownVector());
        }

        initial_unknown = *getUnknownVector();

        floatVector deltaX(getNumUnknowns(), 0);

        callResidualPreIterativeSolve();

        step->damping->reset();

        if (getFailureVerbosityLevel() > 0) {
            addToFailureOutput("Initial Unknown:\n");
            addToFailureOutput(*getUnknownVector());
        }

        while (!checkConvergence() && checkIteration()) {

            addIterationHeader();

            step->incrementSolution();

            // Call residual end of a successful nonlinear step functions
            callResidualSuccessfulIterativeStep();

            // Increment the iteration count
            incrementIteration();

            // Reset the nonlinear step data
            resetNLStepData();

            addIterationFooter();
        }

        if (!checkConvergence()) {
            throw convergence_error("Failure to converge main loop\n");
        }

        callResidualPostIterativeSolve();
    }

    /*!
     * Add the iteration header to the failure output message
     */
    void IterativeSolverBase::addIterationHeader(){

        if (getFailureVerbosityLevel() > 0) {
            addToFailureOutput("\n\n  iteration: ");
            addToFailureOutput(getIteration());
        }

    }

    /*!
     * Add the iteration footer to the failure output message
     */
    void IterativeSolverBase::addIterationFooter(){

        if (getFailureVerbosityLevel() > 0) {
            addToFailureOutput("  final residual: ");
            addToFailureOutput(tardigradeVectorTools::l2norm(*getResidual()));
            addToFailureOutput("\n");
        }
    }

    /*!
     * Signal to the residuals that we are about to start an iterative solve
     */
    void IterativeSolverBase::callResidualPreIterativeSolve() {
        setCurrentResidualIndexMeaningful(true);

        for (auto residual_ptr = std::begin(*getResidualClasses()); residual_ptr != std::end(*getResidualClasses());
             ++residual_ptr) {
            setCurrentResidualIndex(residual_ptr - std::begin(*getResidualClasses()));

            (*residual_ptr)->preIterativeSolve();
        }

        setCurrentResidualIndexMeaningful(false);
    }

    /*!
     * Signal to the residuals that a successful iterative step has been performed
     */
    void IterativeSolverBase::callResidualSuccessfulIterativeStep() {
        setAllowModifyGlobalResidual(true);

        setCurrentResidualIndexMeaningful(true);

        for (auto residual_ptr = std::begin(*getResidualClasses()); residual_ptr != std::end(*getResidualClasses());
             ++residual_ptr) {
            setCurrentResidualIndex(residual_ptr - std::begin(*getResidualClasses()));

            (*residual_ptr)->successfulIterativeStep();
        }

        setCurrentResidualIndexMeaningful(false);

        setAllowModifyGlobalResidual(false);
    }

    /*!
     * Signal to the residuals that we have finished an iterative solve
     */
    void IterativeSolverBase::callResidualPostIterativeSolve() {
        setCurrentResidualIndexMeaningful(true);

        for (auto residual_ptr = std::begin(*getResidualClasses()); residual_ptr != std::end(*getResidualClasses());
             ++residual_ptr) {
            setCurrentResidualIndex(residual_ptr - std::begin(*getResidualClasses()));

            try {
                (*residual_ptr)->postIterativeSolve();

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
    }

    /*!
     * Check the convergence
     */
    bool IterativeSolverBase::checkConvergence() {
        const floatVector *tolerance = getTolerance();

        const floatVector *residual = getResidual();

        TARDIGRADE_ERROR_TOOLS_CHECK(tolerance->size() == residual->size(),
                                     "The residual and tolerance vectors don't have the same size\n  tolerance: " +
                                         std::to_string(tolerance->size()) +
                                         "\n  residual:  " + std::to_string(residual->size()) + "\n");

        for (unsigned int i = 0; i < tolerance->size(); i++) {
            if (std::fabs((*residual)[i]) > (*tolerance)[i]) {
                return false;
            }
        }

        return true;
    }

    /*!
     * Get the tolerance
     */
    const floatVector *IterativeSolverBase::getTolerance() {
        if (!_tolerance.first) {
            TARDIGRADE_ERROR_TOOLS_CATCH(setTolerance());
        }

        return &_tolerance.second;
    }

    /*!
     * Set the tolerance
     *
     * \f$ tol = tolr * ( |R_0| + |X| ) + tola \f$
     */
    void IterativeSolverBase::setTolerance() {
        auto tolerance = get_SetDataStorage_tolerance();

        *tolerance.value = tardigradeVectorTools::abs(*getResidual()) + tardigradeVectorTools::abs(*getUnknownVector());

        *tolerance.value = getRelativeTolerance() * (*tolerance.value) + getAbsoluteTolerance();
    }

    /*!
     * Set the tolerance
     *
     * \param tolerance: The tolerance vector for each value of the residual
     */
    void IterativeSolverBase::setTolerance(const floatVector &tolerance) { setConstantData(tolerance, _tolerance); }

    /*!
     * Return a SetDataStorageConstant setter for the tolerance
     */
    IterativeSolverBase::SetDataStorageConstant<floatVector> IterativeSolverBase::get_SetDataStorage_tolerance() {
        return SetDataStorageConstant<floatVector>(&_tolerance);
    }

    /*!
     * Check if the number of nonlinear iterations has exceeded the allowable count
     */
    bool IterativeSolverBase::checkIteration() { return getIteration() < getMaxIterations(); }

    /*!
     * Set the initialize unknown vector flag
     *
     * \param &value: The value of the flag
     */
    void IterativeSolverBase::setInitializeUnknownVector(const bool &value) { _initializeUnknownVector = value; }

    /*!
     * Set the maximum number of allowable nonlinear iterations
     *
     * \param &value: The maximum number of iterations
     */
    void IterativeSolverBase::setMaxIterations(const unsigned int &value) { _maxIterations = value; }

    /*!
     * Increment the iteration
     */
    void IterativeSolverBase::incrementIteration() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");
        _iteration++;
        step->damping->reset();
    }

    // BEGIN LEVENBERG MARQUARDT FUNCTIONS

    /*!
     * Temporary function used in extraction of Levenberg Marquardt
     */
    void IterativeSolverBase::performLevenbergMarquardtSolve() {
        TARDIGRADE_ERROR_TOOLS_CHECK(step != nullptr, "The step has not been defined");

        setRankDeficientError(false);

        step->setUseLevenbergMarquardt(true);

        // Turn on projection
        step->enableProjection();

        resetIterations();
        updateUnknownVector(initial_unknown);

        try {
            solve();

        } catch (const convergence_error &e) {
            throw;

        } catch (std::exception &e) {
            step->setUseLevenbergMarquardt(false);

            TARDIGRADE_ERROR_TOOLS_CATCH(throw;)
        }
    }
    // END LEVENBERG MARQUARDT FUNCTIONS

}  // namespace tardigradeHydra
