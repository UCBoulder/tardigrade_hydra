/**
 ******************************************************************************
 * \file tardigrade_RelaxedSolver.h
 ******************************************************************************
 * A C++ library for the nonlinear solvers which attempt to relax the problem
 * during its solution
 ******************************************************************************
 */

#ifndef TARDIGRADE_RELAXEDSOLVER_H
#define TARDIGRADE_RELAXEDSOLVER_H

#include "tardigrade_CoreDefinitions.h"
#include "tardigrade_NewtonSolver.h"
#include "tardigrade_SolverBase.h"

namespace tardigradeHydra {

    namespace unit_test {

        class RelaxedSolverTester;  //!< Friend class for RelaxedSolver unit testing

    }

    /*!
     * Class which controls a solve of a problem which may need to be
     * systematically relaxed in order to achieve the solution
     */
    class RelaxedSolver : public IterativeSolverBase {
       public:
        RelaxedSolver();

        RelaxedSolver(hydraBase *_hydra);

        RelaxedSolver(hydraBase *_hydra, SolverBase *_internal_solver_ptr);

        virtual void initialSolveAttempt() override;

        virtual void convergenceErrorFunction() override;

        virtual void reset() override;

        const unsigned int getRelaxedIteration();

        const unsigned int getMaxRelaxedIterations();

        const void setMaxRelaxedIterations(const unsigned int &value);

        bool checkRelaxedConvergence();

        virtual void performRelaxedSolve();

        SolverBase *internal_solver = &_internal_solver;  //!< A pointer to the solver which will be relaxed

       protected:
        NewtonSolver _internal_solver;  //!< The default internal solver

        void setRelaxedIteration(const unsigned int &value);

        void resetRelaxedIteration();

        void incrementRelaxedIteration();

        void initializeResiduals();

        bool attemptInternalSolve();

        virtual bool callResidualRelaxedStepFailure();

        void setupNextRelaxedStep();

       private:
        friend class tardigradeHydra::unit_test::RelaxedSolverTester;  //!< The unit tester for the class

        unsigned int _relaxedIteration = 0;  //!< The current relaxed iteration of the non-linear problem

        unsigned int _maxRelaxedIterations = 5;  //!< The number of allowed relaxed iterations

        void logRelaxedIterationHeader();
    };

}  // namespace tardigradeHydra

#endif
