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

#include"tardigrade_CoreDefinitions.h"
#include"tardigrade_SolverBase.h"

namespace tardigradeHydra{

    /*!
     * Class which controls a solve of a problem which may need to be
     * systematically relaxed in order to achieve the solution
     */
    class RelaxedSolver : public SolverBase{

        public:

            using tardigradeHydra::SolverBase::SolverBase;

            using tardigradeHydra::SolverBase::solve;

            const unsigned int getRelaxedIteration( );

            const unsigned int getMaxRelaxedIterations( );

            const void setMaxRelaxedIterations( const unsigned int &value );

            bool checkRelaxedConvergence( );

            void setInternalSolver( SolverBase *_solver );

            virtual void reset( ) override;

//            virtual void performRelaxedSolve( );

        protected:

            SolverBase *internal_solver = NULL; //!< A pointer to the solver which will be relaxed

            void setRelaxedIteration( const unsigned int &value );

            void resetRelaxedIteration( );

            void incrementRelaxedIteration( );

            void initializeResiduals( );

            bool attemptInternalSolve( );

            virtual bool callResidualRelaxedStepFailure( );

            void setupNextRelaxedStep( );

        private:

            friend class tardigradeHydra::hydraBase; //!< The base class for hydra TEMP

            unsigned int _relaxedIteration = 0; //!< The current relaxed iteration of the non-linear problem

            unsigned int _maxRelaxedIterations = 5; //!< The number of allowed relaxed iterations

            void logRelaxedIterationHeader( );

    };

}

#endif
