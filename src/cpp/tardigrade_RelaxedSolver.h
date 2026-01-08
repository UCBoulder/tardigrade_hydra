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

    class RelaxedSolver : public SolverBase{

        public:

            using tardigradeHydra::SolverBase::SolverBase;

            using tardigradeHydra::SolverBase::solve;

            const unsigned int getRelaxedIteration( );

            bool checkRelaxedConvergence( );

        protected:

            void setRelaxedIteration( const unsigned int &value );

            void resetRelaxedIteration( );

            void incrementRelaxedIteration( );

            void initializeResiduals( );

        private:

            friend class tardigradeHydra::hydraBase; //!< The base class for hydra TEMP

            unsigned int _relaxedIteration = 0; //!< The current relaxed iteration of the non-linear problem

    };

}

#endif
