/**
  ******************************************************************************
  * \file tardigrade_NonlinearStepBase.h
  ******************************************************************************
  * The header file for the base class to determine the nonlinear step
  ******************************************************************************
  */

#ifndef TARDIGRADE_NONLINEARSTEPBASE_H
#define TARDIGRADE_NONLINEARSTEPBASE_H

#include "tardigrade_TrialStepBase.h"
#include "tardigrade_MaxRowPreconditioner.h"
#include "tardigrade_PreconditionerBase.h"

namespace tardigradeHydra{

    /*!
     * The base nonlinear step class
     */
    class NonlinearStepBase : public TrialStepBase{

        public:

            NonlinearStepBase();

            NonlinearStepBase(SolverStepBase *_step);

            NonlinearStepBase(SolverStepBase *_step, PreconditionerBase *_preconditioner);

            virtual void reset() override;

            virtual void computeTrial() override;

            virtual const floatVector *getNonlinearRHS();

            virtual const floatVector *getFlatNonlinearLHS();

            // BEGIN NEWTON SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            void solveNewtonUpdate(floatVector &deltaX_tr);

            // END NEWTON SOLVER FUNCTIONS

            // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            //! Return a flag for whether to use the SQP solver
            const bool getUseSQPSolver() { return _useSQPSolver; }

            virtual void solveConstrainedQP(floatVector &dx, const unsigned int kmax = 100);

            // END SQP SOLVER FUNCTIONS

            MaxRowPreconditioner  _preconditioner;  //!< Default preconditioner
            PreconditionerBase *preconditioner =
                &_preconditioner;  //!< The object that defines the preconditioner TODO: Make this an incoming pointer

        protected:

            // SQP SOLVER FUNCTIONS

            /*!
             * Set whether to use the SQP solver
             *
             * \param &value: The updated value
             */
            void setUseSQPSolver(const unsigned int &value) { _useSQPSolver = value; }

            virtual void initializeActiveConstraints(std::vector<bool> &active_constraints);

            virtual void assembleKKTRHSVector(const floatVector &dx, floatVector &KKTRHSVector,
                                              const std::vector<bool> &active_constraints);

            virtual void assembleKKTMatrix(floatVector &KKTMatrix, const std::vector<bool> &active_constraints);

            virtual void updateKKTMatrix(floatVector &KKTMatrix, const std::vector<bool> &active_constraints);

            // END SQP SOLVER FUNCTIONS

            void addTrialStepOutput();

        private:

            // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            bool _useSQPSolver = false;  //!< The flag for whether to use the SQP solver

            // END SQP SOLVER FUNCTIONS
    };

}

#endif
