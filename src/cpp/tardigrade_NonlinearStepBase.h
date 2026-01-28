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

            virtual const floatVector *getNonlinearRHS();

            virtual const floatVector *getFlatNonlinearLHS();

            MaxRowPreconditioner  _preconditioner;  //!< Default preconditioner
            PreconditionerBase *preconditioner =
                &_preconditioner;  //!< The object that defines the preconditioner TODO: Make this an incoming pointer

        protected:

            void addTrialStepOutput();

        private:

    };

}

#endif
