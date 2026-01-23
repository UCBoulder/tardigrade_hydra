/**
  ******************************************************************************
  * \file tardigrade_NonlinearStepBase.h
  ******************************************************************************
  * The header file for the base class to determine the nonlinear step
  ******************************************************************************
  */

#ifndef TARDIGRADE_NONLINEARSTEPBASE_H
#define TARDIGRADE_NONLINEARSTEPBASE_H

#include"tardigrade_TrialStepBase.h"

namespace tardigradeHydra{

    /*!
     * The base nonlinear step class
     */
    class NonlinearStepBase : public TrialStepBase{

        public:

            using tardigradeHydra::TrialStepBase::TrialStepBase;

            virtual void computeTrial() override;

            // BEGIN NEWTON SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            void solveNewtonUpdate(floatVector &deltaX_tr);

            // END NEWTON SOLVER FUNCTIONS

            // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            virtual void solveConstrainedQP(floatVector &dx, const unsigned int kmax = 100);

            // END SQP SOLVER FUNCTIONS

        protected:

        private:

    };

}

#endif
