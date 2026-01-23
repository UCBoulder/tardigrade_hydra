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

        protected:

        private:

    };

}

#endif
