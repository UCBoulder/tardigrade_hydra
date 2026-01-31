/**
 ******************************************************************************
 * \file tardigrade_NewtonStep.h
 ******************************************************************************
 * A class which defines a Newton-Raphson step
 ******************************************************************************
 */

#ifndef TARDIGRADE_HYDRA_NEWTONSTEP
#define TARDIGRADE_HYDRA_NEWTONSTEP

#include "tardigrade_NonlinearStepBase.h"

namespace tardigradeHydra {

    namespace unit_test {

        class NewtonStepTester;  //!< The test class for NewtonStep

    }

    /*!
     * A class which proposes a Newton-Raphson step to solve a nonlinear
     * problem
     */
    class NewtonStep : public NonlinearStepBase {
       public:
        using tardigradeHydra::NonlinearStepBase::NonlinearStepBase;

        void computeTrial() override;
    };

}  // namespace tardigradeHydra

#endif
