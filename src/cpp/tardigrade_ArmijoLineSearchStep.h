/**
 ******************************************************************************
 * \file tardigrade_ArmijoLineSearchStep.h
 ******************************************************************************
 * A C++ library for the Armijo line search step
 ******************************************************************************
 */

#ifndef TARDIGRADE_ARMIJOLINESEARCHSTEP_H
#define TARDIGRADE_ARMIJOLINESEARCHSTEP_H

#include "tardigrade_CoreDefinitions.h"
#include "tardigrade_SolverBase.h"

namespace tardigradeHydra {

    /*!
     * The Armijo line search step class
     */
    class ArmijoLineSearchStep : virtual public SolverStepBase {
       public:
        using tardigradeHydra::SolverStepBase::SolverStepBase;

       protected:
       private:
    };

}  // namespace tardigradeHydra

#include "tardigrade_ArmijoLineSearchStep.tpp"

#endif
