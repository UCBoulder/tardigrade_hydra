/**
 ******************************************************************************
 * \file tardigrade_SolverStep.h
 ******************************************************************************
 * Core definitions for tardigrade hydra
 ******************************************************************************
 */

#ifndef TARDIGRADE_HYDRA_SOLVERSTEP
#define TARDIGRADE_HYDRA_SOLVERSTEP

#include "tardigrade_ArmijoGradientDamping.h"
#include "tardigrade_NewtonStep.h"
#include "tardigrade_SolverStepBase.h"

namespace tardigradeHydra {

    /*!
     * The default SolverStep
     */
    class SolverStep : public SolverStepBase {
       public:
        /*!
         * Constructor for the default SolverStep class
         */
        SolverStep() : SolverStepBase(NULL) { initializeDefaults(); }

        /*!
         * Constructor for the default SolverStep class
         *
         * \param *_solver: The containing solver class
         */
        SolverStep(SolverBase *_solver) : SolverStepBase(_solver) { initializeDefaults(); }

        /*!
         * Initialize the default classes
         */
        void initializeDefaults() {
            damping       = &_damping;
            damping->step = this;

            trial_step       = &_trial_step;
            trial_step->step = this;
        }

       protected:
        ArmijoGradientDamping _damping;     //!< The default step damping
        NewtonStep            _trial_step;  //!< The default trial step
    };

}  // namespace tardigradeHydra

#include "tardigrade_SolverStep.cpp"
#include "tardigrade_SolverStep.tpp"

#endif
