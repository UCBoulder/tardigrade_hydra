/**
 ******************************************************************************
 * \file tardigrade_SequentialQuadraticProgrammingStep.h
 ******************************************************************************
 * A class which defines a Sequential Quadratic Programming step
 ******************************************************************************
 */

#ifndef TARDIGRADE_HYDRA_NEWTONSTEP
#define TARDIGRADE_HYDRA_NEWTONSTEP

#include "tardigrade_NonlinearStepBase.h"

namespace tardigradeHydra {

    namespace unit_test {

        class SequentialQuadraticProgrammingStepTester; //!< The test class for SequentialQuadraticProgrammingStep

    }

    /*!
     * A class which proposes a Newton-Raphson step to solve a nonlinear
     * problem
     */
    class SequentialQuadraticProgrammingStep : public NonlinearStepBase {

        public:

            using tardigradeHydra::NonlinearStepBase::NonlinearStepBase;

            void computeTrial() override;

            //! Return a flag for whether to use the SQP solver
            const bool getUseSQPSolver() { return _useSQPSolver; } // TODO: remove this

            unsigned int kmax = 100; //!< The maximum number of iterations

        protected:

            /*!
             * Set whether to use the SQP solver
             *
             * \param &value: The updated value
             */
            void setUseSQPSolver(const unsigned int &value) { _useSQPSolver = value; } // TODO: Remove this

            virtual void initializeActiveConstraints(std::vector<bool> &active_constraints);

            virtual void assembleKKTRHSVector(const floatVector &dx, floatVector &KKTRHSVector,
                                              const std::vector<bool> &active_constraints);

            virtual void assembleKKTMatrix(floatVector &KKTMatrix, const std::vector<bool> &active_constraints);

            virtual void updateKKTMatrix(floatVector &KKTMatrix, const std::vector<bool> &active_constraints);

        private:

            bool _useSQPSolver = false;  //!< The flag for whether to use the SQP solver TODO: Remove this

    };

}

#endif
