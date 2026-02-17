/**
 ******************************************************************************
 * \file tardigrade_SequentialQuadraticProgrammingStep.h
 ******************************************************************************
 * A class which defines a Sequential Quadratic Programming step
 ******************************************************************************
 */

#ifndef TARDIGRADE_HYDRA_SEQUENTIALQUADRATICPROGRAMMINGSTEP
#define TARDIGRADE_HYDRA_SEQUENTIALQUADRATICPROGRAMMINGSTEP

#include "tardigrade_NonlinearStepBase.h"

namespace tardigradeHydra {

    namespace unit_test {

        class SequentialQuadraticProgrammingStepTester;  //!< The test class for SequentialQuadraticProgrammingStep

    }

    /*!
     * A class which proposes a Newton-Raphson step to solve a nonlinear
     * problem
     */
    class SequentialQuadraticProgrammingStep : public NonlinearStepBase {
       public:
        using tardigradeHydra::NonlinearStepBase::NonlinearStepBase;

        void computeTrial() override;

        unsigned int kmax = 100;  //!< The maximum number of iterations

       protected:
        virtual void initializeActiveConstraints(std::vector<bool> &active_constraints);

        virtual void assembleKKTRHSVector(const floatVector &dx, floatVector &KKTRHSVector,
                                          const std::vector<bool> &active_constraints);

        virtual void assembleKKTMatrix(floatVector &KKTMatrix, const std::vector<bool> &active_constraints);

        virtual void updateKKTMatrix(floatVector &KKTMatrix, const std::vector<bool> &active_constraints);

       private:
    };

}  // namespace tardigradeHydra

#include "tardigrade_SequentialQuadraticProgrammingStep.tpp"

#endif
