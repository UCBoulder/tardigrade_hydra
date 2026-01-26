/**
 ******************************************************************************
 * \file tardigrade_IterativeSolverBase.h
 ******************************************************************************
 * The base class for iterative solver objects
 ******************************************************************************
 */

#ifndef TARDIGRADE_ITERATIVESOLVERBASE
#define TARDIGRADE_ITERATIVESOLVERBASE

#include "tardigrade_SolverBase.h"

namespace tardigradeHydra {

    namespace unit_test {

        class IterativeSolverBaseTester;

    }

    /*!
     * The base class for step damping operations to improve
     * stability
     */
    class IterativeSolverBase : public SolverBase {
       public:
        using tardigradeHydra::SolverBase::SolverBase;

        virtual void initialSolveAttempt() override;

        virtual bool checkConvergence();

        const floatVector *getTolerance();

        bool checkIteration();

       protected:
        virtual void callResidualPreNLSolve();

        virtual void callResidualSuccessfulNLStep();

        virtual void callResidualPostNLSolve();

        virtual void setTolerance();

        void setTolerance(const floatVector &tolerance);

        virtual tardigradeHydra::SolverBase::SetDataStorageConstant<floatVector> get_SetDataStorage_tolerance();

       private:
        DataStorage<floatVector> _tolerance;  //!< The tolerance vector for the non-linear solve

        friend class tardigradeHydra::unit_test::IterativeSolverBaseTester; //!< The unit test access class
    };

}  // namespace tardigradeHydra

#endif
