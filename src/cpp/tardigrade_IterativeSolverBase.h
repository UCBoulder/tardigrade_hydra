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

        //! Return the maximum number of allowable iterations
        const unsigned int getMaxIterations() { return _maxIterations; }

        void setMaxIterations(const unsigned int &value);

        //! Return the flag which indicates whether hydra should initialize the unknown vector
        const bool getInitializeUnknownVector() { return _initializeUnknownVector; }

       protected:
        virtual void callResidualPreIterativeSolve();

        virtual void callResidualSuccessfulIterativeStep();

        virtual void callResidualPostIterativeSolve();

        void setInitializeUnknownVector(const bool &value);

        virtual void setTolerance();

        void setTolerance(const floatVector &tolerance);

        virtual tardigradeHydra::SolverBase::SetDataStorageConstant<floatVector> get_SetDataStorage_tolerance();

       private:
        friend class tardigradeHydra::unit_test::IterativeSolverBaseTester; //!< The unit test access class

        DataStorage<floatVector> _tolerance;  //!< The tolerance vector for the non-linear solve

        bool _initializeUnknownVector =
            true;  //!< Flag for whether to initialize the unknown vector in the non-linear solve

        unsigned int _maxIterations = 20;  //!< The maximum number of allowable iterations

    };

}  // namespace tardigradeHydra

#endif
