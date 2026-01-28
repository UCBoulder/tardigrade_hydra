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

        virtual void reset() override;

        virtual void initializeSolve() override;

        virtual void initialSolveAttempt() override;

        virtual bool checkConvergence();

        const floatVector *getTolerance();

        bool checkIteration();

        //! Return the maximum number of allowable iterations
        const unsigned int getMaxIterations() { return _maxIterations; }

        void setMaxIterations(const unsigned int &value);

        //! Return the flag which indicates whether hydra should initialize the unknown vector
        const bool getInitializeUnknownVector() { return _initializeUnknownVector; }

        //! Get the current nonlinear iteration number
        const unsigned int getIteration() { return _iteration; }

        //! Reset the number of iterations TODO: Determine if there is another way rather than making this public
        void resetIterations() { _iteration = 0; }

        void addIterationHeader();

        void addIterationFooter();

       protected:
        virtual void callResidualPreIterativeSolve();

        virtual void callResidualSuccessfulIterativeStep();

        virtual void callResidualPostIterativeSolve();

        void setInitializeUnknownVector(const bool &value);

        virtual void setTolerance();

        void setTolerance(const floatVector &tolerance);

        virtual tardigradeHydra::SolverBase::SetDataStorageConstant<floatVector> get_SetDataStorage_tolerance();

        void incrementIteration();

       private:
        friend class tardigradeHydra::unit_test::IterativeSolverBaseTester; //!< The unit test access class

        DataStorage<floatVector> _tolerance;  //!< The tolerance vector for the non-linear solve

        bool _initializeUnknownVector =
            true;  //!< Flag for whether to initialize the unknown vector in the non-linear solve

        unsigned int _maxIterations = 20;  //!< The maximum number of allowable iterations

        unsigned int _iteration = 0;  //!< The current iteration of the non-linear problem

    };

}  // namespace tardigradeHydra

#endif
