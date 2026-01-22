/**
 ******************************************************************************
 * \file tardigrade_NonlinearSolverBase.h
 ******************************************************************************
 * The base class for nonlinear solver objects
 ******************************************************************************
 */

#ifndef TARDIGRADE_NONLINEARSOLVERBASE
#define TARDIGRADE_NONLINEARSOLVERBASE

#include"tardigrade_SolverBase.h"

namespace tardigradeHydra{

    namespace unit_test{

        class NonlinearSolveBaseTester;

    }

    /*!
     * The base class for step damping operations to improve
     * stability
     */
    class NonlinearSolverBase : public SolverBase{

        public:

            using tardigradeHydra::SolverBase::SolverBase;

            virtual void initialSolveAttempt( ) override;

            virtual bool checkConvergence( );

            const floatVector* getTolerance( );

        protected:

            virtual void callResidualPreNLSolve( );

            virtual void callResidualSuccessfulNLStep( );

            virtual void callResidualPostNLSolve( );

            virtual void setTolerance( );

            void setTolerance( const floatVector &tolerance );

            virtual tardigradeHydra::SolverBase::SetDataStorageConstant<floatVector> get_SetDataStorage_tolerance( );

        private:

            DataStorage< floatVector > _tolerance; //!< The tolerance vector for the non-linear solve

    };

}

#endif
