/**
 ******************************************************************************
 * \file tardigrade_TrialStepBase.h
 ******************************************************************************
 * The base class for trial step proposal operations
 ******************************************************************************
 */

#ifndef TARDIGRADE_TRIALSTEPBASE
#define TARDIGRADE_TRIALSTEPBASE

#include "tardigrade_CoreDefinitions.h"
#include "tardigrade_PreconditionerBase.h"
#include "tardigrade_SetDataStorage.h"

namespace tardigradeHydra {

    namespace unit_test {

        class TrialStepBaseTester;

    }

    /*!
     * The base class for step damping operations to improve
     * stability
     */
    class TrialStepBase : public CachingDataBase {
       public:
        TrialStepBase();

        TrialStepBase(SolverStepBase *_step);

        TrialStepBase(SolverStepBase *_step, PreconditionerBase *_preconditioner_ptr);

        SolverStepBase *step;  //!< The containing step class

        virtual void resetCounts();

        virtual void reset();

        virtual void computeTrial();

        // CACHED DATA STORAGE OPERATIONS
        virtual void addIterationData(dataBase *data) override;

        virtual void addNLStepData(dataBase *data) override;
        // END CACHED DATA STORAGE OPERATIONS

        // PASS-THROUGH FUNCTIONS

        const floatType getRelativeTolerance();

        const floatType getAbsoluteTolerance();

        const floatVector *getResidual();

        const unsigned int getNumUnknowns();

        const floatVector *getFlatJacobian();

        const unsigned int getNumConstraints();

        const floatVector *getConstraints();

        const floatVector *getConstraintJacobians();

        bool getRankDeficientError();

        const unsigned int getFailureVerbosityLevel();

        void addToFailureOutput(const std::string &string);

        void addToFailureOutput(const floatVector &value, bool add_endline = true);

        void addToFailureOutput(const std::vector<bool> &value, bool add_endline = true);

        void addToFailureOutput(const floatType &value, bool add_endline = true);

        // END PASS-THROUGH FUNCTIONS

        // NONLINEAR FUNCTIONS (MOVE TO OWN CLASS)

        virtual const floatVector *getNonlinearRHS();

        virtual const floatVector *getFlatNonlinearLHS();

        // END NONLINEAR FUNCTIONS

        // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

        //! Return a flag for whether to use the SQP solver
        const bool getUseSQPSolver() { return _useSQPSolver; }

        // END SQP SOLVER FUNCTIONS

        void performPreconditionedSolve(floatVector &deltaX_tr);  // TEMP REMOVE THIS

        PreconditionerBase  _preconditioner;  //!< Temporary object
        PreconditionerBase *preconditioner =
            &_preconditioner;  //!< The object that defines the preconditioner TODO: Make this an incoming pointer

       protected:
        // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

        /*!
         * Set whether to use the SQP solver
         *
         * \param &value: The updated value
         */

        void setUseSQPSolver(const unsigned int &value) { _useSQPSolver = value; }

        virtual void initializeActiveConstraints(std::vector<bool> &active_constraints);

        virtual void assembleKKTMatrix(floatVector &KKTMatrix, const std::vector<bool> &active_constraints);

        virtual void updateKKTMatrix(floatVector &KKTMatrix, const std::vector<bool> &active_constraints);

        // END SQP SOLVER FUNCTIONS

       private:
        friend class tardigradeHydra::SolverStepBase;                  //!< TEMP REMOVE THIS
        friend class tardigradeHydra::unit_test::TrialStepBaseTester;  //!< The unit tester for the class
        // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

        bool _useSQPSolver = false;  //!< The flag for whether to use the SQP solver

        // END SQP SOLVER FUNCTIONS
    };

}  // namespace tardigradeHydra

#endif
