/**
 ******************************************************************************
 * \file tardigrade_SolverStepBase.h
 ******************************************************************************
 * The base class for solver steps
 ******************************************************************************
 */

#ifndef TARDIGRADE_SOLVERSTEPBASE
#define TARDIGRADE_SOLVERSTEPBASE

#include "tardigrade_CoreDefinitions.h"
#include "tardigrade_SetDataStorage.h"
#include "tardigrade_NonlinearStepBase.h"
#include "tardigrade_StepDampingBase.h"
// Default classes
#include "tardigrade_NewtonStep.h"
#include "tardigrade_ArmijoGradientDamping.h"

namespace tardigradeHydra {

    namespace unit_test {
        class SolverStepBaseTester;      //!< Friend class for SolverStepBase for unit testing
    }  // namespace unit_test

    /*!
     * Base class for Solver Steps
     */
    class SolverStepBase : public CachingDataBase {
       public:
        /*!
         * Constructor for NonlinearStepBase
         */
        SolverStepBase() : solver(NULL) { initializeDefaults(); }

        /*!
         * Constructor for NonlinearStepBase
         *
         * \param *_solver: The containing solver object
         */
        SolverStepBase(SolverBase *_solver) : solver(_solver) { initializeDefaults(); }

        virtual void reset();

        void incrementSolution();

        floatVector X0;  //!< The initial value of the unknown vector

        floatVector deltaX;  //!< The change in the unknown vector

        /*! Set the containing solver object
         * \param *_solver: The containing solver object
         */
        void setSolver(SolverBase *_solver) { solver = _solver; }

        // CACHED DATA STORAGE OPERATIONS
        virtual void addIterationData(dataBase *data) override;

        virtual void addNLStepData(dataBase *data) override;
        // END CACHED DATA STORAGE OPERATIONS

        // PASS-THROUGH functions

        const floatType getRelativeTolerance();

        const floatType getAbsoluteTolerance();

        const floatVector *getResidual();

        const unsigned int getNumUnknowns();

        const floatVector *getUnknownVector();

        void updateUnknownVector(const floatVector &value);

        const floatVector *getFlatJacobian();

        const unsigned int getNumConstraints();

        const floatVector *getConstraints();

        const floatVector *getConstraintJacobians();

        const floatType getToleranceScaleFactor();

        void resetToleranceScaleFactor();

        const unsigned int getFailureVerbosityLevel();

        void addToFailureOutput(const std::string &string);

        void addToFailureOutput(const floatVector &value, bool add_endline = true);

        void addToFailureOutput(const std::vector<bool> &value, bool add_endline = true);

        void addToFailureOutput(const floatType &value, bool add_endline = true);

        void setCurrentResidualIndexMeaningful(const bool &value);

        void setCurrentResidualIndex(const unsigned int &value);

        const std::vector<tardigradeHydra::ResidualBase<> *> *getResidualClasses();

        // END PASS-THROUGH FUNCTIONS

        const bool getRankDeficientError();

        void setRankDeficientError(const bool &value);

        //! Get the number of undamped steps performed
        unsigned int getNumUndamped() { return _NUM_UNDAMPED; }

        TrialStepBase   *trial_step;  //!< The trial step class which proposes a step to reduce the residual
        StepDampingBase *damping;     //!< The damping class which reduces the proposed step to improve stability

       protected:
        SolverBase *solver;  //!< Pointer to the containing SolverBase object

        ArmijoGradientDamping _damping;     //!< The default step damping
        NewtonStep            _trial_step;  //!< The default trial step

        void initializeDefaults();

        //! Reset the number of undamped steps
        void resetNumUndamped() { _NUM_UNDAMPED = 0; }

        //! Increment the number of undamped steps
        void incrementNumUndamped() { _NUM_UNDAMPED++; }

       private:
        friend class tardigradeHydra::hydraBase;                        //!< TEMP REMOVE THIS
        friend class tardigradeHydra::unit_test::SolverStepBaseTester;  //!< The unit tester for the class

        bool _rank_deficient_error = false;  //!< Flag for whether a rank-deficient Jacobian should cause an error

        unsigned int _NUM_UNDAMPED = 0;  //!< The number of undamped steps performed

    };

}  // namespace tardigradeHydra

#endif
