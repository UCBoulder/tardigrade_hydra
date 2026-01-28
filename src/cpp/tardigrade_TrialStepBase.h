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

        StepDampingBase *getDamping( );

        void setCurrentResidualIndexMeaningful(const bool &value);

        void setCurrentResidualIndex(const unsigned int &value);

        const std::vector<tardigradeHydra::ResidualBase<> *> *getResidualClasses();

        // END PASS-THROUGH FUNCTIONS

       protected:

       private:
        friend class tardigradeHydra::SolverStepBase;                  //!< TEMP REMOVE THIS
        friend class tardigradeHydra::unit_test::TrialStepBaseTester;  //!< The unit tester for the class
    };

}  // namespace tardigradeHydra

#endif
