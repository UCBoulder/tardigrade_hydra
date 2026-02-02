/**
 ******************************************************************************
 * \file tardigrade_PreconditionerBase.h
 ******************************************************************************
 * A C++ library for the base classes for solvers
 ******************************************************************************
 */

#ifndef TARDIGRADE_PRECONDITIONERBASE_H
#define TARDIGRADE_PRECONDITIONERBASE_H

#include "tardigrade_CoreDefinitions.h"
#include "tardigrade_SetDataStorage.h"

namespace tardigradeHydra {

    namespace unit_test {

        class PreconditionerBaseTester;  //!< Friend class for PreconditionerBase testing

    }

    /*!
     * The base class for preconditioners to be used in tardigrade hydra solves
     */
    class PreconditionerBase : public CachingDataBase {
       public:
        /*!
         * Initialize the preconditioner object
         */
        PreconditionerBase() { trial_step = nullptr; }

        /*!
         * Constructor for the preconditioner object
         *
         * \param *_trial_step: The containing TrialStepBase object
         */
        PreconditionerBase(NonlinearStepBase *_trial_step) { trial_step = _trial_step; }

        virtual void reset();

        //! Get a pointer to the row-major form of the preconditioner
        const floatVector *getFlatPreconditioner();

        NonlinearStepBase *trial_step;  //!< Pointer to the containing NonlinearStepBase class

        // PASS THROUGH FUNCTIONS
        const unsigned int getNumUnknowns();

        const floatVector *getFlatNonlinearLHS();

        const floatVector *getNonlinearRHS();
        // END PASS THROUGH FUNCTIONS

        virtual void preconditionVector(const floatVector &X, floatVector &Y);

        virtual void preconditionMatrix(const floatVector &A, floatVector &B);

       protected:
        virtual void formPreconditioner();

        // CACHED DATA STORAGE OPERATIONS
        virtual void addIterationData(dataBase *data) override;

        virtual void addNLStepData(dataBase *data) override;

        DataStorage<floatVector>
            _preconditioner;  //!< The pre-conditioner matrix in row-major form for the global solve
                              // END CACHED DATA STORAGE OPERATIONS

       private:
        friend class tardigradeHydra::unit_test::PreconditionerBaseTester;  //!< The unit tester for the class
    };

}  // namespace tardigradeHydra

#endif
