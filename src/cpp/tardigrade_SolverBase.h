/**
 ******************************************************************************
 * \file tardigrade_SolverBase.h
 ******************************************************************************
 * A C++ library for the base classes for solvers
 ******************************************************************************
 */

#ifndef TARDIGRADE_SOLVERBASE_H
#define TARDIGRADE_SOLVERBASE_H

#include "tardigrade_CoreDefinitions.h"
#include "tardigrade_SetDataStorage.h"
#include "tardigrade_SolverStepBase.h"

namespace tardigradeHydra {

    namespace unit_test {
        class SolverBaseTester;          //!< Friend class for SolverBase for unit testing
    }  // namespace unit_test

    /*!
     * Base Solver class
     */
    class SolverBase : public CachingDataBase {
       public:
        SolverBase();

        SolverBase(hydraBase *_hydra);

        SolverBase(hydraBase *_hydra, SolverStepBase *_step_ptr);

        hydraBase *hydra;  //!< Pointer to the containing hydra object

        SolverStepBase  _step;  //!< Temporary object
        SolverStepBase *step =
            &_step;  //!< The object that defines the step to be taken by the solver TODO: Make this an incoming pointer

        floatVector initial_unknown;  //!< The initial unknown vector for the solver

        virtual void solve();

        virtual void initializeSolve();

        virtual void initialSolveAttempt();

        virtual void convergenceErrorFunction();

        virtual void unexpectedErrorFunction();

        virtual void reset();

        const bool getRankDeficientError();

        void setRankDeficientError(const bool &value);

        // CACHED DATA STORAGE OPERATIONS
        virtual void addIterationData(dataBase *data) override;

        virtual void addNLStepData(dataBase *data) override;
        // END CACHED DATA STORAGE OPERATIONS

        // NONLINEAR FUNCTIONS (MOVE TO OWN CLASS)

        void resetNLStepData();

        // END NONLINEAR FUNCTIONS

        // Pass-through functions
        const floatType getRelativeTolerance();  // TODO: Want to allow this to be constexpr

        const floatType getAbsoluteTolerance();  // TODO: Want to allow this to be constexpr

        const unsigned int getNumUnknowns();  // TODO: Want to allow this to be constexpr

        const floatVector *getUnknownVector();  // TODO: Want to generalize this

        void initializeUnknownVector();

        void updateUnknownVector(const floatVector &value);  // TODO: Want to generalize this

        const floatVector *getResidual();  // TODO: Want to generalize this

        const floatVector *getFlatJacobian();  // TODO: Want to generalize this

        const unsigned int getNumConstraints();  // TODO: Want to allow this to be constexpr

        const floatVector *getConstraints();  // TODO: Want to generalize this

        const floatVector *getConstraintJacobians();  // TODO: Want to generalize this

        const unsigned int getFailureVerbosityLevel();

        void addToFailureOutput(const std::string &string);

        void addToFailureOutput(const floatVector &value, bool add_endline = true);

        void addToFailureOutput(const std::vector<bool> &value, bool add_endline = true);

        void addToFailureOutput(const floatType &value, bool add_endline = true);

        const floatType getToleranceScaleFactor();

        void resetToleranceScaleFactor();

        void setCurrentResidualIndexMeaningful(const bool &value);

        void setCurrentResidualIndex(const unsigned int &value);

        const std::vector<tardigradeHydra::ResidualBase<> *> *getResidualClasses();

        void setAllowModifyGlobalResidual(const bool &value);

        /*!
         * Add a general iterable object to the output string
         *
         * \param &v_begin: The starting iterator
         * \param &v_end: The stopping iterator
         * \param add_endline: Whether to add an endline to the string or not
         */
        template <class v_iterator>
        void addToFailureOutput(const v_iterator &v_begin, const v_iterator &v_end, bool add_endline = true) {
            std::stringstream failure_output;

            for (auto v = v_begin; v != v_end; ++v) {
                failure_output << *v << ", ";
            }

            if (add_endline) {
                failure_output << "\n";
            }

            addToFailureOutput(failure_output.str());
        }

        // Levenberg marquardt functions (move to own class)

        void performLevenbergMarquardtSolve();  // TEMP remove this

        // end Levenberg marquard functions
       protected:

       private:
        bool _rank_deficient_error = false;  //!< Flag for whether a rank-deficient Jacobian should cause an error

        friend class tardigradeHydra::hydraBase;                    //!< TEMP REMOVE THIS
        friend class tardigradeHydra::unit_test::SolverBaseTester;  //!< The unit tester for the class

    };

}  // namespace tardigradeHydra

#endif
