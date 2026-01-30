/**
 ******************************************************************************
 * \file tardigrade_SubcyclerSolver.h
 ******************************************************************************
 * A C++ library for the nonlinear solvers which attempt to relax the problem
 * during its solution
 ******************************************************************************
 */

#ifndef TARDIGRADE_SUBCYCLERSOLVER_H
#define TARDIGRADE_SUBCYCLERSOLVER_H

#include "tardigrade_CoreDefinitions.h"
#include "tardigrade_RelaxedSolver.h"//#include "tardigrade_IterativeSolverBase.h"
#include "tardigrade_SolverBase.h"
#include "tardigrade_RelaxedSolver.h"

namespace tardigradeHydra {

    namespace unit_test {

        class SubcyclerSolverTester;  //!< Friend class for SubcyclerSolver unit testing

    }

    /*!
     * Class which controls a solve of a problem which may need to be
     * systematically relaxed in order to achieve the solution
     */
    class SubcyclerSolver : public RelaxedSolver {//IterativeSolverBase {
       public:
        SubcyclerSolver();

        SubcyclerSolver(hydraBase *_hydra);

        SubcyclerSolver(hydraBase *_hydra, SolverBase *_internal_solver_ptr);

//        virtual void initialSolveAttempt() override;
//
//        virtual void convergenceErrorFunction() override;
//
//        virtual void unexpectedErrorFunction() override;
//
//        virtual void reset() override;

        const floatType getCutbackFactor();

        const unsigned int getNumGoodControl();

        const floatType getGrowthFactor();

        const floatType getMinDS();

        void setCutbackFactor(const floatType &value);

        void setNumGoodControl(const unsigned int &value);

        void setGrowthFactor(const floatType &value);

        void setMinDS(const floatType &value);

//        const bool allowStepGrowth();
//
//       protected:
//        RelaxedSolver _internal_solver;  //!< The default internal solver
//
//        SolverBase *internal_solver = &_internal_solver;  //!< A pointer to the solver which will be relaxed
//
//        virtual void performSubcyclerSolve();
//
//        void addSubcyclerHeader();
//
//        void addSubcyclerStepHeader();
//
//        void initializeSubcycler();
//
//        void updatePseudoTimestep();
//
//        void performSubcyclerStep();
//
//        void subcyclerStepSuccess();
//
//        void subcyclerStepFailure();

        virtual void callResidualPreSubcycler();

        virtual void callResidualPostSubcyclerSuccess();

        virtual void callResidualPostSubcyclerFailure();

        floatType sp;  //!< The previous subcycler pseudo-time

        floatType ds;  //!< The subcycler pseudo-timestep

        unsigned int num_good = 0;  //!< The number of good subcycler steps

       private:
//        friend class tardigradeHydra::hydraBase;                       //!< The base class for hydra TEMP
        friend class tardigradeHydra::unit_test::SubcyclerSolverTester;  //!< The unit tester for the class

        floatType _cutback_factor = 0.5;  //!< The factor by which the pseudo-time will be scaled if a solve fails

        floatType _growth_factor =
            1.2;  //!< The factor by which the pseudo-time will be scaled if we can grow the pseudo-timestep

        unsigned int _num_good_control =
            2;  //!< The number of good iterations we need to have before we try and increase the timestep

        floatType _minDS = 1e-2;  //!< The minimum allowable pseudo-timestep

    };

}

#endif
