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
#include "tardigrade_IterativeSolverBase.h"
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
    class SubcyclerSolver : public IterativeSolverBase {
       public:
        SubcyclerSolver();

        SubcyclerSolver(hydraBase *_hydra);

        SubcyclerSolver(hydraBase *_hydra, SolverBase *_internal_solver_ptr);

        virtual void initialSolveAttempt() override;

        virtual void convergenceErrorFunction() override;

        virtual void unexpectedErrorFunction() override;

        const floatType getCutbackFactor() { /*! Get the value of the cutback factor */ return _cutback_factor; }

        const unsigned int
        getNumGoodControl() { /*! Get the number of good iterations we need to have before increasing the timestep */
            return _num_good_control;
        }

        const floatType getGrowthFactor() { /*! Get the growth factor for the timestep increase */
            return _growth_factor;
        }

        const floatType
        getMinDS() { /*! Get the minimum allowable ratio of the total timestep to the cutback timestep */
            return _minDS;
        }

        void setCutbackFactor(const floatType &value) { /*! Get the current value of the cutback factor. \param &value:
                                                           The value of the cutback */
            _cutback_factor = value;
        }

        void setNumGoodControl(
            const unsigned int &value) { /*! Set the number of good iterations that need to happen before the timestep
                                            increases. \param &value: The value of the number of good iterations prior
                                            to increasing the relative timestep */
            _num_good_control = value;
        }

        void setGrowthFactor(const floatType &value) { /*! Set the relative growth factor for the local timestep
                                                          increase \param &value: The new value */
            _growth_factor = value;
        }

        void setMinDS(const floatType &value) { /*! Set the minimum value of the relative cutback timestep \param
                                                   &value: The new value */
            _minDS = value;
        }

        const bool allowStepGrowth(const unsigned int &num_good);

       protected:
        RelaxedSolver _internal_solver;  //!< The default internal solver

        SolverBase *internal_solver = &_internal_solver;  //!< A pointer to the solver which will be relaxed

        virtual void performSubcyclerSolve(); //TEMP

        void addSubcyclerHeader(); //TEMP

        void addSubcyclerStepHeader(); //TEMP

        void initializeSubcycler(); //TEMP

        void updatePseudoTimestep(); //TEMP

        void performSubcyclerStep(); //TEMP

        void subcyclerStepSuccess(); //TEMP

        void subcyclerStepFailure(); //TEMP

        floatType sp;  //!< The previous subcycler pseudo-time

        floatType ds;  //!< The subcycler pseudo-timestep

        unsigned int num_good = 0;  //!< The number of good subcycler steps

       private:
        friend class tardigradeHydra::hydraBase;                       //!< The base class for hydra TEMP
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
