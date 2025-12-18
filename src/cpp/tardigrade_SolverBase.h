/**
  ******************************************************************************
  * \file tardigrade_SolverBase.h
  ******************************************************************************
  * A C++ library for the base classes for solvers
  ******************************************************************************
  */

#ifndef TARDIGRADE_SOLVERBASE_H
#define TARDIGRADE_SOLVERBASE_H

#include"tardigrade_SolverStepBase.h"

namespace tardigradeHydra{

    /*!
     * Base Solver class
     */
    class SolverBase{

        public:

            SolverBase( ) : hydra(NULL){//, step(NULL){ ///!< TODO: Re-enable this

            }

            SolverBase( hydraBase * _hydra ) : hydra( _hydra ){
                /*!
                 * Constructor for NonlinearStepBase
                 *
                 * \param *_hydra: The containing hydraBase object
                 */
            }

            hydraBase *hydra; //!< Pointer to the containing hydra object
            SolverStepBase _step; //!< Temporary object
            SolverStepBase *step = &_step; //!< The object that defines the step to be taken by the solver TODO: Make this an incoming pointer

            virtual void solve( );

            const unsigned int getIteration( );

        protected:

        private:

            friend class tardigradeHydra::hydraBase; //!< TEMP REMOVE THIS
            friend class tardigradeHydra::unit_test::SolverBaseTester; //!< The unit tester for the class

    };

}

#endif
