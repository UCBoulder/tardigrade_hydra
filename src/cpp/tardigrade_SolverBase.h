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
#include"tardigrade_PreconditionerBase.h"

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

            PreconditionerBase _preconditioner; //!< Temporary object
            PreconditionerBase *preconditioner = &_preconditioner; //!< The object that defines the preconditioner TODO: Make this an incoming pointer

            virtual void solve( );

            const unsigned int getIteration( );

            const bool getRankDeficientError( );

            void setRankDeficientError( const bool &value );

        protected:

        private:

            bool _rank_deficient_error = false; //!< Flag for whether a rank-deficient Jacobian should cause an error

            friend class tardigradeHydra::hydraBase; //!< TEMP REMOVE THIS
            friend class tardigradeHydra::unit_test::SolverBaseTester; //!< The unit tester for the class

    };

}

#endif
