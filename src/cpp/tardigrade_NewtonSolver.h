/**
 ******************************************************************************
 * \file tardigrade_NewtonSolver.h
 ******************************************************************************
 * A C++ library for a Newton solver that uses Armijo line search and, if
 * required, gradient descent to solve the non-linear problem
 ******************************************************************************
 */

#ifndef TARDIGRADE_NEWTONSOLVER_H
#define TARDIGRADE_NEWTONSOLVER_H

#include "tardigrade_IterativeSolverBase.h"

namespace tardigradeHydra {

    namespace unit_test {
        class NewtonSolverTester;  //!< Friend class for NewtonSolver for unit testing
    }  // namespace unit_test

    /*!
     * The Newton solver class
     */
    class NewtonSolver : public IterativeSolverBase {
       public:
        /*!
         * Constructor for the Newton solver
         *
         * \param *_hydra: A pointer to the containing hydraBase object
         */
        NewtonSolver(hydraBase *_hydra = nullptr) : IterativeSolverBase(_hydra) {
            step         = &_step;
            step->solver = this;
        }

       protected:
        //! The incrementation object
        tardigradeHydra::SolverStep _step;
    };

}  // namespace tardigradeHydra

#include "tardigrade_NewtonSolver.cpp"
#include "tardigrade_NewtonSolver.tpp"

#endif
