/**
 ******************************************************************************
 * \file tardigrade_RelaxedSolver.h
 ******************************************************************************
 * A RelaxedSolverBase object that specifies a default internal solver
 ******************************************************************************
 */

#ifndef TARDIGRADE_RELAXEDSOLVER_H
#define TARDIGRADE_RELAXEDSOLVER_H

#include "tardigrade_CoreDefinitions.h"
#include "tardigrade_NewtonSolver.h"
#include "tardigrade_RelaxedSolverBase.h"

namespace tardigradeHydra {

    namespace unit_test {

        class RelaxedSolverBaseTester;  //!< Friend class for RelaxedSolverBase unit testing

    }

    /*!
     * Class which provides a RelaxedSolverBase object with default values
     */
    class RelaxedSolver : public RelaxedSolverBase {
       public:
        RelaxedSolver(hydraBase *_hydra = nullptr);

       protected:
        //! The default internal solver
        tardigradeHydra::NewtonSolver _internal_solver;
    };

}  // namespace tardigradeHydra

#endif
