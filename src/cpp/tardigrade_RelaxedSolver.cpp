/**
 ******************************************************************************
 * \file tardigrade_RelaxedSolver.cpp
 ******************************************************************************
 * A RelaxedSolverBase object with a default internal_solver
 ******************************************************************************
 */

#include "tardigrade_RelaxedSolver.h"

namespace tardigradeHydra {

    /*!
     * Constructor for RelaxedSolverBase
     *
     * \param *_hydra: The containing hydra object
     */
    RelaxedSolver::RelaxedSolver(hydraBase *_hydra) : RelaxedSolverBase(_hydra) {
        internal_solver        = &_internal_solver;
        internal_solver->hydra = _hydra;
        step                   = internal_solver->step;
    }

}  // namespace tardigradeHydra
