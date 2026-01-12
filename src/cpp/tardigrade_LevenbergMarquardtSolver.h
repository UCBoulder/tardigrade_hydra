/**
  ******************************************************************************
  * \file tardigrade_LevenbergMarquardtSolver.h
  ******************************************************************************
  * A C++ library for the Levenberg Marquardt solver
  ******************************************************************************
  */

#ifndef TARDIGRADE_LEVENBERGMARQUARDTSOLVER_H
#define TARDIGRADE_LEVENBERGMARQUARDTSOLVER_H

#include"tardigrade_CoreDefinitions.h"
#include"tardigrade_SolverBase.h"

namespace tardigradeHydra{

    /*!
     * The Levenberg Marquardt solver class
     */
    class LevenbergMarquardtSolver : public SolverBase {

        public:

            using tardigradeHydra::SolverBase::SolverBase;

            using tardigradeHydra::SolverBase::solve;

    };

}

#endif
