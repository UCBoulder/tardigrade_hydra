
/**
  ******************************************************************************
  * \file tardigrade_SolverBase.cpp
  ******************************************************************************
  * A C++ library for the base classes for solvers
  ******************************************************************************
  */

#include"tardigrade_PreconditionerBase.h"
#include"tardigrade_hydra.h"

namespace tardigradeHydra{

    void PreconditionerBase::addIterationData( dataBase *data ){
        /*!
         * Add data to the vector of values which will be cleared after each iteration
         * 
         * \param *data: The dataBase object to be cleared
         */

        solver->hydra->addIterationData( data );

    }

    void PreconditionerBase::addNLStepData( dataBase *data ){
        /*!
         * Add data to the vector of values which will be cleared after each nonlinear step
         * 
         * \param *data: The dataBase object to be cleared
         */

        solver->hydra->addNLStepData( data );

    }

}
