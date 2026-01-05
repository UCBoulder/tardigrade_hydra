/**
  ******************************************************************************
  * \file tardigrade_SolverBase.cpp
  ******************************************************************************
  * A C++ library for the base classes for solvers
  ******************************************************************************
  */

#include"tardigrade_SolverBase.h"
#include"tardigrade_hydra.h"

namespace tardigradeHydra{

    /*!
     * Get whether the Jacobian being rank-deficient will throw an error
     */
    const bool SolverBase::getRankDeficientError( ){

        return _rank_deficient_error;

    }

    /*!
     * Set whether the Jacobian being rank-deficient will throw an error
     *
     * \param &value: The incoming value
     */
    void SolverBase::setRankDeficientError( const bool &value ){

        _rank_deficient_error = value;
    }

    /*!
     * Add data to the vector of values which will be cleared after each iteration
     * 
     * \param *data: The dataBase object to be cleared
     */
    void SolverBase::addIterationData( dataBase *data ){

        hydra->addIterationData( data );

    }

    /*!
     * Add data to the vector of values which will be cleared after each nonlinear step
     * 
     * \param *data: The dataBase object to be cleared
     */
    void SolverBase::addNLStepData( dataBase *data ){

        hydra->addNLStepData( data );

    }

}
