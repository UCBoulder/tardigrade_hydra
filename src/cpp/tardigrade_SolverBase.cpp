/**
  ******************************************************************************
  * \file tardigrade_SolverBase.cpp
  ******************************************************************************
  * A C++ library for the base classes for solvers
  ******************************************************************************
  */

#include"tardigrade_SolverBase.h"

namespace tardigradeHydra{

    const bool SolverBase::getRankDeficientError( ){
        /*!
         * Get whether the Jacobian being rank-deficient will throw an error
         */

        return _rank_deficient_error;

    }

    void SolverBase::setRankDeficientError( const bool &value ){
        /*!
         * Set whether the Jacobian being rank-deficient will throw an error
         *
         * \param &value: The incoming value
         */

        _rank_deficient_error = value;
    }

}
