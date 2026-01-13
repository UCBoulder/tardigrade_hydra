/**
 ******************************************************************************
 * \file tardigrade_TrialStepBase.cpp
 ******************************************************************************
 * The base class for step damping operations
 ******************************************************************************
 */

#include"tardigrade_TrialStepBase.h"

namespace tardigradeHydra{

    /*!
     * Reset the counters
     */
    void TrialStepBase::resetCounts( ){

    }

    /*!
     * Reset the trial step class
     */
    void TrialStepBase::reset( ){

        resetCounts( );

    }

    /*!
     * Compute the trial step
     */
    void TrialStepBase::computeTrial( ){

    }

}
