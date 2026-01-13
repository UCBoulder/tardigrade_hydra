/**
 ******************************************************************************
 * \file tardigrade_StepDampingBase.cpp
 ******************************************************************************
 * The base class for step damping operations
 ******************************************************************************
 */

#include"tardigrade_StepDampingBase.h"

namespace tardigradeHydra{

    /*!
     * Apply the damping to the proposed step
     */
    void StepDampingBase::applyDamping( ){

    }


// LINE SEARCH FUNCTIONS


    /*!
     * Set the value of the line-search alpha parameter
     *
     * \param &value: The incoming value of the line-search alpha parameter
     */
    void StepDampingBase::setLSAlpha( const floatType &value ){

        _lsAlpha = value;

    }


// END LINE SEARCH FUNCTIONS

}
