/**
 ******************************************************************************
 * \file tardigrade_StepDampingBase.cpp
 ******************************************************************************
 * The base class for step damping operations
 ******************************************************************************
 */

#include"tardigrade_StepDampingBase.h"
#include"tardigrade_SolverStepBase.h"
#include"tardigrade_vector_tools.h"

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

    /*!
     * Set the maximum number of line-search iterations
     *
     * \param &value: The incoming value
     */
    void StepDampingBase::setMaxLSIterations( const unsigned int &value ){

        _maxLSIterations = value;

    }

    /*!
     * Get the residual vector
     */
    const floatVector *StepDampingBase::getResidual( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        return step->getResidual( );

    }

    /*!
     * Get the residual norm for the line-search convergence criterion
     */
    const floatType* StepDampingBase::getLSResidualNorm( ){

        if ( !_lsResidualNorm.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( resetLSIteration( ) );

        }

        return &_lsResidualNorm.second;

    }

    /*!
     * Set the line-search residual norm
     */
    void StepDampingBase::setLSResidualNorm( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        _lsResidualNorm.second = tardigradeVectorTools::l2norm( *getResidual( ) );

        _lsResidualNorm.first = true;

    }

    /*!
     * Reset the line search iteration
     */
    void StepDampingBase::resetLSIteration( ){

        _LSIteration = 0;

        _lambda = 1.0;

        setLSResidualNorm( );

    };

// END LINE SEARCH FUNCTIONS

}
