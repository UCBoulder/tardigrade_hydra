/**
 ******************************************************************************
 * \file tardigrade_TrialStepBase.cpp
 ******************************************************************************
 * The base class for step damping operations
 ******************************************************************************
 */

#include"tardigrade_TrialStepBase.h"
#include"tardigrade_SolverStepBase.h"

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

    /*!
     * Get the residual vector
     */
    const floatVector *TrialStepBase::getResidual( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        return step->getResidual( );

    }

    /*!
     * Get the number of unknowns
     */
    const unsigned int TrialStepBase::getNumUnknowns( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        return step->getNumUnknowns( );

    }

    /*!
     * Get the Jacobian in row-major format
     */
    const floatVector *TrialStepBase::getFlatJacobian( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        return step->getFlatJacobian( );

    }
    /*!
     * Get the number of constraint equations
     */
    const unsigned int TrialStepBase::getNumConstraints( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        return step->getNumConstraints( );
    }

    /*!
     * Get the current constraint values
     */
    const floatVector *TrialStepBase::getConstraints( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        return step->getConstraints( );
    }

    /*!
     * Get the constraint Jacobians
     */
    const floatVector *TrialStepBase::getConstraintJacobians( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        return step->getConstraintJacobians( );
    }

    // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

    /*!
     * Assemble the right hand side vector for the KKT matrix
     * 
     * \param &dx: The delta vector being solved for
     * \param &KKTRHSVector: The right hand size vector for the KKT matrix
     * \param &active_constraints: The active constraint vector
     */
    void TrialStepBase::assembleKKTRHSVector( const floatVector &dx, floatVector &KKTRHSVector, const std::vector< bool > &active_constraints ){

        const unsigned int numUnknowns = getNumUnknowns( );

        const unsigned int numConstraints = getNumConstraints( );

        KKTRHSVector = floatVector( numUnknowns + numConstraints, 0 );

        Eigen::Map< const Eigen::Vector< floatType, -1 > > _dx( dx.data( ), numUnknowns );

        Eigen::Map< Eigen::Vector< floatType, -1 > > RHS( KKTRHSVector.data( ), ( numUnknowns + numConstraints ), ( numUnknowns + numConstraints ) );

        Eigen::Map< const Eigen::Vector< floatType, -1 > > R( getResidual( )->data( ), numUnknowns );

        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > J( getFlatJacobian( )->data( ), numUnknowns, numUnknowns );

        RHS.head( numUnknowns ) = ( J.transpose( ) * ( R + J * _dx ) + step->damping->getMuk( ) * _dx ).eval( );

        for ( unsigned int i = 0; i < numConstraints; i++ ){

            if ( active_constraints[ i ] ){

                KKTRHSVector[ numUnknowns + i ] = ( *( getConstraints( ) ) )[ i ];

                for ( unsigned int I = 0; I < numUnknowns; I++ ){

                    KKTRHSVector[ numUnknowns + i ] += ( *( getConstraintJacobians( ) ) )[ numUnknowns * i + I ] * dx[ I ];

                }

            }

        }

    }

    // END SQP SOLVER FUNCTIONS

}
