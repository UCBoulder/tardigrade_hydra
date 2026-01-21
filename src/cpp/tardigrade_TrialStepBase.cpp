/**
 ******************************************************************************
 * \file tardigrade_TrialStepBase.cpp
 ******************************************************************************
 * The base class for step damping operations
 ******************************************************************************
 */

#include"tardigrade_TrialStepBase.h"
#define USE_EIGEN
#include"tardigrade_vector_tools.h"
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
        TARDIGRADE_ERROR_TOOLS_CHECK( preconditioner != nullptr, "The preconditioner has not been defined" );
        preconditioner->reset( );

    }

    /*!
     * Compute the trial step
     */
    void TrialStepBase::computeTrial( ){

    }

    /*!
     * Add data to the vector of values which will be cleared after each iteration
     * 
     * \param *data: The dataBase object to be cleared
     */
    void TrialStepBase::addIterationData( dataBase *data ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        step->addIterationData( data );

    }

    /*!
     * Add data to the vector of values which will be cleared after each nonlinear step
     * 
     * \param *data: The dataBase object to be cleared
     */
    void TrialStepBase::addNLStepData( dataBase *data ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        step->addNLStepData( data );

    }

    /*!
     * Get the relative tolerance value
     */
    const floatType TrialStepBase::getRelativeTolerance( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        return step->getRelativeTolerance( );

    }

    /*!
     * Get the absolute tolerance value
     */
    const floatType TrialStepBase::getAbsoluteTolerance( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        return step->getAbsoluteTolerance( );

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

    /*!
     * Get the failure verbosity level
     */
    const unsigned int TrialStepBase::getFailureVerbosityLevel( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        return step->getFailureVerbosityLevel( );

    }

    /*!
     * Add the string to the failure output message
     *
     * \param &string: The string to add to the failure output message
     */
    void TrialStepBase::addToFailureOutput( const std::string &string ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        step->addToFailureOutput( string );

    }

    /*!
     * Add a floatVector to the failure output message
     *
     * \param &value: The floatVector to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void TrialStepBase::addToFailureOutput( const floatVector &value, bool add_endline ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        step->addToFailureOutput( value, add_endline );

    }

    /*!
     * Add a vector of booleans to the failure output message
     *
     * \param &value: The vector of booleans to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void TrialStepBase::addToFailureOutput( const std::vector<bool> &value, bool add_endline ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        step->addToFailureOutput( value, add_endline );

    }

    /*!
     * Add a floatType to the failure output message
     *
     * \param &value: The floatType to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void TrialStepBase::addToFailureOutput( const floatType &value, bool add_endline ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        step->addToFailureOutput( value, add_endline );

    }

// BEGIN NONLINEAR SOLVER FUNCTIONS

    /*!
     * Get the RHS vector for the non-linear problem
     */
    const floatVector* TrialStepBase::getNonlinearRHS( ){

        return getResidual( );

    }

    /*!
     * Get the flat LHS matrix for the non-linear problem
     */
    const floatVector* TrialStepBase::getFlatNonlinearLHS( ){

        return getFlatJacobian( );

    }

// END NONLINEAR SOLVER FUNCTIONS


    // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

    /*!
     * Initialize the active constraint vector
     * 
     * \param &active_constraints: The current constraints that are active
     */
    void TrialStepBase::initializeActiveConstraints( std::vector< bool > &active_constraints ){

        active_constraints = std::vector< bool >( getNumConstraints( ), false );

        for ( auto c = getConstraints( )->begin( ); c != getConstraints( )->end( ); c++ ){

            unsigned int index = ( unsigned int )( c - getConstraints( )->begin( ) );

            active_constraints[ index ] = ( ( *c ) < 0. );

        }

    }

    /*!
     * Assemble the right hand side vector for the KKT matrix
     * 
     * \param &dx: The delta vector being solved for
     * \param &KKTRHSVector: The right hand size vector for the KKT matrix
     * \param &active_constraints: The active constraint vector
     */
    void TrialStepBase::assembleKKTRHSVector( const floatVector &dx, floatVector &KKTRHSVector, const std::vector< bool > &active_constraints ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
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

    /*!
     * Assemble the Karush-Kuhn-Tucker matrix for an inequality constrained Newton-Raphson solve
     * 
     * \param &KKTMatrix: The Karush-Kuhn-Tucker matrix
     * \param &active_constraints: The vector of currently active constraints.
     */
    void TrialStepBase::assembleKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        const unsigned int numUnknowns = getNumUnknowns( );

        const unsigned int numConstraints = getNumConstraints( );

        KKTMatrix = floatVector( ( numUnknowns + numConstraints ) * ( numUnknowns + numConstraints ), 0 );

        Eigen::Map< Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > K( KKTMatrix.data( ), ( numUnknowns + numConstraints ), ( numUnknowns + numConstraints ) );

        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > J( getFlatJacobian( )->data( ), numUnknowns, numUnknowns );

        K.block( 0, 0, numUnknowns, numUnknowns ) = ( J.transpose( ) * J ).eval( );

        for ( unsigned int I = 0; I < numUnknowns; I++ ){

            KKTMatrix[ ( numUnknowns + numConstraints ) * I + I ] += step->damping->getMuk( );

        }

        for ( unsigned int i = 0; i < numConstraints; i++ ){

            if ( active_constraints[ i ] ){

                for ( unsigned int I = 0; I < numUnknowns; I++ ){

                    KKTMatrix[ ( numUnknowns + numConstraints ) * ( I ) + numUnknowns + i ] = ( *getConstraintJacobians( ) )[ numUnknowns * i + I ];
                    KKTMatrix[ ( numUnknowns + numConstraints ) * ( numUnknowns + i ) + I ] = ( *getConstraintJacobians( ) )[ numUnknowns * i + I ];

                }

            }
            else{

                KKTMatrix[ ( numUnknowns + numConstraints ) * ( numUnknowns + i ) + numUnknowns + i ] = 1;

            }

        }

    }

    /*!
     * Update the KKTMatrix if the active constraints have changed
     * 
     * \param &KKTMatrix: The Karush-Kuhn-Tucker matrix
     * \param &active_constraints: The vector of currently active constraints.
     */
    void TrialStepBase::updateKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints ){

        const unsigned int numUnknowns = getNumUnknowns( );

        const unsigned int numConstraints = getNumConstraints( );

        Eigen::Map< Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > K( KKTMatrix.data( ), ( numUnknowns + numConstraints ), ( numUnknowns + numConstraints ) );

        K.block( 0, numUnknowns, numUnknowns, numConstraints ).setZero( );

        K.block( numUnknowns, 0, numConstraints, numUnknowns ).setZero( );

        K.block( numUnknowns, numUnknowns, numConstraints, numConstraints ).setZero( );

        for ( unsigned int i = 0; i < numConstraints; i++ ){

            if ( active_constraints[ i ] ){

                for ( unsigned int I = 0; I < numUnknowns; I++ ){

                    KKTMatrix[ ( numUnknowns + numConstraints ) * ( I ) + numUnknowns + i ] = ( *getConstraintJacobians( ) )[ numUnknowns * i + I ];
                    KKTMatrix[ ( numUnknowns + numConstraints ) * ( numUnknowns + i ) + I ] = ( *getConstraintJacobians( ) )[ numUnknowns * i + I ];

                }

            }
            else{

                KKTMatrix[ ( numUnknowns + numConstraints ) * ( numUnknowns + i ) + numUnknowns + i ] = 1;

            }

        }

    }

    /*!
     * Solve the constrained QP problem to estimate the desired step size
     * 
     * \param &dx: The change in the unknown vector
     * \param kmax: The maximum number of iterations (defaults to 100)
     */
    void TrialStepBase::solveConstrainedQP( floatVector &dx, const unsigned int kmax ){

        const unsigned int numUnknowns = getNumUnknowns( );

        const unsigned int numConstraints = getNumConstraints( );

        floatVector K;

        floatVector RHS;

        std::vector< bool > active_constraints;
        initializeActiveConstraints( active_constraints );

        assembleKKTRHSVector( dx, RHS, active_constraints );

        assembleKKTMatrix( K, active_constraints );

        floatType tol = getRelativeTolerance( ) * ( tardigradeVectorTools::l2norm( RHS ) ) + getAbsoluteTolerance( );

        unsigned int k = 0;

        floatVector y( numUnknowns + numConstraints, 0 );

        floatVector ck = *getConstraints( );

        floatVector ctilde( numConstraints, 0 );

        floatVector negp( numUnknowns, 0 );

        floatVector lambda( numConstraints, 0 );

        floatVector P( numUnknowns + numConstraints, 0 );

        for ( unsigned int i = 0; i < ( numUnknowns + numConstraints ); i++ ){

            P[ i ] = 1 / std::max( std::fabs( *std::max_element( K.begin( ) + ( numUnknowns + numConstraints ) * i,
                                                                 K.begin( ) + ( numUnknowns + numConstraints ) * ( i + 1 ),
                                                                 [ ]( const floatType &a, const floatType &b ){ return std::fabs( a ) < std::fabs( b ); } ) ), 1e-15 );

        }

        Eigen::Map< const Eigen::Vector< floatType, -1 > > _P( P.data( ), numUnknowns + numConstraints );

        tardigradeVectorTools::solverType< floatType > linearSolver;

        while ( k < kmax ){

            Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > _K( K.data( ), numConstraints + numUnknowns, numConstraints + numUnknowns );
            Eigen::Map< const Eigen::Vector< floatType, -1 > > _RHS( RHS.data( ), numConstraints + numUnknowns );
            Eigen::Map< Eigen::Vector< floatType, -1 > > _y( y.data( ), numConstraints + numUnknowns );

            linearSolver = tardigradeVectorTools::solverType< floatType >( _P.asDiagonal( ) * _K );

            _y = linearSolver.solve( _P.asDiagonal( ) * _RHS );

            std::copy( y.begin( ), y.begin( ) + numUnknowns, negp.begin( ) );

            std::copy( y.begin( ) + numUnknowns, y.end( ), lambda.begin( ) );

            if ( tardigradeVectorTools::l2norm( negp ) <= tol ){

                bool negLambda = false;

                floatType minLambda = 1;

                unsigned int imin = 0;

                for ( auto v = std::begin( lambda ); v != std::end( lambda ); v++ ){

                    if ( *v < 0 ){

                        negLambda = true;

                        if ( ( *v ) < minLambda ){

                            imin = ( unsigned int )( v - std::begin( lambda ) );

                            minLambda = *v;

                        }

                    }

                }

                if ( negLambda ){

                    active_constraints[ imin ] = false;

                }
                else{

                    return;

                }

            }
            else{

                ck     = *getConstraints( );
                ctilde = *getConstraints( );
                for ( unsigned int i = 0; i < numConstraints; i++ ){
                    for ( unsigned int j = 0; j < numUnknowns; j++ ){
                        ck[ i ]     += ( *getConstraintJacobians( ) )[ numUnknowns * i + j ] * dx[ j ];
                        ctilde[ i ] += ( *getConstraintJacobians( ) )[ numUnknowns * i + j ] * ( dx[ j ] - negp[ j ] );
                    }
                }

                floatType alpha = 1.0;

                unsigned int iblock = 0;

                bool newBlock = false;

                for ( unsigned int i = 0; i < numConstraints; i++ ){

                    if ( !active_constraints[ i ] ){

                        if ( ctilde[ i ] < -tol ){

                            floatType alpha_trial = -ck[ i ] / ( ctilde[ i ] - ck[ i ] );

                            if ( alpha_trial <= alpha ){

                                iblock = i;

                                alpha = alpha_trial;

                                newBlock = true;

                            }

                        }

                    }

                }

                if ( newBlock ){

                    active_constraints[ iblock ] = true;

                }

                dx -= alpha * negp;

            }

            updateKKTMatrix( K, active_constraints );

            assembleKKTRHSVector(  dx, RHS, active_constraints );

            for ( unsigned int i = 0; i < ( numUnknowns + numConstraints ); i++ ){

                P[ i ] = 1 / std::max( std::fabs( *std::max_element( K.begin( ) + ( numUnknowns + numConstraints ) * i,
                                                                     K.begin( ) + ( numUnknowns + numConstraints ) * ( i + 1 ),
                                                                     [ ]( const floatType &a, const floatType &b ){ return std::fabs( a ) < std::fabs( b ); } ) ), 1e-15 );

            }

            k++;

        }

    }

    // END SQP SOLVER FUNCTIONS

}
