/**
 ******************************************************************************
 * \file tardigrade_SolverStepBase.cpp
 ******************************************************************************
 * The base class for solver steps
 ******************************************************************************
 */

#include"tardigrade_SolverStepBase.h"
#define USE_EIGEN
#include"tardigrade_vector_tools.h"
#include"tardigrade_CustomErrors.h"
#include"tardigrade_SolverBase.h"
#include"tardigrade_StepDampingBase.h"
#include"tardigrade_hydra.h"

namespace tardigradeHydra{

    /*!
     * Initialize the default damping and trial step classes
     */
    void SolverStepBase::initializeDefaults( ){

        damping = &_damping;
        damping->step = this;

        trial_step = &_trial_step;
        trial_step->step = this;
    }

    /*!
     * Reset the step back to an initial state
     */
    void SolverStepBase::reset( ){

        resetNumNewton( );
        damping->resetCounts( );

    }

    /*!
     * Get the relative tolerance value
     */
    const floatType SolverStepBase::getRelativeTolerance( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getRelativeTolerance( );

    }

    /*!
     * Get the absolute tolerance value
     */
    const floatType SolverStepBase::getAbsoluteTolerance( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getAbsoluteTolerance( );

    }

    /*!
     * Add data to the vector of values which will be cleared after each iteration
     * 
     * \param *data: The dataBase object to be cleared
     */
    void SolverStepBase::addIterationData( dataBase *data ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        solver->addIterationData( data );

    }

    /*!
     * Add data to the vector of values which will be cleared after each nonlinear step
     * 
     * \param *data: The dataBase object to be cleared
     */
    void SolverStepBase::addNLStepData( dataBase *data ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        solver->addNLStepData( data );

    }

    /*!
     * Get the residual vector
     */
    const unsigned int SolverStepBase::getIteration( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getIteration( );

    }

    /*!
     * Get the residual vector
     */
    const floatVector *SolverStepBase::getResidual( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getResidual( );

    }

    /*!
     * Update the unknown vector
     *
     * \param &value: The new value of the unknown vector
     */
    void SolverStepBase::updateUnknownVector( const floatVector &value ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        solver->updateUnknownVector( value );

    }

    /*!
     * Get the number of unknowns
     */
    const unsigned int SolverStepBase::getNumUnknowns( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getNumUnknowns( );

    }

    /*!
     * Get the unknown vector
     */
    const floatVector *SolverStepBase::getUnknownVector( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getUnknownVector( );
    }

    /*!
     * Get the Jacobian in row-major format
     */
    const floatVector *SolverStepBase::getFlatJacobian( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getFlatJacobian( );

    }

    /*!
     * Get the number of constraint equations
     */
    const unsigned int SolverStepBase::getNumConstraints( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getNumConstraints( );
    }

    /*!
     * Get the current constraint values
     */
    const floatVector *SolverStepBase::getConstraints( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getConstraints( );
    }

    /*!
     * Get the constraint Jacobians
     */
    const floatVector *SolverStepBase::getConstraintJacobians( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getConstraintJacobians( );
    }

    /*!
     * Get the scale factor for the tolerance
     */
    const floatType SolverStepBase::getToleranceScaleFactor( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getToleranceScaleFactor( );

    }

    /*!
     * Reset the tolerance scale factor
     */
    void SolverStepBase::resetToleranceScaleFactor( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        solver->resetToleranceScaleFactor( );

    }

    /*!
     * Get whether a rank-deficient matrix will throw an error
     */
    bool SolverStepBase::getRankDeficientError( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getRankDeficientError( );

    }

    /*!
     * Get the failure verbosity level
     */
    const unsigned int SolverStepBase::getFailureVerbosityLevel( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getFailureVerbosityLevel( );

    }

    /*!
     * Add the string to the failure output message
     *
     * \param &string: The string to add to the failure output message
     */
    void SolverStepBase::addToFailureOutput( const std::string &string ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        solver->addToFailureOutput( string );

    }

    /*!
     * Add a floatVector to the failure output message
     *
     * \param &value: The floatVector to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void SolverStepBase::addToFailureOutput( const floatVector &value, bool add_endline ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        solver->addToFailureOutput( value, add_endline );

    }

    /*!
     * Add a vector of booleans to the failure output message
     *
     * \param &value: The vector of booleans to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void SolverStepBase::addToFailureOutput( const std::vector<bool> &value, bool add_endline ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        solver->addToFailureOutput( value, add_endline );

    }

    /*!
     * Add a floatType to the failure output message
     *
     * \param &value: The floatType to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void SolverStepBase::addToFailureOutput( const floatType &value, bool add_endline ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        solver->addToFailureOutput( value, add_endline );

    }

    /*!
     * Perform a pre-conditioned solve
     *
     * \param &deltaX_tr: The trial chcange in the unknown vector
     */
    void SolverStepBase::performPreconditionedSolve( floatVector &deltaX_tr ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        TARDIGRADE_ERROR_TOOLS_CHECK( solver->preconditioner != nullptr, "The preconditioner has not been defined" ); //TODO: Move to the trial_step class
        tardigradeVectorTools::solverType< floatType > linearSolver;

        auto dx_map = tardigradeHydra::getDynamicSizeVectorMap( deltaX_tr.data( ), getNumUnknowns( ) );

        auto J_map = tardigradeHydra::getDynamicSizeMatrixMap( getFlatNonlinearLHS( )->data( ), getNumUnknowns( ), getNumUnknowns( ) );

        auto R_map = tardigradeHydra::getDynamicSizeVectorMap( getNonlinearRHS( )->data( ), getNumUnknowns( ) );

        if( solver->preconditioner->getPreconditionerIsDiagonal( ) ){

            auto p_map = tardigradeHydra::getDynamicSizeVectorMap( solver->preconditioner->getFlatPreconditioner( )->data( ), getNumUnknowns( ) );

            linearSolver = tardigradeVectorTools::solverType< floatType >( p_map.asDiagonal( ) * J_map );

            dx_map = -linearSolver.solve( p_map.asDiagonal( ) * R_map );

        }
        else{

            auto p_map = tardigradeHydra::getDynamicSizeMatrixMap( solver->preconditioner->getFlatPreconditioner( )->data( ), getNumUnknowns( ), getNumUnknowns( ) );

            linearSolver = tardigradeVectorTools::solverType< floatType >( p_map * J_map );

            dx_map = -linearSolver.solve( p_map * R_map );

        }

        unsigned int rank = linearSolver.rank( );

        if ( getRankDeficientError( ) && ( rank != getResidual( )->size( ) ) ){

            TARDIGRADE_ERROR_TOOLS_CATCH( throw convergence_error( "The Jacobian is not full rank" ) );

        }

    }

    /*!
     * Increment the solution of the problem
     */
    void SolverStepBase::incrementSolution( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( trial_step != nullptr, "The trial step has not been defined" );
        TARDIGRADE_ERROR_TOOLS_CHECK( damping != nullptr, "The damping has not been defined" );
        if ( getFailureVerbosityLevel( ) > 0 ){
            addToFailureOutput( "\n\n  iteration: " );
            addToFailureOutput( getIteration( ) );
        }

        X0 = *getUnknownVector( );
        deltaX = floatVector( getNumUnknowns( ), 0 );

        if ( getFailureVerbosityLevel( ) > 0 ){
            addToFailureOutput( "  X0:\n" );
            addToFailureOutput( "  " );
            addToFailureOutput( *getUnknownVector( ) );
        }

        damping->setBaseQuantities( );

        trial_step->computeTrial( );

        if ( getUseSQPSolver( ) ){

            solveConstrainedQP( deltaX );

        }
        else{

            solveNewtonUpdate( deltaX );

        }

        if ( getFailureVerbosityLevel( ) > 0 ){
            addToFailureOutput( "  trial deltaX:\n" );
            addToFailureOutput( "  " );
            addToFailureOutput( deltaX );
        }

        updateUnknownVector( X0 + damping->getLambda( ) * deltaX );

        damping->applyDamping( );

        // Refine the estimate if the new point has a higher residual
        if ( !damping->checkLSConvergence( ) ){

            if ( damping->checkDescentDirection( deltaX ) || !damping->getUseGradientDescent( ) ){

                // Perform an Armijo type line search when the search direction is aligned with the gradient
                damping->performArmijoTypeLineSearch( X0, deltaX );

            }
            else{

                // Perform gradient descent if the search direction is not aligned with the gradient
                performGradientStep( X0 );

            }

        }
        else{

            resetToleranceScaleFactor( );

            incrementNumNewton( );

        }

    }
// BEGIN NONLINEAR SOLVER FUNCTIONS

    /*!
     * Get the RHS vector for the non-linear problem
     */
    const floatVector* SolverStepBase::getNonlinearRHS( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return getResidual( );

    }

    /*!
     * Get the flat LHS matrix for the non-linear problem
     */
    const floatVector* SolverStepBase::getFlatNonlinearLHS( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getFlatJacobian( );

    }

// END NONLINEAR SOLVER FUNCTIONS

// BEGIN NEWTON SOLVER FUNCTIONS

    /*!
     * Solve the Newton update returning the trial value of the unknown vector
     *
     * \param &deltaX_tr: The trial change in the unknown vector
     */
    void SolverStepBase::solveNewtonUpdate( floatVector &deltaX_tr ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        TARDIGRADE_ERROR_TOOLS_CHECK( solver->preconditioner != nullptr, "The preconditioner has not been defined" ); //TODO: Move to the trial_step class
        if ( solver->preconditioner->getUsePreconditioner( ) ){

            performPreconditionedSolve( deltaX_tr );

        }
        else{

            auto dx_map = tardigradeHydra::getDynamicSizeVectorMap( deltaX_tr.data( ), getNumUnknowns( ) );

            auto J_map = tardigradeHydra::getDynamicSizeMatrixMap( getFlatNonlinearLHS( )->data( ), getNumUnknowns( ), getNumUnknowns( ) );

            auto R_map = tardigradeHydra::getDynamicSizeVectorMap( getNonlinearRHS( )->data( ), getNumUnknowns( ) );

            tardigradeVectorTools::solverType< floatType > linearSolver( J_map );
            dx_map = -linearSolver.solve( R_map );

            unsigned int rank = linearSolver.rank( );

            if ( getRankDeficientError( ) && ( rank != getResidual( )->size( ) ) ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw convergence_error( "The Jacobian is not full rank" ) );

            }

        }

    }

// END NEWTON SOLVER FUNCTIONS

// BEGIN LM FUNCTIONS

    /*!
     * Set whether to attempt a Levenberg-Marquardt step
     * 
     * \param &value: The value of the parameter
     */
    void SolverStepBase::setUseLevenbergMarquardt( const bool &value ){
    
        damping->setUseGradientDescent( value );
    
        _use_LM_step = value;
    
    }

// END LM FUNCTIONS

// BEGIN SQP SOLVER FUNCTIONS

    /*!
     * Initialize the active constraint vector
     * 
     * \param &active_constraints: The current constraints that are active
     */
    void SolverStepBase::initializeActiveConstraints( std::vector< bool > &active_constraints ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
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
    void SolverStepBase::assembleKKTRHSVector( const floatVector &dx, floatVector &KKTRHSVector, const std::vector< bool > &active_constraints ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        const unsigned int numUnknowns = getNumUnknowns( );

        const unsigned int numConstraints = getNumConstraints( );

        KKTRHSVector = floatVector( numUnknowns + numConstraints, 0 );

        Eigen::Map< const Eigen::Vector< floatType, -1 > > _dx( dx.data( ), numUnknowns );

        Eigen::Map< Eigen::Vector< floatType, -1 > > RHS( KKTRHSVector.data( ), ( numUnknowns + numConstraints ), ( numUnknowns + numConstraints ) );

        Eigen::Map< const Eigen::Vector< floatType, -1 > > R( getResidual( )->data( ), numUnknowns );

        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > J( getFlatJacobian( )->data( ), numUnknowns, numUnknowns );

        RHS.head( numUnknowns ) = ( J.transpose( ) * ( R + J * _dx ) + damping->getMuk( ) * _dx ).eval( );

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
    void SolverStepBase::assembleKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        const unsigned int numUnknowns = getNumUnknowns( );

        const unsigned int numConstraints = getNumConstraints( );

        KKTMatrix = floatVector( ( numUnknowns + numConstraints ) * ( numUnknowns + numConstraints ), 0 );

        Eigen::Map< Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > K( KKTMatrix.data( ), ( numUnknowns + numConstraints ), ( numUnknowns + numConstraints ) );

        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > J( getFlatJacobian( )->data( ), numUnknowns, numUnknowns );

        K.block( 0, 0, numUnknowns, numUnknowns ) = ( J.transpose( ) * J ).eval( );

        for ( unsigned int I = 0; I < numUnknowns; I++ ){

            KKTMatrix[ ( numUnknowns + numConstraints ) * I + I ] += damping->getMuk( );

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
    void SolverStepBase::updateKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
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
    void SolverStepBase::solveConstrainedQP( floatVector &dx, const unsigned int kmax ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
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

// GRADIENT DESCENT FUNCTIONS

    /*!
     * Perform a gradient descent step
     *
     * \param &X0: The base value of the unknown vector
     */
    void SolverStepBase::performGradientStep( const floatVector &X0 ){

        TARDIGRADE_ERROR_TOOLS_CHECK( damping != nullptr, "The damping has not been defined" );
        const floatVector *dResidualNormdX = damping->get_basedResidualNormdX( );

        unsigned int l                     = 0;

        const unsigned int maxiter         = damping->getMaxGradientIterations( );

        while( damping->checkGradientIteration( ) ){

            floatType t = std::pow( damping->getGradientBeta( ), l );

            updateUnknownVector( X0 - t * ( *dResidualNormdX ) );

            if ( damping->checkGradientConvergence( X0 ) ){

                break;

            }

            l++;

            damping->incrementGradientIteration( );

        }

        resetToleranceScaleFactor( );

        if ( l >= maxiter ){

            throw convergence_error( "Failure in gradient step" );

        }

        damping->incrementNumGrad( );

        damping->resetGradientIteration( );

    }

// END GRADIENT DESCENT FUNCTIONS

}
