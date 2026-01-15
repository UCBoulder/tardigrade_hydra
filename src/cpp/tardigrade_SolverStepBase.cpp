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
        damping->resetNumLS( ); //TODO: Replace with general reset function
        damping->resetNumGrad( );

    }

    void SolverStepBase::addIterationData( dataBase *data ){
        /*!
         * Add data to the vector of values which will be cleared after each iteration
         * 
         * \param *data: The dataBase object to be cleared
         */

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        solver->addIterationData( data );

    }

    void SolverStepBase::addNLStepData( dataBase *data ){
        /*!
         * Add data to the vector of values which will be cleared after each nonlinear step
         * 
         * \param *data: The dataBase object to be cleared
         */

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
     * Get the unknown vector
     */
    const floatVector *SolverStepBase::getUnknownVector( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        return solver->getUnknownVector( );
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

    void SolverStepBase::setResidualNorm( ){
        /*!
         * Set the norm of the residual vector
         */

        auto residualNorm = get_SetDataStorage_residualNorm( );

        auto residual = getResidual( );

        using residual_type = std::remove_reference_t<decltype( ( *residual )[ 0 ] )>;

        *residualNorm.value = std::inner_product( std::begin( *residual ), std::end( *residual ), std::begin( *residual ), residual_type( ) );

    }

    void SolverStepBase::setdResidualNormdX( ){
        /*!
         * Set the derivative of the residual norm w.r.t. the unknown vector
         */

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        const unsigned int xsize = solver->getNumUnknowns( );

        auto dResidualNormdX = get_SetDataStorage_dResidualNormdX( );

        dResidualNormdX.zero( xsize );

        const floatVector *residual = getResidual( );

        const floatVector *jacobian = solver->getFlatJacobian( );

        for ( unsigned int i = 0; i < xsize; i++ ){
            for ( unsigned int j = 0; j < xsize; j++ ){
                ( *dResidualNormdX.value )[ j ] += 2 * ( *jacobian )[ xsize * i + j ] * ( *residual )[ i ];
            }
        }

    }

    const floatType *SolverStepBase::get_baseResidualNorm( ){
        /*!
         * Get the base value for the residual norm.
         */
        if ( !_baseResidualNorm.first ){

            throw std::runtime_error( "The base residual norm must be set with set_baseResidualNorm before it can be called" );

        }

        return &_baseResidualNorm.second;

    }

    const floatVector *SolverStepBase::get_basedResidualNormdX( ){
        /*!
         * Get the base value for the derivative of the residual norm w.r.t. the unknown vector
         */

        if ( !_basedResidualNormdX.first ){

            throw std::runtime_error( "The base residual norm must be set with set_dbaseResidualNormdX before it can be called" );

        }

        return &_basedResidualNormdX.second;

    }

    void SolverStepBase::set_baseResidualNorm( const floatType &value ){
        /*! Set the base value of the residual norm
         *
         * \param &value: The new value
         */

        setNLStepData( value, _baseResidualNorm );

    }

    void SolverStepBase::set_basedResidualNormdX( const floatVector &value ){
        /*!
         * Set the base derivative of the residual norm w.r.t. the unknown vector
         *
         * \param &value: The new value
         */

        setNLStepData( value, _basedResidualNormdX );

    }

    /*!
     * Set the base quantities prior to updating the unknown vector
     */
    void SolverStepBase::setBaseQuantities( ){

        // TEMP
        set_baseResidualNorm( *get_residualNorm( ) );

        set_basedResidualNormdX( *get_dResidualNormdX( ) );

        if ( _mu_k < 0 ){

            setMuk( 0.5 * getLMMu( ) * ( *get_baseResidualNorm( ) ) );

        }
        else{

            setMuk( std::fmin( _mu_k, ( *get_baseResidualNorm( ) ) ) );

        }
        // END TEMP

        return;

    }

    void SolverStepBase::performPreconditionedSolve( floatVector &deltaX_tr ){
        /*!
         * Perform a pre-conditioned solve
         *
         * \param &deltaX_tr: The trial chcange in the unknown vector
         */

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        TARDIGRADE_ERROR_TOOLS_CHECK( solver->preconditioner != nullptr, "The preconditioner has not been defined" ); //TODO: Move to the trial_step class
        tardigradeVectorTools::solverType< floatType > linearSolver;

        auto dx_map = tardigradeHydra::getDynamicSizeVectorMap( deltaX_tr.data( ), solver->getNumUnknowns( ) );

        auto J_map = tardigradeHydra::getDynamicSizeMatrixMap( getFlatNonlinearLHS( )->data( ), solver->getNumUnknowns( ), solver->getNumUnknowns( ) );

        auto R_map = tardigradeHydra::getDynamicSizeVectorMap( getNonlinearRHS( )->data( ), solver->getNumUnknowns( ) );

        if( solver->preconditioner->getPreconditionerIsDiagonal( ) ){

            auto p_map = tardigradeHydra::getDynamicSizeVectorMap( solver->preconditioner->getFlatPreconditioner( )->data( ), solver->getNumUnknowns( ) );

            linearSolver = tardigradeVectorTools::solverType< floatType >( p_map.asDiagonal( ) * J_map );

            dx_map = -linearSolver.solve( p_map.asDiagonal( ) * R_map );

        }
        else{

            auto p_map = tardigradeHydra::getDynamicSizeMatrixMap( solver->preconditioner->getFlatPreconditioner( )->data( ), solver->getNumUnknowns( ), solver->getNumUnknowns( ) );

            linearSolver = tardigradeVectorTools::solverType< floatType >( p_map * J_map );

            dx_map = -linearSolver.solve( p_map * R_map );

        }

        unsigned int rank = linearSolver.rank( );

        if ( solver->getRankDeficientError( ) && ( rank != getResidual( )->size( ) ) ){

            TARDIGRADE_ERROR_TOOLS_CATCH( throw convergence_error( "The Jacobian is not full rank" ) );

        }

    }

    /*!
     * Increment the solution of the problem
     */
    void SolverStepBase::incrementSolution( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        TARDIGRADE_ERROR_TOOLS_CHECK( trial_step != nullptr, "The trial step has not been defined" );
        TARDIGRADE_ERROR_TOOLS_CHECK( damping != nullptr, "The damping has not been defined" );
        if ( getFailureVerbosityLevel( ) > 0 ){
            addToFailureOutput( "\n\n  iteration: " );
            addToFailureOutput( getIteration( ) );
        }

        X0 = *getUnknownVector( );
        deltaX = floatVector( solver->getNumUnknowns( ), 0 );

        if ( getFailureVerbosityLevel( ) > 0 ){
            addToFailureOutput( "  X0:\n" );
            addToFailureOutput( "  " );
            addToFailureOutput( *getUnknownVector( ) );
        }

        setBaseQuantities( );

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

        solver->updateUnknownVector( X0 + damping->getLambda( ) * deltaX );

        damping->applyDamping( );

        // Refine the estimate if the new point has a higher residual
        if ( !damping->checkLSConvergence( ) ){

            if ( checkDescentDirection( deltaX ) || !damping->getUseGradientDescent( ) ){

                // Perform an Armijo type line search when the search direction is aligned with the gradient
                damping->performArmijoTypeLineSearch( X0, deltaX );

            }
            else{

                // Perform gradient descent if the search direction is not aligned with the gradient
                performGradientStep( X0 );

            }

        }
        else{

            solver->resetToleranceScaleFactor( );

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

    void SolverStepBase::solveNewtonUpdate( floatVector &deltaX_tr ){
        /*!
         * Solve the Newton update returning the trial value of the unknown vector
         *
         * \param &deltaX_tr: The trial change in the unknown vector
         */

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        TARDIGRADE_ERROR_TOOLS_CHECK( solver->preconditioner != nullptr, "The preconditioner has not been defined" ); //TODO: Move to the trial_step class
        if ( solver->preconditioner->getUsePreconditioner( ) ){

            performPreconditionedSolve( deltaX_tr );

        }
        else{

            auto dx_map = tardigradeHydra::getDynamicSizeVectorMap( deltaX_tr.data( ), solver->getNumUnknowns( ) );

            auto J_map = tardigradeHydra::getDynamicSizeMatrixMap( getFlatNonlinearLHS( )->data( ), solver->getNumUnknowns( ), solver->getNumUnknowns( ) );

            auto R_map = tardigradeHydra::getDynamicSizeVectorMap( getNonlinearRHS( )->data( ), solver->getNumUnknowns( ) );

            tardigradeVectorTools::solverType< floatType > linearSolver( J_map );
            dx_map = -linearSolver.solve( R_map );

            unsigned int rank = linearSolver.rank( );

            if ( solver->getRankDeficientError( ) && ( rank != getResidual( )->size( ) ) ){

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

    /*!
     * Set the value of the mu_k parameter for Levenberg-Marquardt steps
     *
     * \param &value: The value of the parameter
     */
    void SolverStepBase::setMuk( const floatType &value ){
 
        _mu_k = value;

    }

    /*!
     * Set the value of the mu parameter for Levenberg-Marquardt steps
     *
     * \param &value: The value of the parameter
     */
    void SolverStepBase::setLMMu( const floatType &value ){
 
        _lm_mu = value;

    }

// END LM FUNCTIONS

// BEGIN SQP SOLVER FUNCTIONS

    void SolverStepBase::initializeActiveConstraints( std::vector< bool > &active_constraints ){
        /*!
         * Initialize the active constraint vector
         * 
         * \param &active_constraints: The current constraints that are active
         */

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        active_constraints = std::vector< bool >( solver->getNumConstraints( ), false );

        for ( auto c = solver->getConstraints( )->begin( ); c != solver->getConstraints( )->end( ); c++ ){

            unsigned int index = ( unsigned int )( c - solver->getConstraints( )->begin( ) );

            active_constraints[ index ] = ( ( *c ) < 0. );

        }

    }

    void SolverStepBase::assembleKKTRHSVector( const floatVector &dx, floatVector &KKTRHSVector, const std::vector< bool > &active_constraints ){
        /*!
         * Assemble the right hand side vector for the KKT matrix
         * 
         * \param &dx: The delta vector being solved for
         * \param &KKTRHSVector: The right hand size vector for the KKT matrix
         * \param &active_constraints: The active constraint vector
         */

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        const unsigned int numUnknowns = solver->getNumUnknowns( );

        const unsigned int numConstraints = solver->getNumConstraints( );

        KKTRHSVector = floatVector( numUnknowns + numConstraints, 0 );

        Eigen::Map< const Eigen::Vector< floatType, -1 > > _dx( dx.data( ), numUnknowns );

        Eigen::Map< Eigen::Vector< floatType, -1 > > RHS( KKTRHSVector.data( ), ( numUnknowns + numConstraints ), ( numUnknowns + numConstraints ) );

        Eigen::Map< const Eigen::Vector< floatType, -1 > > R( getResidual( )->data( ), numUnknowns );

        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > J( solver->getFlatJacobian( )->data( ), numUnknowns, numUnknowns );

        RHS.head( numUnknowns ) = ( J.transpose( ) * ( R + J * _dx ) + getMuk( ) * _dx ).eval( );

        for ( unsigned int i = 0; i < numConstraints; i++ ){

            if ( active_constraints[ i ] ){

                KKTRHSVector[ numUnknowns + i ] = ( *( solver->getConstraints( ) ) )[ i ];

                for ( unsigned int I = 0; I < numUnknowns; I++ ){

                    KKTRHSVector[ numUnknowns + i ] += ( *( solver->getConstraintJacobians( ) ) )[ numUnknowns * i + I ] * dx[ I ];

                }

            }

        }

    }

    void SolverStepBase::assembleKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints ){
        /*!
         * Assemble the Karush-Kuhn-Tucker matrix for an inequality constrained Newton-Raphson solve
         * 
         * \param &KKTMatrix: The Karush-Kuhn-Tucker matrix
         * \param &active_constraints: The vector of currently active constraints.
         */

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        const unsigned int numUnknowns = solver->getNumUnknowns( );

        const unsigned int numConstraints = solver->getNumConstraints( );

        KKTMatrix = floatVector( ( numUnknowns + numConstraints ) * ( numUnknowns + numConstraints ), 0 );

        Eigen::Map< Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > K( KKTMatrix.data( ), ( numUnknowns + numConstraints ), ( numUnknowns + numConstraints ) );

        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > J( solver->getFlatJacobian( )->data( ), numUnknowns, numUnknowns );

        K.block( 0, 0, numUnknowns, numUnknowns ) = ( J.transpose( ) * J ).eval( );

        for ( unsigned int I = 0; I < numUnknowns; I++ ){

            KKTMatrix[ ( numUnknowns + numConstraints ) * I + I ] += getMuk( );

        }

        for ( unsigned int i = 0; i < numConstraints; i++ ){

            if ( active_constraints[ i ] ){

                for ( unsigned int I = 0; I < numUnknowns; I++ ){

                    KKTMatrix[ ( numUnknowns + numConstraints ) * ( I ) + numUnknowns + i ] = ( *solver->getConstraintJacobians( ) )[ numUnknowns * i + I ];
                    KKTMatrix[ ( numUnknowns + numConstraints ) * ( numUnknowns + i ) + I ] = ( *solver->getConstraintJacobians( ) )[ numUnknowns * i + I ];

                }

            }
            else{

                KKTMatrix[ ( numUnknowns + numConstraints ) * ( numUnknowns + i ) + numUnknowns + i ] = 1;

            }

        }

    }

    void SolverStepBase::updateKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints ){
        /*!
         * Update the KKTMatrix if the active constraints have changed
         * 
         * \param &KKTMatrix: The Karush-Kuhn-Tucker matrix
         * \param &active_constraints: The vector of currently active constraints.
         */

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        const unsigned int numUnknowns = solver->getNumUnknowns( );

        const unsigned int numConstraints = solver->getNumConstraints( );

        Eigen::Map< Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > K( KKTMatrix.data( ), ( numUnknowns + numConstraints ), ( numUnknowns + numConstraints ) );

        K.block( 0, numUnknowns, numUnknowns, numConstraints ).setZero( );

        K.block( numUnknowns, 0, numConstraints, numUnknowns ).setZero( );

        K.block( numUnknowns, numUnknowns, numConstraints, numConstraints ).setZero( );

        for ( unsigned int i = 0; i < numConstraints; i++ ){

            if ( active_constraints[ i ] ){

                for ( unsigned int I = 0; I < numUnknowns; I++ ){

                    KKTMatrix[ ( numUnknowns + numConstraints ) * ( I ) + numUnknowns + i ] = ( *solver->getConstraintJacobians( ) )[ numUnknowns * i + I ];
                    KKTMatrix[ ( numUnknowns + numConstraints ) * ( numUnknowns + i ) + I ] = ( *solver->getConstraintJacobians( ) )[ numUnknowns * i + I ];

                }

            }
            else{

                KKTMatrix[ ( numUnknowns + numConstraints ) * ( numUnknowns + i ) + numUnknowns + i ] = 1;

            }

        }

    }

    void SolverStepBase::solveConstrainedQP( floatVector &dx, const unsigned int kmax ){
        /*!
         * Solve the constrained QP problem to estimate the desired step size
         * 
         * \param &dx: The change in the unknown vector
         * \param kmax: The maximum number of iterations (defaults to 100)
         */

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        const unsigned int numUnknowns = solver->getNumUnknowns( );

        const unsigned int numConstraints = solver->getNumConstraints( );

        floatVector K;

        floatVector RHS;

        std::vector< bool > active_constraints;
        initializeActiveConstraints( active_constraints );

        assembleKKTRHSVector( dx, RHS, active_constraints );

        assembleKKTMatrix( K, active_constraints );

        floatType tol = solver->getRelativeTolerance( ) * ( tardigradeVectorTools::l2norm( RHS ) ) + solver->getAbsoluteTolerance( );

        unsigned int k = 0;

        floatVector y( numUnknowns + numConstraints, 0 );

        floatVector ck = *solver->getConstraints( );

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

                ck     = *solver->getConstraints( );
                ctilde = *solver->getConstraints( );
                for ( unsigned int i = 0; i < numConstraints; i++ ){
                    for ( unsigned int j = 0; j < numUnknowns; j++ ){
                        ck[ i ]     += ( *solver->getConstraintJacobians( ) )[ numUnknowns * i + j ] * dx[ j ];
                        ctilde[ i ] += ( *solver->getConstraintJacobians( ) )[ numUnknowns * i + j ] * ( dx[ j ] - negp[ j ] );
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
     * Check if the search direction is a descent direction of the Jacobian
     * 
     * \param &dx: The proposed change in x
     */
    bool SolverStepBase::checkDescentDirection( const floatVector &dx ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        const unsigned int xsize = solver->getNumUnknowns( );

        const floatType RHS = -damping->getGradientRho( ) * std::pow( tardigradeVectorTools::l2norm( dx ), damping->getGradientP( ) );

        floatType LHS = 0;

        const floatVector *dResidualNormdX = get_basedResidualNormdX( );

        for ( unsigned int i = 0; i < xsize; i++ ){

            LHS += ( *dResidualNormdX )[ i ] * dx[ i ];

        }

        return LHS <= RHS;

    }

    /*!
     * Check the convergence of a gradient step
     *
     * \param &X0: The initial value of the unknown vector
     */
    bool SolverStepBase::checkGradientConvergence( const floatVector &X0 ){

        TARDIGRADE_ERROR_TOOLS_CHECK( solver != nullptr, "The solver has not been defined" );
        const unsigned int xsize = solver->getNumUnknowns( );

        floatVector dx = ( *getUnknownVector( ) ) - X0;

        floatType RHS = *get_baseResidualNorm( );

        for ( unsigned int i = 0; i < xsize; ++i ){

            RHS += damping->getGradientSigma( ) * ( *( get_basedResidualNormdX( ) ) )[ i ] * dx[ i ];

        }

        return ( *get_residualNorm( ) ) < solver->getToleranceScaleFactor( ) * RHS;

    }

    /*!
     * Perform a gradient descent step
     *
     * \param &X0: The base value of the unknown vector
     */
    void SolverStepBase::performGradientStep( const floatVector &X0 ){

        TARDIGRADE_ERROR_TOOLS_CHECK( damping != nullptr, "The damping has not been defined" );
        const floatVector *dResidualNormdX = get_basedResidualNormdX( );

        unsigned int l                     = 0;

        const unsigned int maxiter         = damping->getMaxGradientIterations( );

        while( damping->checkGradientIteration( ) ){

            floatType t = std::pow( damping->getGradientBeta( ), l );

            updateUnknownVector( X0 - t * ( *dResidualNormdX ) );

            if ( checkGradientConvergence( X0 ) ){

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
