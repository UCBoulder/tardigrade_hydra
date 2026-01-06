/**
 ******************************************************************************
 * \file tardigrade_SolverStepBase.cpp
 ******************************************************************************
 * The base class for solver steps
 ******************************************************************************
 */

#include"tardigrade_SolverStepBase.h"
#include"tardigrade_hydra.h"

namespace tardigradeHydra{

    void SolverStepBase::addIterationData( dataBase *data ){
        /*!
         * Add data to the vector of values which will be cleared after each iteration
         * 
         * \param *data: The dataBase object to be cleared
         */

        solver->addIterationData( data );

    }

    void SolverStepBase::addNLStepData( dataBase *data ){
        /*!
         * Add data to the vector of values which will be cleared after each nonlinear step
         * 
         * \param *data: The dataBase object to be cleared
         */

        solver->addNLStepData( data );

    }

    void SolverStepBase::setResidualNorm( ){
        /*!
         * Set the norm of the residual vector
         */

        auto residualNorm = get_SetDataStorage_residualNorm( );

        auto residual = solver->getResidual( );

        using residual_type = std::remove_reference_t<decltype( ( *residual )[ 0 ] )>;

        *residualNorm.value = std::inner_product( std::begin( *residual ), std::end( *residual ), std::begin( *residual ), residual_type( ) );

    }

    void SolverStepBase::setdResidualNormdX( ){
        /*!
         * Set the derivative of the residual norm w.r.t. the unknown vector
         */

        const unsigned int xsize = solver->getNumUnknowns( );

        auto dResidualNormdX = get_SetDataStorage_dResidualNormdX( );

        dResidualNormdX.zero( xsize );

        const floatVector *residual = solver->getResidual( );

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

    void SolverStepBase::setBaseQuantities( ){
        /*!
         * Set the base quantities required for gradient steps
         */

        set_baseResidualNorm( *get_residualNorm( ) );

        set_basedResidualNormdX( *get_dResidualNormdX( ) );

        if ( _mu_k < 0 ){

            setMuk( 0.5 * getLMMu( ) * ( *get_baseResidualNorm( ) ) );

        }
        else{

            setMuk( std::fmin( _mu_k, ( *get_baseResidualNorm( ) ) ) );

        }

    }

// BEGIN NONLINEAR SOLVER FUNCTIONS

    const floatVector* SolverStepBase::getNonlinearRHS( ){
        /*!
         * Get the RHS vector for the non-linear problem
         */

        if ( solver->hydra->getUseLevenbergMarquardt( ) ){

            if ( !_nonlinearRHS.first ){

                const unsigned int xsize = solver->getNumUnknowns( );

                _nonlinearRHS.first = true;

                _nonlinearRHS.second = floatVector( xsize, 0 );

                const floatVector *residual = solver->getResidual( );

                const floatVector *jacobian = solver->getFlatJacobian( );

                auto Jmap = tardigradeHydra::getDynamicSizeMatrixMap( jacobian->data( ), xsize, xsize );

                auto Rmap = tardigradeHydra::getDynamicSizeVectorMap( residual->data( ), xsize );

                auto RHSmap = tardigradeHydra::getDynamicSizeVectorMap( _nonlinearRHS.second.data( ), xsize );

                RHSmap = ( Jmap.transpose( ) * Rmap ).eval( );

                addIterationData( &_nonlinearRHS );

            }

            return &_nonlinearRHS.second;

        }

        return solver->getResidual( );

    }

    const floatVector* SolverStepBase::getFlatNonlinearLHS( ){
        /*!
         * Get the flat LHS matrix for the non-linear problem
         */

        if ( solver->hydra->getUseLevenbergMarquardt( ) ){

            if ( !_flatNonlinearLHS.first ){

                const unsigned int xsize = solver->getNumUnknowns( );

                _flatNonlinearLHS.first = true;

                _flatNonlinearLHS.second = floatVector( xsize * xsize, 0 );

                const floatVector * jacobian = solver->getFlatJacobian( );

                auto Jmap = tardigradeHydra::getDynamicSizeMatrixMap( jacobian->data( ), xsize, xsize );

                auto LHSmap = tardigradeHydra::getDynamicSizeMatrixMap( _flatNonlinearLHS.second.data( ), xsize, xsize );

                LHSmap = ( Jmap.transpose( ) * Jmap ).eval( );

                for ( unsigned int i = 0; i < xsize; i++ ){
                    _flatNonlinearLHS.second[ xsize * i + i ] += getMuk( );
                }

                addIterationData( &_flatNonlinearLHS );

            }

            return &_flatNonlinearLHS.second;

        }

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

        if ( solver->preconditioner->getUsePreconditioner( ) ){

            solver->hydra->performPreconditionedSolve( deltaX_tr );

        }
        else{

            auto dx_map = tardigradeHydra::getDynamicSizeVectorMap( deltaX_tr.data( ), solver->getNumUnknowns( ) );

            auto J_map = tardigradeHydra::getDynamicSizeMatrixMap( getFlatNonlinearLHS( )->data( ), solver->getNumUnknowns( ), solver->getNumUnknowns( ) );

            auto R_map = tardigradeHydra::getDynamicSizeVectorMap( getNonlinearRHS( )->data( ), solver->getNumUnknowns( ) );

            tardigradeVectorTools::solverType< floatType > linearSolver( J_map );
            dx_map = -linearSolver.solve( R_map );

            unsigned int rank = linearSolver.rank( );

            if ( solver->getRankDeficientError( ) && ( rank != solver->getResidual( )->size( ) ) ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw convergence_error( "The Jacobian is not full rank" ) );

            }

        }

    }

// END NEWTON SOLVER FUNCTIONS

// BEGIN SQP SOLVER FUNCTIONS

    void SolverStepBase::initializeActiveConstraints( std::vector< bool > &active_constraints ){
        /*!
         * Initialize the active constraint vector
         * 
         * \param &active_constraints: The current constraints that are active
         */

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

        const unsigned int numUnknowns = solver->getNumUnknowns( );

        const unsigned int numConstraints = solver->getNumConstraints( );

        KKTRHSVector = floatVector( numUnknowns + numConstraints, 0 );

        Eigen::Map< const Eigen::Vector< floatType, -1 > > _dx( dx.data( ), numUnknowns );

        Eigen::Map< Eigen::Vector< floatType, -1 > > RHS( KKTRHSVector.data( ), ( numUnknowns + numConstraints ), ( numUnknowns + numConstraints ) );

        Eigen::Map< const Eigen::Vector< floatType, -1 > > R( solver->getResidual( )->data( ), numUnknowns );

        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > J( solver->getFlatJacobian( )->data( ), numUnknowns, numUnknowns );

        RHS.head( numUnknowns ) = ( J.transpose( ) * ( R + J * _dx ) + solver->step->getMuk( ) * _dx ).eval( );

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

}
