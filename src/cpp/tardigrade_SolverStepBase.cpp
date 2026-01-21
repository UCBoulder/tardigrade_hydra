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

        resetNumUndamped( );
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

        TARDIGRADE_ERROR_TOOLS_CHECK( trial_step != nullptr, "The trial step has not been defined" );
        TARDIGRADE_ERROR_TOOLS_CHECK( trial_step->preconditioner != nullptr, "The preconditioner has not been defined" ); //TODO: Move to the trial_step class
                                                                                                              TARDIGRADE_ERROR_TOOLS_CHECK( trial_step != nullptr, "The trial step has not been defined" );
        tardigradeVectorTools::solverType< floatType > linearSolver;

        auto dx_map = tardigradeHydra::getDynamicSizeVectorMap( deltaX_tr.data( ), getNumUnknowns( ) );

        auto J_map = tardigradeHydra::getDynamicSizeMatrixMap( trial_step->getFlatNonlinearLHS( )->data( ), getNumUnknowns( ), getNumUnknowns( ) );

        auto R_map = tardigradeHydra::getDynamicSizeVectorMap( trial_step->getNonlinearRHS( )->data( ), getNumUnknowns( ) );

        if( trial_step->preconditioner->getPreconditionerIsDiagonal( ) ){

            auto p_map = tardigradeHydra::getDynamicSizeVectorMap( trial_step->preconditioner->getFlatPreconditioner( )->data( ), getNumUnknowns( ) );

            linearSolver = tardigradeVectorTools::solverType< floatType >( p_map.asDiagonal( ) * J_map );

            dx_map = -linearSolver.solve( p_map.asDiagonal( ) * R_map );

        }
        else{

            auto p_map = tardigradeHydra::getDynamicSizeMatrixMap( trial_step->preconditioner->getFlatPreconditioner( )->data( ), getNumUnknowns( ), getNumUnknowns( ) );

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

        if ( trial_step->getUseSQPSolver( ) ){

            trial_step->solveConstrainedQP( deltaX );

        }
        else{

            solveNewtonUpdate( deltaX );

        }

        if ( getFailureVerbosityLevel( ) > 0 ){
            addToFailureOutput( "  trial deltaX:\n" );
            addToFailureOutput( "  " );
            addToFailureOutput( deltaX );
        }

        if( !damping->applyDamping( ) ){ incrementNumUndamped( ); }

    }

// BEGIN NEWTON SOLVER FUNCTIONS

    /*!
     * Solve the Newton update returning the trial value of the unknown vector
     *
     * \param &deltaX_tr: The trial change in the unknown vector
     */
    void SolverStepBase::solveNewtonUpdate( floatVector &deltaX_tr ){

        TARDIGRADE_ERROR_TOOLS_CHECK( trial_step != nullptr, "The trial step has not been defined" );
        TARDIGRADE_ERROR_TOOLS_CHECK( trial_step->preconditioner != nullptr, "The preconditioner has not been defined" ); //TODO: Move to the trial_step class
                                                                                                              TARDIGRADE_ERROR_TOOLS_CHECK( trial_step != nullptr, "The trial step has not been defined" );
        if ( trial_step->preconditioner->getUsePreconditioner( ) ){

            performPreconditionedSolve( deltaX_tr );

        }
        else{

            auto dx_map = tardigradeHydra::getDynamicSizeVectorMap( deltaX_tr.data( ), getNumUnknowns( ) );

            auto J_map = tardigradeHydra::getDynamicSizeMatrixMap( trial_step->getFlatNonlinearLHS( )->data( ), getNumUnknowns( ), getNumUnknowns( ) );

            auto R_map = tardigradeHydra::getDynamicSizeVectorMap( trial_step->getNonlinearRHS( )->data( ), getNumUnknowns( ) );

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

}
