/**
 ******************************************************************************
 * \file tardigrade_NonlinearSolverBase.cpp
 ******************************************************************************
 * The base class for nonlinear solver classes
 ******************************************************************************
 */

#include"tardigrade_NonlinearSolverBase.h"
#include"tardigrade_hydra.h"

namespace tardigradeHydra{

    /*!
     * The function that is called when first attempting to
     * solve the problem
     */
    void NonlinearSolverBase::initialSolveAttempt( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );

        // Reset the internal steps
        step->reset( );

        setRankDeficientError( false );

        // Form the initial unknown vector
        if ( getInitializeUnknownVector( ) ){
            TARDIGRADE_ERROR_TOOLS_CATCH( initializeUnknownVector( ) );
        }

        initial_unknown = *getUnknownVector( );

        floatVector deltaX( getNumUnknowns( ), 0 );

        callResidualPreNLSolve( );

        step->damping->reset( );

        if ( getFailureVerbosityLevel( ) > 0 ){
            addToFailureOutput( "Initial Unknown:\n" );
            addToFailureOutput( *getUnknownVector( ) );
        }

        while( !checkConvergence( ) && checkIteration( ) ){

            step->incrementSolution( );

            // Call residual end of a successful nonlinear step functions
            callResidualSuccessfulNLStep( );

            // Increment the iteration count
            incrementIteration( );

            // Reset the nonlinear step data
            resetNLStepData( );

            if ( getFailureVerbosityLevel( ) > 0 ){
                addToFailureOutput( "  final residual: " );
                addToFailureOutput( tardigradeVectorTools::l2norm( *getResidual( ) ) );
                addToFailureOutput( "\n" );
            }

        }

        if ( !checkConvergence( ) ){

            throw convergence_error( "Failure to converge main loop\n" );

        }

        callResidualPostNLSolve( );

    }

    /*!
     * Signal to the residuals that we are about to start a nonlinear solve
     */
    void NonlinearSolverBase::callResidualPreNLSolve( ){

        setCurrentResidualIndexMeaningful( true );

        for ( auto residual_ptr = std::begin( *getResidualClasses( ) ); residual_ptr != std::end( *getResidualClasses( ) ); ++residual_ptr ){

            setCurrentResidualIndex( residual_ptr - std::begin( *getResidualClasses( ) ) );

            ( *residual_ptr )->preNLSolve( );

        }

        setCurrentResidualIndexMeaningful( false );

    }

    /*!
     * Signal to the residuals that a successful nonlinear step has been performed
     */
    void NonlinearSolverBase::callResidualSuccessfulNLStep( ){

        setAllowModifyGlobalResidual( true );

        setCurrentResidualIndexMeaningful( true );

        for ( auto residual_ptr = std::begin( *getResidualClasses( ) ); residual_ptr != std::end( *getResidualClasses( ) ); ++residual_ptr ){

            setCurrentResidualIndex( residual_ptr - std::begin( *getResidualClasses( ) ) );

            ( *residual_ptr )->successfulNLStep( );

        }

        setCurrentResidualIndexMeaningful( false );

        setAllowModifyGlobalResidual( false );

    }

    /*!
     * Signal to the residuals that we have finished a nonlinear solve
     */
    void NonlinearSolverBase::callResidualPostNLSolve( ){

        setCurrentResidualIndexMeaningful( true );

        for ( auto residual_ptr = std::begin( *getResidualClasses( ) ); residual_ptr != std::end( *getResidualClasses( ) ); ++residual_ptr ){

            setCurrentResidualIndex( residual_ptr - std::begin( *getResidualClasses( ) ) );

            try{

                ( *residual_ptr )->postNLSolve( );

            }
            catch( std::exception &e ){

                if ( getFailureVerbosityLevel( ) > 0 ){

                    addToFailureOutput( "Failure in residual " + std::to_string( residual_ptr - std::begin( *getResidualClasses( ) ) ) + "\n" );
                    std::string message;
                    tardigradeErrorTools::captureNestedExceptions(e, message);
                    addToFailureOutput( message );

                }

                throw;

            }

        }

        setCurrentResidualIndexMeaningful( false );

    }

}
