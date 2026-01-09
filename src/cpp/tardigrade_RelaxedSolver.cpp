/**
  ******************************************************************************
  * \file tardigrade_RelaxedSolver.cpp
  ******************************************************************************
  * A C++ library for the nonlinear solvers which attempt to relax the problem
  * during its solution
  ******************************************************************************
  */

#include"tardigrade_RelaxedSolver.h"
#include"tardigrade_hydra.h"

namespace tardigradeHydra{

    /*!
     * Reset the solver
     */
    void RelaxedSolver::reset( ){

        resetRelaxedIteration( );
        internal_solver->reset( );
        tardigradeHydra::SolverBase::reset( );

    }

    /*!
     * Get the current relaxed iteration
     */
    const unsigned int RelaxedSolver::getRelaxedIteration( ){

        return _relaxedIteration;
    }

    /*!
     * Get the maximum number of relaxed iterations
     */
    const unsigned int RelaxedSolver::getMaxRelaxedIterations( ){

        return _maxRelaxedIterations;

    }

    /*!
     * Set the maximum allowable number of relaxed iterations
     * 
     * \param &value: The number of relaxed iterations
     */
    const void RelaxedSolver::setMaxRelaxedIterations( const unsigned int &value ){

        _maxRelaxedIterations = value;

    }

    /*!
     * Set the relaxed iteration number
     *
     * \param &value: The incoming value
     */
    void RelaxedSolver::setRelaxedIteration( const unsigned int &value ){

        _relaxedIteration = value;

    }

    /*!
     * Reset the relaxed iteration number
     */
    void RelaxedSolver::resetRelaxedIteration( ){

        setRelaxedIteration( 0 );

    }

    /*!
     * Increment the relaxed iteration number
     */
    void RelaxedSolver::incrementRelaxedIteration( ){

        _relaxedIteration++;

    }

    /*!
     * Initialize the residuals for a relaxed solve
     */
    void RelaxedSolver::initializeResiduals( ){

        hydra->setCurrentResidualIndexMeaningful( true );

        for ( auto residual = std::begin( *( hydra->getResidualClasses( ) ) ); residual != std::end( *( hydra->getResidualClasses( ) ) ); ++residual ){
            hydra->setCurrentResidualIndex( residual - std::begin( *( hydra->getResidualClasses( ) ) ) );

            // Prepare the residuals to take a relaxed step
            ( *residual )->setupRelaxedStep( getRelaxedIteration( ) );

        }

        hydra->setCurrentResidualIndexMeaningful( false );
    }

    /*!
     * Check if the relaxation iterations have converged
     */
    bool RelaxedSolver::checkRelaxedConvergence( ){

        bool relaxedConverged = true;

        hydra->setCurrentResidualIndexMeaningful( true );

        for ( auto residual = std::begin( *( hydra->getResidualClasses( ) ) ); residual != std::end( *( hydra->getResidualClasses( ) ) ); ++residual ){
            hydra->setCurrentResidualIndex( residual - std::begin( *( hydra->getResidualClasses( ) ) ) );

            if ( !( *residual )->checkRelaxedConvergence( ) ){

                relaxedConverged = false;
                break;

            }

        }

        hydra->setCurrentResidualIndexMeaningful( false );

        return relaxedConverged;

    }

    /*!
     * Signal to the residuals that we have a failed relaxed solve step and
     * determine if a new relaxed step should be taken
     */
    bool RelaxedSolver::callResidualRelaxedStepFailure( ){

        bool attempt_relaxed_step = false;

        hydra->setCurrentResidualIndexMeaningful( true );

        for ( auto residual_ptr = std::begin( *( hydra->getResidualClasses( ) ) ); residual_ptr != std::end( *( hydra->getResidualClasses( ) ) ); ++residual_ptr ){

            hydra->setCurrentResidualIndex( residual_ptr - std::begin( *( hydra->getResidualClasses( ) ) ) );

            try{

                auto val = ( *residual_ptr )->relaxedStepFailure( );

                attempt_relaxed_step = attempt_relaxed_step || val;

            }
            catch( std::exception &e ){

                if ( getFailureVerbosityLevel( ) > 0 ){

                    addToFailureOutput( "Failure in residual " + std::to_string( residual_ptr - std::begin( *( hydra->getResidualClasses( ) ) ) ) + "\n" );
                    std::string message;
                    tardigradeErrorTools::captureNestedExceptions( e, message );
                    addToFailureOutput( message );

                }

                throw;

            }

        }

        hydra->setCurrentResidualIndexMeaningful( false );

        return attempt_relaxed_step;

    }

    /*!
     * Set the internal solver which will attempt to solve the relaxed problem
     *
     * \param *_solver: The solver to relax
     */
    void RelaxedSolver::setInternalSolver( SolverBase *_solver ){

        internal_solver = _solver;

    }

    /*!
     * Attempt to perform a solve of the non-linear problem
     */
    bool RelaxedSolver::attemptInternalSolve( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( internal_solver != nullptr, "The solver which is to be relaxed (i.e., the internal solver) has not been defined" );

        try{

            internal_solver->solve( );

            // Exit if the relaxed solver has converged
            return checkRelaxedConvergence( );

        }
        catch( convergence_error &e ){

            if ( !callResidualRelaxedStepFailure( ) ){

                throw;

            }

            return false;

        }
        catch( std::exception &e ){

            throw;

        }

    }

}
