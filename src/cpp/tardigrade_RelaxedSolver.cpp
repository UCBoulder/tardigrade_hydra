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
     * Get the current relaxed iteration
     */
    const unsigned int RelaxedSolver::getRelaxedIteration( ){

        return _relaxedIteration;
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

}
