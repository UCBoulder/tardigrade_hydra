
/**
  ******************************************************************************
  * \file tardigrade_SolverBase.cpp
  ******************************************************************************
  * A C++ library for the base classes for solvers
  ******************************************************************************
  */

#include"tardigrade_PreconditionerBase.h"
#include"tardigrade_hydra.h"

namespace tardigradeHydra{

    /*!
     * Add data to the vector of values which will be cleared after each iteration
     * 
     * \param *data: The dataBase object to be cleared
     */
    void PreconditionerBase::addIterationData( dataBase *data ){

        solver->addIterationData( data );

    }

    /*!
     * Add data to the vector of values which will be cleared after each nonlinear step
     * 
     * \param *data: The dataBase object to be cleared
     */
    void PreconditionerBase::addNLStepData( dataBase *data ){

        solver->addNLStepData( data );

    }

    /*!
     * Get the flattened row-major preconditioner for the non-linear problem
     */
    const floatVector* PreconditionerBase::getFlatPreconditioner( ){

        if ( !_preconditioner.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( formPreconditioner( ) );

        }
        
        return &_preconditioner.second;

    }

    /*!
     * Form the preconditioner matrix
     */
    void PreconditionerBase::formPreconditioner( ){

        if ( _preconditioner_type == 0 ){

            formMaxRowPreconditioner( );

        }
        else{

            throw std::runtime_error( "Preconditioner type not recognized" );

        }

        _preconditioner.first = true;

        addIterationData( &_preconditioner );

    }

    /*!
     * Form a left preconditioner comprised of the inverse of the maximum value of each row
     */
    void PreconditionerBase::formMaxRowPreconditioner( ){

        const unsigned int problem_size = solver->getNumUnknowns( );

        _preconditioner.second = floatVector( problem_size, 0 );

        // Find the absolute maximum value in each row
        for ( unsigned int i = 0; i < problem_size; i++ ){

            _preconditioner.second[ i ] = 1 / std::max( std::fabs( *std::max_element( solver->hydra->getFlatNonlinearLHS( )->begin( ) + problem_size * i,
                                                                                      solver->hydra->getFlatNonlinearLHS( )->begin( ) + problem_size * ( i + 1 ),
                                                                                      [ ]( const floatType &a, const floatType &b ){ return std::fabs( a ) < std::fabs( b ); } ) ), 1e-15 );

        }

    }

}
