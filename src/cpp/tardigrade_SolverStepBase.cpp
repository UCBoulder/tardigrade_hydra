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

        solver->hydra->addIterationData( data );

    }

    void SolverStepBase::addNLStepData( dataBase *data ){
        /*!
         * Add data to the vector of values which will be cleared after each nonlinear step
         * 
         * \param *data: The dataBase object to be cleared
         */

        solver->hydra->addNLStepData( data );

    }

    void SolverStepBase::setResidualNorm( ){
        /*!
         * Set the norm of the residual vector
         */

        auto residualNorm = get_SetDataStorage_residualNorm( );

        auto residual = solver->hydra->getResidual( );

        using residual_type = std::remove_reference_t<decltype( ( *residual )[ 0 ] )>;

        *residualNorm.value = std::inner_product( std::begin( *residual ), std::end( *residual ), std::begin( *residual ), residual_type( ) );

    }

    void SolverStepBase::setdResidualNormdX( ){
        /*!
         * Set the derivative of the residual norm w.r.t. the unknown vector
         */

        const unsigned int xsize = solver->hydra->getNumUnknowns( );

        auto dResidualNormdX = get_SetDataStorage_dResidualNormdX( );

        dResidualNormdX.zero( xsize );

        const floatVector *residual = solver->hydra->getResidual( );

        const floatVector *jacobian = solver->hydra->getFlatJacobian( );

        for ( unsigned int i = 0; i < xsize; i++ ){
            for ( unsigned int j = 0; j < xsize; j++ ){
                ( *dResidualNormdX.value )[ j ] += 2 * ( *jacobian )[ xsize * i + j ] * ( *residual )[ i ];
            }
        }

    }

}
