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

}
