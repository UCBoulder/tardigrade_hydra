/**
 ******************************************************************************
 * \file tardigrade_StepDampingBase.cpp
 ******************************************************************************
 * The base class for step damping operations
 ******************************************************************************
 */

#include"tardigrade_StepDampingBase.h"
#include"tardigrade_SolverStepBase.h"
#include"tardigrade_vector_tools.h"
#include"tardigrade_CustomErrors.h"

namespace tardigradeHydra{

    /*!
     * Apply the damping to the proposed step
     */
    void StepDampingBase::applyDamping( ){

    }

    /*!
     * Reset the internal count
     */
    void StepDampingBase::resetCount( ){

        resetNumLS( );
        resetNumGrad( );
    }

    /*!
     * Set the base quantities prior to updating the unknown vector
     */
    void StepDampingBase::setBaseQuantities( ){

        set_baseResidualNorm( *get_residualNorm( ) );

        set_basedResidualNormdX( *get_dResidualNormdX( ) );

        if ( getMuk( ) < 0 ){

            setMuk( 0.5 * getLMMu( ) * ( *get_baseResidualNorm( ) ) );

        }
        else{

            setMuk( std::fmin( getMuk( ), ( *get_baseResidualNorm( ) ) ) );

        }

        return;

    }

    /*!
     * Add data to the vector of values which will be cleared after each iteration
     * 
     * \param *data: The dataBase object to be cleared
     */
    void StepDampingBase::addIterationData( dataBase *data ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        step->addIterationData( data );

    }

    /*!
     * Add data to the vector of values which will be cleared after each nonlinear step
     * 
     * \param *data: The dataBase object to be cleared
     */
    void StepDampingBase::addNLStepData( dataBase *data ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        step->addNLStepData( data );

    }

    /*!
     * Get the residual vector
     */
    const floatVector *StepDampingBase::getResidual( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        return step->getResidual( );

    }

    /*!
     * Update the unknown vector
     *
     * \param &value: The new value of the unknown vector
     */
    void StepDampingBase::updateUnknownVector( const floatVector &value ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The solver has not been defined" );
        step->updateUnknownVector( value );

    }

    /*!
     * Get the number of unknowns
     */
    const unsigned int StepDampingBase::getNumUnknowns( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The solver has not been defined" );
        return step->getNumUnknowns( );

    }

    /*!
     * Get the unknown vector
     */
    const floatVector *StepDampingBase::getUnknownVector( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The solver has not been defined" );
        return step->getUnknownVector( );
    }

    /*!
     * Get the Jacobian in row-major format
     */
    const floatVector *StepDampingBase::getFlatJacobian( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The solver has not been defined" );
        return step->getFlatJacobian( );

    }

// LINE SEARCH FUNCTIONS


    /*!
     * Set the value of the line-search alpha parameter
     *
     * \param &value: The incoming value of the line-search alpha parameter
     */
    void StepDampingBase::setLSAlpha( const floatType &value ){

        _lsAlpha = value;

    }

    /*!
     * Set the maximum number of line-search iterations
     *
     * \param &value: The incoming value
     */
    void StepDampingBase::setMaxLSIterations( const unsigned int &value ){

        _maxLSIterations = value;

    }

    /*!
     * Get the residual norm for the line-search convergence criterion
     */
    const floatType* StepDampingBase::getLSResidualNorm( ){

        if ( !_lsResidualNorm.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( resetLSIteration( ) );

        }

        return &_lsResidualNorm.second;

    }

    /*!
     * Set the line-search residual norm
     */
    void StepDampingBase::setLSResidualNorm( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        _lsResidualNorm.second = tardigradeVectorTools::l2norm( *getResidual( ) );

        _lsResidualNorm.first = true;

    }

    /*!
     * Reset the line search iteration
     */
    void StepDampingBase::resetLSIteration( ){

        _LSIteration = 0;

        _lambda = 1.0;

        setLSResidualNorm( );

    };

    /*!
     * Check the line-search convergence
     */
    bool StepDampingBase::checkLSConvergence( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        if ( tardigradeVectorTools::l2norm( *getResidual( ) ) < step->getToleranceScaleFactor( ) * ( 1 - getLSAlpha( ) ) * ( *getLSResidualNorm( ) ) ){

            return true;

        }

        return false;

    }

    /*!
     * Check the current line search iteration
     */
    bool StepDampingBase::checkLSIteration( ){
        return getLSIteration( ) < getMaxLSIterations( );
    }

    /*!
     * Perform an Armijo-type line search
     *
     * \param &X0: The base value of the unknown vector
     * \param &deltaX: The proposed change in X
     */
    void StepDampingBase::performArmijoTypeLineSearch( const floatVector &X0, const floatVector &deltaX ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The solver has not been defined" );
        while ( !checkLSConvergence( ) && checkLSIteration( ) ){

            if ( step->getFailureVerbosityLevel( ) > 0 ){
                step->addToFailureOutput( "    lambda, |R|: " );
                step->addToFailureOutput( getLambda( ), false );
                step->addToFailureOutput( ", " );
                step->addToFailureOutput( tardigradeVectorTools::l2norm( *getResidual( ) ) );
            }

            updateLambda( );

            incrementLSIteration( );

            step->updateUnknownVector( X0 + getLambda( ) * deltaX );

        }

        if ( step->getFailureVerbosityLevel( ) > 0 ){
            step->addToFailureOutput( "    lambda, |R|: " );
            step->addToFailureOutput( getLambda( ), false );
            step->addToFailureOutput( ", " );
            step->addToFailureOutput( tardigradeVectorTools::l2norm( *getResidual( ) ) );
        }

        if ( !checkLSConvergence( ) ){

            step->resetToleranceScaleFactor( );

            throw convergence_error( "Failure in line search\n" );

        }

        step->resetToleranceScaleFactor( );

        incrementNumLS( );

        resetLSIteration( );

    }
// END LINE SEARCH FUNCTIONS

// BEGIN GRADIENT FUNCTIONS

    void StepDampingBase::setUseGradientDescent( const bool &value ){
        /*!
         * Set whether to attempt a gradient descent step
         * 
         * \param &value: The value of the parameter
         */

        _use_gradient_descent = value;

    }

    /*!
     * Set the value of the rho parameter for gradient descent steps
     *
     * \param &value: The value of the parameter
     */
    void StepDampingBase::setGradientRho( const floatType &value ){
 
        _gradientRho = value;

    }

    /*!
     * Set the value of the p parameter for gradient descent steps
     *
     * \param &value: The value of the parameter
     */
    void StepDampingBase::setGradientP( const floatType &value ){
 
        _gradientP = value;

    }

    /*!
     * Set the value of the beta parameter for gradient descent steps
     *
     * \param &value: The value of the parameter
     */
    void StepDampingBase::setGradientBeta( const floatType &value ){
 
        _gradientBeta = value;

    }

    /*!
     * Set the value of the sigma parameter for gradient descent steps
     *
     * \param &value: The value of the parameter
     */
    void StepDampingBase::setGradientSigma( const floatType &value ){

        _gradientSigma = value;

    }

    /*!
     * Set the value of the maximum number of iterations for gradient descent steps
     *
     * \param &value: The value of the parameter
     */
    void StepDampingBase::setMaxGradientIterations( const unsigned int &value ){

        _maxGradientIterations = value;

    }

    /*!
     * Check if the gradient hasn't exceeded the number of allowed iterations
     */
    bool StepDampingBase::checkGradientIteration( ){

        return getGradientIteration( ) < getMaxGradientIterations( );

    }

    /*!
     * Set the norm of the residual vector
     */
    void StepDampingBase::setResidualNorm( ){

        auto residualNorm = get_SetDataStorage_residualNorm( );

        auto residual = getResidual( );

        using residual_type = std::remove_reference_t<decltype( ( *residual )[ 0 ] )>;

        *residualNorm.value = std::inner_product( std::begin( *residual ), std::end( *residual ), std::begin( *residual ), residual_type( ) );

    }

    /*!
     * Set the derivative of the residual norm w.r.t. the unknown vector
     */
    void StepDampingBase::setdResidualNormdX( ){

        const unsigned int xsize = getNumUnknowns( );

        auto dResidualNormdX = get_SetDataStorage_dResidualNormdX( );

        dResidualNormdX.zero( xsize );

        const floatVector *residual = getResidual( );

        const floatVector *jacobian = getFlatJacobian( );

        for ( unsigned int i = 0; i < xsize; i++ ){
            for ( unsigned int j = 0; j < xsize; j++ ){
                ( *dResidualNormdX.value )[ j ] += 2 * ( *jacobian )[ xsize * i + j ] * ( *residual )[ i ];
            }
        }

    }

    /*!
     * Set the value of the mu_k parameter for gradient descent steps
     *
     * \param &value: The value of the parameter
     */
    void StepDampingBase::setMuk( const floatType &value ){
 
        _mu_k = value;

    }

    /*!
     * Set the initializing scaling value of the mu parameter
     *
     * \param &value: The value of the parameter
     */
    void StepDampingBase::setLMMu( const floatType &value ){
 
        _lm_mu = value;

    }

    /*!
     * Get the base value for the residual norm.
     */
    const floatType *StepDampingBase::get_baseResidualNorm( ){

        if ( !_baseResidualNorm.first ){

            throw std::runtime_error( "The base residual norm must be set with set_baseResidualNorm before it can be called" );

        }

        return &_baseResidualNorm.second;

    }

    /*!
     * Get the base value for the derivative of the residual norm w.r.t. the unknown vector
     */
    const floatVector *StepDampingBase::get_basedResidualNormdX( ){

        if ( !_basedResidualNormdX.first ){

            throw std::runtime_error( "The base residual norm must be set with set_dbaseResidualNormdX before it can be called" );

        }

        return &_basedResidualNormdX.second;

    }

    /*! Set the base value of the residual norm
     *
     * \param &value: The new value
     */
    void StepDampingBase::set_baseResidualNorm( const floatType &value ){

        setNLStepData( value, _baseResidualNorm );

    }

    /*!
     * Set the base derivative of the residual norm w.r.t. the unknown vector
     *
     * \param &value: The new value
     */
    void StepDampingBase::set_basedResidualNormdX( const floatVector &value ){

        setNLStepData( value, _basedResidualNormdX );

    }

// END GRADIENT FUNCTIONS

}
