/**
 ******************************************************************************
 * \file tardigrade_ArmijoGradientDamping.cpp
 ******************************************************************************
 * The base class for a combined Armijo - gradient descent damping
 ******************************************************************************
 */

#include"tardigrade_ArmijoGradientDamping.h"
#include"tardigrade_vector_tools.h"
#include"tardigrade_CustomErrors.h"
#include"tardigrade_SolverStepBase.h"

namespace tardigradeHydra{

    /*!
     * Apply the damping to the proposed step
     *
     * This function must call the updateUnknownVector function
     *
     * Returns true if any damping was applied
     */
    const bool ArmijoGradientDamping::applyDamping( ){

        updateUnknownVector( step->X0 + getLambda( ) * step->deltaX );

        // Refine the estimate if the new point has a higher residual
        if ( !checkLSConvergence( ) ){

            if ( checkDescentDirection( step->deltaX ) || !getUseGradientDescent( ) ){

                // Perform an Armijo type line search when the search direction is aligned with the gradient
                performArmijoTypeLineSearch( step->X0, step->deltaX );

            }
            else{

                // Perform gradient descent if the search direction is not aligned with the gradient
                performGradientStep( step->X0 );

            }

            return true;

        }
        else{

            resetToleranceScaleFactor( );

            return false;

        }

    }

    /*!
     * Reset the internal count
     */
    void ArmijoGradientDamping::resetCounts( ){

        resetNumLS( );
        resetNumGrad( );
    }

    /*!
     * Set the base quantities prior to updating the unknown vector
     */
    void ArmijoGradientDamping::setBaseQuantities( ){

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

//// LINE SEARCH FUNCTIONS
//
//
//    /*!
//     * Set the value of the line-search alpha parameter
//     *
//     * \param &value: The incoming value of the line-search alpha parameter
//     */
//    void ArmijoGradientDamping::setLSAlpha( const floatType &value ){
//
//        _lsAlpha = value;
//
//    }
//
//    /*!
//     * Set the maximum number of line-search iterations
//     *
//     * \param &value: The incoming value
//     */
//    void ArmijoGradientDamping::setMaxLSIterations( const unsigned int &value ){
//
//        _maxLSIterations = value;
//
//    }
//
//    /*!
//     * Get the residual norm for the line-search convergence criterion
//     */
//    const floatType* ArmijoGradientDamping::getLSResidualNorm( ){
//
//        if ( !_lsResidualNorm.first ){
//
//            TARDIGRADE_ERROR_TOOLS_CATCH( resetLSIteration( ) );
//
//        }
//
//        return &_lsResidualNorm.second;
//
//    }
//
//    /*!
//     * Set the line-search residual norm
//     */
//    void ArmijoGradientDamping::setLSResidualNorm( ){
//
//        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
//        _lsResidualNorm.second = tardigradeVectorTools::l2norm( *getResidual( ) );
//
//        _lsResidualNorm.first = true;
//
//    }
//
//    /*!
//     * Reset the line search iteration
//     */
//    void ArmijoGradientDamping::resetLSIteration( ){
//
//        _LSIteration = 0;
//
//        _lambda = 1.0;
//
//        setLSResidualNorm( );
//
//    };
//
//    /*!
//     * Check the line-search convergence
//     */
//    bool ArmijoGradientDamping::checkLSConvergence( ){
//
//        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
//        if ( tardigradeVectorTools::l2norm( *getResidual( ) ) < step->getToleranceScaleFactor( ) * ( 1 - getLSAlpha( ) ) * ( *getLSResidualNorm( ) ) ){
//
//            return true;
//
//        }
//
//        return false;
//
//    }
//
//    /*!
//     * Check the current line search iteration
//     */
//    bool ArmijoGradientDamping::checkLSIteration( ){
//        return getLSIteration( ) < getMaxLSIterations( );
//    }
//
//    /*!
//     * Perform an Armijo-type line search
//     *
//     * \param &X0: The base value of the unknown vector
//     * \param &deltaX: The proposed change in X
//     */
//    void ArmijoGradientDamping::performArmijoTypeLineSearch( const floatVector &X0, const floatVector &deltaX ){
//
//        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The solver has not been defined" );
//        while ( !checkLSConvergence( ) && checkLSIteration( ) ){
//
//            if ( step->getFailureVerbosityLevel( ) > 0 ){
//                step->addToFailureOutput( "    lambda, |R|: " );
//                step->addToFailureOutput( getLambda( ), false );
//                step->addToFailureOutput( ", " );
//                step->addToFailureOutput( tardigradeVectorTools::l2norm( *getResidual( ) ) );
//            }
//
//            updateLambda( );
//
//            incrementLSIteration( );
//
//            step->updateUnknownVector( X0 + getLambda( ) * deltaX );
//
//        }
//
//        if ( step->getFailureVerbosityLevel( ) > 0 ){
//            step->addToFailureOutput( "    lambda, |R|: " );
//            step->addToFailureOutput( getLambda( ), false );
//            step->addToFailureOutput( ", " );
//            step->addToFailureOutput( tardigradeVectorTools::l2norm( *getResidual( ) ) );
//        }
//
//        if ( !checkLSConvergence( ) ){
//
//            step->resetToleranceScaleFactor( );
//
//            throw convergence_error( "Failure in line search\n" );
//
//        }
//
//        step->resetToleranceScaleFactor( );
//
//        incrementNumLS( );
//
//        resetLSIteration( );
//
//    }
//// END LINE SEARCH FUNCTIONS
//
//// BEGIN GRADIENT FUNCTIONS
//
//    void ArmijoGradientDamping::setUseGradientDescent( const bool &value ){
//        /*!
//         * Set whether to attempt a gradient descent step
//         * 
//         * \param &value: The value of the parameter
//         */
//
//        _use_gradient_descent = value;
//
//    }
//
//    /*!
//     * Set the value of the rho parameter for gradient descent steps
//     *
//     * \param &value: The value of the parameter
//     */
//    void ArmijoGradientDamping::setGradientRho( const floatType &value ){
// 
//        _gradientRho = value;
//
//    }
//
//    /*!
//     * Set the value of the p parameter for gradient descent steps
//     *
//     * \param &value: The value of the parameter
//     */
//    void ArmijoGradientDamping::setGradientP( const floatType &value ){
// 
//        _gradientP = value;
//
//    }
//
//    /*!
//     * Set the value of the beta parameter for gradient descent steps
//     *
//     * \param &value: The value of the parameter
//     */
//    void ArmijoGradientDamping::setGradientBeta( const floatType &value ){
// 
//        _gradientBeta = value;
//
//    }
//
//    /*!
//     * Set the value of the sigma parameter for gradient descent steps
//     *
//     * \param &value: The value of the parameter
//     */
//    void ArmijoGradientDamping::setGradientSigma( const floatType &value ){
//
//        _gradientSigma = value;
//
//    }
//
//    /*!
//     * Set the value of the maximum number of iterations for gradient descent steps
//     *
//     * \param &value: The value of the parameter
//     */
//    void ArmijoGradientDamping::setMaxGradientIterations( const unsigned int &value ){
//
//        _maxGradientIterations = value;
//
//    }
//
//    /*!
//     * Check if the gradient hasn't exceeded the number of allowed iterations
//     */
//    bool ArmijoGradientDamping::checkGradientIteration( ){
//
//        return getGradientIteration( ) < getMaxGradientIterations( );
//
//    }
//
//    /*!
//     * Set the norm of the residual vector
//     */
//    void ArmijoGradientDamping::setResidualNorm( ){
//
//        auto residualNorm = get_SetDataStorage_residualNorm( );
//
//        auto residual = getResidual( );
//
//        using residual_type = std::remove_reference_t<decltype( ( *residual )[ 0 ] )>;
//
//        *residualNorm.value = std::inner_product( std::begin( *residual ), std::end( *residual ), std::begin( *residual ), residual_type( ) );
//
//    }
//
//    /*!
//     * Set the derivative of the residual norm w.r.t. the unknown vector
//     */
//    void ArmijoGradientDamping::setdResidualNormdX( ){
//
//        const unsigned int xsize = getNumUnknowns( );
//
//        auto dResidualNormdX = get_SetDataStorage_dResidualNormdX( );
//
//        dResidualNormdX.zero( xsize );
//
//        const floatVector *residual = getResidual( );
//
//        const floatVector *jacobian = getFlatJacobian( );
//
//        for ( unsigned int i = 0; i < xsize; i++ ){
//            for ( unsigned int j = 0; j < xsize; j++ ){
//                ( *dResidualNormdX.value )[ j ] += 2 * ( *jacobian )[ xsize * i + j ] * ( *residual )[ i ];
//            }
//        }
//
//    }
//
//    /*!
//     * Set the value of the mu_k parameter for gradient descent steps
//     *
//     * \param &value: The value of the parameter
//     */
//    void ArmijoGradientDamping::setMuk( const floatType &value ){
// 
//        _mu_k = value;
//
//    }
//
//    /*!
//     * Set the initializing scaling value of the mu parameter
//     *
//     * \param &value: The value of the parameter
//     */
//    void ArmijoGradientDamping::setLMMu( const floatType &value ){
// 
//        _lm_mu = value;
//
//    }
//
//    /*!
//     * Get the base value for the residual norm.
//     */
//    const floatType *ArmijoGradientDamping::get_baseResidualNorm( ){
//
//        if ( !_baseResidualNorm.first ){
//
//            throw std::runtime_error( "The base residual norm must be set with set_baseResidualNorm before it can be called" );
//
//        }
//
//        return &_baseResidualNorm.second;
//
//    }
//
//    /*!
//     * Get the base value for the derivative of the residual norm w.r.t. the unknown vector
//     */
//    const floatVector *ArmijoGradientDamping::get_basedResidualNormdX( ){
//
//        if ( !_basedResidualNormdX.first ){
//
//            throw std::runtime_error( "The base residual norm must be set with set_dbaseResidualNormdX before it can be called" );
//
//        }
//
//        return &_basedResidualNormdX.second;
//
//    }
//
//    /*! Set the base value of the residual norm
//     *
//     * \param &value: The new value
//     */
//    void ArmijoGradientDamping::set_baseResidualNorm( const floatType &value ){
//
//        setNLStepData( value, _baseResidualNorm );
//
//    }
//
//    /*!
//     * Set the base derivative of the residual norm w.r.t. the unknown vector
//     *
//     * \param &value: The new value
//     */
//    void ArmijoGradientDamping::set_basedResidualNormdX( const floatVector &value ){
//
//        setNLStepData( value, _basedResidualNormdX );
//
//    }
//
//    /*!
//     * Check if the search direction is a descent direction of the Jacobian
//     * 
//     * \param &dx: The proposed change in x
//     */
//    bool ArmijoGradientDamping::checkDescentDirection( const floatVector &dx ){
//
//        const unsigned int xsize = getNumUnknowns( );
//
//        const floatType RHS = -getGradientRho( ) * std::pow( tardigradeVectorTools::l2norm( dx ), getGradientP( ) );
//
//        floatType LHS = 0;
//
//        const floatVector *dResidualNormdX = get_basedResidualNormdX( );
//
//        for ( unsigned int i = 0; i < xsize; i++ ){
//
//            LHS += ( *dResidualNormdX )[ i ] * dx[ i ];
//
//        }
//
//        return LHS <= RHS;
//
//    }
//
//    /*!
//     * Check the convergence of a gradient step
//     *
//     * \param &X0: The initial value of the unknown vector
//     */
//    bool ArmijoGradientDamping::checkGradientConvergence( const floatVector &X0 ){
//
//        const unsigned int xsize = getNumUnknowns( );
//
//        floatVector dx = *getUnknownVector( ) - X0;
//
//        floatType RHS = *get_baseResidualNorm( );
//
//        for ( unsigned int i = 0; i < xsize; ++i ){
//
//            RHS += getGradientSigma( ) * ( *get_basedResidualNormdX( ) )[ i ] * dx[ i ];
//
//        }
//
//        return ( *get_residualNorm( ) ) < getToleranceScaleFactor( ) * RHS;
//
//    }
//
//    /*!
//     * Perform a gradient descent step
//     *
//     * \param &X0: The base value of the unknown vector
//     */
//    void ArmijoGradientDamping::performGradientStep( const floatVector &X0 ){
//
//        const floatVector *dResidualNormdX = get_basedResidualNormdX( );
//
//        unsigned int l                     = 0;
//
//        const unsigned int maxiter         = getMaxGradientIterations( );
//
//        while( checkGradientIteration( ) ){
//
//            floatType t = std::pow( getGradientBeta( ), l );
//
//            updateUnknownVector( X0 - t * ( *dResidualNormdX ) );
//
//            if ( checkGradientConvergence( X0 ) ){
//
//                break;
//
//            }
//
//            l++;
//
//            incrementGradientIteration( );
//
//        }
//
//        resetToleranceScaleFactor( );
//
//        if ( l >= maxiter ){
//
//            throw convergence_error( "Failure in gradient step" );
//
//        }
//
//        incrementNumGrad( );
//
//        resetGradientIteration( );
//
//    }
//
//// END GRADIENT FUNCTIONS

}
