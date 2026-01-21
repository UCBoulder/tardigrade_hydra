/**
  ******************************************************************************
  * \file tardigrade_SolverBase.cpp
  ******************************************************************************
  * A C++ library for the base classes for solvers
  ******************************************************************************
  */

#include"tardigrade_SolverBase.h"
#include"tardigrade_hydra.h"

namespace tardigradeHydra{

    /*!
     * Get whether the Jacobian being rank-deficient will throw an error
     */
    const bool SolverBase::getRankDeficientError( ){

        return _rank_deficient_error;

    }

    /*!
     * Set whether the Jacobian being rank-deficient will throw an error
     *
     * \param &value: The incoming value
     */
    void SolverBase::setRankDeficientError( const bool &value ){

        _rank_deficient_error = value;
    }

    /*!
     * Add data to the vector of values which will be cleared after each iteration
     * 
     * \param *data: The dataBase object to be cleared
     */
    void SolverBase::addIterationData( dataBase *data ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        hydra->addIterationData( data );

    }

    /*!
     * Add data to the vector of values which will be cleared after each nonlinear step
     * 
     * \param *data: The dataBase object to be cleared
     */
    void SolverBase::addNLStepData( dataBase *data ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        hydra->addNLStepData( data );

    }

    /*!
     * Get the relative tolerance from hydra
     */
    const floatType SolverBase::getRelativeTolerance( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        return hydra->getRelativeTolerance( );

    }

    /*!
     * Get the absolute tolerance from hydra
     */
    const floatType SolverBase::getAbsoluteTolerance( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        return hydra->getAbsoluteTolerance( );

    }

    /*!
     * Get the number of unknowns from hydra
     */
    const unsigned int SolverBase::getNumUnknowns( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        return hydra->getNumUnknowns( );

    }

    /*!
     * Get the unknown vector from hydra
     */
    const floatVector *SolverBase::getUnknownVector( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        return hydra->getUnknownVector( );

    }

    /*!
     * Update the unknown vector
     *
     * \param &value: The new value of the unknown vector
     */
    void SolverBase::updateUnknownVector( const floatVector &value ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        hydra->updateUnknownVector( value );

    }

    /*!
     * Get the residual from hydra
     */
    const floatVector *SolverBase::getResidual( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        return hydra->getResidual( );

    }

    /*!
     * Get the flat Jacobian from hydra
     */
    const floatVector *SolverBase::getFlatJacobian( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        return hydra->getFlatJacobian( );

    }

    /*!
     * Get the number of constraints
     */
    const unsigned int SolverBase::getNumConstraints( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        return hydra->getNumConstraints( );

    }

    /*!
     * Get the constraints
     */
    const floatVector *SolverBase::getConstraints( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        return hydra->getConstraints( );

    }

    /*!
     * Get the constraint Jacobians
     */
    const floatVector *SolverBase::getConstraintJacobians( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        return hydra->getConstraintJacobians( );

    }

    /*!
     * Get the verbosity level for failure messages
     */
    const unsigned int SolverBase::getFailureVerbosityLevel( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        return hydra->getFailureVerbosityLevel( );

    }

    /*!
     * Add the string to the failure output message
     *
     * \param &string: The string to add to the failure output message
     */
    void SolverBase::addToFailureOutput( const std::string &string ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        hydra->addToFailureOutput( string );

    }

    /*!
     * Add a floatVector to the failure output message
     *
     * \param &value: The floatVector to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void SolverBase::addToFailureOutput( const floatVector &value, bool add_endline ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        hydra->addToFailureOutput( value, add_endline );

    }

    /*!
     * Add a vector of booleans to the failure output message
     *
     * \param &value: The vector of booleans to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void SolverBase::addToFailureOutput( const std::vector<bool> &value, bool add_endline ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        hydra->addToFailureOutput( value, add_endline );

    }

    /*!
     * Add a floatType to the failure output message
     *
     * \param &value: The floatType to add to the failure output message
     * \param add_endline: Whether to add an endline after the value or not
     */
    void SolverBase::addToFailureOutput( const floatType &value, bool add_endline ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        hydra->addToFailureOutput( value, add_endline );

    }

    /*!
     * Get the scale factor for the tolerance
     */
    const floatType SolverBase::getToleranceScaleFactor( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        return hydra->getToleranceScaleFactor( );

    }

    /*!
     * Reset the tolerance scale factor
     */
    void SolverBase::resetToleranceScaleFactor( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        hydra->resetToleranceScaleFactor( );

    }

    /*!
     * Reset the solver
     */
    void SolverBase::reset( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        resetIterations( );
        step->reset( );

    }

    /*!
     * Solve the problem
     */
    void SolverBase::solve( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        // Form the initial unknown vector
        if ( getInitializeUnknownVector( ) ){
            TARDIGRADE_ERROR_TOOLS_CATCH( hydra->initializeUnknownVector( ) );
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

// NONLINEAR FUNCTIONS

    /*!
     * Set the initialize unknown vector flag
     * 
     * \param &value: The value of the flag
     */
    void SolverBase::setInitializeUnknownVector( const bool &value ){

        _initializeUnknownVector = value;

    }

    /*!
     * Set the maximum number of allowable nonlinear iterations
     *
     * \param &value: The maximum number of iterations
     */
    void SolverBase::setMaxIterations( const unsigned int &value ){

        _maxIterations = value;

    }

    /*!
     * Increment the iteration
     */
    void SolverBase::incrementIteration( ){
        TARDIGRADE_ERROR_TOOLS_CHECK( step != nullptr, "The step has not been defined" );
        _iteration++; step->damping->reset( );
    }

    /*!
     * Reset all nonlinear step data
     */
    void SolverBase::resetNLStepData( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        hydra->resetNLStepData( );

    }

    /*!
     * Signal to the residuals that a successful nonlinear step has been performed
     */
    void SolverBase::callResidualSuccessfulNLStep( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        hydra->setAllowModifyGlobalResidual( true );

        hydra->setCurrentResidualIndexMeaningful( true );

        for ( auto residual_ptr = std::begin( *( hydra->getResidualClasses( ) ) ); residual_ptr != std::end( *( hydra->getResidualClasses( ) ) ); ++residual_ptr ){

            hydra->setCurrentResidualIndex( residual_ptr - std::begin( *( hydra->getResidualClasses( ) ) ) );

            ( *residual_ptr )->successfulNLStep( );

        }

        hydra->setCurrentResidualIndexMeaningful( false );

        hydra->setAllowModifyGlobalResidual( false );

    }

    /*!
     * Signal to the residuals that we are about to start a nonlinear solve
     */
    void SolverBase::callResidualPreNLSolve( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        hydra->setCurrentResidualIndexMeaningful( true );

        for ( auto residual_ptr = std::begin( *( hydra->getResidualClasses( ) ) ); residual_ptr != std::end( *( hydra->getResidualClasses( ) ) ); ++residual_ptr ){

            hydra->setCurrentResidualIndex( residual_ptr - std::begin( *( hydra->getResidualClasses( ) ) ) );

            ( *residual_ptr )->preNLSolve( );

        }

        hydra->setCurrentResidualIndexMeaningful( false );

    }


    /*!
     * Signal to the residuals that we have finisehd a nonlinear solve
     */
    void SolverBase::callResidualPostNLSolve( ){

        TARDIGRADE_ERROR_TOOLS_CHECK( hydra != nullptr, "Hydra has not been defined" );
        hydra->setCurrentResidualIndexMeaningful( true );

        for ( auto residual_ptr = std::begin( *( hydra->getResidualClasses( ) ) ); residual_ptr != std::end( *( hydra->getResidualClasses( ) ) ); ++residual_ptr ){

            hydra->setCurrentResidualIndex( residual_ptr - std::begin( *( hydra->getResidualClasses( ) ) ) );

            try{

                ( *residual_ptr )->postNLSolve( );

            }
            catch( std::exception &e ){

                if ( getFailureVerbosityLevel( ) > 0 ){

                    addToFailureOutput( "Failure in residual " + std::to_string( residual_ptr - std::begin( *( hydra->getResidualClasses( ) ) ) ) + "\n" );
                    std::string message;
                    tardigradeErrorTools::captureNestedExceptions(e, message);
                    addToFailureOutput( message );

                }

                throw;

            }

        }

        hydra->setCurrentResidualIndexMeaningful( false );

    }

    /*!
     * Check the convergence
     */
    bool SolverBase::checkConvergence( ){

        const floatVector *tolerance = getTolerance( );

        const floatVector *residual = getResidual( );

        TARDIGRADE_ERROR_TOOLS_CHECK(
            tolerance->size( ) == residual->size( ),
            "The residual and tolerance vectors don't have the same size\n  tolerance: " + std::to_string( tolerance->size( ) ) +
            "\n  residual:  " + std::to_string( residual->size( ) ) + "\n" 
        );

        for ( unsigned int i = 0; i < tolerance->size( ); i++ ){

            if ( std::fabs( ( *residual )[ i ] ) > ( *tolerance )[ i ] ){

                return false;

            }

        }

        return true;

    }

    /*!
     * Set the tolerance
     * 
     * \f$ tol = tolr * ( |R_0| + |X| ) + tola \f$
     */
    void SolverBase::setTolerance( ){

        auto tolerance = get_SetDataStorage_tolerance( );

        *tolerance.value = tardigradeVectorTools::abs( *getResidual( ) ) + tardigradeVectorTools::abs( *getUnknownVector( ) );

        *tolerance.value = getRelativeTolerance( ) * ( *tolerance.value ) + getAbsoluteTolerance( );

    }

    /*!
     * Set the tolerance
     *
     * \param tolerance: The tolerance vector for each value of the residual
     */
    void SolverBase::setTolerance( const floatVector &tolerance ){

        setConstantData( tolerance, _tolerance );

    }

    /*!
     * Return a SetDataStorageConstant setter for the tolerance
     */
    SolverBase::SetDataStorageConstant<floatVector> SolverBase::get_SetDataStorage_tolerance( ){

        return SetDataStorageConstant<floatVector>( &_tolerance );

    }

    /*!
     * Get the tolerance
     */
    const floatVector* SolverBase::getTolerance( ){

        if ( !_tolerance.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setTolerance( ) );

        }

        return &_tolerance.second;

    }

// END NONLINEAR FUNCTIONS

}
