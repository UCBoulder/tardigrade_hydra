/**
  ******************************************************************************
  * \file tardigrade_hydra.cpp
  ******************************************************************************
  * A C++ library for defining frameworks to solve finite deformation material
  * models.
  ******************************************************************************
  */

#include"tardigrade_hydra.h"
#include"tardigrade_SolverStepBase.h"

namespace tardigradeHydra{

    //Define hydra global constants in a place that Doxygen can pick up for documentation
    /** \brief Define the expected number of tensor spatial dimensions. */
    const int spatialDimensions = 3;

    hydraBase::hydraBase( const floatType &time, const floatType &deltaTime,
                          const floatType &temperature, const floatType &previousTemperature,
                          const secondOrderTensor &deformationGradient, const secondOrderTensor &previousDeformationGradient,
                          const floatVector &additionalDOF, const floatVector &previousAdditionalDOF,
                          const floatVector &previousStateVariables, const floatVector &parameters,
                          const unsigned int numConfigurations, const unsigned int numNonLinearSolveStateVariables,
                          const unsigned int dimension, const unsigned int configuration_unknown_count, const floatType tolr, const floatType tola, const unsigned int maxIterations,
                          const unsigned int maxLSIterations, const floatType lsAlpha,
                          const bool use_preconditioner, const unsigned int preconditioner_type ) : _dimension( dimension ),
                                                           _configuration_unknown_count( configuration_unknown_count ),
                                                           _stress_size( configuration_unknown_count ),
                                                           _time( time ), _deltaTime( deltaTime ),
                                                           _temperature( temperature ), _previousTemperature( previousTemperature ),
                                                           _deformationGradient( deformationGradient ),
                                                           _previousDeformationGradient( previousDeformationGradient ),
                                                           _additionalDOF( additionalDOF ),
                                                           _previousAdditionalDOF( previousAdditionalDOF ),
                                                           _previousStateVariables( previousStateVariables ),
                                                           _parameters( parameters ),
                                                           _numConfigurations( numConfigurations ),
                                                           _numNonLinearSolveStateVariables( numNonLinearSolveStateVariables ),
                                                           _tolr( tolr ), _tola( tola ){
        /*!
         * The main constructor for the hydra base class. Inputs are all the required values for most solves.
         * 
         * \param &time: The current time
         * \param &deltaTime: The change in time
         * \param &temperature: The current temperature
         * \param &previousTemperature: The previous temperature
         * \param &deformationGradient: The current deformation gradient
         * \param &additionalDOF: Any additional degrees of freedom required for the model
         * \param &previousAdditionalDOF: Any previous additional degrees of freedom required for the model
         * \param &previousDeformationGradient The previous deformation gradient
         * \param &previousStateVariables: The previous state variables
         * \param &parameters: The model parameters
         * \param &numConfigurations: The number of configurations
         * \param &numNonLinearSolveStateVariables: The number of state variables which will contribute terms to the non-linear solve's residual
         * \param &dimension: The dimension of the problem (defaults to 3)
         * \param &configuration_unknown_count: The number of unknowns in each configuration (defaults to 9)
         * \param &tolr: The relative tolerance (defaults to 1e-9)
         * \param &tola: The absolute tolerance (defaults to 1e-9)
         * \param &maxIterations: The maximum number of non-linear iterations (defaults to 20)
         * \param &maxLSIterations: The maximum number of line-search iterations (defaults to 5)
         * \param &lsAlpha: The alpha term for the line search (defaults to 1e-4)
         * \param &use_preconditioner: A flag for whether to pre-condition the Jacobian (can help with scaling issues)
         * \param &preconditioner_type: The type of pre-conditioner to use. Options are
         *     0. A diagonal pre-conditioner populate by the inverse of the absolute largest entries of the Jacobian's rows
         */

        // Initialize the scaled-quantities
        setScaledQuantities( );

        // Decompose the state variable vector initializing all of the configurations
        decomposeStateVariableVector( );

        // Set the residual classes
        setResidualClasses( );

        // TEMP
        _solver.hydra = this;
        solver->setMaxIterations( maxIterations );
        solver->step->setSolver( solver );
        solver->step->setLSAlpha( lsAlpha );
        solver->step->setMaxLSIterations( maxLSIterations );
        solver->preconditioner->setSolver( solver );
        solver->preconditioner->_use_preconditioner = use_preconditioner;
        solver->preconditioner->_preconditioner_type = preconditioner_type;
        // END TEMP

    }

    void hydraBase::setStress( const floatVector &stress ){
        /*!
         * Set the value of the stress
         * 
         * \param &stress: The stress in row-major form
         */

        setIterationData( stress, _stress );

    }

    hydraBase::SetDataStorageIteration<secondOrderTensor> hydraBase::get_SetDataStorage_stress( ){
        /*!
         * Get a SetDataStorage object for the stress
         */

        return hydraBase::SetDataStorageIteration<secondOrderTensor>( &_stress, this );

    }

    void hydraBase::extractStress( ){
        /*!
         * Extract the stresses out of the unknown vector
         */

        const floatVector *unknownVector = getUnknownVector( );

        auto stress = get_SetDataStorage_stress( );

        stress.zero( getConfigurationUnknownCount( ) );

        std::copy( std::begin( *unknownVector ), std::begin( *unknownVector ) + getConfigurationUnknownCount( ), std::begin( *stress.value ) );

    }

    void hydraBase::computeConfigurations( const floatVector *data_vector, const unsigned int start_index,
                                           const floatVector &total_transformation,
                                           floatVector &configurations, floatVector &inverseConfigurations,
                                           const bool add_eye ){
        /*!
         * Compute the configurations from the provided vector. Each configuration is assumed to have a dimension
         * of dimension x dimension
         *
         * \param *data_vector: A pointer to the vector of data which contains the configurations and other information
         * \param &start_index: The starting index for the vector
         * \param &total_transformation: The total transformation from the reference to the current configuration
         * \param &configurations: The resulting collection of configurations
         * \param &inverseConfigurations: The resulting inverse configurations
         * \param add_eye: A flag for whether to add the identity matrix to each of the configurations except for the
         *     total transformation. Defaults to false.
         */

        auto dim = getDimension( );
        auto sot_dim = getSOTDimension( );

        auto num_configs = getNumConfigurations( );

        // Set the configurations
        configurations = floatVector( num_configs * sot_dim, 0 );

        inverseConfigurations = floatVector( num_configs * sot_dim, 0 );

        auto mat = tardigradeHydra::getFixedSizeMatrixMap<floatType, 3, 3>( inverseConfigurations.data() );
#ifdef TARDIGRADE_HYDRA_USE_LLXSMM
        kernel_type kernel(LIBXSMM_GEMM_FLAG_NONE, dim, dim, dim, 1, 0 );

        // Initialize the first configuration with the total deformation gradient
        secondOrderTensor temp( sot_dim, 0 );
#else
        auto mat2 = tardigradeHydra::getFixedSizeMatrixMap<floatType, 3, 3>( configurations.data() );
#endif

        std::copy( total_transformation.begin( ), total_transformation.end( ), configurations.begin( ) );

        for ( int i = num_configs - 2; i >= 0; i-- ){

            // Set the current configuration as being equal to the previous
            std::copy( data_vector->begin( ) + i * sot_dim + start_index,
                       data_vector->begin( ) + ( i + 1 ) * sot_dim + start_index,
                       configurations.begin( ) + sot_dim * ( i + 1 ) );

            if ( add_eye ){

                for ( unsigned int j = 0; j < dim; j++ ){ configurations[ sot_dim * ( i + 1 ) + dim * j + j ] += 1; }

            }

            // Compute the inverse of the current configuration and store it
            std::copy( configurations.begin( ) + sot_dim * ( i + 1 ),
                       configurations.begin( ) + sot_dim * ( i + 2 ),
                       inverseConfigurations.begin( ) + sot_dim * ( i + 1 ) );
            new (&mat) Eigen::Map< Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> >( inverseConfigurations.data() + sot_dim * ( i + 1 ), 3, 3 );
            mat = mat.inverse( ).eval( );

#ifdef TARDIGRADE_HYDRA_USE_LLXSMM
            std::copy( configurations.begin( ),
                       configurations.begin( ) + sot_dim,
                       temp.begin( ) );

            kernel( &inverseConfigurations[ sot_dim * ( i + 1 ) ], &temp[ 0 ], &configurations[ 0 ] );
#else
            // Add contribution of deformation gradient to the first configuration

            new (&mat2) Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor>>( configurations.data( ), 3, 3 );

            mat2 *= mat;
#endif

        }

        std::copy( configurations.begin( ),
                   configurations.begin( ) + sot_dim,
                   inverseConfigurations.begin( ) );

        new (&mat) Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> >( inverseConfigurations.data(), 3, 3 );
        mat = mat.inverse( ).eval( );

        return;

    }

    void hydraBase::updateConfigurationsFromUnknownVector( ){
        /*!
         * Update the configurations from the unknown vector
         */

        const floatVector *unknownVector = getUnknownVector( );

        // Set the configurations
        auto configurations = get_SetDataStorage_configurations( );

        auto inverseConfigurations = get_SetDataStorage_inverseConfigurations( );

        computeConfigurations( unknownVector, getStressSize( ), *getDeformationGradient( ), *configurations.value, *inverseConfigurations.value );

        // Extract the remaining state variables required for the non-linear solve
        auto nonLinearSolveStateVariables = get_SetDataStorage_nonLinearSolveStateVariables( );

        auto nNLISV = getNumNonLinearSolveStateVariables( );

        nonLinearSolveStateVariables.zero( nNLISV );

        std::copy( std::begin( *unknownVector ) + getNumConfigurations( ) * getConfigurationUnknownCount( ),
                   std::end(   *unknownVector ),
                   std::begin( *nonLinearSolveStateVariables.value ) );

    }

    void hydraBase::decomposeUnknownVector( ){
        /*!
         * Decompose the unknown vector into the cauchy stress, configurations, and state variables used for the non-linear solve
         */

        // Set the stress
        extractStress( );

        updateConfigurationsFromUnknownVector( );

    }

    void hydraBase::decomposeStateVariableVector( ){
        /*!
         * Decompose the incoming state variable vector setting the different configurations along the way
         * 
         * The state variable vector is assumed to be of the form:
         * 
         * \f$ \text{ISV} = \left\{\bf{F}^2 - \bf{I}, \bf{F}^3 - \bf{I}, \cdots, \bf{F}^n - \bf{I}, \xi^1, \xi^2, \cdots, \xi^m, \eta^1, \cdots\right\} \f$
         * 
         * where the \f$\bf{F}^x\f$ are the different configurations, \f$\xi^y\f$ are the other variables to be solved during the non-linear
         * solve and \f$\eta^z\f$ are other state variables. Note that we decompose the deformation gradient as
         * 
         * \f$\bf{F} = \bf{F}^1 \bf{F}^2 \cdots \bf{F}^n\f$
         * 
         * and so because \f$\bf{F}\f$ is provided we can solve for \f$\bf{F}^1\f$. Typically, this configuration would be the elastic
         * configuration (i.e., the configuration that generates the stress) though we do not insist that users follow convention.
         */

        auto nConfig = getNumConfigurations( );

        auto nNLISV  = getNumNonLinearSolveStateVariables( );

        // Extract the previous configurations
        if ( getPreviousStateVariables( )->size( ) < ( ( nConfig - 1 ) * getConfigurationUnknownCount( ) + nNLISV ) ){

            std::string message = "The number of state variables is less than required for the configurations and ";
            message            += "non-linear state variables\n";
            message            += "  # previousStateVariables                               : " + std::to_string( getPreviousStateVariables( )->size( ) ) + "\n";
            message            += "  # ( configurations - 1 ) * configuration_unknown_count : " + std::to_string( ( nConfig - 1 ) * getConfigurationUnknownCount( ) ) + "\n";
            message            += "  # non-linear solve ISVs                                : " + std::to_string( nNLISV ) + "\n";
            message            += "  # minimum required ISVs                                : " + std::to_string( ( nConfig - 1 ) * getConfigurationUnknownCount( ) + nNLISV );

            TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

        }

        auto configurations = get_SetDataStorage_configurations( );

        auto previousConfigurations = get_SetDataStorage_previousConfigurations( );

        auto inverseConfigurations = get_SetDataStorage_inverseConfigurations( );

        auto previousInverseConfigurations = get_SetDataStorage_previousInverseConfigurations( );

        // Compute the configurations
        computeConfigurations( getPreviousStateVariables( ), 0, *getDeformationGradient( ), *configurations.value, *inverseConfigurations.value, true );

        computeConfigurations( getPreviousStateVariables( ), 0, *getPreviousDeformationGradient( ), *previousConfigurations.value, *previousInverseConfigurations.value, true );

        // Extract the remaining state variables required for the non-linear solve
        auto nonLinearSolveStateVariables         = get_SetDataStorage_nonLinearSolveStateVariables( );
        auto previousNonLinearSolveStateVariables = get_SetDataStorage_previousNonLinearSolveStateVariables( );

        previousNonLinearSolveStateVariables.zero( nNLISV );

        std::copy( std::begin( *getPreviousStateVariables( ) ) + ( nConfig - 1 ) * getConfigurationUnknownCount( ),
                   std::begin( *getPreviousStateVariables( ) ) + ( nConfig - 1 ) * getConfigurationUnknownCount( ) + nNLISV,
                   std::begin( *previousNonLinearSolveStateVariables.value ) );

        *nonLinearSolveStateVariables.value = *previousNonLinearSolveStateVariables.value;

        // Extract the additional state variables
        auto additionalStateVariables         = get_SetDataStorage_additionalStateVariables( );
        auto previousAdditionalStateVariables = get_SetDataStorage_previousAdditionalStateVariables( );

        unsigned int nAISV = ( unsigned int )( std::end( *getPreviousStateVariables( ) ) - ( std::begin( *getPreviousStateVariables( ) ) + ( nConfig - 1 ) * getConfigurationUnknownCount( ) + nNLISV ) );

        previousAdditionalStateVariables.zero( nAISV );

        std::copy( std::begin( *getPreviousStateVariables( ) ) + ( nConfig - 1 ) * getConfigurationUnknownCount( ) + nNLISV,
                   std::end(   *getPreviousStateVariables( ) ),
                   std::begin( *previousAdditionalStateVariables.value ) );

        *additionalStateVariables.value = *previousAdditionalStateVariables.value;

    }

    std::string hydraBase::build_upper_index_out_of_range_error_string( const unsigned int upperIndex, const unsigned int num_configurations ){
        /*!
         * Build an error message for when the upper index is larger than the number of configurations
         *
         * \param upperIndex: The upper index
         * \param num_configurations: The number of configurations
         */

        std::string message = "The upper index must be less than or equal to the total number of configurations\n";
        message            += "  upperIndex      : " + std::to_string( upperIndex ) + "\n";
        message            += "  # configurations: " + std::to_string( num_configurations );

        return message;
    }

    std::string hydraBase::build_lower_index_out_of_range_error_string( const unsigned int lowerIndex, const unsigned int upperIndex ){
        /*!
         * Build an error message for when the lower index is larger than the upper index
         *
         * \param lowerIndex: The lower configuration index
         * \param upperIndex: The upper configuration index
         */

        std::string message = "The upper index must be greater than or equal to the lower index\n";
        message            += "  lowerIndex: " + std::to_string( lowerIndex ) + "\n";
        message            += "  upperIndex: " + std::to_string( upperIndex ) + "\n";

        return message;

    }

    floatVector hydraBase::getSubConfiguration( const floatVector &configurations, const unsigned int &lowerIndex,
                                                const unsigned int &upperIndex ){
        /*!
         * Get a sub-configuration \f$\bf{F}^{sc}\f$ defined as
         *
         * \f$ F^{sc}_{iI} = F^{\text{lowerIndex}}_{i\hat{I}} F^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots F^{\text{upperIndex-1}}_{\bar{I}I} \f$
         * \param &configurations: The configurations to operate on
         * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
         * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
         *   Note, the configuration indicated by the index is NOT included in the sub-configuration
         */

        const unsigned int dim = getDimension( );
        const unsigned int sot_dim = getSOTDimension( );
        const unsigned int local_num_configurations = configurations.size( ) / sot_dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( configurations.size( ) % sot_dim == 0, "The configurations vector must be a multiple of the size of a second order tensor" )

        TARDIGRADE_ERROR_TOOLS_CHECK( upperIndex <= local_num_configurations, build_upper_index_out_of_range_error_string( upperIndex, local_num_configurations ) )

        TARDIGRADE_ERROR_TOOLS_CHECK( lowerIndex <= upperIndex, build_lower_index_out_of_range_error_string( lowerIndex, upperIndex ) )

        secondOrderTensor Fsc( sot_dim, 0 );
        for ( unsigned int i = 0; i < 3; i++ ){ Fsc[ dim * i + i ] = 1.; }

#ifdef TARDIGRADE_HYDRA_USE_LLXSMM
        floatVector temp;
        kernel_type kernel(LIBXSMM_GEMM_FLAG_NONE, dim, dim, dim, 1, 0 );
#else
        auto Fsc_mat = tardigradeHydra::getFixedSizeMatrixMap<floatType, 3, 3>( Fsc.data( ) );
        auto mat = tardigradeHydra::getFixedSizeMatrixMap<floatType, 3, 3>( configurations.data() );
#endif

        for ( unsigned int i = lowerIndex; i < upperIndex; i++ ){

#ifdef TARDIGRADE_HYDRA_USE_LLXSMM
            temp = Fsc;

            kernel( &configurations[ sot_dim * i ], &temp[ 0 ], &Fsc[ 0 ] );
#else
            new (&mat) Eigen::Map<const Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor>>( configurations.data( ) + sot_dim * i, 3, 3 );
            Fsc_mat *= mat;
#endif

        }

        return Fsc;

    }

    floatVector hydraBase::getSubConfigurationJacobian( const floatVector &configurations, const unsigned int &lowerIndex,
                                                        const unsigned int &upperIndex ){
        /*!
         * Get the jacobian of a sub-configuration \f$\bf{F}^{sc}\f$ defined as
         *
         * \f$ F^{sc}_{iI} = F^{\text{lowerIndex}}_{i\hat{I}} F^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots F^{\text{upperIndex-1}}_{\bar{I}I} \f$
         * 
         * with respect to all of the configurations. The returned matrix will be of size ( dimensions**2, configurations.size( ) * dimensions**2 )
         * 
         * \param &configurations: The configurations to operate on
         * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
         * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
         *   Note, the configuration indicated by the index is NOT included in the sub-configuration
         */

        const unsigned int dim = getDimension( );
        const unsigned int sot_dim = getSOTDimension( );
        const unsigned int num_incoming_configs = configurations.size( ) / sot_dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( configurations.size( ) % sot_dim == 0, "The configurations vector must be a scalar multiple of the second order tensor size (9 for 3D)" )

        floatVector gradient( sot_dim * sot_dim * num_incoming_configs, 0 );

        secondOrderTensor Fm, FpT;

        auto map = tardigradeHydra::getFixedSizeMatrixMap<floatType,3,3>(FpT.data());

        for ( unsigned int index = lowerIndex; index < upperIndex; index++ ){

            TARDIGRADE_ERROR_TOOLS_CATCH( Fm = getSubConfiguration( configurations, lowerIndex, index ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( FpT = getSubConfiguration( configurations, index + 1, upperIndex ) );
            new (&map) Eigen::Map< Eigen::Matrix< floatType, 3, 3 > >( FpT.data( ), 3, 3 );
            map = map.transpose( ).eval( );

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int I = 0; I < dim; I++ ){

                    for ( unsigned int a = 0; a < dim; a++ ){

                        for ( unsigned int A = 0; A < dim; A++ ){

                            gradient[ dim * num_incoming_configs * sot_dim * i + num_incoming_configs * sot_dim * I + sot_dim * index + dim * a + A ] = Fm[ dim * i + a ] * FpT[ dim * I + A ];

                        }

                    }

                }

            } 

        }

        return gradient;

    }

    secondOrderTensor hydraBase::getSubConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex ){
        /*!
         * Get a sub-configuration \f$\bf{F}^{sc}\f$ defined as
         *
         * \f$ F^{sc}_{iI} = F^{\text{lowerIndex}}_{i\hat{I}} F^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots F^{\text{upperIndex-1}}_{\bar{I}I} \f$
         * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
         * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
         *   Note, the configuration indicated by the index is NOT included in the sub-configuration
         */

        return getSubConfiguration( *get_configurations( ), lowerIndex, upperIndex );

    }

    secondOrderTensor hydraBase::getPrecedingConfiguration( const unsigned int &index ){
        /*!
         * Get the sub-configuration preceding but not including the index
         * 
         * \param &index: The index of the configuration immediately following the sub-configuration
         */

        return getSubConfiguration( 0, index );

    }

    secondOrderTensor hydraBase::getFollowingConfiguration( const unsigned int &index ){
        /*!
         * Get the sub-configuration following but not including the index
         * 
         * \param &index: The index of the current configuration immediately before the sub-configuration
         */

        return getSubConfiguration( index + 1, getNumConfigurations( ) );

    }

    secondOrderTensor hydraBase::getConfiguration( const unsigned int &index ){
        /*!
         * Get the configuration indicated by the provided index
         * 
         * \param &index: The index of the current configuration to be extracted
         */

        return getSubConfiguration( index, index + 1 );

    }

    secondOrderTensor hydraBase::getPreviousSubConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex ){
        /*!
         * Get a previous sub-configuration \f$\bf{F}^{sc}\f$ defined as
         *
         * \f$ F^{sc}_{iI} = F^{\text{lowerIndex}}_{i\hat{I}} F^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots F^{\text{upperIndex-1}}_{\bar{I}I} \f$
         * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
         * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
         *   Note, the configuration indicated by the index is NOT included in the sub-configuration
         */

        return getSubConfiguration( *get_previousConfigurations( ), lowerIndex, upperIndex );

    }

    secondOrderTensor hydraBase::getPreviousPrecedingConfiguration( const unsigned int &index ){
        /*!
         * Get the previous sub-configuration preceding but not including the index
         * 
         * \param &index: The index of the configuration immediately following the sub-configuration
         */

        return getPreviousSubConfiguration( 0, index );

    }

    secondOrderTensor hydraBase::getPreviousFollowingConfiguration( const unsigned int &index ){
        /*!
         * Get the previous sub-configuration following but not including the index
         * 
         * \param &index: The index of the current configuration immediately before the sub-configuration
         */

        return getPreviousSubConfiguration( index + 1, getNumConfigurations( ) );

    }

    secondOrderTensor hydraBase::getPreviousConfiguration( const unsigned int &index ){
        /*!
         * Get the previous configuration indicated by the provided index
         * 
         * \param &index: The index of the current configuration to be extracted
         */

        return getPreviousSubConfiguration( index, index + 1 );

    }

    floatVector hydraBase::getSubConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex ){
        /*!
         * Get the jacobian of a sub-configuration \f$\bf{F}^{sc}\f$ defined as
         *
         * \f$ F^{sc}_{iI} = F^{\text{lowerIndex}}_{i\hat{I}} F^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots F^{\text{upperIndex-1}}_{\bar{I}I} \f$
         * 
         * with respect to the current configurations.
         *
         * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
         * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
         *   Note, the configuration indicated by the index is NOT included in the sub-configuration
         */

        return getSubConfigurationJacobian( *get_configurations( ), lowerIndex, upperIndex );

    }

    floatVector hydraBase::getPrecedingConfigurationJacobian( const unsigned int &index ){
        /*!
         * Get the jacobian of the sub-configuration preceding but not including the index with respect to the current configurations.
         * 
         * \param &index: The index of the configuration immediately following the sub-configuration
         */

        return getSubConfigurationJacobian( 0, index );

    }

    floatVector hydraBase::getFollowingConfigurationJacobian( const unsigned int &index ){
        /*!
         * Get the jacobian of the sub-configuration following but not including the index with respect to the current configurations.
         * 
         * \param &index: The index of the current configuration immediately before the sub-configuration
         */

        return getSubConfigurationJacobian( index + 1, getNumConfigurations( ) );

    }

    floatVector hydraBase::getPreviousSubConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex ){
        /*!
         * Get the jacobian of a previous sub-configuration \f$\bf{F}^{sc}\f$ defined as
         *
         * \f$ F^{sc}_{iI} = F^{\text{lowerIndex}}_{i\hat{I}} F^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots F^{\text{upperIndex-1}}_{\bar{I}I} \f$
         * 
         * with respect to the previous configurations.
         *
         * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
         * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
         *   Note, the configuration indicated by the index is NOT included in the sub-configuration
         */

        return getSubConfigurationJacobian( *get_previousConfigurations( ), lowerIndex, upperIndex );

    }

    floatVector hydraBase::getPreviousPrecedingConfigurationJacobian( const unsigned int &index ){
        /*!
         * Get the jacobian of the previous sub-configuration preceding but not including the index with
         * respect to the previous configurations.
         * 
         * \param &index: The index of the configuration immediately following the sub-configuration
         */

        return getPreviousSubConfigurationJacobian( 0, index );

    }

    floatVector hydraBase::getPreviousFollowingConfigurationJacobian( const unsigned int &index ){
        /*!
         * Get the jacobian of the previous sub-configuration following but not including the index with
         * respect to the previous configurations
         * 
         * \param &index: The index of the current configuration immediately before the sub-configuration
         */

        return getPreviousSubConfigurationJacobian( index + 1, getNumConfigurations( ) );

    }

    void hydraBase::calculateFirstConfigurationJacobians( const floatVector &configurations, fourthOrderTensor &dC1dC, floatVector &dC1dCn ){
        /*!
         * Get the Jacobian of the first configuration w.r.t. the total mapping and the remaining configurations.
         * 
         * \param &configurations: The configurations which describe the mapping from the current to the reference configuration
         * \param &dC1dC: The Jacobian of the first entry w.r.t. the total
         * \param &dC1dCn: The Jacobian of the first entry w.r.t. the remaining terms
         * 
         * whre \f$C^n = C^2, C^3, \cdots \f$
         */

        constexpr unsigned int dim = 3;
        constexpr unsigned int sot_dim = dim * dim;
        auto num_configs = getNumConfigurations( );

        dC1dC  = secondOrderTensor( sot_dim * sot_dim, 0 );

        dC1dCn = floatVector( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

        secondOrderTensor fullConfiguration = getSubConfiguration( configurations, 0, num_configs );

        secondOrderTensor invCsc = getSubConfiguration( configurations, 1, num_configs );
        auto mat = tardigradeHydra::getFixedSizeMatrixMap<floatType, 3, 3>( invCsc.data( ) );
        mat = mat.inverse( ).eval( );

        fourthOrderTensor dInvCscdCsc = tardigradeVectorTools::computeFlatDInvADA( invCsc, dim, dim );
        auto map_dInvCscdCsc = tardigradeHydra::getFixedSizeMatrixMap<floatType, sot_dim, sot_dim>( dInvCscdCsc.data( ) );

        floatVector dCscdCs = getSubConfigurationJacobian( configurations, 1, num_configs );
        auto map_dCscdCs = tardigradeHydra::getDynamicSizeMatrixMap( dCscdCs.data( ), sot_dim, num_configs * sot_dim );

        floatVector dInvCscdCs( sot_dim * num_configs * sot_dim, 0 );
        auto map_dInvCscdCs = tardigradeHydra::getDynamicSizeMatrixMap( dInvCscdCs.data( ), sot_dim, num_configs * sot_dim ); 

        map_dInvCscdCs = ( map_dInvCscdCsc * map_dCscdCs ).eval( );

        // Compute the gradients
        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int barI = 0; barI < dim; barI++ ){

                for ( unsigned int A = 0; A < dim; A++ ){

                    dC1dC[ dim * sot_dim * i + sot_dim * barI + dim * i + A ] += invCsc[ dim * A + barI ];

                }
            }
        }
        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int J = 0; J < dim; J++ ){

                for ( unsigned int barI = 0; barI < dim; barI++ ){

                    for ( unsigned int indexaA = 0; indexaA < ( num_configs - 1 ) * sot_dim; indexaA++ ){

                        dC1dCn[ dim * ( num_configs - 1 ) * sot_dim * i + ( num_configs - 1 ) * sot_dim * barI + indexaA ]
                            += fullConfiguration[ dim * i + J ]
                             * dInvCscdCs[ dim * num_configs * sot_dim * J + num_configs * sot_dim * barI + indexaA + sot_dim ];

                    }

                }

            }

        }

    }

    void hydraBase::setFirstConfigurationJacobians( ){
        /*!
         * Set the Jacobians of the first configuration w.r.t. the total configuration and the remaining sub-configurations
         */

        auto dF1dF = get_SetDataStorage_dF1dF( );

        auto dF1dFn = get_SetDataStorage_dF1dFn( );

        calculateFirstConfigurationJacobians( *get_configurations( ), *dF1dF.value, *dF1dFn.value );

    }

    void hydraBase::setPreviousFirstConfigurationJacobians( ){
        /*!
         * Set the Jacobians of the previous first configuration w.r.t. the total configuration and the remaining sub-configurations
         */

        auto dF1dF = get_SetDataStorage_previousdF1dF( );

        auto dF1dFn = get_SetDataStorage_previousdF1dFn( );

        calculateFirstConfigurationJacobians( *get_previousConfigurations( ), *dF1dF.value, *dF1dFn.value );

    }

    void hydraBase::resetIterationData( ){
        /*!
         * Reset the iteration data to the new base state
         */

        for ( auto d = _iterationData.begin( ); d != _iterationData.end( ); d++ ){

            ( *d )->clear( );

        }

        _iterationData.clear( );

    }

    void hydraBase::resetNLStepData( ){
        /*!
         * Reset the nonlinear step data to the new base state
         */

        for ( auto d = _nlStepData.begin( ); d != _nlStepData.end( ); d++ ){

            ( *d )->clear( );

        }

        _nlStepData.clear( );

    }

    void hydraBase::setResidualClasses( ){
        /*!
         * Set the vectors for the residuals.
         * 
         * The expected form of the residual vector is cauchy stress, configurations,
         * state variables though only the residual on the cauchy stress must come
         * in this specific order if the deconstructSolutionVector function is redefined.
         * 
         * The user should define a vector of ResidualBase objects and use the
         * setResidualClasses( std::vector< ResidualBase > & ) function here.
         * 
         * The resulting residual should have the form
         * 
         * residual = { cauchyResidual, F2residual, ... Fnresidual, xiresidual1, xiresidual2, ... }
         * 
         * and can be formed by any number of residual classes. The first residual class must also
         * have the method `void getStress( )` defined which will return the current value
         * of the stress.
         */

    }

    void hydraBase::setResidualClasses( std::vector< ResidualBase<hydraBase>* > &residualClasses ){
        /*!
         * Set the residual classes
         * 
         * \param &residualClasses: A vector of residual classes which will be used to
         *     populate the residual and jacobian matrices for the non-linear solve
         */

        unsigned int numEquations = 0;

        _residualClasses.second = std::vector< ResidualBase<hydraBase>* >( residualClasses.size( ) );

        for ( auto c = residualClasses.begin( ); c != residualClasses.end( ); c++ ){

            numEquations += ( *c )->getNumEquations( );

            _residualClasses.second[ c - residualClasses.begin( ) ] = *c;

        }

        if ( numEquations != ( getNumConfigurations( ) * getConfigurationUnknownCount( ) + getNumNonLinearSolveStateVariables( ) ) ){

            std::string message = "The number of equations for the non-linear solve is not equal to the number of equations defined\n";
            message            += "  expected number of equations: " + std::to_string( getNumConfigurations( ) * getConfigurationUnknownCount( ) + getNumNonLinearSolveStateVariables( ) ) + "\n";
            message            += "  number of defined equations:  " + std::to_string( numEquations ) + "\n";

            TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

        }

        _residualClasses.first = true;

    }

    std::vector< ResidualBase<hydraBase>* >* hydraBase::getResidualClasses( ){
        /*!
         * Get a pointer to the vector of residual class pointers
         */

        if ( !_residualClasses.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setResidualClasses( ) );

        }

        return &_residualClasses.second;

    }

    void hydraBase::formNonLinearResidual( ){
        /*!
         * Form the residual
         */

        auto configurationUnknownCount = getConfigurationUnknownCount( );

        auto residualSize = getNumConfigurations( ) * configurationUnknownCount + getNumNonLinearSolveStateVariables( );

        auto numAdditionalDOF = getAdditionalDOF( )->size( );

        _residual.second = floatVector( residualSize, 0 );

        _jacobian.second = floatVector( residualSize * residualSize, 0 );

        _dRdF.second = floatVector( residualSize * configurationUnknownCount, 0 );

        _dRdT.second = floatVector( residualSize, 0 );

        _dRdAdditionalDOF.second = floatVector( residualSize * numAdditionalDOF, 0 );

        _additionalDerivatives.second.clear( );

        unsigned int offset = 0;

        setCurrentResidualIndexMeaningful( true );
        for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); ++residual_ptr ){
            setCurrentResidualIndex( residual_ptr - getResidualClasses( )->begin() );

            // Extract the terms

            const floatVector* localResidual;
            TARDIGRADE_ERROR_TOOLS_CATCH( localResidual = ( *residual_ptr )->getResidual( ) );

            // Check the contributions to make sure they are consistent sizes

            TARDIGRADE_ERROR_TOOLS_CHECK( localResidual->size( ) == ( *residual_ptr )->getNumEquations( ),
                  "The residual for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n"
                + "  expected: " + std::to_string( ( *residual_ptr )->getNumEquations( ) ) + "\n"
                + "  actual:   " + std::to_string( localResidual->size( ) ) + "\n"
            )

            // Store the values in the global quantities

            // Copy over the values of the local vector to the global structures
            std::copy( localResidual->begin( ), localResidual->end( ), _residual.second.begin( ) + offset );

            offset += ( *residual_ptr )->getNumEquations( );

        }
        setCurrentResidualIndexMeaningful( false );

        // Allow the residuals to modify the global residual if needed
        setAllowModifyGlobalResidual( true );
        setCurrentResidualIndexMeaningful( true );
        for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); ++residual_ptr ){
            setCurrentResidualIndex( residual_ptr - getResidualClasses( )->begin( ) );

            ( *residual_ptr )->modifyGlobalResidual( );

        }
        setCurrentResidualIndexMeaningful( false );
        setAllowModifyGlobalResidual( false );

        _residual.first = true;

        addIterationData( &_residual );

    }

    void hydraBase::formNonLinearDerivatives( ){
        /*!
         * Form the jacobian and gradient matrices
         */

        auto configurationUnknownCount = getConfigurationUnknownCount( );

        auto residualSize = getNumConfigurations( ) * configurationUnknownCount + getNumNonLinearSolveStateVariables( );

        auto numAdditionalDOF = getAdditionalDOF( )->size( );

        _jacobian.second = floatVector( residualSize * residualSize, 0 );

        _dRdF.second = floatVector( residualSize * configurationUnknownCount, 0 );

        _dRdT.second = floatVector( residualSize, 0 );

        _dRdAdditionalDOF.second = floatVector( residualSize * numAdditionalDOF, 0 );

        _additionalDerivatives.second.clear( );

        unsigned int offset = 0;

        unsigned int numAdditionalDerivatives = 0;

        setCurrentResidualIndexMeaningful( true );
        for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); ++residual_ptr ){
            setCurrentResidualIndex( residual_ptr - getResidualClasses( )->begin( ) );

            // Extract the terms

            const floatVector* localJacobian;
            TARDIGRADE_ERROR_TOOLS_CATCH( localJacobian = ( *residual_ptr )->getJacobian( ) );

            const floatVector* localdRdF;
            TARDIGRADE_ERROR_TOOLS_CATCH( localdRdF = ( *residual_ptr )->getdRdF( ) );

            const floatVector* localdRdT;
            TARDIGRADE_ERROR_TOOLS_CATCH( localdRdT = ( *residual_ptr )->getdRdT( ) );

            const floatVector* localdRdAdditionalDOF;
            TARDIGRADE_ERROR_TOOLS_CATCH( localdRdAdditionalDOF = ( *residual_ptr )->getdRdAdditionalDOF( ) );

            const floatVector* localAdditionalDerivatives;
            TARDIGRADE_ERROR_TOOLS_CATCH( localAdditionalDerivatives = ( *residual_ptr )->getAdditionalDerivatives( ) );

            // Check the contributions to make sure they are consistent sizes

            TARDIGRADE_ERROR_TOOLS_CHECK( localJacobian->size( ) == ( *residual_ptr )->getNumEquations( ) * residualSize,
                  "The jacobian for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n"
                + "  expected: " + std::to_string( ( *residual_ptr )->getNumEquations( ) * residualSize ) + "\n"
                + "  actual:   " + std::to_string( localJacobian->size( ) ) + "\n"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK( localdRdF->size( ) == ( *residual_ptr )->getNumEquations( ) * configurationUnknownCount,
                  "dRdF for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n"
                + "  expected: " + std::to_string( ( *residual_ptr )->getNumEquations( ) * configurationUnknownCount ) + "\n"
                + "  actual:   " + std::to_string( localdRdF->size( ) ) + "\n"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK( localdRdT->size( ) == ( *residual_ptr )->getNumEquations( ),
                  "dRdT for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n"
                + "  expected: " + std::to_string( ( *residual_ptr )->getNumEquations( ) ) + "\n"
                + "  actual:   " + std::to_string( localdRdT->size( ) ) + "\n"
            )

            if ( localdRdAdditionalDOF->size( ) != 0 ){

                TARDIGRADE_ERROR_TOOLS_CHECK( localdRdAdditionalDOF->size( ) == ( ( *residual_ptr )->getNumEquations( ) * numAdditionalDOF ),
                                              "dRdAdditionalDOF for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n"
                                            + "  expected: " + std::to_string( ( *residual_ptr )->getNumEquations( ) * numAdditionalDOF ) + "\n"
                                            + "  actual  : " + std::to_string( localdRdAdditionalDOF->size( ) ) + "\n"
                )

                std::copy ( localdRdAdditionalDOF->begin( ), localdRdAdditionalDOF->end( ), _dRdAdditionalDOF.second.begin( ) + numAdditionalDOF * offset );

            }

            if ( localAdditionalDerivatives->size( ) != 0 ){

                if ( ( *localAdditionalDerivatives ).size( ) != ( *residual_ptr )->getNumEquations( ) * numAdditionalDerivatives ){
    
                    if ( numAdditionalDerivatives == 0 ){
    
                        numAdditionalDerivatives = ( *localAdditionalDerivatives ).size( ) / ( *residual_ptr )->getNumEquations( );
    
                        _additionalDerivatives.second = floatVector( residualSize * numAdditionalDerivatives, 0 );
    
                    }
                    else{
    
                        std::string message = "The additional derivatives for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " are not the expected length as determined from the first residual\n";
                        message            += "  expected: " + std::to_string( numAdditionalDerivatives ) + "\n";
                        message            += "  actual:   " + std::to_string( ( *localAdditionalDerivatives ).size( ) ) + "\n";
    
                        TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );
    
                    }
    
                }

            }

            // Store the values in the global quantities

            // Copy over the values of the local vector to the global structures
            std::copy( localJacobian->begin( ), localJacobian->end( ), _jacobian.second.begin( ) + residualSize * offset );

            std::copy( localdRdF->begin( ), localdRdF->end( ), _dRdF.second.begin( ) + configurationUnknownCount * offset );

            std::copy( localdRdT->begin( ), localdRdT->end( ), _dRdT.second.begin( ) + offset );

            std::copy( localAdditionalDerivatives->begin( ), localAdditionalDerivatives->end( ), _additionalDerivatives.second.begin( ) + numAdditionalDerivatives * offset );

            offset += ( *residual_ptr )->getNumEquations( );

        }
        setCurrentResidualIndexMeaningful( false );

        setAllowModifyGlobalJacobian( true );
        setAllowModifyGlobaldRdT( true );
        setAllowModifyGlobaldRdF( true );
        setAllowModifyGlobaldRdAdditionalDOF( true );
        setCurrentResidualIndexMeaningful( true );
        for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); residual_ptr++ ){

            setCurrentResidualIndexMeaningful( residual_ptr - getResidualClasses( )->begin( ) );

            ( *residual_ptr )->modifyGlobalJacobian( );
            ( *residual_ptr )->modifyGlobaldRdT( );
            ( *residual_ptr )->modifyGlobaldRdF( );
            ( *residual_ptr )->modifyGlobaldRdAdditionalDOF( );

        }
        setCurrentResidualIndexMeaningful( false );
        setAllowModifyGlobalJacobian( false );
        setAllowModifyGlobaldRdT( false );
        setAllowModifyGlobaldRdF( false );
        setAllowModifyGlobaldRdAdditionalDOF( false );

        _jacobian.first = true;

        _dRdF.first = true;

        _dRdT.first = true;

        _dRdAdditionalDOF.first = true;

        _additionalDerivatives.first = true;

        addIterationData( &_jacobian );

        addIterationData( &_dRdF );

        addIterationData( &_dRdT );

        addIterationData( &_dRdAdditionalDOF );

        addIterationData( &_additionalDerivatives );

    }

    const floatVector* hydraBase::getResidual( ){
        /*!
         * Get the residual vector for the non-linear problem
         */

        if ( !_residual.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( formNonLinearResidual( ) );

        }

        return &_residual.second;

    }

    const floatVector* hydraBase::getFlatJacobian( ){
        /*!
         * Get the flattened row-major jacobian for the non-linear problem
         */

        if ( !_jacobian.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( formNonLinearDerivatives( ) );

        }

        return &_jacobian.second;

    }

    floatMatrix hydraBase::getJacobian( ){
        /*!
         * Get the jacobian for the non-linear problem
         */

        return tardigradeVectorTools::inflate( *getFlatJacobian( ), getResidual( )->size( ), getResidual( )->size( ) );

    }

    const floatVector* hydraBase::getFlatdRdF( ){
        /*!
         * Get the flattened row-major dRdF for the non-linear problem
         */

        if ( !_dRdF.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( formNonLinearDerivatives( ) );

        }

        return &_dRdF.second;

    }

    floatMatrix hydraBase::getdRdF( ){
        /*!
         * Get dRdF for the non-linear problem
         */

        return tardigradeVectorTools::inflate( *getFlatdRdF( ), getResidual( )->size( ), getSOTDimension( ) );
    }

    const floatVector* hydraBase::getFlatdRdAdditionalDOF( ){
        /*!
         * Get the flattened row-major dRdAdditional for the non-linear problem
         */

        if ( !_dRdAdditionalDOF.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( formNonLinearDerivatives( ) );

        }

        return &_dRdAdditionalDOF.second;

    }

    floatMatrix hydraBase::getdRdAdditionalDOF( ){
        /*!
         * Get dRdAdditionalDOF for the non-linear problem
         */

        return tardigradeVectorTools::inflate( *getFlatdRdAdditionalDOF( ), getResidual( )->size( ), getAdditionalDOF( )->size( ) );
    }

    const floatVector* hydraBase::getdRdT( ){
        /*!
         * Get dRdT for the non-linear problem
         */

        if ( !_dRdT.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( formNonLinearDerivatives( ) );

        }

        return &_dRdT.second;

    }

    const floatVector* hydraBase::getFlatAdditionalDerivatives( ){
        /*!
         * Get the flattened row-major additional derivatives for the non-linear problem
         */

        if ( !_additionalDerivatives.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( formNonLinearDerivatives( ) );

        }

        return &_additionalDerivatives.second;

    }

    floatMatrix hydraBase::getAdditionalDerivatives( ){
        /*!
         * Get the additional derivatives for the non-linear problem
         */

        if ( getFlatAdditionalDerivatives( )->size( ) > 0 ){

            return tardigradeVectorTools::inflate( *getFlatAdditionalDerivatives( ), getResidual( )->size( ), getFlatAdditionalDerivatives( )->size( ) / getResidual( )->size( ) );

        }

        return floatMatrix( 0, floatVector( 0, 0 ) );

    }

    const floatVector* hydraBase::getStress( ){
        /*!
         * Get the stress
         */

        if ( !_stress.first ){

            if ( getResidualClasses( )->size( ) == 0 ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "No residual classes are defined." ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( _stress.second = *( *getResidualClasses( ) )[ 0 ]->getStress( ) );

            if ( getViscoplasticDampingSet( ) ){

                auto previouslyConvergedStress = getPreviouslyConvergedStress( );

                for ( auto v = std::begin( *previouslyConvergedStress ); v != std::end( *previouslyConvergedStress ); ++v ){

                    _stress.second[ v - std::begin( *previouslyConvergedStress ) ] -= getViscoplasticDamping( ) * ( *v );

                }

            }

            _stress.first = true;

            addIterationData( &_stress );

        }

        return &_stress.second;

    }

    const floatVector* hydraBase::getPreviousStress( ){
        /*!
         * Get the previous value of the stress
         */

        if ( !_previousStress.first ){

            if ( getResidualClasses( )->size( ) == 0 ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "No residual classes are defined." ) );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH( _previousStress.second = *( *getResidualClasses( ) )[ 0 ]->getPreviousStress( ) );

            _previousStress.first = true;

        }

        return &_previousStress.second;

    }

    const floatVector* hydraBase::getPreviouslyConvergedStress( ){
        /*!
         * Get the previously converged stress value
         */

        if ( !_previouslyConvergedStress.first ){

            setPreviouslyConvergedStress( *getPreviousStress( ) );

        }

        return &_previouslyConvergedStress.second;

    }

    void hydraBase::setViscoplasticDamping( const floatType &factor ){
        /*!
         * Use viscoplastic damping to reduce the current stress levels
         *
         * This should be used with care and likely only in the context of a
         * relaxed solve where it will be removed as the solution is obtained.
         *
         * \param &factor: The fraction of the difference between the last converged
         *     stress (which may be the previous stress) and the trial stress which
         *     will be suppressed
         */

        _viscoplastic_damping_factor = factor;
        _viscoplastic_damping_set    = true;

    }

    void hydraBase::clearViscoplasticDamping( ){
        /*!
         * Clear the viscoplastic damping
         */

        _viscoplastic_damping_factor = 0.;
        _viscoplastic_damping_set    = false;
    }

    const bool hydraBase::getViscoplasticDampingSet( ){
        /*!
         * Get whether the viscoplastic damping has been set or not
         */

        return _viscoplastic_damping_set;

    }

    void hydraBase::resetProblem( ){
        /*!
         * Reset the problem to the initial state
         */

        decomposeStateVariableVector( );

        initializeUnknownVector( );

    }

    void hydraBase::initializeUnknownVector( ){
        /*!
         * Initialize the unknown vector for the non-linear solve.
         * 
         * \f$X = \left\{ \bf{\sigma}, \bf{F}^2, \bf{F}^3, ..., \bf{F}n, \xi^1, \xi^2, ..., \xi^m \right\} \f$
         * 
         * It is assumed that the first residual calculation also has a method `void getStress( )`
         * which returns a pointer to the current value of the stress.
         */

        const unsigned int sot_dim = getSOTDimension( );

        const floatVector *cauchyStress;
        TARDIGRADE_ERROR_TOOLS_CATCH( cauchyStress = getStress( ) );

        const floatVector *configurations = get_configurations( );

        const unsigned int num_local_configs = configurations->size( ) / sot_dim;

        TARDIGRADE_ERROR_TOOLS_CHECK( configurations->size( ) % num_local_configs == 0, "The size of the configurations vector must be a scalar multiple of the second order tensor size" )

        const floatVector *nonLinearSolveStateVariables = get_nonLinearSolveStateVariables( );

        floatVector X( getNumUnknowns( ), 0 );

        std::copy( std::begin( *cauchyStress ), std::end( *cauchyStress ), std::begin( X ) );

        std::copy( std::begin( *configurations ) + sot_dim, std::end( *configurations ), std::begin( X ) + sot_dim );

        std::copy( std::begin( *nonLinearSolveStateVariables ), std::end( *nonLinearSolveStateVariables ), std::begin( X ) + num_local_configs * sot_dim );

        bool resetRequired = false;

        setCurrentResidualIndexMeaningful( true );
        for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); residual_ptr++ ){
            setCurrentResidualIndex( residual_ptr - getResidualClasses( )->begin( ) );

            std::vector< unsigned int > indices;

            std::vector< floatType > values;

            ( *residual_ptr )->suggestInitialIterateValues( indices, values );

            if ( indices.size( ) > 0 ){
                resetRequired = true;
            }

            for ( auto i = indices.begin( ); i != indices.end( ); i++ ){

                X[ *i ] = values[ ( unsigned int )( i - indices.begin( ) ) ];

            }

        }
        setCurrentResidualIndexMeaningful( false );

        if ( resetRequired ){

            updateUnknownVector( X );

        }
        else{

            setX( X );

        }

    }

    const floatVector* hydraBase::getUnknownVector( ){
        /*!
         * Get the unknown vector
         */

        if ( !_X.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( initializeUnknownVector( ) );

        }

        return &_X.second;

    }

    void hydraBase::updateUnknownVector( const floatVector &newUnknownVector ){
        /*!
         * Update the unknown vector
         * 
         * \param &newUnknownVector: The new unknown vector
         */

        // Project the trial unknown vector to the allowable space
        floatVector trialX = newUnknownVector;
        floatVector Xp;

        if ( !_X.first ){

            Xp = trialX;

        }
        else{

            Xp = *getUnknownVector( );

        }

        setCurrentResidualIndexMeaningful( true );
        for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); residual_ptr++ ){
            setCurrentResidualIndex( residual_ptr - getResidualClasses( )->begin( ) );
            if ( ( *residual_ptr )->getUseProjection( ) ){
                ( *residual_ptr )->projectSuggestedX( trialX, Xp );
            }
        }
        setCurrentResidualIndexMeaningful( false );

        // Reset all of the iteration data
        resetIterationData( );

        // Set the unknown vector
        setX( trialX );

        // Decompose the unknown vector and update the state
        TARDIGRADE_ERROR_TOOLS_CATCH( decomposeUnknownVector( ) );

    }

    void hydraBase::callResidualPreSubcycler( ){
        /*!
         * Signal to the residuals that we are entering the subcycler
         */

        setCurrentResidualIndexMeaningful( true );

        for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); ++residual_ptr ){

            setCurrentResidualIndex( residual_ptr - getResidualClasses( )->begin( ) );

            try{

                ( *residual_ptr )->preSubcycler( );

            }
            catch( std::exception &e ){

                if ( getFailureVerbosityLevel( ) > 0 ){

                    addToFailureOutput( "Failure in residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + "\n" );
                    std::string message;
                    tardigradeErrorTools::captureNestedExceptions( e, message );
                    addToFailureOutput( message );

                }

                throw;

            }

        }

        setCurrentResidualIndexMeaningful( false );

    }

    void hydraBase::callResidualPostSubcyclerSuccess( ){
        /*!
         * Signal to the residuals that we have a successful subcycle increment
         */

        setCurrentResidualIndexMeaningful( true );

        for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); ++residual_ptr ){

            setCurrentResidualIndex( residual_ptr - getResidualClasses( )->begin( ) );

            try{

                ( *residual_ptr )->postSubcyclerSuccess( );

            }
            catch( std::exception &e ){

                if ( getFailureVerbosityLevel( ) > 0 ){

                    addToFailureOutput( "Failure in residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + "\n" );
                    std::string message;
                    tardigradeErrorTools::captureNestedExceptions( e, message );
                    addToFailureOutput( message );

                }

                throw;

            }

        }

        setCurrentResidualIndexMeaningful( false );

    }

    void hydraBase::callResidualPostSubcyclerFailure( ){
        /*!
         * Signal to the residuals that we have a failed subcycle increment
         */

        setCurrentResidualIndexMeaningful( true );

        for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); ++residual_ptr ){

            setCurrentResidualIndex( residual_ptr - getResidualClasses( )->begin( ) );

            try{

                ( *residual_ptr )->postSubcyclerFailure( );

            }
            catch( std::exception &e ){

                if ( getFailureVerbosityLevel( ) > 0 ){

                    addToFailureOutput( "Failure in residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + "\n" );
                    std::string message;
                    tardigradeErrorTools::captureNestedExceptions( e, message );
                    addToFailureOutput( message );

                }

                throw;

            }

        }

        setCurrentResidualIndexMeaningful( false );

    }

    bool hydraBase::callResidualRelaxedStepFailure( ){
        /*!
         * Signal to the residuals that we have a failed relaxed solve step and
         * determine if a new relaxed step should be taken
         */

        bool attempt_relaxed_step = false;

        setCurrentResidualIndexMeaningful( true );

        for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); ++residual_ptr ){

            setCurrentResidualIndex( residual_ptr - getResidualClasses( )->begin( ) );

            try{

                auto val = ( *residual_ptr )->relaxedStepFailure( );

                attempt_relaxed_step = attempt_relaxed_step || val;

            }
            catch( std::exception &e ){

                if ( getFailureVerbosityLevel( ) > 0 ){

                    addToFailureOutput( "Failure in residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + "\n" );
                    std::string message;
                    tardigradeErrorTools::captureNestedExceptions( e, message );
                    addToFailureOutput( message );

                }

                throw;

            }

        }

        setCurrentResidualIndexMeaningful( false );

        return attempt_relaxed_step;

    }

    /*!
     * Solve the non-linear problem by relaxing difficult sub-problems
     * to achieve a series of solutions.
     */
    void hydraBase::performRelaxedSolve( ){

        // TEMP: Remove when we extract this to RelaxedSolver.solve
        auto local_solver = dynamic_cast<RelaxedSolver*>(solver);
        TARDIGRADE_ERROR_TOOLS_CHECK(local_solver,"The solver must be a relaxed solver");
        // END TEMP

        local_solver->resetRelaxedIteration( );

        // Initialize the residuals
        local_solver->hydra->setCurrentResidualIndexMeaningful( true );
        for ( auto residual = std::begin( *( local_solver->hydra->getResidualClasses( ) ) ); residual != std::end( *( local_solver->hydra->getResidualClasses( ) ) ); ++residual ){
            setCurrentResidualIndexMeaningful( residual - std::begin( *( local_solver->hydra->getResidualClasses( ) ) ) );

            // Prepare the residuals to take a relaxed step
            ( *residual )->setupRelaxedStep( local_solver->getRelaxedIteration( ) );

        }
        local_solver->hydra->setCurrentResidualIndexMeaningful( false );

        while ( local_solver->getRelaxedIteration( ) < getMaxRelaxedIterations( ) ){

            if ( getFailureVerbosityLevel( ) > 0 ){
                addToFailureOutput( "\n\n###  relaxed iteration: " );
                addToFailureOutput( local_solver->getRelaxedIteration( ) );
                addToFailureOutput( "\n\n" );
            }
            // Solve the non-linear problem
            try{

                solver->solve( );

                // Check if the relaxation has converged
                bool relaxedConverged = true;
                setCurrentResidualIndexMeaningful( true );
                for ( auto residual = getResidualClasses( )->begin( ); residual != getResidualClasses( )->end( ); residual++ ){
                    setCurrentResidualIndex( residual - getResidualClasses( )->begin( ) );

                    if ( !( *residual )->checkRelaxedConvergence( ) ){

                        relaxedConverged = false;
                        break;

                    }

                }
                setCurrentResidualIndexMeaningful( false );

                if ( relaxedConverged ){

                    return;

                }

            }
            catch( convergence_error &e ){

                if ( !callResidualRelaxedStepFailure( ) ){

                    throw;

                }

            }
            catch( std::exception &e ){

                throw;

            }

            // Use the current unknown vector as the initial estimate
            if ( solver->getInitializeUnknownVector( ) ){

                solver->setInitializeUnknownVector( false );

            }

            local_solver->incrementRelaxedIteration( );

            // Initialize the residuals
            setCurrentResidualIndexMeaningful( true );
            for ( auto residual = getResidualClasses( )->begin( ); residual != getResidualClasses( )->end( ); residual++ ){
                setCurrentResidualIndex( residual - getResidualClasses( )->begin( ) );

                // Prepare the residuals to take a relaxed step
                ( *residual )->setupRelaxedStep( local_solver->getRelaxedIteration( ) );

            }
            setCurrentResidualIndexMeaningful( false );

            // Reset hydra
            updateUnknownVector( *getUnknownVector( ) ); //This allows for the relaxed to change the projection and adjust the decomposition
            solver->resetIterations( );

        }

        if ( local_solver->getRelaxedIteration( ) >= getMaxRelaxedIterations( ) ){

            throw convergence_error( "Failure in relaxed solve:\n  scale_factor: " + std::to_string( getScaleFactor( ) ) );

        }

    }

    void hydraBase::setScaledQuantities( ){
        /*!
         * Set the scaled quantities
         */

        _scaled_time = ( _scale_factor - 1 ) * _deltaTime + _time;

        _scaled_deltaTime = _scale_factor * _deltaTime;

        _scaled_temperature = _scale_factor * ( _temperature - _previousTemperature ) + _previousTemperature;

        _scaled_deformationGradient = _scale_factor * ( _deformationGradient - _previousDeformationGradient ) + _previousDeformationGradient;

        _scaled_additionalDOF = _scale_factor * ( _additionalDOF - _previousAdditionalDOF ) + _previousAdditionalDOF;

    }

    void hydraBase::setScaleFactor( const floatType &value ){
        /*!
         * Set the value of the scale factor. Will automatically re-calculate the deformation and trial stresses
         * 
         * \param &value: The value of the scale factor
         */

        // Update the scale factor
        _scale_factor = value;

        // Update the scaled quantities
        setScaledQuantities( );

        // Copy the current unknown vector
        floatVector unknownVector = *getUnknownVector( );

        // Reset the iteration data
        resetIterationData( );

        // Set the unknown vector
        setX( unknownVector );

        // Update the deformation quantities
        updateConfigurationsFromUnknownVector( );

        // Compute the new trial stress
        std::copy( getStress( )->begin( ), getStress( )->end( ), unknownVector.begin( ) );

        // Re-set the unknown vector
        setX( unknownVector );

        // Extract the stress
        extractStress( );

    }

    const bool hydraBase::allowStepGrowth( const unsigned int &num_good ){
        /*!
         * Function to determine if we can increase the step-size for the sub-cycler
         * 
         * \param &num_good: The number of good increments since the last failure
         */

        if ( num_good >= getNumGoodControl( ) ){

            return true;

        }

        return false;

    }

    void hydraBase::evaluate( const bool &use_subcycler ){
        /*!
         * Solver the non-linear problem and update the variables
         * 
         * \param &use_subcycler: Flag for if the subcycler should be used for difficult analyses (defaults to false)
         */

        try{

            evaluateInternal( );

            return;

        }
        catch( std::exception &e ){

            if ( !use_subcycler ){

                throw;

            }

            if ( getFailureVerbosityLevel( ) > 0 ){
                addToFailureOutput( "\n\n" );
                addToFailureOutput( "#########################################\n" );
                addToFailureOutput( "###        ENTERING SUB-CYCLER        ###\n" );
                addToFailureOutput( "#########################################\n" );
                addToFailureOutput( "\n\n" );
            }

            floatType sp = 0.0;

            floatType ds = getCutbackFactor( );

            unsigned int num_good = 0;

            callResidualPreSubcycler( );

            resetProblem( );

            while ( sp < 1.0 ){

                try{

                    if ( getFailureVerbosityLevel( ) > 0 ){
                        addToFailureOutput( "\n\n" );
                        addToFailureOutput( "######### PSEUDO-TIME INCREMENT #########\n" );
                        addToFailureOutput( "\n\n    sp, ds: " + std::to_string( sp ) + ", " + std::to_string( ds ) );
                        addToFailureOutput( "\n" );
                    }

                    setScaleFactor( sp + ds ); // Update the scaling factor

                    solver->resetIterations( ); // Reset the non-linear iteration count

                    evaluateInternal( ); // Try to solve the non-linear problem

                    setPreviouslyConvergedStress( *getStress( ) ); // Set the previously converged stress

                    callResidualPostSubcyclerSuccess( ); // Let the residuals know the subcycle step was successful

                    sp += ds; // Update the pseudo-time

                    num_good++; // Update the number of good iterations

                    // Grow the step if possible
                    if ( allowStepGrowth( num_good ) ){

                        ds *= getGrowthFactor( );

                    }

                    // Make sure s will be less than or equal to 1
                    if ( sp + ds > 1.0 ){

                        ds = 1.0 - sp;

                    }

                }
                catch( std::exception &e ){

                    callResidualPostSubcyclerFailure( ); // Let the residuals know the subcycle step failed

                    // Reduce the time-step and try again
                    num_good = 0;

                    ds *= getCutbackFactor( );

                    if ( getUseRelaxedSolve( ) ){

                        _initialX = _prerelaxed_initialX;

                    }

                    setX( _initialX ); // Reset X to the last good point

                    if ( ds < getMinDS( ) ){

                        throw;

                    }

                }

            }

        }

    }

    void hydraBase::evaluateInternal( ){
        /*!
         * Solve the non-linear problem with the current scaling and update the variables
         */

        // Reset the counters for the number of steps being performed
        solver->step->resetNumNewton( );

        solver->step->resetNumLS( );

        solver->step->resetNumGrad( );

        solver->setRankDeficientError( false );

        try{

            solver->solve( );

        }
        catch( const convergence_error &e ){

            if ( getUseRelaxedSolve( ) ){

                if ( getFailureVerbosityLevel( ) > 0 ){
                    addToFailureOutput( "Failure in conventional solve. Starting relaxed solve.\n" );
                }

                try{

                    solver->resetIterations( );
                    _prerelaxed_initialX = _initialX;
                    updateUnknownVector( _initialX );
                    performRelaxedSolve( );

                }
                catch( const convergence_error &e ){

                    throw;

                }
                catch( std::exception *e ){

                    TARDIGRADE_ERROR_TOOLS_CATCH( throw; )

                }

            }
            else{
                //Try a Levenberg-Marquardt solve if there is a convergence error
                solver->setRankDeficientError( false );
    
                solver->step->setUseLevenbergMarquardt( true );
    
                // Turn on projection
                setCurrentResidualIndexMeaningful( true );
                for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); residual_ptr++ ){
                    setCurrentResidualIndex( residual_ptr - getResidualClasses( )->begin( ) );
                    ( *residual_ptr )->setUseProjection( true );
                }
                setCurrentResidualIndexMeaningful( false );
    
                solver->resetIterations( );
                updateUnknownVector( _initialX );
    
                try{
    
                    solver->solve( );
    
                }
                catch( const convergence_error &e ){
    
                    throw;
    
                }
                catch( std::exception &e ){
    
                    solver->step->setUseLevenbergMarquardt( false );
    
                    TARDIGRADE_ERROR_TOOLS_CATCH( throw; )
    
                }

            }

        }
        catch( std::exception &e ){

            TARDIGRADE_ERROR_TOOLS_CATCH( throw; )

        }

    }

    void hydraBase::computeTangents( ){
        /*!
         * Compute the values of the consistent tangents
         */

        //Form the solver based on the current value of the jacobian
        auto Amat = tardigradeHydra::getDynamicSizeMatrixMap( getFlatJacobian( )->data( ), getResidual( )->size( ), getResidual( )->size( ) );

        // Form the maps for dXdF
        auto dRdFmat = tardigradeHydra::getDynamicSizeMatrixMap( getFlatdRdF( )->data( ), getResidual( )->size( ), getConfigurationUnknownCount( ) );

        _flatdXdF.second = floatVector( getNumUnknowns( ) * getConfigurationUnknownCount( ) );
        auto dXdFmat = tardigradeHydra::getDynamicSizeMatrixMap( _flatdXdF.second.data( ), getNumUnknowns( ), getConfigurationUnknownCount( ) );

        // Form the maps for dXdT
        auto dRdTmat = tardigradeHydra::getDynamicSizeVectorMap( getdRdT( )->data( ), getResidual( )->size( ) );

        _flatdXdT.second = floatVector( getNumUnknowns( ) );
        auto dXdTmat = tardigradeHydra::getDynamicSizeVectorMap( _flatdXdT.second.data( ), getNumUnknowns( ) );

        // Solve
        tardigradeVectorTools::solverType< floatType > linear_solver;

        if ( solver->preconditioner->getUsePreconditioner( ) ){

            if( solver->preconditioner->getPreconditionerIsDiagonal( ) ){

                auto p_map = tardigradeHydra::getDynamicSizeVectorMap( solver->preconditioner->getFlatPreconditioner( )->data( ), getResidual( )->size( ) );

                linear_solver = tardigradeVectorTools::solverType< floatType >( p_map.asDiagonal( ) * Amat );

                dXdFmat = -linear_solver.solve( p_map.asDiagonal( ) * dRdFmat );

                dXdTmat = -linear_solver.solve( p_map.asDiagonal( ) * dRdTmat );

            }
            else{

                auto p_map = tardigradeHydra::getDynamicSizeMatrixMap( solver->preconditioner->getFlatPreconditioner( )->data( ), getResidual( )->size( ), getResidual( )->size( ) );

                linear_solver = tardigradeVectorTools::solverType< floatType >( p_map * Amat );

                dXdFmat = -linear_solver.solve( p_map * dRdFmat );

                dXdTmat = -linear_solver.solve( p_map * dRdTmat );

            }

        }
        else{

            linear_solver = tardigradeVectorTools::solverType< floatType >( Amat );

            dXdFmat = -linear_solver.solve( dRdFmat );

            dXdTmat = -linear_solver.solve( dRdTmat );

        }

        unsigned int rank = linear_solver.rank( );

        TARDIGRADE_ERROR_TOOLS_CATCH(

            if ( solver->getRankDeficientError( ) && ( rank != getResidual( )->size( ) ) ){

                throw convergence_error( "The Jacobian is not full rank" );

            }

        )

        _flatdXdF.first = true;

        _flatdXdT.first = true;

    }

    void hydraBase::computedXdAdditionalDOF( ){
        /*!
         * Compute the consistent tangent w.r.t. the additional dof
         */

        //Form the solver based on the current value of the jacobian
        auto Amat = tardigradeHydra::getDynamicSizeMatrixMap( getFlatJacobian( )->data( ), getResidual( )->size( ), getResidual( )->size( ) );

        auto dRdAdditionalDOF = tardigradeHydra::getDynamicSizeMatrixMap( getFlatdRdAdditionalDOF( )->data( ), getResidual( )->size( ), getAdditionalDOF( )->size( ) );

        // Form the map for dXdF
        _flatdXdAdditionalDOF.second = floatVector( getNumUnknowns( ) * getAdditionalDOF( )->size( ), 0 );
        auto dXdAdditionalDOF = tardigradeHydra::getDynamicSizeMatrixMap( _flatdXdAdditionalDOF.second.data( ), getNumUnknowns( ), getAdditionalDOF( )->size( ) );

        // Solve
        tardigradeVectorTools::solverType< floatType > linear_solver;

        if ( solver->preconditioner->getUsePreconditioner( ) ){

            if( solver->preconditioner->getPreconditionerIsDiagonal( ) ){

                auto p_map = tardigradeHydra::getDynamicSizeVectorMap( solver->preconditioner->getFlatPreconditioner( )->data( ), getResidual( )->size( ) );

                linear_solver = tardigradeVectorTools::solverType< floatType >( p_map.asDiagonal( ) * Amat );

                dXdAdditionalDOF = -linear_solver.solve( p_map.asDiagonal( ) * dRdAdditionalDOF );

            }
            else{

                auto p_map = tardigradeHydra::getDynamicSizeMatrixMap( solver->preconditioner->getFlatPreconditioner( )->data( ), getResidual( )->size( ), getResidual( )->size( ) );

                linear_solver = tardigradeVectorTools::solverType< floatType >( p_map * Amat );

                dXdAdditionalDOF = -linear_solver.solve( p_map * dRdAdditionalDOF );

            }

        }
        else{

            linear_solver = tardigradeVectorTools::solverType< floatType >( Amat );

            dXdAdditionalDOF = -linear_solver.solve( dRdAdditionalDOF );

        }

        _flatdXdAdditionalDOF.first = true;

    }

    const floatVector *hydraBase::getFlatdXdF( ){
        /*!
         * Get the total derivative of X w.r.t. the deformation in row-major format
         */

        if ( !_flatdXdF.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( computeTangents( ) )

        }

        return &_flatdXdF.second;

    }

    const floatVector *hydraBase::getFlatdXdT( ){
        /*!
         * Get the total derivative of X w.r.t. the temperature in row-major format
         */

        if ( !_flatdXdT.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( computeTangents( ) )

        }

        return &_flatdXdT.second;

    }

    const floatVector *hydraBase::getFlatdXdAdditionalDOF( ){
        /*!
         * Get the total derivative of X w.r.t. the additional degrees of freedom
         */

        if ( !_flatdXdAdditionalDOF.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( computedXdAdditionalDOF( ) );

        }

        return &_flatdXdAdditionalDOF.second;

    }

    void hydraBase::setConstraints( ){
        /*!
         * Set the constraint values
         */

        const unsigned int numConstraints = getNumConstraints( );

        auto constraints = get_SetDataStorage_constraints( );

        constraints.zero( numConstraints );

        unsigned int offset = 0;

        setCurrentResidualIndexMeaningful( true );
        for ( auto v = getResidualClasses( )->begin( ); v != getResidualClasses( )->end( ); v++ ){
            setCurrentResidualIndexMeaningful( v - getResidualClasses( )->begin( ) );

            if ( ( *v )->getNumConstraints( ) > 0 ){

                std::copy( ( *v )->getConstraints( )->begin( ), ( *v )->getConstraints( )->end( ),  constraints.value->begin( ) + offset );

                offset += ( *v )->getConstraints( )->size( );

            }

        }
        setCurrentResidualIndexMeaningful( false );

    }

    void hydraBase::setConstraintJacobians( ){
        /*!
         * Set the constraint Jacobians values
         */

        const unsigned int numUnknowns    = getNumUnknowns( );

        const unsigned int numConstraints = getNumConstraints( );

        auto constraintJacobians = get_SetDataStorage_constraintJacobians( );

        constraintJacobians.zero( numConstraints * numUnknowns );

        unsigned int offset = 0;

        setCurrentResidualIndexMeaningful( true );
        for ( auto v = getResidualClasses( )->begin( ); v != getResidualClasses( )->end( ); v++ ){
            setCurrentResidualIndex( v - getResidualClasses( )->begin( ) );

            if ( ( *v )->getNumConstraints( ) > 0 ){

                std::copy( ( *v )->getConstraintJacobians( )->begin( ), ( *v )->getConstraintJacobians( )->end( ),  constraintJacobians.value->begin( ) + offset );

                offset += ( *v )->getConstraintJacobians( )->size( );

            }

        }
        setCurrentResidualIndexMeaningful( false );

    }

    const unsigned int hydraBase::getCurrentResidualOffset( ){
        /*!
         * Get the offset of the current residual
         */
        unsigned int offset = 0;
        for ( auto v = getResidualClasses( )->begin( ); v != getResidualClasses( )->begin( ) + getCurrentResidualIndex( ); ++v ){
            offset += ( *v )->getNumEquations( );
        }
        return offset;
    }

    const unsigned int hydraBase::getNumConstraints( ){
        /*!
         * Get the value of the number of constraint equations
         */

        unsigned int value = 0;

        for ( auto v = getResidualClasses( )->begin( ); v != getResidualClasses( )->end( ); v++ ){

            value += ( *v )->getNumConstraints( );

        }

        return value;

    }

    std::string hydraBase::getResidualParameterizationInfo( ){
        /*!
         * Get the parameterization information of the residual classes
         */

        std::string parameterization_info = "########################################\n# RESIDUAL PARAMETERIZATION INFORMATION#\n########################################\n\n";

        for ( auto v = std::begin( *getResidualClasses( ) ); v != std::end( *getResidualClasses( ) ); ++v ){

            parameterization_info += "RESIDUAL CLASS:";
            parameterization_info += " " + std::to_string( ( unsigned int )( v - std::begin( *getResidualClasses( ) ) ) ) + "\n\n";

            ( *v )->addParameterizationInfo( parameterization_info );

            parameterization_info += "\n";

        }

        return parameterization_info;

    }


}
