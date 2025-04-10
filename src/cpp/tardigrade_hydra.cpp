/**
  ******************************************************************************
  * \file tardigrade_hydra.cpp
  ******************************************************************************
  * A C++ library for defining frameworks to solve finite deformation material
  * models.
  ******************************************************************************
  */

#include<tardigrade_hydra.h>

#include<tardigrade_abaqus_tools.h>

namespace tardigradeHydra{

    //Define hydra global constants in a place that Doxygen can pick up for documentation
    /** \brief Define the expected number of tensor spatial dimensions for the Abaqus interface. */
    const int spatialDimensions = 3;

    /** \brief Define required number of Abaqus material constants for the Abaqus interface. */
    const int nStateVariables = 2;

    /** \brief Define required number of Abaqus material constants for the Abaqus interface. */
    const int nMaterialParameters = 2;

    void residualBase::addIterationData( dataBase *data ){
        /*!
         * Add data to the vector of values which will be cleared after each iteration
         * 
         * \param *data: The dataBase object to be cleared
         */

        hydra->addIterationData( data );

    }

    void residualBase::addNLStepData( dataBase *data ){
        /*!
         * Add data to the vector of values which will be cleared after each nonlinear step
         * 
         * \param *data: The dataBase object to be cleared
         */

        hydra->addNLStepData( data );

    }

    void residualBase::setupRelaxedStep( const unsigned int &relaxedStep ){
        /*!
         * When performing a relaxed iteration this function is called prior to the solution of the non-linear
         * problem. Users can use this function to dynamically adjust parameters or perform other tuning tasks.
         *
         * \param &relaxedStep: The current relaxed step.
         */
    }

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
                                                           _tolr( tolr ), _tola( tola ),
                                                           _maxIterations( maxIterations ), _maxLSIterations( maxLSIterations ),
                                                           _lsAlpha( lsAlpha ),
                                                           _use_preconditioner( use_preconditioner ), _preconditioner_type( preconditioner_type ){
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

        // Initialize the preconditioner if required
        if ( _use_preconditioner ){

            initializePreconditioner( );

        }

    }

    void hydraBase::setStress( const floatVector &stress ){
        /*!
         * Set the value of the stress
         * 
         * \param &stress: The stress in row-major form
         */

        setIterationData( stress, _stress );

    }

    hydraBase::setDataStorageIteration<secondOrderTensor> hydraBase::get_setDataStorage_stress( ){
        /*!
         * Get a setDataStorage object for the stress
         */

        return hydraBase::setDataStorageIteration<secondOrderTensor>( &_stress, this );

    }

    void hydraBase::extractStress( ){
        /*!
         * Extract the stresses out of the unknown vector
         */

        const floatVector *unknownVector = getUnknownVector( );

        auto stress = get_setDataStorage_stress( );

        stress.zero( *getConfigurationUnknownCount( ) );

        std::copy( std::begin( *unknownVector ), std::begin( *unknownVector ) + *getConfigurationUnknownCount( ), std::begin( *stress.value ) );

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

        const unsigned int dim = getDimension( );
        const unsigned int sot_dim = getSOTDimension( );

        const unsigned int num_configs = *getNumConfigurations( );

        // Set the configurations
        configurations = floatVector( num_configs * sot_dim, 0 );

        inverseConfigurations = floatVector( num_configs * sot_dim, 0 );

        Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> > mat( NULL, 3, 3 );
#ifdef TARDIGRADE_HYDRA_USE_LLXSMM
        kernel_type kernel(LIBXSMM_GEMM_FLAG_NONE, dim, dim, dim, 1, 0 );

        // Initialize the first configuration with the total deformation gradient
        secondOrderTensor temp( sot_dim, 0 );
#else
        Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> > mat2( NULL, 3, 3 );
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
        auto configurations = get_setDataStorage_configurations( );

        auto inverseConfigurations = get_setDataStorage_inverseConfigurations( );

        computeConfigurations( unknownVector, *getStressSize( ), *getDeformationGradient( ), *configurations.value, *inverseConfigurations.value );

        // Extract the remaining state variables required for the non-linear solve
        auto nonLinearSolveStateVariables = get_setDataStorage_nonLinearSolveStateVariables( );

        const unsigned int *nNLISV = getNumNonLinearSolveStateVariables( );

        nonLinearSolveStateVariables.zero( *nNLISV );

        std::copy( std::begin( *unknownVector ) + ( *getNumConfigurations( ) ) * ( *getConfigurationUnknownCount( ) ),
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

        const unsigned int* nConfig = getNumConfigurations( );

        const unsigned int* nNLISV  = getNumNonLinearSolveStateVariables( );

        // Extract the previous configurations
        if ( getPreviousStateVariables( )->size( ) < ( ( ( *nConfig ) - 1 ) * ( *getConfigurationUnknownCount( ) ) + ( *nNLISV ) ) ){

            std::string message = "The number of state variables is less than required for the configurations and ";
            message            += "non-linear state variables\n";
            message            += "  # previousStateVariables                               : " + std::to_string( getPreviousStateVariables( )->size( ) ) + "\n";
            message            += "  # ( configurations - 1 ) * configuration_unknown_count : " + std::to_string( ( ( *nConfig ) - 1 ) * ( *getConfigurationUnknownCount( ) ) ) + "\n";
            message            += "  # non-linear solve ISVs                                : " + std::to_string( ( *nNLISV ) ) + "\n";
            message            += "  # minimum required ISVs                                : " + std::to_string( ( ( *nConfig ) - 1 ) * ( *getConfigurationUnknownCount( ) ) + ( *nNLISV ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

        }

        auto configurations = get_setDataStorage_configurations( );

        auto previousConfigurations = get_setDataStorage_previousConfigurations( );

        auto inverseConfigurations = get_setDataStorage_inverseConfigurations( );

        auto previousInverseConfigurations = get_setDataStorage_previousInverseConfigurations( );

        // Compute the configurations
        computeConfigurations( getPreviousStateVariables( ), 0, *getDeformationGradient( ), *configurations.value, *inverseConfigurations.value, true );

        computeConfigurations( getPreviousStateVariables( ), 0, *getPreviousDeformationGradient( ), *previousConfigurations.value, *previousInverseConfigurations.value, true );

        // Extract the remaining state variables required for the non-linear solve
        auto nonLinearSolveStateVariables         = get_setDataStorage_nonLinearSolveStateVariables( );
        auto previousNonLinearSolveStateVariables = get_setDataStorage_previousNonLinearSolveStateVariables( );

        previousNonLinearSolveStateVariables.zero( *nNLISV );

        std::copy( std::begin( *getPreviousStateVariables( ) ) + ( ( *nConfig ) - 1 ) * ( *getConfigurationUnknownCount( ) ),
                   std::begin( *getPreviousStateVariables( ) ) + ( ( *nConfig ) - 1 ) * ( *getConfigurationUnknownCount( ) ) + *nNLISV,
                   std::begin( *previousNonLinearSolveStateVariables.value ) );

        *nonLinearSolveStateVariables.value = *previousNonLinearSolveStateVariables.value;

        // Extract the additional state variables
        auto additionalStateVariables         = get_setDataStorage_additionalStateVariables( );
        auto previousAdditionalStateVariables = get_setDataStorage_previousAdditionalStateVariables( );

        unsigned int nAISV = ( unsigned int )( std::end( *getPreviousStateVariables( ) ) - ( std::begin( *getPreviousStateVariables( ) ) + ( ( *nConfig ) - 1 ) * ( *getConfigurationUnknownCount( ) ) + *nNLISV ) );

        previousAdditionalStateVariables.zero( nAISV );

        std::copy( std::begin( *getPreviousStateVariables( ) ) + ( ( *nConfig ) - 1 ) * ( *getConfigurationUnknownCount( ) ) + *nNLISV,
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
        Eigen::Map<Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor>> Fsc_mat( Fsc.data( ), 3, 3 );
        Eigen::Map<const Eigen::Matrix<floatType, 3, 3, Eigen::RowMajor>> mat( NULL, 3, 3 );
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

        Eigen::Map< Eigen::Matrix< floatType, 3, 3 > > map( NULL, dim, dim );

        for ( unsigned int index = lowerIndex; index < upperIndex; index++ ){

            secondOrderTensor Fm, FpT;

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

        return getSubConfiguration( index + 1, *getNumConfigurations( ) );

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

        return getPreviousSubConfiguration( index + 1, *getNumConfigurations( ) );

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

        return getSubConfigurationJacobian( index + 1, *getNumConfigurations( ) );

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

        return getPreviousSubConfigurationJacobian( index + 1, *getNumConfigurations( ) );

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
        const unsigned int num_configs = *getNumConfigurations( );

        dC1dC  = secondOrderTensor( sot_dim * sot_dim, 0 );

        dC1dCn = floatVector( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

        secondOrderTensor fullConfiguration = getSubConfiguration( configurations, 0, num_configs );

        secondOrderTensor invCsc = getSubConfiguration( configurations, 1, num_configs );
        Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor > > mat( invCsc.data( ), dim, dim );
        mat = mat.inverse( ).eval( );

        fourthOrderTensor dInvCscdCsc = tardigradeVectorTools::computeFlatDInvADA( invCsc, dim, dim );
        Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dInvCscdCsc( dInvCscdCsc.data( ), sot_dim, sot_dim );

        floatVector dCscdCs = getSubConfigurationJacobian( configurations, 1, num_configs );
        Eigen::Map< const Eigen::Matrix< floatType, sot_dim, -1, Eigen::RowMajor > > map_dCscdCs( dCscdCs.data( ), sot_dim, num_configs * sot_dim );

        floatVector dInvCscdCs( sot_dim * num_configs * sot_dim, 0 );
        Eigen::Map< Eigen::Matrix< floatType, sot_dim, -1, Eigen::RowMajor > > map_dInvCscdCs( dInvCscdCs.data( ), sot_dim, num_configs * sot_dim ); 

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

        auto dF1dF = get_setDataStorage_dF1dF( );

        auto dF1dFn = get_setDataStorage_dF1dFn( );

        calculateFirstConfigurationJacobians( *get_configurations( ), *dF1dF.value, *dF1dFn.value );

    }

    void hydraBase::setPreviousFirstConfigurationJacobians( ){
        /*!
         * Set the Jacobians of the previous first configuration w.r.t. the total configuration and the remaining sub-configurations
         */

        auto dF1dF = get_setDataStorage_previousdF1dF( );

        auto dF1dFn = get_setDataStorage_previousdF1dFn( );

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
         * The user should define a vector of residualBase objects and use the
         * setResidualClasses( std::vector< residualBase > & ) function here.
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

    void hydraBase::setResidualClasses( std::vector< residualBase* > &residualClasses ){
        /*!
         * Set the residual classes
         * 
         * \param &residualClasses: A vector of residual classes which will be used to
         *     populate the residual and jacobian matrices for the non-linear solve
         */

        unsigned int numEquations = 0;

        _residualClasses.second = std::vector< residualBase* >( residualClasses.size( ) );

        for ( auto c = residualClasses.begin( ); c != residualClasses.end( ); c++ ){

            numEquations += *( *c )->getNumEquations( );

            _residualClasses.second[ c - residualClasses.begin( ) ] = *c;

        }

        if ( numEquations != ( *getNumConfigurations( ) * ( *getConfigurationUnknownCount( ) ) + *getNumNonLinearSolveStateVariables( ) ) ){

            std::string message = "The number of equations for the non-linear solve is not equal to the number of equations defined\n";
            message            += "  expected number of equations: " + std::to_string( ( *getNumConfigurations( ) ) * ( *getConfigurationUnknownCount( ) ) + *getNumNonLinearSolveStateVariables( ) ) + "\n";
            message            += "  number of defined equations:  " + std::to_string( numEquations ) + "\n";

            TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

        }

        _residualClasses.first = true;

    }

    std::vector< residualBase* >* hydraBase::getResidualClasses( ){
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

        const unsigned int configurationUnknownCount = *getConfigurationUnknownCount( );

        const unsigned int residualSize = ( *getNumConfigurations( ) ) * configurationUnknownCount + *getNumNonLinearSolveStateVariables( );

        const unsigned int numAdditionalDOF = getAdditionalDOF( )->size( );

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

            TARDIGRADE_ERROR_TOOLS_CHECK( localResidual->size( ) == *( *residual_ptr )->getNumEquations( ),
                  "The residual for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n"
                + "  expected: " + std::to_string( *( *residual_ptr )->getNumEquations( ) ) + "\n"
                + "  actual:   " + std::to_string( localResidual->size( ) ) + "\n"
            )

            // Store the values in the global quantities

            // Copy over the values of the local vector to the global structures
            std::copy( localResidual->begin( ), localResidual->end( ), _residual.second.begin( ) + offset );

            offset += *( *residual_ptr )->getNumEquations( );

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

        const unsigned int configurationUnknownCount = *getConfigurationUnknownCount( );

        const unsigned int residualSize = ( *getNumConfigurations( ) ) * configurationUnknownCount + *getNumNonLinearSolveStateVariables( );

        const unsigned int numAdditionalDOF = getAdditionalDOF( )->size( );

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

            TARDIGRADE_ERROR_TOOLS_CHECK( localJacobian->size( ) == *( *residual_ptr )->getNumEquations( ) * residualSize,
                  "The jacobian for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n"
                + "  expected: " + std::to_string( *( *residual_ptr )->getNumEquations( ) * residualSize ) + "\n"
                + "  actual:   " + std::to_string( localJacobian->size( ) ) + "\n"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK( localdRdF->size( ) == *( *residual_ptr )->getNumEquations( ) * configurationUnknownCount,
                  "dRdF for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n"
                + "  expected: " + std::to_string( *( *residual_ptr )->getNumEquations( ) * configurationUnknownCount ) + "\n"
                + "  actual:   " + std::to_string( localdRdF->size( ) ) + "\n"
            )

            TARDIGRADE_ERROR_TOOLS_CHECK( localdRdT->size( ) == *( *residual_ptr )->getNumEquations( ),
                  "dRdT for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n"
                + "  expected: " + std::to_string( *( *residual_ptr )->getNumEquations( ) ) + "\n"
                + "  actual:   " + std::to_string( localdRdT->size( ) ) + "\n"
            )

            if ( localdRdAdditionalDOF->size( ) != 0 ){

                TARDIGRADE_ERROR_TOOLS_CHECK( localdRdAdditionalDOF->size( ) == ( ( *( *residual_ptr )->getNumEquations( ) ) * numAdditionalDOF ),
                                              "dRdAdditionalDOF for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n"
                                            + "  expected: " + std::to_string( ( *( *residual_ptr )->getNumEquations( ) ) * numAdditionalDOF ) + "\n"
                                            + "  actual  : " + std::to_string( localdRdAdditionalDOF->size( ) ) + "\n"
                )

                std::copy ( localdRdAdditionalDOF->begin( ), localdRdAdditionalDOF->end( ), _dRdAdditionalDOF.second.begin( ) + numAdditionalDOF * offset );

            }

            if ( localAdditionalDerivatives->size( ) != 0 ){

                if ( ( *localAdditionalDerivatives ).size( ) != *( *residual_ptr )->getNumEquations( ) * numAdditionalDerivatives ){
    
                    if ( numAdditionalDerivatives == 0 ){
    
                        numAdditionalDerivatives = ( *localAdditionalDerivatives ).size( ) / ( *( *residual_ptr )->getNumEquations( ) );
    
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

            offset += *( *residual_ptr )->getNumEquations( );

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

    void hydraBase::formPreconditioner( ){
        /*!
         * Form the preconditioner matrix
         */

        if ( _preconditioner_type == 0 ){

            formMaxRowPreconditioner( );

        }
        else{

            throw std::runtime_error( "Preconditioner type not recognized" );

        }

        _preconditioner.first = true;

        addIterationData( &_preconditioner );

    }

    void hydraBase::formMaxRowPreconditioner( ){
        /*!
         * Form a left preconditioner comprised of the inverse of the maximum value of each row
         */

        const unsigned int problem_size = getNumUnknowns( );

        _preconditioner.second = floatVector( problem_size, 0 );

        // Find the absolute maximum value in each row
        for ( unsigned int i = 0; i < problem_size; i++ ){

            _preconditioner.second[ i ] = 1 / std::max( std::fabs( *std::max_element( getFlatNonlinearLHS( )->begin( ) + problem_size * i,
                                                                                      getFlatNonlinearLHS( )->begin( ) + problem_size * ( i + 1 ),
                                                                                      [ ]( const floatType &a, const floatType &b ){ return std::fabs( a ) < std::fabs( b ); } ) ), 1e-15 );

        }

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

    const floatVector* hydraBase::getNonlinearRHS( ){
        /*!
         * Get the RHS vector for the non-linear problem
         */

        if ( _use_LM_step ){

            if ( !_nonlinearRHS.first ){

                const unsigned int xsize = getNumUnknowns( );

                _nonlinearRHS.first = true;

                _nonlinearRHS.second = floatVector( xsize, 0 );

                const floatVector *residual = getResidual( );

                const floatVector *jacobian = getFlatJacobian( );

                Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > Jmap( jacobian->data( ), xsize, xsize );

                Eigen::Map< const Eigen::Matrix< floatType, -1,  1 > > Rmap( residual->data( ), xsize );

                Eigen::Map< Eigen::Matrix< floatType, -1,  1 > > RHSmap( _nonlinearRHS.second.data( ), xsize );

                RHSmap = ( Jmap.transpose( ) * Rmap ).eval( );

                addIterationData( &_nonlinearRHS );

            }

            return &_nonlinearRHS.second;

        }

        return getResidual( );

    }

    const floatVector* hydraBase::getFlatNonlinearLHS( ){
        /*!
         * Get the flat LHS matrix for the non-linear problem
         */

        if ( _use_LM_step ){

            if ( !_flatNonlinearLHS.first ){

                const unsigned int xsize = getNumUnknowns( );

                _flatNonlinearLHS.first = true;

                _flatNonlinearLHS.second = floatVector( xsize * xsize, 0 );

                const floatVector * jacobian = getFlatJacobian( );

                Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > Jmap( jacobian->data( ), xsize, xsize );

                Eigen::Map< Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > LHSmap( _flatNonlinearLHS.second.data( ), xsize, xsize );

                LHSmap = ( Jmap.transpose( ) * Jmap ).eval( );

                for ( unsigned int i = 0; i < xsize; i++ ){
                    _flatNonlinearLHS.second[ xsize * i + i ] += *getMuk( );
                }

                addIterationData( &_flatNonlinearLHS );

            }

            return &_flatNonlinearLHS.second;

        }

        return getFlatJacobian( );

    }

    const floatVector* hydraBase::getFlatPreconditioner( ){
        /*!
         * Get the flattened row-major preconditioner for the non-linear problem
         */

        if ( !_preconditioner.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( formPreconditioner( ) );

        }
        
        return &_preconditioner.second;

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

    /// Say hello
    /// @param message The message to print
    void sayHello( std::string message ) {
        TARDIGRADE_ERROR_TOOLS_CHECK( ( message.compare( "George" ) != 0 ), "ERROR: George is a wolf in sheep's clothing!");
        std::cout << "Hello " << message << std::endl;
        return;
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

    void hydraBase::setTolerance( ){
        /*!
         * Set the tolerance
         * 
         * \f$ tol = tolr * ( |R_0| + |X| ) + tola \f$
         */

        auto tolerance = get_setDataStorage_tolerance( );

        *tolerance.value = tardigradeVectorTools::abs( *getResidual( ) ) + tardigradeVectorTools::abs( *getUnknownVector( ) );

        *tolerance.value = *getRelativeTolerance( ) * ( *tolerance.value ) + *getAbsoluteTolerance( );

    }

    void hydraBase::setTolerance( const floatVector &tolerance ){
        /*!
         * Set the tolerance
         *
         * \param tolerance: The tolerance vector for each value of the residual
         */

        setConstantData( tolerance, _tolerance );

    }

    hydraBase::setDataStorageConstant<floatVector> hydraBase::get_setDataStorage_tolerance( ){
        /*!
         * Return a setDataStorageConstant setter for the tolerance
         */

        return setDataStorageConstant<floatVector>( &_tolerance );

    }

    const floatVector* hydraBase::getTolerance( ){
        /*!
         * Get the tolerance
         */

        if ( !_tolerance.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setTolerance( ) );

        }

        return &_tolerance.second;

    }

    bool hydraBase::checkConvergence( ){
        /*!
         * Check the convergence
         */

        const floatVector *tolerance = getTolerance( );

        const floatVector *residual = getResidual( );

        if ( tolerance->size( ) != residual->size( ) ){

            std::string message = "The residual and tolerance vectors don't have the same size\n";
            message            += "  tolerance: " + std::to_string( tolerance->size( ) ) + "\n";
            message            += "  residual:  " + std::to_string( residual->size( ) ) + "\n";

            TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

        }

        for ( unsigned int i = 0; i < tolerance->size( ); i++ ){

            if ( std::fabs( ( *residual )[ i ] ) > ( *tolerance )[ i ] ){

                return false;

            }

        }

        return true;

    }

    const floatType* hydraBase::getLSResidualNorm( ){
        /*!
         * Get the residual norm for the line-search convergence criterion
         */

        if ( !_lsResidualNorm.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( resetLSIteration( ) );

        }

        return &_lsResidualNorm.second;

    }

    bool hydraBase::checkLSConvergence( ){
        /*!
         * Check the line-search convergence
         */

        if ( tardigradeVectorTools::l2norm( *getResidual( ) ) < ( 1 - *getLSAlpha( ) ) * ( *getLSResidualNorm( ) ) ){

            return true;

        }

        return false;

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
            if ( *( *residual_ptr )->getUseProjection( ) ){
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

    void hydraBase::performPreconditionedSolve( floatVector &deltaX_tr ){
        /*!
         * Perform a pre-conditioned solve
         *
         * \param &deltaX_tr: The trial chcange in the unknown vector
         */

        tardigradeVectorTools::solverType< floatType > linearSolver;

        Eigen::Map< Eigen::Vector< floatType, -1 > > dx_map( deltaX_tr.data( ), getNumUnknowns( ) );

        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > J_map( getFlatNonlinearLHS( )->data( ), getNumUnknowns( ), getNumUnknowns( ) );

        Eigen::Map< const Eigen::Vector< floatType, -1 > > R_map( getNonlinearRHS( )->data( ), getNumUnknowns( ) );

        if( *getPreconditionerIsDiagonal( ) ){

            Eigen::Map< const Eigen::Vector< floatType, -1 > > p_map( getFlatPreconditioner( )->data( ), getNumUnknowns( ) );

            linearSolver = tardigradeVectorTools::solverType< floatType >( p_map.asDiagonal( ) * J_map );

            dx_map = -linearSolver.solve( p_map.asDiagonal( ) * R_map );

        }
        else{

            Eigen::Map< const Eigen::Matrix< floatType, -1, -1 > > p_map( getFlatPreconditioner( )->data( ), getNumUnknowns( ), getNumUnknowns( ) );

            linearSolver = tardigradeVectorTools::solverType< floatType >( p_map * J_map );

            dx_map = -linearSolver.solve( p_map * R_map );

        }

        unsigned int rank = linearSolver.rank( );

        if ( *getRankDeficientError( ) && ( rank != getResidual( )->size( ) ) ){

            TARDIGRADE_ERROR_TOOLS_CATCH( throw convergence_error( "The Jacobian is not full rank" ) );

        }

    }

    void hydraBase::solveNewtonUpdate( floatVector &deltaX_tr ){
        /*!
         * Solve the Newton update returning the trial value of the unknown vector
         *
         * \param &deltaX_tr: The trial change in the unknown vector
         */

        if ( *getUsePreconditioner( ) ){

            performPreconditionedSolve( deltaX_tr );

        }
        else{

            Eigen::Map< Eigen::Vector< floatType, -1 > > dx_map( deltaX_tr.data( ), getNumUnknowns( ) );

            Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > J_map( getFlatNonlinearLHS( )->data( ), getNumUnknowns( ), getNumUnknowns( ) );

            Eigen::Map< const Eigen::Vector< floatType, -1 > > R_map( getNonlinearRHS( )->data( ), getNumUnknowns( ) );

            tardigradeVectorTools::solverType< floatType > linearSolver( J_map );
            dx_map = -linearSolver.solve( R_map );

            unsigned int rank = linearSolver.rank( );

            if ( *getRankDeficientError( ) && ( rank != getResidual( )->size( ) ) ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw convergence_error( "The Jacobian is not full rank" ) );

            }

        }

    }

    void hydraBase::performArmijoTypeLineSearch( const floatVector &X0, const floatVector &deltaX ){
        /*!
         * Perform an Armijo-type line search
         *
         * \param &X0: The base value of the unknown vector
         * \param &deltaX: The proposed change in X
         */

        while ( !checkLSConvergence( ) && checkLSIteration( ) ){

            if ( ( *getFailureVerbosityLevel( ) ) > 0 ){
                addToFailureOutput( "    lambda, |R|: " );
                addToFailureOutput( *getLambda( ) );
                addToFailureOutput( ", " );
                addToFailureOutput( tardigradeVectorTools::l2norm( *getResidual( ) ) );
                addToFailureOutput( "\n" );
            }

            updateLambda( );

            incrementLSIteration( );

            updateUnknownVector( X0 + *getLambda( ) * deltaX );

        }

        if ( !checkLSConvergence( ) ){

            throw convergence_error( "Failure in line search:\n  scale factor: " + std::to_string( *getScaleFactor( ) ) + "\n" );

        }

        incrementNumLS( );

        resetLSIteration( );

    }

    const floatType *hydraBase::get_baseResidualNorm( ){
        /*!
         * Get the base value for the residual norm.
         */

        if ( !_baseResidualNorm.first ){

            throw std::runtime_error( "The base residual norm must be set with set_baseResidualNorm before it can be called" );

        }

        return &_baseResidualNorm.second;

    }

    const floatVector *hydraBase::get_basedResidualNormdX( ){
        /*!
         * Get the base value for the derivative of the residual norm w.r.t. the unknown vector
         */

        if ( !_basedResidualNormdX.first ){

            throw std::runtime_error( "The base residual norm must be set with set_dbaseResidualNormdX before it can be called" );

        }

        return &_basedResidualNormdX.second;

    }

    bool hydraBase::checkGradientConvergence( const floatVector &X0 ){
        /*!
         * Check the convergence of a gradient step
         *
         * \param &X0: The initial value of the unknown vector
         */

        const unsigned int xsize = getNumUnknowns( );

        floatVector dx = ( *getUnknownVector( ) ) - X0;

        floatType RHS = *get_baseResidualNorm( );

        for ( unsigned int i = 0; i < xsize; i++ ){

            RHS += ( *getGradientSigma( ) ) * ( *get_basedResidualNormdX( ) )[ i ] * dx[ i ];

        }

        return ( *get_residualNorm( ) ) < RHS;

    }

    void hydraBase::performGradientStep( const floatVector &X0 ){
        /*!
         * Perform a gradient descent step
         *
         * \param &X0: The base value of the unknown vector
         */

        const floatVector *dResidualNormdX = get_basedResidualNormdX( );

        unsigned int l                     = 0;

        const unsigned int maxiter         = *getMaxGradientIterations( );

        while( checkGradientIteration( ) ){

            floatType t = std::pow( *getGradientBeta( ), l );

            updateUnknownVector( X0 - t * ( *dResidualNormdX ) );

            if ( checkGradientConvergence( X0 ) ){

                break;

            }

            l++;

            incrementGradientIteration( );

        }

        if ( l >= maxiter ){

            throw convergence_error( "Failure in gradient step:\n  scale_factor: " + std::to_string( *getScaleFactor( ) ) );

        }

        incrementNumGrad( );

        resetGradientIteration( );

    }

    bool hydraBase::checkDescentDirection( const floatVector &dx ){
        /*!
         * Check if the search direction is a descent direction of the Jacobian
         * 
         * \param &dx: The proposed change in x
         */

        const unsigned int xsize = getNumUnknowns( );

        const floatType RHS = -( *getGradientRho( ) ) * std::pow( tardigradeVectorTools::l2norm( dx ), *getGradientP( ) );

        floatType LHS = 0;

        const floatVector *dResidualNormdX = get_basedResidualNormdX( );

        for ( unsigned int i = 0; i < xsize; i++ ){

            LHS += ( *dResidualNormdX )[ i ] * dx[ i ];

        }

        return LHS <= RHS;

    }

    void hydraBase::setBaseQuantities( ){
        /*!
         * Set the base quantities required for gradient steps
         */

        set_baseResidualNorm( *get_residualNorm( ) );

        set_basedResidualNormdX( *get_dResidualNormdX( ) );

        if ( _mu_k < 0 ){

            setMuk( 0.5 * ( *getLMMu( ) ) * ( *get_baseResidualNorm( ) ) );

        }
        else{

            setMuk( std::fmin( _mu_k, ( *get_baseResidualNorm( ) ) ) );

        }

    }

    void hydraBase::callResidualSuccessfulNLStep( ){
        /*!
         * Signal to the residuals that a successful nonlinear step has been performed
         */

        setAllowModifyGlobalResidual( true );
        setCurrentResidualIndexMeaningful( true );
        for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); residual_ptr++ ){
            setCurrentResidualIndex( residual_ptr - getResidualClasses( )->begin( ) );

            ( *residual_ptr )->successfulNLStep( );

        }
        setCurrentResidualIndexMeaningful( false );
        setAllowModifyGlobalResidual( false );

    }

    void hydraBase::callResidualPreNLSolve( ){
        /*!
         * Signal to the residuals that we are about to start a nonlinear solve
         */

        setCurrentResidualIndexMeaningful( true );
        for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); residual_ptr++ ){
            setCurrentResidualIndex( residual_ptr - getResidualClasses( )->begin( ) );

            ( *residual_ptr )->preNLSolve( );

        }
        setCurrentResidualIndexMeaningful( false );

    }

    void hydraBase::solveNonLinearProblem( ){
        /*!
         * Solve the non-linear problem
         */

        // Form the initial unknown vector
        if ( getInitializeUnknownVector( ) ){
            TARDIGRADE_ERROR_TOOLS_CATCH( initializeUnknownVector( ) );
        }

        _initialX = *getUnknownVector( );

        floatVector deltaX( getNumUnknowns( ), 0 );

        resetLSIteration( );

        resetGradientIteration( );

        if ( ( *getFailureVerbosityLevel( ) ) > 0 ){
            addToFailureOutput( "Initial Unknown:\n" );
            addToFailureOutput( *getUnknownVector( ) );
        }

        callResidualPreNLSolve( );

        while( !checkConvergence( ) && checkIteration( ) ){

            if ( ( *getFailureVerbosityLevel( ) ) > 0 ){
                addToFailureOutput( "\n\n  iteration: " );
                addToFailureOutput( _iteration );
            }
            floatVector X0 = *getUnknownVector( );

            if ( ( *getFailureVerbosityLevel( ) ) > 0 ){
                addToFailureOutput( "  X0:\n" );
                addToFailureOutput( "  " );
                addToFailureOutput( *getUnknownVector( ) );
            }
            setBaseQuantities( );

            if ( *getUseSQPSolver( ) ){

                std::fill( deltaX.begin( ), deltaX.end( ), 0 );

                solveConstrainedQP( deltaX );

            }
            else{

                solveNewtonUpdate( deltaX );

            }
            if ( ( *getFailureVerbosityLevel( ) ) > 0 ){
                addToFailureOutput( "  deltaX:\n" );
                addToFailureOutput( "  " );
                addToFailureOutput( deltaX );
            }

            updateUnknownVector( X0 + *getLambda( ) * deltaX );

            // Refine the estimate if the new point has a higher residual
            if ( !checkLSConvergence( ) ){

                if ( checkDescentDirection( deltaX ) || !( *getUseGradientDescent( ) ) ){

                    // Perform an Armijo type line search when the search direction is aligned with the gradient
                    performArmijoTypeLineSearch( X0, deltaX );

                }
                else{

                    // Perform gradient descent if the search direction is not aligned with the gradient
                    performGradientStep( X0 );

                }

            }
            else{

                incrementNumNewton( );

            }

            // Increment the iteration count
            incrementIteration( );

            // Reset the nonlinear step data
            resetNLStepData( );

            // Call residual end of a successful nonlinear step functions
            callResidualSuccessfulNLStep( );

            if ( ( *getFailureVerbosityLevel( ) ) > 0 ){
                addToFailureOutput( "  final residual: " );
                addToFailureOutput( tardigradeVectorTools::l2norm( *getResidual( ) ) );
                addToFailureOutput( "\n" );
            }

        }

        if ( !checkConvergence( ) ){

            throw convergence_error( "Failure to converge main loop:\n  scale_factor: " + std::to_string( *getScaleFactor( ) ) );

        }

    }

    void hydraBase::performRelaxedSolve( ){
        /*!
         * Solve the non-linear problem by relaxing difficult sub-problems
         * to achieve a series of solutions.
         */

        unsigned int relaxedIteration = 0;

        // Initialize the residuals
        setCurrentResidualIndexMeaningful( true );
        for ( auto residual = getResidualClasses( )->begin( ); residual != getResidualClasses( )->end( ); residual++ ){
            setCurrentResidualIndexMeaningful( residual - getResidualClasses( )->begin( ) );

            // Prepare the residuals to take a relaxed step
            ( *residual )->setupRelaxedStep( relaxedIteration );

        }
        setCurrentResidualIndexMeaningful( false );

        while ( relaxedIteration < *getMaxRelaxedIterations( ) ){

            if ( ( *getFailureVerbosityLevel( ) ) > 0 ){
                addToFailureOutput( "\n\n###  relaxed iteration: " );
                addToFailureOutput( relaxedIteration );
                addToFailureOutput( "\n\n" );
            }
            // Solve the non-linear problem
            TARDIGRADE_ERROR_TOOLS_CATCH( solveNonLinearProblem( ) );

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

            // Use the current unknown vector as the initial estimate
            if ( *getInitializeUnknownVector( ) ){

                setInitializeUnknownVector( false );

            }

            relaxedIteration++;

            // Initialize the residuals
            setCurrentResidualIndexMeaningful( true );
            for ( auto residual = getResidualClasses( )->begin( ); residual != getResidualClasses( )->end( ); residual++ ){
                setCurrentResidualIndex( residual - getResidualClasses( )->begin( ) );

                // Prepare the residuals to take a relaxed step
                ( *residual )->setupRelaxedStep( relaxedIteration );

            }
            setCurrentResidualIndexMeaningful( false );

            // Reset hydra
            updateUnknownVector( *getUnknownVector( ) ); //This allows for the relaxed to change the projection and adjust the decomposition
            resetIterations( );

        }

        if ( relaxedIteration >= *getMaxRelaxedIterations( ) ){

            throw convergence_error( "Failure in relaxed solve:\n  scale_factor: " + std::to_string( *getScaleFactor( ) ) );

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

        if ( num_good >= ( *getNumGoodControl( ) ) ){

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

            if ( ( *getFailureVerbosityLevel( ) ) > 0 ){
                addToFailureOutput( "\n\n" );
                addToFailureOutput( "#########################################\n" );
                addToFailureOutput( "###        ENTERING SUB-CYCLER        ###\n" );
                addToFailureOutput( "#########################################\n" );
                addToFailureOutput( "\n\n" );
            }

            floatType sp = 0.0;

            floatType ds = ( *getCutbackFactor( ) );

            unsigned int num_good = 0;

            // Set the unknown vector to the initial unknown. We're using setX because we call setScaleFactor right away which will update the unknown vector
            setX( _initialX );

            while ( sp < 1.0 ){

                try{

                    if ( ( *getFailureVerbosityLevel( ) ) > 0 ){
                        addToFailureOutput( "\n\n" );
                        addToFailureOutput( "######### PSEUDO-TIME INCREMENT #########\n" );
                        addToFailureOutput( "\n\n    sp, ds: " + std::to_string( sp ) + ", " + std::to_string( ds ) );
                        addToFailureOutput( "\n" );
                    }

                    setScaleFactor( sp + ds ); // Update the scaling factor

                    resetIterations( ); // Reset the non-linear iteration count

                    evaluateInternal( ); // Try to solve the non-linear problem

                    sp += ds; // Update the pseudo-time

                    num_good++; // Update the number of good iterations

                    // Grow the step if possible
                    if ( allowStepGrowth( num_good ) ){

                        ds *= ( *getGrowthFactor( ) );

                    }

                    // Make sure s will be less than or equal to 1
                    if ( sp + ds > 1.0 ){

                        ds = 1.0 - sp;

                    }

                }
                catch( std::exception &e ){

                    // Reduce the time-step and try again
                    num_good = 0;

                    ds *= ( *getCutbackFactor( ) );

                    setX( _initialX ); // Reset X to the last good point

                    if ( ds < *getMinDS( ) ){

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
        resetNumNewton( );

        resetNumLS( );

        resetNumGrad( );

        setRankDeficientError( false );

        try{

            solveNonLinearProblem( );

        }
        catch( const convergence_error &e ){

            if ( *getUseRelaxedSolve( ) ){

                try{

                    resetIterations( );
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
                setRankDeficientError( false );
    
                setUseLevenbergMarquardt( true );
    
                // Turn on projection
                setCurrentResidualIndexMeaningful( true );
                for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); residual_ptr++ ){
                    setCurrentResidualIndex( residual_ptr - getResidualClasses( )->begin( ) );
                    ( *residual_ptr )->setUseProjection( true );
                }
                setCurrentResidualIndexMeaningful( false );
    
                resetIterations( );
                updateUnknownVector( _initialX );
    
                try{
    
                    solveNonLinearProblem( );
    
                }
                catch( const convergence_error &e ){
    
                    throw;
    
                }
                catch( std::exception &e ){
    
                    setUseLevenbergMarquardt( false );
    
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
        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > Amat( getFlatJacobian( )->data( ), getResidual( )->size( ), getResidual( )->size( ) );

        // Form the maps for dXdF
        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > dRdFmat( getFlatdRdF( )->data( ), getResidual( )->size( ), *getConfigurationUnknownCount( ) );

        _flatdXdF.second = floatVector( getNumUnknowns( ) * ( *getConfigurationUnknownCount( ) ) );
        Eigen::Map< Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > dXdFmat( _flatdXdF.second.data( ), getNumUnknowns( ), ( *getConfigurationUnknownCount( ) ) );

        // Form the maps for dXdT
        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > dRdTmat( getdRdT( )->data( ), getResidual( )->size( ), 1 );

        _flatdXdT.second = floatVector( getNumUnknowns( ) );
        Eigen::Map< Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > dXdTmat( _flatdXdT.second.data( ), getNumUnknowns( ), 1 );

        // Solve
        tardigradeVectorTools::solverType< floatType > solver;

        if ( *getUsePreconditioner( ) ){

            if( *getPreconditionerIsDiagonal( ) ){

                Eigen::Map< const Eigen::Vector< floatType, -1 > > p_map( getFlatPreconditioner( )->data( ), getResidual( )->size( ) );

                solver = tardigradeVectorTools::solverType< floatType >( p_map.asDiagonal( ) * Amat );

                dXdFmat = -solver.solve( p_map.asDiagonal( ) * dRdFmat );

                dXdTmat = -solver.solve( p_map.asDiagonal( ) * dRdTmat );

            }
            else{

                Eigen::Map< const Eigen::Matrix< floatType, -1, -1 > > p_map( getFlatPreconditioner( )->data( ), getResidual( )->size( ), getResidual( )->size( ) );

                solver = tardigradeVectorTools::solverType< floatType >( p_map * Amat );

                dXdFmat = -solver.solve( p_map * dRdFmat );

                dXdTmat = -solver.solve( p_map * dRdTmat );

            }

        }
        else{

            solver = tardigradeVectorTools::solverType< floatType >( Amat );

            dXdFmat = -solver.solve( dRdFmat );

            dXdTmat = -solver.solve( dRdTmat );

        }

        unsigned int rank = solver.rank( );

        TARDIGRADE_ERROR_TOOLS_CATCH(

            if ( *getRankDeficientError( ) && ( rank != getResidual( )->size( ) ) ){

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
        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > Amat( getFlatJacobian( )->data( ), getResidual( )->size( ), getResidual( )->size( ) );

        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > dRdAdditionalDOF( getFlatdRdAdditionalDOF( )->data( ), getResidual( )->size( ), getAdditionalDOF( )->size( ) );

        // Form the map for dXdF
        _flatdXdAdditionalDOF.second = floatVector( getNumUnknowns( ) * getAdditionalDOF( )->size( ), 0 );
        Eigen::Map< Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > dXdAdditionalDOF( _flatdXdAdditionalDOF.second.data( ), getNumUnknowns( ), getAdditionalDOF( )->size( ) );

        // Solve
        tardigradeVectorTools::solverType< floatType > solver;

        if ( *getUsePreconditioner( ) ){

            if( *getPreconditionerIsDiagonal( ) ){

                Eigen::Map< const Eigen::Vector< floatType, -1 > > p_map( getFlatPreconditioner( )->data( ), getResidual( )->size( ) );

                solver = tardigradeVectorTools::solverType< floatType >( p_map.asDiagonal( ) * Amat );

                dXdAdditionalDOF = -solver.solve( p_map.asDiagonal( ) * dRdAdditionalDOF );

            }
            else{

                Eigen::Map< const Eigen::Matrix< floatType, -1, -1 > > p_map( getFlatPreconditioner( )->data( ), getResidual( )->size( ), getResidual( )->size( ) );

                solver = tardigradeVectorTools::solverType< floatType >( p_map * Amat );

                dXdAdditionalDOF = -solver.solve( p_map * dRdAdditionalDOF );

            }

        }
        else{

            solver = tardigradeVectorTools::solverType< floatType >( Amat );

            dXdAdditionalDOF = -solver.solve( dRdAdditionalDOF );

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

    void hydraBase::initializePreconditioner( ){
        /*!
         * Initialize the preconditioner
         */

        _preconditioner_is_diagonal = false;

        if ( _preconditioner_type == 0 ){

            _preconditioner_is_diagonal = true;

        }
        else{

            throw std::runtime_error( "Preconditioner type " + std::to_string( _preconditioner_type ) + " is not recognized." );

        }

    }

    void hydraBase::setResidualNorm( ){
        /*!
         * Set the norm of the residual vector
         */

        auto residualNorm = get_setDataStorage_residualNorm( );

        residualNorm.zero( );

        const unsigned int xsize = getNumUnknowns( );

        const floatVector *residual = getResidual( );

        for ( unsigned int i = 0; i < xsize; i++ ){

            *residualNorm.value += ( *residual )[ i ] * ( *residual )[ i ];

        }

    }

    void hydraBase::setdResidualNormdX( ){
        /*!
         * Set the derivative of the residual norm w.r.t. the unknown vector
         */

        const unsigned int xsize = getNumUnknowns( );

        auto dResidualNormdX = get_setDataStorage_dResidualNormdX( );

        dResidualNormdX.zero( xsize );

        const floatVector *residual = getResidual( );

        const floatVector *jacobian = getFlatJacobian( );

        for ( unsigned int i = 0; i < xsize; i++ ){
            for ( unsigned int j = 0; j < xsize; j++ ){
                ( *dResidualNormdX.value )[ j ] += 2 * ( *jacobian )[ xsize * i + j ] * ( *residual )[ i ];
            }
        }

    }

    void hydraBase::assembleKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints ){
        /*!
         * Assemble the Karush-Kuhn-Tucker matrix for an inequality constrained Newton-Raphson solve
         * 
         * \param &KKTMatrix: The Karush-Kuhn-Tucker matrix
         * \param &active_constraints: The vector of currently active constraints.
         */

        const unsigned int numUnknowns = getNumUnknowns( );

        const unsigned int numConstraints = getNumConstraints( );

        KKTMatrix = floatVector( ( numUnknowns + numConstraints ) * ( numUnknowns + numConstraints ), 0 );


        Eigen::Map< Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > K( KKTMatrix.data( ), ( numUnknowns + numConstraints ), ( numUnknowns + numConstraints ) );

        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > J( getFlatJacobian( )->data( ), numUnknowns, numUnknowns );

        K.block( 0, 0, numUnknowns, numUnknowns ) = ( J.transpose( ) * J ).eval( );

        for ( unsigned int I = 0; I < numUnknowns; I++ ){

            KKTMatrix[ ( numUnknowns + numConstraints ) * I + I ] += ( *getMuk( ) );

        }

        for ( unsigned int i = 0; i < numConstraints; i++ ){

            if ( active_constraints[ i ] ){

                for ( unsigned int I = 0; I < numUnknowns; I++ ){

                    KKTMatrix[ ( numUnknowns + numConstraints ) * ( I ) + numUnknowns + i ] = ( *getConstraintJacobians( ) )[ numUnknowns * i + I ];
                    KKTMatrix[ ( numUnknowns + numConstraints ) * ( numUnknowns + i ) + I ] = ( *getConstraintJacobians( ) )[ numUnknowns * i + I ];

                }

            }
            else{

                KKTMatrix[ ( numUnknowns + numConstraints ) * ( numUnknowns + i ) + numUnknowns + i ] = 1;

            }

        }

    }

    void hydraBase::updateKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints ){
        /*!
         * Update the KKTMatrix if the active constraints have changed
         * 
         * \param &KKTMatrix: The Karush-Kuhn-Tucker matrix
         * \param &active_constraints: The vector of currently active constraints.
         */

        const unsigned int numUnknowns = getNumUnknowns( );

        const unsigned int numConstraints = getNumConstraints( );

        Eigen::Map< Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > K( KKTMatrix.data( ), ( numUnknowns + numConstraints ), ( numUnknowns + numConstraints ) );

        K.block( 0, numUnknowns, numUnknowns, numConstraints ).setZero( );

        K.block( numUnknowns, 0, numConstraints, numUnknowns ).setZero( );

        K.block( numUnknowns, numUnknowns, numConstraints, numConstraints ).setZero( );

        for ( unsigned int i = 0; i < numConstraints; i++ ){

            if ( active_constraints[ i ] ){

                for ( unsigned int I = 0; I < numUnknowns; I++ ){

                    KKTMatrix[ ( numUnknowns + numConstraints ) * ( I ) + numUnknowns + i ] = ( *getConstraintJacobians( ) )[ numUnknowns * i + I ];
                    KKTMatrix[ ( numUnknowns + numConstraints ) * ( numUnknowns + i ) + I ] = ( *getConstraintJacobians( ) )[ numUnknowns * i + I ];

                }

            }
            else{

                KKTMatrix[ ( numUnknowns + numConstraints ) * ( numUnknowns + i ) + numUnknowns + i ] = 1;

            }

        }

    }

    void hydraBase::assembleKKTRHSVector( const floatVector &dx, floatVector &KKTRHSVector, const std::vector< bool > &active_constraints ){
        /*!
         * Assemble the right hand side vector for the KKT matrix
         * 
         * \param &dx: The delta vector being solved for
         * \param &KKTRHSVector: The right hand size vector for the KKT matrix
         * \param &active_constraints: The active constraint vector
         */

        const unsigned int numUnknowns = getNumUnknowns( );

        const unsigned int numConstraints = getNumConstraints( );

        KKTRHSVector = floatVector( numUnknowns + numConstraints, 0 );

        Eigen::Map< const Eigen::Vector< floatType, -1 > > _dx( dx.data( ), numUnknowns );

        Eigen::Map< Eigen::Vector< floatType, -1 > > RHS( KKTRHSVector.data( ), ( numUnknowns + numConstraints ), ( numUnknowns + numConstraints ) );

        Eigen::Map< const Eigen::Vector< floatType, -1 > > R( getResidual( )->data( ), numUnknowns );

        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > J( getFlatJacobian( )->data( ), numUnknowns, numUnknowns );

        RHS.head( numUnknowns ) = ( J.transpose( ) * ( R + J * _dx ) + ( *getMuk( ) ) * _dx ).eval( );

        for ( unsigned int i = 0; i < numConstraints; i++ ){

            if ( active_constraints[ i ] ){

                KKTRHSVector[ numUnknowns + i ] = ( *getConstraints( ) )[ i ];

                for ( unsigned int I = 0; I < numUnknowns; I++ ){

                    KKTRHSVector[ numUnknowns + i ] += ( *getConstraintJacobians( ) )[ numUnknowns * i + I ] * dx[ I ];

                }

            }

        }

    }

    void hydraBase::initializeActiveConstraints( std::vector< bool > &active_constraints ){
        /*!
         * Initialize the active constraint vector
         * 
         * \param &active_constraints: The current constraints that are active
         */

        active_constraints = std::vector< bool >( getNumConstraints( ), false );

        for ( auto c = getConstraints( )->begin( ); c != getConstraints( )->end( ); c++ ){

            unsigned int index = ( unsigned int )( c - getConstraints( )->begin( ) );

            active_constraints[ index ] = ( ( *c ) < 0. );

        }

    }

    void hydraBase::solveConstrainedQP( floatVector &dx, const unsigned int kmax ){
        /*!
         * Solve the constrained QP problem to estimate the desired step size
         * 
         * \param &dx: The change in the unknown vector
         * \param kmax: The maximum number of iterations (defaults to 100)
         */

        const unsigned int numUnknowns = getNumUnknowns( );

        const unsigned int numConstraints = getNumConstraints( );

        floatVector K;

        floatVector RHS;

        std::vector< bool > active_constraints;
        initializeActiveConstraints( active_constraints );

        assembleKKTRHSVector( dx, RHS, active_constraints );

        assembleKKTMatrix( K, active_constraints );

        floatType tol = ( *getRelativeTolerance( ) ) * ( tardigradeVectorTools::l2norm( RHS ) ) + ( *getAbsoluteTolerance( ) );

        unsigned int k = 0;

        floatVector y( numUnknowns + numConstraints, 0 );

        floatVector ck = *getConstraints( );

        floatVector ctilde( numConstraints, 0 );

        floatVector negp( numUnknowns, 0 );

        floatVector lambda( numConstraints, 0 );

        floatVector P( numUnknowns + numConstraints, 0 );

        for ( unsigned int i = 0; i < ( numUnknowns + numConstraints ); i++ ){

            P[ i ] = 1 / std::max( std::fabs( *std::max_element( K.begin( ) + ( numUnknowns + numConstraints ) * i,
                                                                 K.begin( ) + ( numUnknowns + numConstraints ) * ( i + 1 ),
                                                                 [ ]( const floatType &a, const floatType &b ){ return std::fabs( a ) < std::fabs( b ); } ) ), 1e-15 );

        }

        Eigen::Map< const Eigen::Vector< floatType, -1 > > _P( P.data( ), numUnknowns + numConstraints );

        tardigradeVectorTools::solverType< floatType > linearSolver;

        while ( k < kmax ){

            Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > _K( K.data( ), numConstraints + numUnknowns, numConstraints + numUnknowns );
            Eigen::Map< const Eigen::Vector< floatType, -1 > > _RHS( RHS.data( ), numConstraints + numUnknowns );
            Eigen::Map< Eigen::Vector< floatType, -1 > > _y( y.data( ), numConstraints + numUnknowns );

            linearSolver = tardigradeVectorTools::solverType< floatType >( _P.asDiagonal( ) * _K );

            _y = linearSolver.solve( _P.asDiagonal( ) * _RHS );

            std::copy( y.begin( ), y.begin( ) + numUnknowns, negp.begin( ) );

            std::copy( y.begin( ) + numUnknowns, y.end( ), lambda.begin( ) );

            if ( tardigradeVectorTools::l2norm( negp ) <= tol ){

                bool negLambda = false;

                floatType minLambda = 1;

                unsigned int imin = 0;

                for ( auto v = std::begin( lambda ); v != std::end( lambda ); v++ ){

                    if ( *v < 0 ){

                        negLambda = true;

                        if ( ( *v ) < minLambda ){

                            imin = ( unsigned int )( v - std::begin( lambda ) );

                            minLambda = *v;

                        }

                    }

                }

                if ( negLambda ){

                    active_constraints[ imin ] = false;

                }
                else{

                    return;

                }

            }
            else{

                ck     = *getConstraints( );
                ctilde = *getConstraints( );
                for ( unsigned int i = 0; i < numConstraints; i++ ){
                    for ( unsigned int j = 0; j < numUnknowns; j++ ){
                        ck[ i ]     += ( *getConstraintJacobians( ) )[ numUnknowns * i + j ] * dx[ j ];
                        ctilde[ i ] += ( *getConstraintJacobians( ) )[ numUnknowns * i + j ] * ( dx[ j ] - negp[ j ] );
                    }
                }

                floatType alpha = 1.0;

                unsigned int iblock = 0;

                bool newBlock = false;

                for ( unsigned int i = 0; i < numConstraints; i++ ){

                    if ( !active_constraints[ i ] ){

                        if ( ctilde[ i ] < -tol ){

                            floatType alpha_trial = -ck[ i ] / ( ctilde[ i ] - ck[ i ] );

                            if ( alpha_trial <= alpha ){

                                iblock = i;

                                alpha = alpha_trial;

                                newBlock = true;

                            }

                        }

                    }

                }

                if ( newBlock ){

                    active_constraints[ iblock ] = true;

                }

                dx -= alpha * negp;

            }

            updateKKTMatrix( K, active_constraints );

            assembleKKTRHSVector(  dx, RHS, active_constraints );

            for ( unsigned int i = 0; i < ( numUnknowns + numConstraints ); i++ ){

                P[ i ] = 1 / std::max( std::fabs( *std::max_element( K.begin( ) + ( numUnknowns + numConstraints ) * i,
                                                                     K.begin( ) + ( numUnknowns + numConstraints ) * ( i + 1 ),
                                                                     [ ]( const floatType &a, const floatType &b ){ return std::fabs( a ) < std::fabs( b ); } ) ), 1e-15 );

            }

            k++;

        }

    }

    void hydraBase::setConstraints( ){
        /*!
         * Set the constraint values
         */

        const unsigned int numConstraints = getNumConstraints( );

        auto constraints = get_setDataStorage_constraints( );

        constraints.zero( numConstraints );

        unsigned int offset = 0;

        setCurrentResidualIndexMeaningful( true );
        for ( auto v = getResidualClasses( )->begin( ); v != getResidualClasses( )->end( ); v++ ){
            setCurrentResidualIndexMeaningful( v - getResidualClasses( )->begin( ) );

            if ( ( *( *v )->getNumConstraints( ) ) > 0 ){

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

        auto constraintJacobians = get_setDataStorage_constraintJacobians( );

        constraintJacobians.zero( numConstraints * numUnknowns );

        unsigned int offset = 0;

        setCurrentResidualIndexMeaningful( true );
        for ( auto v = getResidualClasses( )->begin( ); v != getResidualClasses( )->end( ); v++ ){
            setCurrentResidualIndex( v - getResidualClasses( )->begin( ) );

            if ( ( *( *v )->getNumConstraints( ) ) > 0 ){

                std::copy( ( *v )->getConstraintJacobians( )->begin( ), ( *v )->getConstraintJacobians( )->end( ),  constraintJacobians.value->begin( ) + offset );

                offset += ( *v )->getConstraintJacobians( )->size( );

            }

        }
        setCurrentResidualIndexMeaningful( false );

    }

    void dummyMaterialModel( floatVector &stress,             floatVector &statev,        floatMatrix &ddsdde,       floatType &SSE,            floatType &SPD,
                             floatType &SCD,                  floatType &RPL,             floatVector &ddsddt,       floatVector &drplde,       floatType &DRPLDT,
                             const floatVector &strain,       const floatVector &dstrain, const floatVector &time,   const floatType &DTIME,    const floatType &TEMP,
                             const floatType &DTEMP,          const floatVector &predef,  const floatVector &dpred,  const std::string &cmname, const int &NDI,
                             const int &NSHR,                 const int &NTENS,           const int &NSTATV,         const floatVector &props,  const int &NPROPS,
                             const floatVector &coords,       const floatMatrix &drot,    floatType &PNEWDT,         const floatType &CELENT,   const floatMatrix &dfgrd0,
                             const floatMatrix &dfgrd1,       const int &NOEL,            const int &NPT,            const int &LAYER,          const int &KSPT,
                             const std::vector< int > &jstep, const int &KINC ){
        /*!
         * A template Abaqus c++ UMAT using c++ STL types. Variables in all caps reference ABAQUS FORTRAN
         * memory directly. Variables in lower case are native c++ type conversions stored separately from the original
         * ABAQUS FORTRAN memory.
         */

        //Call functions of constitutive model to do things
        TARDIGRADE_ERROR_TOOLS_CATCH( sayHello( "Abaqus" ) );

        return;
    }

    void abaqusInterface( double *STRESS,       double *STATEV,       double *DDSDDE,       double &SSE,          double &SPD,
                          double &SCD,          double &RPL,          double *DDSDDT,       double *DRPLDE,       double &DRPLDT,
                          const double *STRAN,  const double *DSTRAN, const double *TIME,   const double &DTIME,  const double &TEMP,
                          const double &DTEMP,  const double *PREDEF, const double *DPRED,  const char *CMNAME,   const int &NDI,
                          const int &NSHR,      const int &NTENS,     const int &NSTATV,    const double *PROPS,  const int &NPROPS,
                          const double *COORDS, const double *DROT,   double &PNEWDT,       const double &CELENT, const double *DFGRD0,
                          const double *DFGRD1, const int &NOEL,      const int &NPT,       const int &LAYER,     const int &KSPT,
                          const int *JSTEP,     const int &KINC ){
        /*!
         * A template Abaqus UMAT c++ interface that performs Fortran to C++ type conversions, calculates the material
         * model's expected input, handles tensor shape changes, and calls a c++ material model.
         */

        //Provide a variable string message for error nodes
        std::ostringstream message;

        //Map FORTRAN UMAT variables to C++ types as necessary. Use case sensitivity to distinguish.
        //TODO: Decide if case sensitive variable names is a terrible idea or not
        //Vectors can be created directly with pointer arithmetic
        std::vector< double > stress( STRESS, STRESS + NTENS );
        std::vector< double > statev( STATEV, STATEV + NSTATV );
        std::vector< double > ddsddt( DDSDDT, DDSDDT + NTENS );
        std::vector< double > drplde( DRPLDE, DRPLDE + NTENS );
        const std::vector< double > strain( STRAN, STRAN + NTENS );
        const std::vector< double > dstrain( DSTRAN, DSTRAN + NTENS );
        const std::vector< double > time( TIME, TIME + 2 );
        const std::vector< double > predef( PREDEF, PREDEF + 1 );
        const std::vector< double > dpred( DPRED, DPRED + 1 );
        const std::string cmname( tardigradeAbaqusTools::FtoCString( 80, CMNAME ) );
        const std::vector< double > props( PROPS, PROPS + NPROPS );
        const std::vector< double > coords( COORDS, COORDS + spatialDimensions );
        const std::vector< int > jstep( JSTEP, JSTEP + 4 );
        //Fortran two-dimensional arrays require careful column to row major conversions to c++ types
        std::vector< std::vector< double > > ddsdde = tardigradeAbaqusTools::columnToRowMajor( DDSDDE, NTENS, NTENS );
        const std::vector< std::vector< double > > drot = tardigradeAbaqusTools::columnToRowMajor( DROT, spatialDimensions, spatialDimensions );
        const std::vector< std::vector< double > > dfgrd0 = tardigradeAbaqusTools::columnToRowMajor( DFGRD0, spatialDimensions, spatialDimensions );
        const std::vector< std::vector< double > > dfgrd1 = tardigradeAbaqusTools::columnToRowMajor( DFGRD1, spatialDimensions, spatialDimensions );

        //Verify number of state variables against hydra expectations
        if ( statev.size( ) != nStateVariables ){
            message.clear();
            message << "ERROR:" << __FILENAME__ << "." << __func__ << ": The hydra Abaqus interface requires exactly "
                << nStateVariables << " state variables. Found " << statev.size( ) << ".";
            throw std::runtime_error( message.str( ) );
        }

        //Verify number of material parameters against hydra expectations
        if ( props.size( ) != nMaterialParameters ){
            message.clear();
            message << "ERROR:" << __FILENAME__ << "." << __func__ << ": The hydra Abaqus interface requires exactly "
                << nMaterialParameters << " material constants. Found " << props.size( ) << ".";
            throw std::runtime_error( message.str( ) );
        }

        //Call the constitutive model c++ interface
        if ( KINC == 1 && NOEL == 1 && NPT == 1 ){

            try{

                dummyMaterialModel( stress, statev,  ddsdde, SSE,    SPD,
                                            SCD,    RPL,     ddsddt, drplde, DRPLDT,
                                            strain, dstrain, time,   DTIME,  TEMP,
                                            DTEMP,  predef,  dpred,  cmname, NDI,
                                            NSHR,   NTENS,   NSTATV, props,  NPROPS,
                                            coords, drot,    PNEWDT, CELENT, dfgrd0,
                                            dfgrd1, NOEL,    NPT,    LAYER,  KSPT,
                                            jstep,  KINC );

            }
            catch(std::exception &e){

                if ( PNEWDT >= 1. ){

                    throw e;

                }

            }
        }

        //Re-pack C++ objects into FORTRAN memory to return values to Abaqus
        //Scalars were passed by reference and will update correctly
        //Vectors don't require row/column major considerations, but do require re-packing to the Fortran pointer
        tardigradeAbaqusTools::rowToColumnMajor( STRESS, stress, 1, NTENS );
        tardigradeAbaqusTools::rowToColumnMajor( DDSDDT, ddsddt, 1, NTENS );
        tardigradeAbaqusTools::rowToColumnMajor( DRPLDE, drplde, 1, NTENS );
        tardigradeAbaqusTools::rowToColumnMajor( STATEV, statev, 1, NSTATV );
        //Arrays require vector of vector to column major conversion
        tardigradeAbaqusTools::rowToColumnMajor( DDSDDE, ddsdde, NTENS, NTENS );

    }

}
