/**
  ******************************************************************************
  * \file tardigrade_hydra.cpp
  ******************************************************************************
  * A C++ library for defining frameworks to solve finite deformation material
  * models.
  ******************************************************************************
  */

#include<tardigrade_hydra.h>

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

    hydraBase::hydraBase( const floatType &time, const floatType &deltaTime,
                          const floatType &temperature, const floatType &previousTemperature,
                          const floatVector &deformationGradient, const floatVector &previousDeformationGradient,
                          const floatVector &previousStateVariables, const floatVector &parameters,
                          const unsigned int numConfigurations, const unsigned int numNonLinearSolveStateVariables,
                          const unsigned int dimension, const unsigned int configuration_unknown_count, const floatType tolr, const floatType tola, const unsigned int maxIterations,
                          const unsigned int maxLSIterations, const floatType lsAlpha ) : _dimension( dimension ),
                                                           _configuration_unknown_count( configuration_unknown_count ),
                                                           _time( time ), _deltaTime( deltaTime ),
                                                           _temperature( temperature ), _previousTemperature( previousTemperature ),
                                                           _deformationGradient( deformationGradient ),
                                                           _previousDeformationGradient( previousDeformationGradient ),
                                                           _previousStateVariables( previousStateVariables ),
                                                           _parameters( parameters ),
                                                           _numConfigurations( numConfigurations ),
                                                           _numNonLinearSolveStateVariables( numNonLinearSolveStateVariables ),
                                                           _tolr( tolr ), _tola( tola ),
                                                           _maxIterations( maxIterations ), _maxLSIterations( maxLSIterations ),
                                                           _lsAlpha( lsAlpha ){
        /*!
         * The main constructor for the hydra base class. Inputs are all the required values for most solves.
         * 
         * \param &time: The current time
         * \param &deltaTime: The change in time
         * \param &temperature: The current temperature
         * \param &previousTemperature: The previous temperature
         * \param &deformationGradient: The current deformation gradient
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
         */

        // Decompose the state variable vector initializing all of the configurations
        decomposeStateVariableVector( );

        // Set the residual classes
        setResidualClasses( );

    }

    void hydraBase::setStress( const floatVector &stress ){
        /*!
         * Set the value of the stress
         * 
         * \param &stress: The stress in row-major form
         */

        setIterationData( stress, _stress );

    }

    void hydraBase::extractStress( ){
        /*!
         * Extract the stresses out of the unknown vector
         */

        const floatVector *unknownVector = getUnknownVector( );

        setStress( floatVector( unknownVector->begin( ),
                                unknownVector->begin( ) + *getConfigurationUnknownCount( ) ) );

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

        floatVector eye( sot_dim );
        tardigradeVectorTools::eye( eye );

        // Set the configurations
        configurations = floatVector( num_configs * sot_dim, 0 );

        inverseConfigurations = floatVector( num_configs * sot_dim, 0 );

        // Initialize the first configuration with the total deformation gradient
        floatMatrix local_configurations( num_configs );
        floatMatrix local_inverseConfigurations( num_configs );

        local_configurations[ 0 ] = total_transformation;

        for ( int i = num_configs - 2; i >= 0; i-- ){

            // Set the current configuration as being equal to the previous
            local_configurations[ i + 1 ] = floatVector( data_vector->begin( ) + i * sot_dim + start_index,
                                                         data_vector->begin( ) + ( i + 1 ) * sot_dim + start_index );

            if ( add_eye ){

                local_configurations[ i + 1 ] += eye;

            }

            // Compute the inverse of the current configuration and store it
            local_inverseConfigurations[ i + 1 ] = local_configurations[ i + 1 ];
            Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> > mat( local_inverseConfigurations[ i + 1 ].data(), 3, 3 );
            mat = mat.inverse( ).eval( );

            // Add contribution of deformation gradient to the first configuration
            local_configurations[ 0 ] = tardigradeVectorTools::matrixMultiply( local_configurations[ 0 ], local_inverseConfigurations[ i + 1 ],
                                                                               dim, dim, dim, dim );

        }

        local_inverseConfigurations[ 0 ] = local_configurations[ 0 ];
        Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> > mat( local_inverseConfigurations[ 0 ].data(), 3, 3 );
        mat = mat.inverse( ).eval( );

        configurations        = tardigradeVectorTools::appendVectors( local_configurations );
        inverseConfigurations = tardigradeVectorTools::appendVectors( local_inverseConfigurations );

        return;

    }

    void hydraBase::decomposeUnknownVector( ){
        /*!
         * Decompose the unknown vector into the cauchy stress, configurations, and state variables used for the non-linear solve
         */

        const floatVector *unknownVector = getUnknownVector( );

        // Set the stress
        extractStress( );

        // Set the configurations
        floatVector configurations;

        floatVector inverseConfigurations;

        computeConfigurations( unknownVector, getStress( )->size( ), *getDeformationGradient( ), configurations, inverseConfigurations );

        set_configurations( configurations );

        set_inverseConfigurations( inverseConfigurations );

        // Extract the remaining state variables required for the non-linear solve
        set_nonLinearSolveStateVariables( floatVector( unknownVector->begin( ) + ( *getNumConfigurations( ) ) * ( *getConfigurationUnknownCount( ) ),
                                                       unknownVector->end( ) ) );

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
        floatVector eye( getSOTDimension( ), 0 );
        tardigradeVectorTools::eye( eye );

        if ( getPreviousStateVariables( )->size( ) < ( ( ( *nConfig ) - 1 ) * ( *getConfigurationUnknownCount( ) ) + ( *nNLISV ) ) ){

            std::string message = "The number of state variables is less than required for the configurations and ";
            message            += "non-linear state variables\n";
            message            += "  # previousStateVariables                               : " + std::to_string( getPreviousStateVariables( )->size( ) ) + "\n";
            message            += "  # ( configurations - 1 ) * configuration_unknown_count : " + std::to_string( ( ( *nConfig ) - 1 ) * ( *getConfigurationUnknownCount( ) ) ) + "\n";
            message            += "  # non-linear solve ISVs                                : " + std::to_string( ( *nNLISV ) ) + "\n";
            message            += "  # minimum required ISVs                                : " + std::to_string( ( ( *nConfig ) - 1 ) * ( *getConfigurationUnknownCount( ) ) + ( *nNLISV ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

        }

        floatVector configurations;

        floatVector previousConfigurations;

        floatVector inverseConfigurations;

        floatVector previousInverseConfigurations;

        // Compute the configurations
        computeConfigurations( getPreviousStateVariables( ), 0, *getDeformationGradient( ), configurations, inverseConfigurations, true );

        computeConfigurations( getPreviousStateVariables( ), 0, *getPreviousDeformationGradient( ), previousConfigurations, previousInverseConfigurations, true );

        // Set the configurations

        set_configurations( configurations );

        set_inverseConfigurations( inverseConfigurations );

        set_previousConfigurations( previousConfigurations );

        set_previousInverseConfigurations( previousInverseConfigurations );

        // Extract the remaining state variables required for the non-linear solve
        set_previousNonLinearSolveStateVariables( floatVector( getPreviousStateVariables( )->begin( ) + ( ( *nConfig ) - 1 ) * ( *getConfigurationUnknownCount( ) ),
                                                               getPreviousStateVariables( )->begin( ) + ( ( *nConfig ) - 1 ) * ( *getConfigurationUnknownCount( ) ) + *nNLISV ) );

        set_nonLinearSolveStateVariables( *get_previousNonLinearSolveStateVariables( ) );

        // Extract the additional state variables
        set_previousAdditionalStateVariables( floatVector( getPreviousStateVariables( )->begin( ) + ( ( *nConfig ) - 1 ) * ( *getConfigurationUnknownCount( ) ) + *nNLISV,
                                                           getPreviousStateVariables( )->end( ) ) );

        set_additionalStateVariables( *get_previousAdditionalStateVariables( ) );

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

        floatVector Fsc( sot_dim, 0 );
        for ( unsigned int i = 0; i < 3; i++ ){ Fsc[ dim * i + i ] = 1.; }

        for ( unsigned int i = lowerIndex; i < upperIndex; i++ ){

            Fsc = tardigradeVectorTools::matrixMultiply( Fsc, floatVector( configurations.begin( ) + sot_dim * i, configurations.begin( ) + sot_dim * ( i + 1 ) ), dim, dim, dim, dim );

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

            floatVector Fm, FpT;

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

    floatVector hydraBase::getSubConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex ){
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

    floatVector hydraBase::getPrecedingConfiguration( const unsigned int &index ){
        /*!
         * Get the sub-configuration preceding but not including the index
         * 
         * \param &index: The index of the configuration immediately following the sub-configuration
         */

        return getSubConfiguration( 0, index );

    }

    floatVector hydraBase::getFollowingConfiguration( const unsigned int &index ){
        /*!
         * Get the sub-configuration following but not including the index
         * 
         * \param &index: The index of the current configuration immediately before the sub-configuration
         */

        return getSubConfiguration( index + 1, *getNumConfigurations( ) );

    }

    floatVector hydraBase::getConfiguration( const unsigned int &index ){
        /*!
         * Get the configuration indicated by the provided index
         * 
         * \param &index: The index of the current configuration to be extracted
         */

        return getSubConfiguration( index, index + 1 );

    }

    floatVector hydraBase::getPreviousSubConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex ){
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

    floatVector hydraBase::getPreviousPrecedingConfiguration( const unsigned int &index ){
        /*!
         * Get the previous sub-configuration preceding but not including the index
         * 
         * \param &index: The index of the configuration immediately following the sub-configuration
         */

        return getPreviousSubConfiguration( 0, index );

    }

    floatVector hydraBase::getPreviousFollowingConfiguration( const unsigned int &index ){
        /*!
         * Get the previous sub-configuration following but not including the index
         * 
         * \param &index: The index of the current configuration immediately before the sub-configuration
         */

        return getPreviousSubConfiguration( index + 1, *getNumConfigurations( ) );

    }

    floatVector hydraBase::getPreviousConfiguration( const unsigned int &index ){
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

    void hydraBase::calculateFirstConfigurationJacobians( const floatVector &configurations, floatVector &dC1dC, floatVector &dC1dCn ){
        /*!
         * Get the Jacobian of the first configuration w.r.t. the total mapping and the remaining configurations.
         * 
         * \param &configurations: The configurations which describe the mapping from the current to the reference configuration
         * \param &dC1dC: The Jacobian of the first entry w.r.t. the total
         * \param &dC1dCn: The Jacobian of the first entry w.r.t. the remaining terms
         * 
         * whre \f$C^n = C^2, C^3, \cdots \f$
         */

        const unsigned int dim = getDimension( );
        const unsigned int sot_dim = getSOTDimension( );
        const unsigned int num_configs = *getNumConfigurations( );

        dC1dC  = floatVector( sot_dim * sot_dim, 0 );

        dC1dCn = floatVector( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

        floatVector eye( sot_dim );
        tardigradeVectorTools::eye( eye );

        floatVector fullConfiguration = getSubConfiguration( configurations, 0, num_configs );

        floatVector invCsc = getSubConfiguration( configurations, 1, num_configs );
        Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor > > mat( invCsc.data( ), dim, dim );
        mat = mat.inverse( ).eval( );

        floatVector dInvCscdCsc = tardigradeVectorTools::computeFlatDInvADA( invCsc, dim, dim );

        floatVector dInvCscdCs = tardigradeVectorTools::matrixMultiply( dInvCscdCsc, getSubConfigurationJacobian( configurations, 1, num_configs ), sot_dim, sot_dim, sot_dim, num_configs * sot_dim );

        // Compute the gradients
        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int barI = 0; barI < dim; barI++ ){

                for ( unsigned int a = 0; a < dim; a++ ){

                    for ( unsigned int A = 0; A < dim; A++ ){

                        dC1dC[ dim * sot_dim * i + sot_dim * barI + dim * a + A ] += eye[ dim * i + a ] * invCsc[ dim * A + barI ];

                        for ( unsigned int index = 0; index < num_configs - 1; index++ ){

                            for ( unsigned int J = 0; J < dim; J++ ){

                                dC1dCn[ dim * ( num_configs - 1 ) * sot_dim * i + ( num_configs - 1 ) * sot_dim * barI + sot_dim * index + dim * a + A ]
                                    += fullConfiguration[ dim * i + J ]
                                     * dInvCscdCs[ dim * num_configs * sot_dim * J + num_configs * sot_dim * barI + sot_dim * ( index + 1 ) + dim * a + A ];

                            }

                        }

                    }

                }

            }

        }

    }

    void hydraBase::setFirstConfigurationJacobians( ){
        /*!
         * Set the Jacobians of the first configuration w.r.t. the total configuration and the remaining sub-configurations
         */

        floatVector dF1dF;

        floatVector dF1dFn;

        calculateFirstConfigurationJacobians( *get_configurations( ), dF1dF, dF1dFn );

        set_dF1dF( dF1dF );

        set_dF1dFn( dF1dFn );

    }

    void hydraBase::setPreviousFirstConfigurationJacobians( ){
        /*!
         * Set the Jacobians of the previous first configuration w.r.t. the total configuration and the remaining sub-configurations
         */

        floatVector dF1dF;

        floatVector dF1dFn;

        calculateFirstConfigurationJacobians( *get_previousConfigurations( ), dF1dF, dF1dFn );

        set_previousdF1dF( dF1dF );

        set_previousdF1dFn( dF1dFn );

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

    void hydraBase::formNonLinearProblem( ){
        /*!
         * Form the residual, jacobian, and gradient matrices
         */

        unsigned int residualSize = ( *getNumConfigurations( ) ) * ( *getConfigurationUnknownCount( ) ) + *getNumNonLinearSolveStateVariables( );

        _residual.second = floatVector( residualSize, 0 );

        _jacobian.second = floatVector( residualSize * residualSize, 0 );

        _dRdF.second = floatVector( residualSize * ( *getConfigurationUnknownCount( ) ), 0 );

        _dRdT.second = floatVector( residualSize, 0 );

        _additionalDerivatives.second.clear( );

        unsigned int offset = 0;

        unsigned int numAdditionalDerivatives = 0;

        for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); residual_ptr++ ){

            residualBase *residual = ( *residual_ptr );

            // Extract the terms

            const floatVector* localResidual;
            TARDIGRADE_ERROR_TOOLS_CATCH( localResidual = residual->getResidual( ) );

            const floatMatrix* localJacobian;
            TARDIGRADE_ERROR_TOOLS_CATCH( localJacobian = residual->getJacobian( ) );

            const floatMatrix* localdRdF;
            TARDIGRADE_ERROR_TOOLS_CATCH( localdRdF = residual->getdRdF( ) );

            const floatVector* localdRdT;
            TARDIGRADE_ERROR_TOOLS_CATCH( localdRdT = residual->getdRdT( ) );

            const floatMatrix* localAdditionalDerivatives;
            TARDIGRADE_ERROR_TOOLS_CATCH( localAdditionalDerivatives = residual->getAdditionalDerivatives( ) );

            // Check the contributions to make sure they are consistent sizes

            if ( localResidual->size( ) != *residual->getNumEquations( ) ){

                std::string message = "The residual for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n";
                message            += "  expected: " + std::to_string( *residual->getNumEquations( ) ) + "\n";
                message            += "  actual:   " + std::to_string( localResidual->size( ) ) + "\n";

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

            }

            if ( localJacobian->size( ) != *residual->getNumEquations( ) ){

                std::string message = "The jacobian for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n";
                message            += "  expected: " + std::to_string( *residual->getNumEquations( ) ) + "\n";
                message            += "  actual:   " + std::to_string( localJacobian->size( ) ) + "\n";

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

            }

            if ( localdRdF->size( ) != *residual->getNumEquations( ) ){

                std::string message = "dRdF for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n";
                message            += "  expected: " + std::to_string( *residual->getNumEquations( ) ) + "\n";
                message            += "  actual:   " + std::to_string( localdRdF->size( ) ) + "\n";

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

            }

            if ( localdRdT->size( ) != *residual->getNumEquations( ) ){

                std::string message = "dRdT for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n";
                message            += "  expected: " + std::to_string( *residual->getNumEquations( ) ) + "\n";
                message            += "  actual:   " + std::to_string( localdRdT->size( ) ) + "\n";

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

            }

            if ( localAdditionalDerivatives->size( ) != 0 ){

                if ( localAdditionalDerivatives->size( ) != *residual->getNumEquations( ) ){

                    std::string message = "additionalDerivatives for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n";
                    message            += "  expected: " + std::to_string( *residual->getNumEquations( ) ) + "\n";
                    message            += "  actual:   " + std::to_string( localAdditionalDerivatives->size( ) ) + "\n";
    
                    TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

                }

                if ( ( *localAdditionalDerivatives )[ 0 ].size( ) != numAdditionalDerivatives ){
    
                    if ( ( residual_ptr - getResidualClasses( )->begin( ) ) == 0 ){
    
                        numAdditionalDerivatives = ( *localAdditionalDerivatives )[ 0 ].size( );
    
                        _additionalDerivatives.second = floatVector( residualSize * numAdditionalDerivatives, 0 );
    
                    }
                    else{
    
                        std::string message = "The additional derivatives for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " are not the expected length as determined from the first residual\n";
                        message            += "  expected: " + std::to_string( numAdditionalDerivatives ) + "\n";
                        message            += "  actual:   " + std::to_string( ( *localAdditionalDerivatives )[ 0 ].size( ) ) + "\n";
    
                        TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );
    
                    }
    
                }

            }

            // Store the values in the global quantities

            for ( unsigned int row = 0; row < *residual->getNumEquations( ); row++ ){

                _residual.second[ row + offset ] = ( *localResidual )[ row ];

                if ( ( *localJacobian )[ row ].size( ) != residualSize ){

                    std::string message = "Row " + std::to_string( row ) + " of the jacobian for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n";
                    message            += "  expected: " + std::to_string( residualSize ) + "\n";
                    message            += "  actual:   " + std::to_string( ( *localJacobian )[ row ].size( ) ) + "\n";

                    TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

                }

                for ( unsigned int col = 0; col < residualSize; col++ ){
                
                    _jacobian.second[ residualSize * ( row + offset ) + col ] = ( *localJacobian )[ row ][ col ];

                }

                if ( ( *localdRdF )[ row ].size( ) != *getConfigurationUnknownCount( ) ){

                    std::string message = "Row " + std::to_string( row ) + " of dRdF for residual " + std::to_string( residual_ptr - getResidualClasses( )->begin( ) ) + " is not the expected length\n";
                    message            += "  expected: " + std::to_string( ( *getConfigurationUnknownCount( ) ) ) + "\n";
                    message            += "  actual:   " + std::to_string( ( *localJacobian )[ row ].size( ) ) + "\n";

                    TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

                }

                for ( unsigned int col = 0; col < ( *getConfigurationUnknownCount( ) ); col++ ){

                    _dRdF.second[ ( *getConfigurationUnknownCount( ) ) * ( row + offset ) + col ] = ( *localdRdF )[ row ][ col ];

                }

                _dRdT.second[ row + offset ] = ( *localdRdT )[ row ];

                for ( unsigned int col = 0; col < numAdditionalDerivatives; col++ ){

                    _additionalDerivatives.second[ numAdditionalDerivatives * ( row + offset ) + col ] = ( *localAdditionalDerivatives )[ row ][ col ];

                }

            }

            offset += *residual->getNumEquations( );

        }

        _residual.first = true;

        _jacobian.first = true;

        _dRdF.first = true;

        _dRdT.first = true;

        _additionalDerivatives.first = true;

        addIterationData( &_residual );

        addIterationData( &_jacobian );

        addIterationData( &_dRdF );

        addIterationData( &_dRdT );

        addIterationData( &_additionalDerivatives );

    }

    const floatVector* hydraBase::getResidual( ){
        /*!
         * Get the residual vector for the non-linear problem
         */

        if ( !_residual.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( formNonLinearProblem( ) );

        }

        return &_residual.second;

    }

    const floatVector* hydraBase::getFlatJacobian( ){
        /*!
         * Get the flattened row-major jacobian for the non-linear problem
         */

        if ( !_jacobian.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( formNonLinearProblem( ) );

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

            TARDIGRADE_ERROR_TOOLS_CATCH( formNonLinearProblem( ) );

        }

        return &_dRdF.second;

    }

    floatMatrix hydraBase::getdRdF( ){
        /*!
         * Get dRdF for the non-linear problem
         */

        return tardigradeVectorTools::inflate( *getFlatdRdF( ), getResidual( )->size( ), getSOTDimension( ) );
    }

    const floatVector* hydraBase::getdRdT( ){
        /*!
         * Get dRdT for the non-linear problem
         */

        if ( !_dRdT.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( formNonLinearProblem( ) );

        }

        return &_dRdT.second;

    }

    const floatVector* hydraBase::getFlatAdditionalDerivatives( ){
        /*!
         * Get the flattened row-major additional derivatives for the non-linear problem
         */

        if ( !_additionalDerivatives.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( formNonLinearProblem( ) );

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
    errorOut sayHello( std::string message ) {
        if ( message.compare( "George" ) == 0 ){
            errorOut result = new errorNode( __func__, "ERROR: George is a wolf in sheep's clothing!");
            return result;
        }
        std::cout << "Hello " << message << std::endl;
        return NULL;
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

        floatMatrix Xmat( 1 + num_local_configs );

        Xmat[ 0 ] = *cauchyStress;

        for ( unsigned int i = 1; i < num_local_configs; i++ ){

            Xmat[ i ] = floatVector( configurations->begin( ) + sot_dim * i, configurations->begin( ) + sot_dim * ( i + 1 ) );

        }

        Xmat[ Xmat.size( ) - 1 ] = *nonLinearSolveStateVariables;

        setX( tardigradeVectorTools::appendVectors( Xmat ) );

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

        floatVector tolerance = tardigradeVectorTools::abs( *getResidual( ) ) + tardigradeVectorTools::abs( *getUnknownVector( ) );

        tolerance = *getRelativeTolerance( ) * tolerance + *getAbsoluteTolerance( );

        setTolerance( tolerance );

    }

    void hydraBase::setTolerance( const floatVector &tolerance ){
        /*!
         * Set the tolerance
         *
         * \param tolerance: The tolerance vector for each value of the residual
         */

        setConstantData( tolerance, _tolerance );

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

        // Reset all of the iteration data
        resetIterationData( );

        // Set the unknown vector
        setX( newUnknownVector );

        // Decompose the unknown vector and update the state
        TARDIGRADE_ERROR_TOOLS_CATCH( decomposeUnknownVector( ) );

    }

    void hydraBase::solveNonLinearProblem( ){
        /*!
         * Solve the non-linear problem
         */

        // Form the initial unknown vector
        TARDIGRADE_ERROR_TOOLS_CATCH( initializeUnknownVector( ) );

        unsigned int rank;

        floatVector deltaX;

        resetLSIteration( );

        while( !checkConvergence( ) && checkIteration( ) ){

            floatVector X0 = *getUnknownVector( );

            TARDIGRADE_ERROR_TOOLS_CATCH( deltaX = -tardigradeVectorTools::solveLinearSystem( *getFlatJacobian( ), *getResidual( ),
                                                                         getResidual( )->size( ), getResidual( )->size( ), rank ) );

            if ( rank != getResidual( )->size( ) ){

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "The Jacobian is not full rank" ) );

            }

            updateUnknownVector( X0 + *getLambda( ) * deltaX );

            while ( !checkLSConvergence( ) && checkLSIteration( ) ){

                updateLambda( );

                incrementLSIteration( );

                updateUnknownVector( X0 + *getLambda( ) * deltaX );

            }

            if ( !checkLSConvergence( ) ){

                throw convergence_error( "Failure in line search" );

            }

            resetLSIteration( );

            // Increment the iteration count
            incrementIteration( );

        }

        if ( !checkConvergence( ) ){

            throw convergence_error( "Failure to converge main loop" );

        }

    }

    void hydraBase::evaluate( ){
        /*!
         * Solve the non-linear problem and update the variables
         */

        try{

            solveNonLinearProblem( );

        }
        catch( const convergence_error &e ){

            throw;

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

        tardigradeVectorTools::solverType< floatType > solver( Amat );

        unsigned int rank = solver.rank( );

        TARDIGRADE_ERROR_TOOLS_CATCH(
            if ( rank != getResidual( )->size( ) ){
                throw std::runtime_error( "The Jacobian is not full rank" );
            }
        )

        // Solve for dXdF
        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > dRdFmat( getFlatdRdF( )->data( ), getResidual( )->size( ), *getConfigurationUnknownCount( ) );

        _flatdXdF.second = floatVector( getUnknownVector( )->size( ) * ( *getConfigurationUnknownCount( ) ) );
        Eigen::Map< Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > dXdFmat( _flatdXdF.second.data( ), getUnknownVector( )->size( ), ( *getConfigurationUnknownCount( ) ) );

        dXdFmat = -solver.solve( dRdFmat );

        _flatdXdF.first = true;

        // Solve for dXdT
        Eigen::Map< const Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > dRdTmat( getdRdT( )->data( ), getResidual( )->size( ), 1 );

        _flatdXdT.second = floatVector( getUnknownVector( )->size( ) );
        Eigen::Map< Eigen::Matrix< floatType, -1, -1, Eigen::RowMajor > > dXdTmat( _flatdXdT.second.data( ), getUnknownVector( )->size( ), 1 );

        dXdTmat = -solver.solve( dRdTmat );

        _flatdXdT.first = true;

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

    errorOut dummyMaterialModel( floatVector &stress,             floatVector &statev,        floatMatrix &ddsdde,       floatType &SSE,            floatType &SPD,
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
        errorOut error = sayHello( "Abaqus" );

        //Error handling
        if ( error ){
            errorOut result = new errorNode( __func__, "Error when calling sayHello" );
            result->addNext( error );
            return result;
        }

        return NULL;
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

        //Initialize error return codes
        errorOut error = NULL;

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
            error = dummyMaterialModel( stress, statev,  ddsdde, SSE,    SPD,
                                        SCD,    RPL,     ddsddt, drplde, DRPLDT,
                                        strain, dstrain, time,   DTIME,  TEMP,
                                        DTEMP,  predef,  dpred,  cmname, NDI,
                                        NSHR,   NTENS,   NSTATV, props,  NPROPS,
                                        coords, drot,    PNEWDT, CELENT, dfgrd0,
                                        dfgrd1, NOEL,    NPT,    LAYER,  KSPT,
                                        jstep,  KINC );
        }

        //Error handling
        if ( error ){
            message.clear();
            message << "ERROR:" << __FILENAME__ << "." << __func__ << ": Error when calling dummyMaterialModel.";
            errorOut result = new errorNode( __func__, message.str( ) );
            result->addNext( error );
            error->print( true );
            //If an error was thrown, but the ratio of new/current time increment is not updated, it was a fatal error.
            if ( PNEWDT >= 1. ){
                throw std::runtime_error( message.str( ) );
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
