/**
  ******************************************************************************
  * \file tardigrade_hydraMicromorphic.cpp
  ******************************************************************************
  * A C++ utility for constructing finite deformation micromorphic constitutive
  * models.
  ******************************************************************************
  */

#include<tardigrade_hydraMicromorphic.h>

namespace tardigradeHydra{

    hydraBaseMicromorphic::hydraBaseMicromorphic( const floatType &time, const floatType &deltaTime,
                                                  const floatType &temperature, const floatType &previousTemperature,
                                                  const secondOrderTensor &deformationGradient, const secondOrderTensor &previousDeformationGradient,
                                                  const secondOrderTensor &microDeformation,    const secondOrderTensor &previousMicroDeformation,
                                                  const thirdOrderTensor  &gradientMicroDeformation,  const thirdOrderTensor  &previousGradientMicroDeformation,
                                                  const floatVector &additionalDOF, const floatVector &previousAdditionalDOF,
                                                  const floatVector &previousStateVariables, const floatVector &parameters,
                                                  const unsigned int numConfigurations, const unsigned int numNonLinearSolveStateVariables,
                                                  const unsigned int dimension, const unsigned int configuration_unknown_count,
                                                  const floatType tolr, const floatType tola, const unsigned int maxIterations,
                                                  const unsigned int maxLSIterations, const floatType lsAlpha,
                                                  const bool use_preconditioner, const unsigned int preconditioner_type ) :
                                                  hydraBase( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                                             additionalDOF, previousAdditionalDOF,
                                                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables,
                                                             dimension, configuration_unknown_count, tolr, tola, maxIterations, maxLSIterations, lsAlpha,
                                                             use_preconditioner, preconditioner_type ),
                                                  _microDeformation( microDeformation ), _previousMicroDeformation( previousMicroDeformation ),
                                                  _gradientMicroDeformation( gradientMicroDeformation ), _previousGradientMicroDeformation( previousGradientMicroDeformation ){

        /*!
         * The main constructor for the micromorphic hydra base class. Inputs are all the required values for most solves.
         * 
         * \param &time: The current time
         * \param &deltaTime: The change in time
         * \param &temperature: The current temperature
         * \param &previousTemperature: The previous temperature
         * \param &deformationGradient: The current deformation gradient
         * \param &previousDeformationGradient The previous deformation gradient
         * \param &microDeformation: The current micro-deformation \f$ \chi \f$
         * \param &previousMicroDeformation: The previous micro-deformation \f$ \chi \f$
         * \param &gradientMicroDeformation: The current reference spatial gradient of the micro-deformation \f$ \frac{\partial}{\partial X} \chi \f$
         * \param &previousGradientMicroDeformation: The previous reference spatial gradient of the micro-deformation \f$ \frac{\partial}{\partial X} \chi \f$
         * \param &previousStateVariables: The previous state variables
         * \param &parameters: The model parameters
         * \param &numConfigurations: The number of configurations
         * \param &numNonLinearSolveStateVariables: The number of state variables which will contribute terms to the non-linear solve's residual
         * \param &dimension: The dimension of the problem (defaults to 3)
         * \param &configuration_unknown_count: The number of unknowns in each configuration (defaults to 27)
         * \param &tolr: The relative tolerance (defaults to 1e-9)
         * \param &tola: The absolute tolerance (defaults to 1e-9)
         * \param &maxIterations: The maximum number of non-linear iterations (defaults to 20)
         * \param &maxLSIterations: The maximum number of line-search iterations (defaults to 5)
         * \param &lsAlpha: The alpha term for the line search (defaults to 1e-4)
         */

        setScaledQuantities( );

        decomposeStateVariableVectorMicroConfigurations( );

    }

    void hydraBaseMicromorphic::setScaledQuantities( ){
        /*!
         * Scale the current values by the scale factor
         */

        hydraBase::hydraBase::setScaledQuantities( );

        _scaled_microDeformation = ( *getScaleFactor( ) ) * ( _microDeformation - _previousMicroDeformation ) + _previousMicroDeformation;

        _scaled_gradientMicroDeformation = ( *getScaleFactor( ) ) * ( _gradientMicroDeformation - _previousGradientMicroDeformation ) + _previousGradientMicroDeformation;

    }

    void hydraBaseMicromorphic::initializeUnknownVector( ){
        /*!
         * Initialize the unknown vector for the non-linear solve.
         * 
         * \f$X = \left\{ \bf{\sigma}, \bf{F}^2, \bf{F}^3, ..., \bf{F}n, \bf{\chi}^2, \bf{\chi}^3, ..., \frac{\partial}{\partial \bm{X}^2} \bf{\chi}^2, \frac{\partial}{\partial \bm{X}^3} \bf{\chi}^3, ..., \xi^1, \xi^2, ..., \xi^m \right\} \f$
         * 
         * It is assumed that the first residual calculation also has a method `void getStress( )`
         * which returns a pointer to the current value of the stress.
         */

        const unsigned int sot_dim = getSOTDimension( );
        const unsigned int tot_dim = getTOTDimension( );

        const floatVector *stress;
        TARDIGRADE_ERROR_TOOLS_CATCH( stress = getStress( ) );

        const floatVector *configurations = get_configurations( );

        const floatVector *microConfigurations = get_microConfigurations( );

        const floatVector *gradientMicroConfigurations = get_gradientMicroConfigurations( );

        const floatVector *nonLinearSolveStateVariables = get_nonLinearSolveStateVariables( );

        floatVector X( getNumUnknowns( ), 0 );

        unsigned int offset = 0;

        std::copy( std::begin( *stress ), std::end( *stress ), std::begin( X ) );

        offset += stress->size( );

        // Set the initial values of the macro configurations
        std::copy( std::begin( *configurations ) + sot_dim, std::end( *configurations ), std::begin( X ) + offset );

        offset += configurations->size( ) - sot_dim;

        // Set the values of the micro configurations
        std::copy( std::begin( *microConfigurations ) + sot_dim, std::end( *microConfigurations ), std::begin( X ) + offset );

        offset += microConfigurations->size( ) - sot_dim;

        // Set the values of the micro-gradient configurations
        std::copy( std::begin( *gradientMicroConfigurations ) + tot_dim, std::end( *gradientMicroConfigurations ), std::begin( X ) + offset );

        offset += gradientMicroConfigurations->size( ) - tot_dim;

        std::copy( std::begin( *nonLinearSolveStateVariables ), std::end( *nonLinearSolveStateVariables ), std::begin( X ) + offset );

        bool resetRequired = false;

        for ( auto residual_ptr = getResidualClasses( )->begin( ); residual_ptr != getResidualClasses( )->end( ); residual_ptr++ ){

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

        if ( resetRequired ){

            updateUnknownVector( X );

        }
        else{

            setX( X );

        }

    }

    void hydraBaseMicromorphic::decomposeStateVariableVector( ){
        /*!
         * Decompose the incoming state variable vector setting the different configurations along the way
         * 
         * The state variable vector is assumed to be of the form:
         * 
         * \f$ \text{ISV} = \left\{\bf{F}^2 - \bf{I}, \bf{F}^3 - \bf{I}, \cdots, \bf{F}^n - \bf{I}, \bf{\chi}^2 - \bf{I}, \bf{\chi}^3 - \bf{I}, \cdots, \bf{\chi}^n - \bf{I} \frac{\partial}{\partial \bf{X}} \bf{\chi}^2, \frac{\partial}{\partial \bf{X}} \bf{\chi}^3, \cdots, \frac{\partial}{\partial \bf{X}} \bf{\chi}^n, \xi^1, \xi^2, \cdots, \xi^m, \eta^1, \cdots\right\} \f$
         * 
         * where the \f$\bf{F}\f$ are the different deformation gradients (configurations), \f$\bf{\chi}\f$ are the micro-deformations,
         * \f$\xi^y\f$ are the other variables to be solved during the non-linear solve, and \f$\eta^z\f$ are other state variables. Note
         * that we decompose the deformation gradient and micro-deformation as
         * 
         * \f$\bf{F} = \bf{F}^1 \bf{F}^2 \cdots \bf{F}^n\f$
         * 
         * \f$\bf{\chi} = \bf{\chi}^1 \bf{\chi}^2 \cdots \bf{\chi}^n\f$
         * 
         * and so because \f$\bf{F}\f$ and \f$\bf{\chi}\f$ are provided we can solve for \f$\bf{F}^1\f$ and \f$\bf{\chi}\f$. Typically,
         * this configuration would be the elastic configuration (i.e., the configuration that generates the stress) though we do not insist that users follow convention.
         * 
         * NOTE: Though we overload the decomposeStateVariableVector in hydraBase this function will not be called in hydraBase's constructor because
         *       virtual functions do not come into being during the construction of parent classes constructors. We could work around this but instead
         *       we will overload and define a local method to do the decomposition of the micro-deformation tensors.
         */

        // Call the parent class decomposition
        hydraBase::decomposeStateVariableVector( );

        // Decompose the micro-deformation
        decomposeStateVariableVectorMicroConfigurations( );

    }

    void hydraBaseMicromorphic::decomposeUnknownVector( ){
        /*!
         * Decompose the incoming unknown vector setting the different configurations along the way
         * 
         * The state variable vector is assumed to be of the form:
         * 
         * \f$ \text{ISV} = \left\{\bf{S}^2, \bf{\Sigma}, \bf{M}, \bf{F}^2, \bf{F}^3, \cdots, \bf{F}^n, \bf{\chi}^2, \bf{\chi}^3, \cdots, \bf{\chi}^n \frac{\partial}{\partial \bf{X}} \bf{\chi}^2, \frac{\partial}{\partial \bf{X}} \bf{\chi}^3, \cdots, \frac{\partial}{\partial \bf{X}} \bf{\chi}^n, \xi^1, \xi^2, \cdots, \xi^m\right\} \f$
         * 
         * where \f$\bf{S}^2\f$ is the second Piola Kirchhoff stress, \f$\bf{\Sigma}\f$ is the reference symmetric micro stress, and \f$\bf{M}\f$
         * is the reference higher order stress the \f$\bf{F}\f$ are the different deformation gradients (configurations), \f$\bf{\chi}\f$ are the micro-deformations,
         * \f$\xi^y\f$ are the other variables to be solved during the non-linear solve, and \f$\eta^z\f$ are other state variables. Note
         * that we decompose the deformation gradient and micro-deformation as
         * 
         * \f$\bf{F} = \bf{F}^1 \bf{F}^2 \cdots \bf{F}^n\f$
         * 
         * \f$\bf{\chi} = \bf{\chi}^1 \bf{\chi}^2 \cdots \bf{\chi}^n\f$
         * 
         * and so because \f$\bf{F}\f$ and \f$\bf{\chi}\f$ are provided we can solve for \f$\bf{F}^1\f$ and \f$\bf{\chi}\f$. Typically,
         * this configuration would be the elastic configuration (i.e., the configuration that generates the stress) though we do not insist that users follow convention.
         * 
         * NOTE: Though we overload the decomposeStateVariableVector in hydraBase this function will not be called in hydraBase's constructor because
         *       virtual functions do not come into being during the construction of parent classes constructors. We could work around this but instead
         *       we will overload and define a local method to do the decomposition of the micro-deformation tensors.
         */

        // Call the parent class decomposition
        hydraBase::decomposeUnknownVector( );

        // Decompose the micro-deformation
        decomposeUnknownVectorMicroConfigurations( );

    }

    void hydraBaseMicromorphic::computeGradientMicroConfigurations( const floatVector *data_vector, unsigned int start_index,
                                                                    const floatVector &configurations,             const floatVector &microConfigurations,
                                                                    const floatVector &gradientMicroConfiguration, floatVector &gradientMicroConfigurations ){
        /*!
         * Compute the gradient of the micro-configurations in their reference configurations
         * 
         * \param *data_vector: The vector of data
         * \param &start_index: The index at which to start reading the data from data_vector
         * \param &configurations: The macro-scale configurations
         * \param &microConfigurations: The micro-configurations
         * \param &gradientMicroConfiguration: The gradient of the total micro-configuration w.r.t. the reference configuration
         * \param &gradientMicroConfigurations: The resulting gradients of the micro configurations
         */

        const unsigned int tot_dim = getTOTDimension( );
        const unsigned int num_configs = *getNumConfigurations( );

        gradientMicroConfigurations = tardigradeVectorTools::appendVectors( { thirdOrderTensor( tot_dim, 0 ),
                                                                              floatVector( data_vector->begin( ) + start_index,
                                                                              data_vector->begin( ) + start_index + ( num_configs - 1 ) * tot_dim ) } );

        calculateFirstConfigurationGradChi( configurations, microConfigurations, gradientMicroConfiguration, gradientMicroConfigurations );

    }

    void hydraBaseMicromorphic::decomposeUnknownVectorMicroConfigurations( ){
        /*!
         * Decompose the micro-deformation parts of the unknown vector
         */

        const unsigned int sot_dim = getSOTDimension( );
        const unsigned int num_configs = *getNumConfigurations( );

        unsigned int start_index = getStress( )->size( ) + ( num_configs - 1 ) * sot_dim;

        auto microConfigurations         = get_setDataStorage_microConfigurations( );

        auto inverseMicroConfigurations  = get_setDataStorage_inverseMicroConfigurations( );

        auto gradientMicroConfigurations = get_setDataStorage_gradientMicroConfigurations( );

        // Compute the micro-configurations

        computeConfigurations( getUnknownVector( ), start_index, *getMicroDeformation( ), *microConfigurations.value, *inverseMicroConfigurations.value );

        start_index += ( num_configs - 1 ) * sot_dim;

        computeGradientMicroConfigurations( getUnknownVector( ), start_index, *get_configurations( ), *microConfigurations.value,
                                            *getGradientMicroDeformation( ), *gradientMicroConfigurations.value );

    }

    void hydraBaseMicromorphic::decomposeStateVariableVectorMicroConfigurations( ){
        /*!
         * Decompose the micro-deformation parts of the state variable vector
         */

        const unsigned int sot_dim = getSOTDimension( );
        const unsigned int num_configs = *getNumConfigurations( );

        unsigned int start_index = ( num_configs - 1 ) * sot_dim;

        auto microConfigurations = get_setDataStorage_microConfigurations( );

        auto inverseMicroConfigurations  = get_setDataStorage_inverseMicroConfigurations( );

        auto gradientMicroConfigurations = get_setDataStorage_gradientMicroConfigurations( );

        auto previousMicroConfigurations = get_setDataStorage_previousMicroConfigurations( );

        auto previousInverseMicroConfigurations  = get_setDataStorage_previousInverseMicroConfigurations( );

        auto previousGradientMicroConfigurations = get_setDataStorage_previousGradientMicroConfigurations( );

        // Compute the micro-configurations

        computeConfigurations( getPreviousStateVariables( ), start_index, *getMicroDeformation( ), *microConfigurations.value, *inverseMicroConfigurations.value, true );

        computeConfigurations( getPreviousStateVariables( ), start_index, *getPreviousMicroDeformation( ), *previousMicroConfigurations.value, *previousInverseMicroConfigurations.value, true );

        start_index += ( num_configs - 1 ) * sot_dim;

        computeGradientMicroConfigurations( getPreviousStateVariables( ), start_index, *get_configurations( ), *microConfigurations.value,
                                            *getGradientMicroDeformation( ), *gradientMicroConfigurations.value );

        computeGradientMicroConfigurations( getPreviousStateVariables( ), start_index, *get_previousConfigurations( ), *previousMicroConfigurations.value,
                                            *getPreviousGradientMicroDeformation( ), *previousGradientMicroConfigurations.value );

    }

    secondOrderTensor hydraBaseMicromorphic::getSubMicroConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex ){
        /*!
         * Get a sub-micro configuration \f$\bf{\chi}^{sc}\f$ defined as
         *
         * \f$ \chi^{sc}_{iI} = \chi^{\text{lowerIndex}}_{i\hat{I}} \chi^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots \chi^{\text{upperIndex-1}}_{\bar{I}I} \f$
         * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
         * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
         *   Note, the configuration indicated by the index is NOT included in the sub-configuration
         */

        return getSubConfiguration( *get_microConfigurations( ), lowerIndex, upperIndex );

    }

    secondOrderTensor hydraBaseMicromorphic::getPrecedingMicroConfiguration( const unsigned int &index ){
        /*!
         * Get the sub-micro configuration preceding but not including the index
         * 
         * \param &index: The index of the configuration immediately following the sub-micro configuration
         */

        return getSubMicroConfiguration( 0, index );

    }

    secondOrderTensor hydraBaseMicromorphic::getFollowingMicroConfiguration( const unsigned int &index ){
        /*!
         * Get the sub-micro configuration following but not including the index
         * 
         * \param &index: The index of the current configuration immediately before the sub-micro configuration
         */

        return getSubMicroConfiguration( index + 1, *getNumConfigurations( ) );

    }

    secondOrderTensor hydraBaseMicromorphic::getMicroConfiguration( const unsigned int &index ){
        /*!
         * Get the micro configuration indicated by the provided index
         * 
         * \param &index: The index of the current configuration to be extracted
         */

        return getSubMicroConfiguration( index, index + 1 );

    }

    secondOrderTensor hydraBaseMicromorphic::getPreviousSubMicroConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex ){
        /*!
         * Get a previous sub-micro configuration \f$\bf{\chi}^{sc}\f$ defined as
         *
         * \f$ \chi^{sc}_{iI} = \chi^{\text{lowerIndex}}_{i\hat{I}} \chi^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots \chi^{\text{upperIndex-1}}_{\bar{I}I} \f$
         * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
         * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
         *   Note, the configuration indicated by the index is NOT included in the sub-configuration
         */

        return getSubConfiguration( *get_previousMicroConfigurations( ), lowerIndex, upperIndex );

    }

    secondOrderTensor hydraBaseMicromorphic::getPreviousPrecedingMicroConfiguration( const unsigned int &index ){
        /*!
         * Get the previous sub-micro configuration preceding but not including the index
         * 
         * \param &index: The index of the configuration immediately following the sub-micro configuration
         */

        return getPreviousSubMicroConfiguration( 0, index );

    }

    secondOrderTensor hydraBaseMicromorphic::getPreviousFollowingMicroConfiguration( const unsigned int &index ){
        /*!
         * Get the previous sub-micro configuration following but not including the index
         * 
         * \param &index: The index of the current configuration immediately before the sub-micro configuration
         */

        return getPreviousSubMicroConfiguration( index + 1, *getNumConfigurations( ) );

    }

    secondOrderTensor hydraBaseMicromorphic::getPreviousMicroConfiguration( const unsigned int &index ){
        /*!
         * Get the previous micro configuration indicated by the provided index
         * 
         * \param &index: The index of the current configuration to be extracted
         */

        return getPreviousSubMicroConfiguration( index, index + 1 );

    }

    floatVector hydraBaseMicromorphic::getSubMicroConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex ){
        /*!
         * Get the jacobian of a sub-micro configuration \f$\bf{\chi}^{sc}\f$ defined as
         *
         * \f$ \chi^{sc}_{iI} = \chi^{\text{lowerIndex}}_{i\hat{I}} \chi^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots \chi^{\text{upperIndex-1}}_{\bar{I}I} \f$
         * 
         * with respect to the current configurations.
         *
         * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
         * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
         *   Note, the configuration indicated by the index is NOT included in the sub-configuration
         */

        return getSubConfigurationJacobian( *get_microConfigurations( ), lowerIndex, upperIndex );

    }

    floatVector hydraBaseMicromorphic::getPrecedingMicroConfigurationJacobian( const unsigned int &index ){
        /*!
         * Get the jacobian of the sub-micro configuration preceding but not including the index with respect to the current configurations.
         * 
         * \param &index: The index of the configuration immediately following the sub-configuration
         */

        return getSubMicroConfigurationJacobian( 0, index );

    }

    floatVector hydraBaseMicromorphic::getFollowingMicroConfigurationJacobian( const unsigned int &index ){
        /*!
         * Get the jacobian of the sub-micro configuration following but not including the index with respect to the current configurations.
         * 
         * \param &index: The index of the current configuration immediately before the sub-configuration
         */

        return getSubMicroConfigurationJacobian( index + 1, *getNumConfigurations( ) );

    }

    floatVector hydraBaseMicromorphic::getPreviousSubMicroConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex ){
        /*!
         * Get the jacobian of a previous sub-micro configuration \f$\bf{\chi}^{sc}\f$ defined as
         *
         * \f$ \chi^{sc}_{iI} = \chi^{\text{lowerIndex}}_{i\hat{I}} \chi^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots \chi^{\text{upperIndex-1}}_{\bar{I}I} \f$
         * 
         * with respect to the current configurations.
         *
         * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
         * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
         *   Note, the configuration indicated by the index is NOT included in the sub-configuration
         */

        return getSubConfigurationJacobian( *get_previousMicroConfigurations( ), lowerIndex, upperIndex );

    }

    floatVector hydraBaseMicromorphic::getPreviousPrecedingMicroConfigurationJacobian( const unsigned int &index ){
        /*!
         * Get the jacobian of the previous sub-micro configuration preceding but not including the index with respect to the current configurations.
         * 
         * \param &index: The index of the configuration immediately following the sub-configuration
         */

        return getPreviousSubMicroConfigurationJacobian( 0, index );

    }

    floatVector hydraBaseMicromorphic::getPreviousFollowingMicroConfigurationJacobian( const unsigned int &index ){
        /*!
         * Get the jacobian of the previous sub-micro configuration following but not including the index with respect to the current configurations.
         * 
         * \param &index: The index of the current configuration immediately before the sub-configuration
         */

        return getPreviousSubMicroConfigurationJacobian( index + 1, *getNumConfigurations( ) );

    }

    void hydraBaseMicromorphic::setFirstMicroConfigurationJacobians( ){
        /*!
         * Set the Jacobians of the first micro configuration w.r.t. the total micro configuration and the remaining sub-micro configurations
         */

        auto dChi1dChi  = get_setDataStorage_dChi1dChi( ); 

        auto dChi1dChin = get_setDataStorage_dChi1dChin( );

        calculateFirstConfigurationJacobians( *get_microConfigurations( ), *dChi1dChi.value, *dChi1dChin.value );

    }

    void hydraBaseMicromorphic::setPreviousFirstMicroConfigurationJacobians( ){
        /*!
         * Set the Jacobians of the previous first micro configuration w.r.t. the total micro configuration and the remaining sub-micro configurations
         */

        auto previousdChi1dChi  = get_setDataStorage_previousdChi1dChi( );

        auto previousdChi1dChin = get_setDataStorage_previousdChi1dChin( );

        calculateFirstConfigurationJacobians( *get_previousMicroConfigurations( ), *previousdChi1dChi.value, *previousdChi1dChin.value );

    }

    void hydraBaseMicromorphic::setFirstGradientMicroConfigurationJacobians( ){
        /*!
         * Set the Jacobians of the gradient of the first micro configuration w.r.t. the total micro configuration and the remaining sub-micro configurations
         */

        auto dGradChi1dFn       = get_setDataStorage_dGradChi1dFn( );

        auto dGradChi1dChi      = get_setDataStorage_dGradChi1dChi( );

        auto dGradChi1dChin     = get_setDataStorage_dGradChi1dChin( );

        auto dGradChi1dGradChi  = get_setDataStorage_dGradChi1dGradChi( );

        auto dGradChi1dGradChin = get_setDataStorage_dGradChi1dGradChin( );

        calculateFirstConfigurationGradChiJacobian( *get_configurations( ), *get_microConfigurations( ),
                                                    *getGradientMicroDeformation( ), *get_gradientMicroConfigurations( ),
                                                    *get_dChi1dChi( ), *get_dChi1dChin( ),
                                                    *dGradChi1dFn.value, *dGradChi1dChi.value, *dGradChi1dChin.value,
                                                    *dGradChi1dGradChi.value, *dGradChi1dGradChin.value );

    }

    void hydraBaseMicromorphic::setPreviousFirstGradientMicroConfigurationJacobians( ){
        /*!
         * Set the Jacobians of the previous gradient of the first micro configuration w.r.t. the total micro configuration and the remaining sub-micro configurations
         */

        auto previousdGradChi1dFn       = get_setDataStorage_previousdGradChi1dFn( );

        auto previousdGradChi1dChi      = get_setDataStorage_previousdGradChi1dChi( );

        auto previousdGradChi1dChin     = get_setDataStorage_previousdGradChi1dChin( );

        auto previousdGradChi1dGradChi  = get_setDataStorage_previousdGradChi1dGradChi( );

        auto previousdGradChi1dGradChin = get_setDataStorage_previousdGradChi1dGradChin( );

        calculateFirstConfigurationGradChiJacobian( *get_previousConfigurations( ), *get_previousMicroConfigurations( ),
                                                    *getPreviousGradientMicroDeformation( ), *get_previousGradientMicroConfigurations( ),
                                                    *get_previousdChi1dChi( ), *get_previousdChi1dChin( ),
                                                    *previousdGradChi1dFn.value, *previousdGradChi1dChi.value, *previousdGradChi1dChin.value,
                                                    *previousdGradChi1dGradChi.value, *previousdGradChi1dGradChin.value );

    }

    void hydraBaseMicromorphic::calculateFirstConfigurationGradChi( const floatVector &configurations, const floatVector &microConfigurations, const thirdOrderTensor &gradientMicroConfiguration, floatVector &gradientMicroConfigurations ){
        /*!
         * Calculate the value of the gradient of the first micro-configuration given all of the configurations, the micro-configurations,
         * the spatial gradient of the micro deformation in the reference configuration, and the gradients of the micro-configurations
         * other than the first in their own reference configurations.
         * 
         * \param &configurations: The configuration matrix
         * \param &microConfigurations: The micro-configuration matrix
         * \param &gradientMicroConfiguration: The gradient of the micro-deformation in the reference configuration
         * \param &gradientMicroConfigurations: The matrix of gradients of the micro-configurations in their reference configurations
         */

        const unsigned int dim = getDimension( );
        const unsigned int sot_dim = getSOTDimension( );
        const unsigned int tot_dim = getTOTDimension( );
        const unsigned int num_configs = *getNumConfigurations( );

        // Compute the gradient in the reference configuration
        thirdOrderTensor gradientChi1Reference( tot_dim, 0 );// = gradientMicroConfiguration; // Initialize to the total gradient in the reference configuration

        thirdOrderTensor  temp_tot1( tot_dim, 0 );

        secondOrderTensor chiPrecede( sot_dim, 0. );
        for ( unsigned int i = 0; i < dim; i++ ){ chiPrecede[ dim * i + i ] = 1.; }

        Eigen::Map< Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor > > chip_map( chiPrecede.data( ), dim, dim );
        Eigen::Map< const Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor > > map1( NULL, dim, dim );

        for ( unsigned int index = 1; index < num_configs; index++ ){

            new (&map1) Eigen::Map< const Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor > >( microConfigurations.data( ) + ( index - 1 ) * sot_dim, 3, 3 );
            chip_map *= map1;

            // Add the contribution of the term
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < dim; j++ ){

                    for ( unsigned int I = 0; I < dim; I++ ){

                        for ( unsigned int J = 0; J < dim; J++ ){

                                    gradientChi1Reference[ dim * dim * i + dim * I + J ]
                                        -= chiPrecede[ dim * i + j ] * gradientMicroConfigurations[ tot_dim * index + dim * dim * j + dim * I + J ];

                        }

                    }

                }

            }

            if ( index == ( num_configs - 1 ) ){ break; }

//            std::copy( gradientChi1Reference.begin( ), gradientChi1Reference.end( ), temp_tot1.begin( ) );
            std::fill( temp_tot1.begin( ), temp_tot1.end( ), 0 );

            secondOrderTensor chiFollow  = getSubConfiguration( microConfigurations, index + 1, *getNumConfigurations( ) );

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int I = 0; I < dim; I++ ){

                    for ( unsigned int J = 0; J < dim; J++ ){

                        for ( unsigned int k = 0; k < dim; k++ ){

                                    temp_tot1[ dim * dim * i + dim * I + J ]
                                        += chiFollow[ dim * k + I ] * gradientChi1Reference[ dim * dim * i + dim * k + J ];

                        }

                    }

                }

            }

//            std::copy( gradientChi1Reference.begin( ), gradientChi1Reference.end( ), temp_tot1.begin( ) );
            std::fill( gradientChi1Reference.begin( ), gradientChi1Reference.end( ), 0 );

            secondOrderTensor FFollow = getSubConfiguration( configurations, index + 1, *getNumConfigurations( ) );

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int I = 0; I < dim; I++ ){

                    for ( unsigned int l = 0; l < dim; l++ ){

                        for ( unsigned int J = 0; J < dim; J++ ){

                                    gradientChi1Reference[ dim * dim * i + dim * I + J ]
                                        += FFollow[ dim * l + J ] * temp_tot1[ dim * dim * i + dim * I + l ];

                        }

                    }

                }

            }

        }

        gradientChi1Reference += gradientMicroConfiguration;

        // Map the gradient of the micro-configuration to the reference of the first configuration
        secondOrderTensor invChiFollow = getSubConfiguration( microConfigurations, 1, num_configs );
        Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> > mat( invChiFollow.data(), 3, 3 );
        mat = mat.inverse( ).eval( );

        secondOrderTensor invFFollow = getSubConfiguration( configurations, 1, num_configs );
        new (&mat) Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> >( invFFollow.data(), 3, 3 );
        mat = mat.inverse( ).eval( );

        std::fill( gradientMicroConfigurations.begin( ), gradientMicroConfigurations.begin() + tot_dim, 0. );

        std::fill( temp_tot1.begin( ), temp_tot1.end( ), 0 );

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int a = 0; a < dim; a++ ){

                for ( unsigned int I = 0; I < dim; I++ ){

                    for ( unsigned int J = 0; J < dim; J++ ){

                        temp_tot1[ dim * dim * i + dim * I + J ] += gradientChi1Reference[ sot_dim * i + dim * a + J ] * invChiFollow[ dim * a + I ];

                    }

                }

            }

        }

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int I = 0; I < dim; I++ ){

                for ( unsigned int b = 0; b < dim; b++ ){

                    for ( unsigned int J = 0; J < dim; J++ ){

                        gradientMicroConfigurations[ dim * dim * i + dim * I + J ] += temp_tot1[ sot_dim * i + dim * I + b ] * invFFollow[ dim * b + J ];

                    }

                }

            }

        }

    }

    void hydraBaseMicromorphic::calculateFirstConfigurationGradChiJacobian( const floatVector &configurations,                  const floatVector &microConfigurations,
                                                                            const thirdOrderTensor &gradientMicroConfiguration, const floatVector &gradientMicroConfigurations,
                                                                            const fourthOrderTensor &dChi1dChi,                 const floatVector &dChi1dChin,
                                                                            floatVector &dGradChi1dCn,
                                                                            fifthOrderTensor &dGradChi1dChi,                    floatVector &dGradChi1dChin,
                                                                            sixthOrderTensor &dGradChi1dGradChi,                floatVector &dGradChi1dGradChin ){
        /*!
         * Calculate the value of the jacobian of the gradient of the first micro-configuration given all of the configurations, the micro-configurations,
         * the spatial gradient of the micro deformation in the reference configuration, and the gradients of the micro-configurations
         * other than the first in their own reference configurations.
         * 
         * \param &configurations: The configuration matrix
         * \param &microConfigurations: The micro-configuration matrix
         * \param &gradientMicroConfiguration: The gradient of the micro-deformation in the reference configuration
         * \param &gradientMicroConfigurations: The gradient of the micro-deformations in their reference configurations
         * \param &dChi1dChi: The gradient of the first micro sub-configuration w.r.t. the total micro deformation
         * \param &dChi1dChin: The gradient of the first micro sub-configuration w.r.t. the remaining sub-micro configurations
         * \param &dGradChi1dCn: The Jacobian of the gradient of the first micro-configuration w.r.t. the remaining configurations
         * \param &dGradChi1dChi: The Jacobian of the gradient of the first micro-configuration w.r.t. the total micro-configuration
         * \param &dGradChi1dChin: The Jacobian of the gradient of the first micro-configuration w.r.t. the remaining micro-configurations
         * \param &dGradChi1dGradChi: The Jacobian of the gradient of the first micro-configuration w.r.t. the gradient of the total micro-configuration
         * \param &dGradChi1dGradChin: The Jacobian of the gradient of the first micro-configuration w.r.t. the gradient of the remaining sub micro-configurations
         */

        const unsigned int dim = getDimension( );
        const unsigned int sot_dim = getSOTDimension( );
        const unsigned int tot_dim = getTOTDimension( );
        const unsigned int num_configs = *getNumConfigurations( );

        // Compute the gradient in the reference configuration
        thirdOrderTensor gradientChi1Reference( tot_dim, 0 );

        floatVector dGradientChi1ReferencedCn( tot_dim * ( num_configs - 1 ) * sot_dim, 0 );

        fifthOrderTensor dGradientChi1ReferencedChi( tot_dim * sot_dim, 0 );

        floatVector dGradientChi1ReferencedChin( tot_dim * ( num_configs - 1 ) * sot_dim, 0 );

        floatVector dGradientChi1ReferencedGradChin( tot_dim * ( num_configs - 1 ) * tot_dim, 0 );

        thirdOrderTensor temp_tot1( tot_dim, 0 );

        thirdOrderTensor temp_tot2( tot_dim, 0 );

        thirdOrderTensor temp_tot2a( tot_dim, 0 );

        thirdOrderTensor temp_tot3( tot_dim, 0 );

        thirdOrderTensor temp_tot3a( tot_dim, 0 );

        secondOrderTensor chiPrecede( sot_dim, 0 ), chiFollow( sot_dim, 0 ), FFollow( sot_dim, 0 );

        for ( unsigned int index = 1; index < num_configs; index++ ){

            chiPrecede = getSubConfiguration( microConfigurations, 0, index );

            // Set the Jacobians of the mapping terms
            floatVector dFFollowdCs = getSubConfigurationJacobian( configurations, index + 1, num_configs );

            floatVector dChiPrecededChis = getSubConfigurationJacobian( microConfigurations, 0, index );

            fourthOrderTensor dChiPrecededChi( sot_dim * sot_dim, 0 );

            floatVector dChiPrecededChin( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int k = 0; k < sot_dim; k++ ){

                    for ( unsigned int j = 0; j < sot_dim; j++ ){

                        dChiPrecededChi[ sot_dim * i + j ] += dChiPrecededChis[ num_configs * sot_dim * i + k ] * dChi1dChi[ sot_dim * k + j ];

                    }

                }
                for ( unsigned int j = 0; j < sot_dim * ( num_configs - 1 ); j++ ){

                    dChiPrecededChin[ ( num_configs - 1 ) * sot_dim * i + j ] += dChiPrecededChis[ num_configs * sot_dim * i + j + sot_dim ];

                    for ( unsigned int k = 0; k < sot_dim; k++ ){

                        dChiPrecededChin[ ( num_configs - 1 ) * sot_dim * i + j ] += dChiPrecededChis[ num_configs * sot_dim * i + k ] * dChi1dChin[ ( num_configs - 1 ) * sot_dim * k + j ];

                    }

                }

            }

            floatVector dChiFollowdChis = getSubConfigurationJacobian( microConfigurations, index + 1, num_configs );

            // Add the contribution of the term
            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int j = 0; j < dim; j++ ){

                    for ( unsigned int I = 0; I < dim; I++ ){

                        for ( unsigned int J = 0; J < dim; J++ ){

                                    gradientChi1Reference[ dim * dim * i + dim * I + J ]
                                        -= chiPrecede[ dim * i + j ] * gradientMicroConfigurations[ tot_dim * index + dim * dim * j + dim * I + J ];

                        }

                    }

                }

            }

            std::copy( gradientChi1Reference.begin( ), gradientChi1Reference.end( ), temp_tot3.begin( ) );

            if ( index != ( num_configs - 1 ) ){

                std::copy( gradientChi1Reference.begin( ), gradientChi1Reference.end( ), temp_tot1.begin( ) );
                std::fill( gradientChi1Reference.begin( ), gradientChi1Reference.end( ), 0 );
                std::fill( temp_tot2.begin( ), temp_tot2.end( ), 0 );

                chiFollow  = getSubConfiguration( microConfigurations, index + 1, *getNumConfigurations( ) );

                for ( unsigned int i = 0; i < dim; i++ ){

                    for ( unsigned int I = 0; I < dim; I++ ){

                        for ( unsigned int J = 0; J < dim; J++ ){

                            for ( unsigned int k = 0; k < dim; k++ ){

                                gradientChi1Reference[ dim * dim * i + dim * I + J ]
                                    += chiFollow[ dim * k + I ] * temp_tot1[ dim * dim * i + dim * k + J ];
                                temp_tot2[ dim * dim * i + dim * I + J ]
                                    -= chiFollow[ dim * k + I ] * gradientMicroConfigurations[ tot_dim * index + dim * dim * i + dim * k + J ];

                            }

                        }

                    }

                }

                std::copy( gradientChi1Reference.begin( ), gradientChi1Reference.end( ), temp_tot1.begin( ) );
                std::fill( gradientChi1Reference.begin( ), gradientChi1Reference.end( ), 0 );
                std::fill( temp_tot2a.begin( ), temp_tot2a.end( ), 0 );
                std::fill( temp_tot3a.begin( ), temp_tot3a.end( ), 0 );

                FFollow = getSubConfiguration( configurations, index + 1, *getNumConfigurations( ) );

                for ( unsigned int i = 0; i < dim; i++ ){

                    for ( unsigned int I = 0; I < dim; I++ ){

                        for ( unsigned int l = 0; l < dim; l++ ){

                            for ( unsigned int J = 0; J < dim; J++ ){

                                gradientChi1Reference[ dim * dim * i + dim * I + J ]
                                    += FFollow[ dim * l + J ] * temp_tot1[ dim * dim * i + dim * I + l ];

                                temp_tot2a[ dim * dim * i + dim * I + J ] += FFollow[ dim * l + J ] * temp_tot2[ dim * dim * i + dim * I + l ];

                                temp_tot3a[ dim * dim * i + dim * I +  J ] += FFollow[ dim * l + J ] * temp_tot3[ dim * dim * i + dim * I + l ];

                                for ( unsigned int A = 0; A < ( num_configs - 1 ) * sot_dim; A++ ){

                                    dGradientChi1ReferencedCn[ dim * dim * ( num_configs - 1 ) * sot_dim * i + dim * ( num_configs - 1 ) * sot_dim * I + ( num_configs - 1 ) * sot_dim * J + A ]
                                         += dFFollowdCs[ dim * num_configs * sot_dim * l + num_configs * sot_dim * J + sot_dim + A ]
                                         * temp_tot1[ dim * dim * i + dim * I + l ];

                                }

                            }

                        }

                    }

                }

                std::copy( temp_tot2a.begin( ), temp_tot2a.end( ), temp_tot2.begin( ) );
                std::copy( temp_tot3a.begin( ), temp_tot3a.end( ), temp_tot3.begin( ) );

            }
            else{

                std::fill( chiFollow.begin( ), chiFollow.end( ), 0. );

                std::fill( FFollow.begin( ), FFollow.end( ), 0. );

                std::transform( gradientMicroConfigurations.begin( ) + tot_dim * index,
                                gradientMicroConfigurations.begin( ) + tot_dim * ( index + 1 ),
                                temp_tot2.begin( ), std::negate<floatType>( ) );

                for ( unsigned int i = 0; i < dim; i++ ){ chiFollow[ dim * i + i ] = 1.; FFollow[ dim * i + i ] = 1.; }

                for ( unsigned int i = 0; i < dim; i++ ){

                    for ( unsigned int I = 0; I < dim; I++ ){

                        for ( unsigned int l = 0; l < dim; l++ ){

                            for ( unsigned int J = 0; J < dim; J++ ){

                                for ( unsigned int A = 0; A < ( num_configs - 1 ) * sot_dim; A++ ){
    
                                    dGradientChi1ReferencedCn[ dim * dim * ( num_configs - 1 ) * sot_dim * i + dim * ( num_configs - 1 ) * sot_dim * I + ( num_configs - 1 ) * sot_dim * J + A ]
                                         += dFFollowdCs[ dim * num_configs * sot_dim * l + num_configs * sot_dim * J + sot_dim + A ]
                                         * gradientChi1Reference[ dim * dim * i + dim * I + l ];
    
                                }

                            }

                        }

                    }

                }

            }

            for ( unsigned int i = 0; i < dim; i++ ){

                for ( unsigned int I = 0; I < dim; I++ ){

                    for ( unsigned int J = 0; J < dim; J++ ){

                        for ( unsigned int j = 0; j < dim; j++ ){

                            for ( unsigned int A = 0; A < sot_dim; A++ ){

                                dGradientChi1ReferencedChi[ dim * dim * sot_dim * i + dim * sot_dim * I + sot_dim * J + A ]
                                    += dChiPrecededChi[ dim * sot_dim * i + sot_dim * j + A ] * temp_tot2[ dim * dim * j + dim * I + J ];

                            }

                            for ( unsigned int A = 0; A < ( num_configs - 1 ) * sot_dim; A++ ){

                                        dGradientChi1ReferencedChin[ dim * dim * ( num_configs - 1 ) * sot_dim * i + dim * ( num_configs - 1 ) * sot_dim * I + ( num_configs - 1 ) * sot_dim * J + A ]
                                            += dChiPrecededChin[ dim * ( num_configs - 1 ) * sot_dim * i + ( num_configs - 1 ) * sot_dim * j + A ] * temp_tot2[ dim * dim * j + dim * I + J ]
                                             + dChiFollowdChis[ dim * num_configs * sot_dim * j + num_configs * sot_dim * I + A + sot_dim ] * temp_tot3[ dim * dim * i + dim * j + J ];

                            }

                            for ( unsigned int k = 0; k < dim; k++ ){

                                for ( unsigned int l = 0; l < dim; l++ ){

                                    dGradientChi1ReferencedGradChin[ dim * dim * ( num_configs - 1 ) * tot_dim * i + dim * ( num_configs - 1 ) * tot_dim * I + ( num_configs - 1 ) * tot_dim * J + tot_dim * ( index - 1 ) + dim * dim * j + dim * k + l ]
                                        -= chiPrecede[ dim * i + j ] * chiFollow[ dim * k + I ] * FFollow[ dim * l + J ];

                                }

                            }

                        }

                    }

                }

            }

        }

        gradientChi1Reference += gradientMicroConfiguration;

        // Map the gradient of the micro-configuration to the reference of the first configuration
        chiFollow = getSubConfiguration( microConfigurations, 1, num_configs );

        FFollow = getSubConfiguration( configurations, 1, num_configs );

        secondOrderTensor invChiFollow = chiFollow;
        Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> > mat( invChiFollow.data(), 3, 3 );
        mat = mat.inverse( ).eval( );
        secondOrderTensor invChiFollow_T = invChiFollow;
        new (&mat) Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> >( invChiFollow_T.data(), 3, 3 );
        mat = mat.transpose( ).eval( );

        secondOrderTensor invFFollow = FFollow;
        new (&mat) Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> >( invFFollow.data(), 3, 3 );
        mat = mat.inverse( ).eval( );
        secondOrderTensor invFFollow_T = invFFollow;
        new (&mat) Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> >( invFFollow_T.data(), 3, 3 );
        mat = mat.transpose( ).eval( );

        floatVector dChiFollowdChis = getSubConfigurationJacobian( microConfigurations, 1, num_configs );

        floatVector dFFollowdFs = getSubConfigurationJacobian( configurations, 1, num_configs );

        fourthOrderTensor dInvChiFollowdChiFollow = tardigradeVectorTools::computeFlatDInvADA( invChiFollow, dim, dim );

        fourthOrderTensor dInvFFollowdFFollow = tardigradeVectorTools::computeFlatDInvADA( invFFollow, dim, dim );

        floatVector dInvChiFollowdChin( sot_dim * sot_dim * ( num_configs - 1 ), 0 );

        floatVector dInvFFollowdFn( sot_dim * sot_dim * ( num_configs - 1 ), 0 );

        for ( unsigned int i = 0; i < sot_dim; i++ ){

            for ( unsigned int k = 0; k < sot_dim; k++ ){

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){

                    dInvChiFollowdChin[ ( num_configs - 1 ) * sot_dim * i + j ] += dInvChiFollowdChiFollow[ sot_dim * i + k ] * dChiFollowdChis[ num_configs * sot_dim * k + j + sot_dim ];

                    dInvFFollowdFn[ ( num_configs - 1 ) * sot_dim * i + j ] += dInvFFollowdFFollow[ sot_dim * i + k ] * dFFollowdFs[ num_configs * sot_dim * k + j + sot_dim ];

                }

            }

        }

        dGradChi1dChi      = fifthOrderTensor( tot_dim * sot_dim, 0 );

        dGradChi1dChin     = floatVector( tot_dim * ( num_configs - 1 ) * sot_dim, 0 );

        dGradChi1dGradChi  = sixthOrderTensor( tot_dim * tot_dim, 0 );

        dGradChi1dCn       = floatVector( tot_dim * ( num_configs - 1 ) * sot_dim, 0 );

        dGradChi1dGradChin = floatVector( tot_dim * ( num_configs - 1 ) * tot_dim, 0 );

        std::fill( temp_tot1.begin( ), temp_tot1.end( ), 0. );
        std::fill( temp_tot2.begin( ), temp_tot2.end( ), 0. );

        fifthOrderTensor temp_fiot( tot_dim * sot_dim, 0 );

        floatVector temp_siot1( tot_dim * ( num_configs - 1 ) * sot_dim, 0 );
        floatVector temp_siot2( tot_dim * ( num_configs - 1 ) * sot_dim, 0 );

        floatVector temp_seot( tot_dim * ( num_configs - 1 ) * tot_dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int I = 0; I < dim; I++ ){

                for ( unsigned int J = 0; J < dim; J++ ){

                    for ( unsigned int a = 0; a < dim; a++ ){

                        temp_tot1[ dim * dim * i + dim * I + J ] += gradientChi1Reference[ dim * dim * i + dim * a + J ] * invChiFollow_T[ dim * I + a ];

                        temp_tot2[ dim * dim * i + dim * I + J ] += gradientChi1Reference[ dim * dim * i + dim * I + a ] * invFFollow_T[ dim * J + a ];

                        for ( unsigned int b = 0; b < dim; b++ ){

                            dGradChi1dGradChi[ dim * dim * tot_dim * i + dim * tot_dim * I + tot_dim * J + dim * dim * i + dim * a + b ]
                                += invChiFollow_T[ dim * I + a ] * invFFollow_T[ dim * J + b ];

                            for ( unsigned int k = 0; k < dim; k++ ){

                                    temp_fiot[ dim * dim * sot_dim * i + dim * sot_dim * I + sot_dim * J + dim * a + b ]
                                        += dGradientChi1ReferencedChi[ dim * dim * sot_dim * i + dim * sot_dim * k + sot_dim * J + dim * a + b ] * invChiFollow_T[ dim * I + k ];

                            }

                        }

                    }

                }

            }

        }
        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int I = 0; I < dim; I++ ){

                for ( unsigned int J = 0; J < dim; J++ ){

                    for ( unsigned int k = 0; k < dim; k++ ){

                        for ( unsigned int a = 0; a < dim; a++ ){

                            for ( unsigned int b = 0; b < dim; b++ ){

                                    dGradChi1dChi[ dim * dim * sot_dim * i + dim * sot_dim * I + sot_dim * J + dim * a + b ]
                                        += temp_fiot[ dim * dim * sot_dim * i + dim * sot_dim * I + sot_dim * k + dim * a + b ] * invFFollow_T[ dim * J + k ];

                            }

                        }

                    }

                }

            }

        }

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int I = 0; I < dim; I++ ){

                for ( unsigned int J = 0; J < dim; J++ ){

                    for ( unsigned int index = 1; index < num_configs; index++ ){

                        for ( unsigned int a = 0; a < dim; a++ ){
    
                            for ( unsigned int b = 0; b < dim; b++ ){

                                for ( unsigned int k = 0; k < dim; k++ ){

                                        temp_siot1[ dim * dim * ( num_configs - 1 ) * sot_dim * i + dim * ( num_configs - 1 ) * sot_dim * I + ( num_configs - 1 ) * sot_dim * J + dim * dim * ( index - 1 ) + dim * a + b ]
                                            += dGradientChi1ReferencedChin[ dim * dim * ( num_configs - 1 ) * sot_dim * i + dim * ( num_configs - 1 ) * sot_dim * k + ( num_configs - 1 ) * sot_dim * J + dim * dim * ( index - 1 ) + dim * a + b ] * invChiFollow_T[ dim * I + k ];

                                        temp_siot2[ dim * dim * ( num_configs - 1 ) * sot_dim * i + dim * ( num_configs -1 ) * sot_dim * I + ( num_configs - 1 ) * sot_dim * J + dim * dim * ( index - 1 ) + dim * a + b ]
                                            += dGradientChi1ReferencedCn[ dim * dim * ( num_configs - 1 ) * sot_dim * i + dim * ( num_configs - 1 ) * sot_dim * k + ( num_configs - 1 ) * sot_dim * J + dim * dim * ( index - 1 ) + dim * a + b ] * invChiFollow_T[ dim * I + k ];
                                }

                            }

                        }

                    }

                }

            }

        }

        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int I = 0; I < dim; I++ ){

                for ( unsigned int c = 0; c < dim; c++ ){

                    for ( unsigned int J = 0; J < dim; J++ ){

                        for ( unsigned int index = 1; index < num_configs; index++ ){

                            for ( unsigned int a = 0; a < dim; a++ ){

                                for ( unsigned int b = 0; b < dim; b++ ){

                                    dGradChi1dCn[ dim * dim * ( num_configs - 1 ) * sot_dim * i + dim * ( num_configs -1 ) * sot_dim * I + ( num_configs - 1 ) * sot_dim * J + dim * dim * ( index - 1 ) + dim * a + b ]
                                        += temp_tot1[ dim * dim * i + dim * I + c ] * dInvFFollowdFn[ dim * ( num_configs - 1 ) * sot_dim  * c + ( num_configs - 1 ) * sot_dim * J + dim * dim * ( index - 1 ) + dim * a + b ];

                                    dGradChi1dChin[ dim * dim * ( num_configs - 1 ) * sot_dim * i + dim * ( num_configs - 1 ) * sot_dim * I + ( num_configs - 1 ) * sot_dim * J + dim * dim * ( index - 1 ) + dim * a + b ]
                                        += temp_tot2[ dim * dim * i + dim * c + J ] * dInvChiFollowdChin[ dim * ( num_configs - 1 ) * sot_dim * c + ( num_configs - 1 ) * sot_dim * I + dim * dim * ( index - 1 ) + dim * a + b ];
                                    dGradChi1dChin[ dim * dim * ( num_configs - 1 ) * sot_dim * i + dim * ( num_configs - 1 ) * sot_dim * I + ( num_configs - 1 ) * sot_dim * J + dim * dim * ( index - 1 ) + dim * a + b ]
                                        += temp_siot1[ dim * dim * ( num_configs - 1 ) * sot_dim * i + dim * ( num_configs - 1 ) * sot_dim * I + ( num_configs - 1 ) * sot_dim * c + dim * dim * ( index - 1 ) + dim * a + b ] * invFFollow_T[ dim * J + c ];

                                    dGradChi1dCn[ dim * dim * ( num_configs - 1 ) * sot_dim * i + dim * ( num_configs -1 ) * sot_dim * I + ( num_configs - 1 ) * sot_dim * J + dim * dim * ( index - 1 ) + dim * a + b ]
                                        += temp_siot2[ dim * dim * ( num_configs - 1 ) * sot_dim * i + dim * ( num_configs - 1 ) * sot_dim * I + ( num_configs - 1 ) * sot_dim * c + dim * dim * ( index - 1 ) + dim * a + b ] * invFFollow[ dim * c + J ];

                                    for ( unsigned int k = 0; k < dim; k++ ){

                                            temp_seot[ dim * dim * ( num_configs - 1 ) * tot_dim * i + dim * ( num_configs - 1 ) * tot_dim * I + ( num_configs - 1 ) * tot_dim * J + tot_dim * ( index - 1 ) + dim * dim * a + dim * b + k ]
                                                += dGradientChi1ReferencedGradChin[ dim * dim * ( num_configs - 1 ) * tot_dim * i + dim * ( num_configs - 1 ) * tot_dim * c + ( num_configs - 1 ) * tot_dim * J + tot_dim * ( index - 1 ) + dim * dim * a + dim * b + k ] * invChiFollow_T[ dim * I + c ];

                                    }

                                }

                            }

                        }

                    }

                }

            }

        }






        for ( unsigned int i = 0; i < dim; i++ ){

            for ( unsigned int I = 0; I < dim; I++ ){

                for ( unsigned int J = 0; J < dim; J++ ){

                    for ( unsigned int c = 0; c < dim; c++ ){

                        for ( unsigned int index = 1; index < num_configs; index++ ){

                            for ( unsigned int a = 0; a < dim; a++ ){

                                for ( unsigned int b = 0; b < dim; b++ ){

                                    for ( unsigned int k = 0; k < dim; k++ ){

                                        dGradChi1dGradChin[ dim * dim * ( num_configs - 1 ) * tot_dim * i + dim * ( num_configs - 1 ) * tot_dim * I + ( num_configs - 1 ) * tot_dim * J + tot_dim * ( index - 1 ) + dim * dim * a + dim * b + k ]
                                            += temp_seot[ dim * dim * ( num_configs - 1 ) * tot_dim * i + dim * ( num_configs - 1 ) * tot_dim * I + ( num_configs - 1 ) * tot_dim * c + tot_dim * ( index - 1 ) + dim * dim * a + dim * b + k ] * invFFollow_T[ dim * J + c ];

                                    }

                                }

                            }

                        }

                    }

                }

            }

        }

    }

}
