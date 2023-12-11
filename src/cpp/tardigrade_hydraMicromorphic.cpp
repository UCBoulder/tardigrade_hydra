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
                                                  const floatVector &deformationGradient, const floatVector &previousDeformationGradient,
                                                  const floatVector &microDeformation, const floatVector &previousMicroDeformation,
                                                  const floatVector &gradientMicroDeformation, const floatVector &previousGradientMicroDeformation,
                                                  const floatVector &previousStateVariables, const floatVector &parameters,
                                                  const unsigned int numConfigurations, const unsigned int numNonLinearSolveStateVariables,
                                                  const unsigned int dimension, const unsigned int configuration_unknown_count,
                                                  const floatType tolr, const floatType tola, const unsigned int maxIterations,
                                                  const unsigned int maxLSIterations, const floatType lsAlpha ) :
                                                  hydraBase( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables,
                                                             dimension, configuration_unknown_count, tolr, tola, maxIterations, maxLSIterations, lsAlpha ),
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

        decomposeStateVariableVectorMicroConfigurations( );

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

    void hydraBaseMicromorphic::computeGradientMicroConfigurations( const floatVector *data_vector, unsigned int start_index,
                                                                    const floatMatrix &configurations, const floatMatrix &microConfigurations,
                                                                    const floatVector &gradientMicroConfiguration, floatMatrix &gradientMicroConfigurations ){
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

        const unsigned int *dim = getDimension( );

        gradientMicroConfigurations = floatMatrix( *getNumConfigurations( ), floatVector( ( *dim ) * ( *dim ) * ( *dim ), 0 ) );

        for ( unsigned int i = 1; i < *getNumConfigurations( ); i++ ){

            gradientMicroConfigurations[ i ] = floatVector( data_vector->begin( ) + ( i - 1 ) * ( *dim ) * ( *dim ) * ( *dim ) + start_index,
                                                            data_vector->begin( ) + i * ( *dim ) * ( *dim ) * ( *dim ) + start_index );

        }

        calculateFirstConfigurationGradChi( configurations, microConfigurations, gradientMicroConfiguration, gradientMicroConfigurations );

    }

    void hydraBaseMicromorphic::decomposeStateVariableVectorMicroConfigurations( ){
        /*!
         * Decompose the micro-deformation parts of the state variable vector
         */

        unsigned int start_index = ( ( *getNumConfigurations( ) ) - 1 ) * ( *getDimension( ) ) * ( *getDimension( ) );

        floatMatrix microConfigurations;

        floatMatrix inverseMicroConfigurations;

        floatMatrix gradientMicroConfigurations;

        floatMatrix previousMicroConfigurations;

        floatMatrix previousInverseMicroConfigurations;

        floatMatrix previousGradientMicroConfigurations;

        // Compute the micro-configurations

        computeConfigurations( getPreviousStateVariables( ), start_index, *getMicroDeformation( ), microConfigurations, inverseMicroConfigurations, true );

        computeConfigurations( getPreviousStateVariables( ), start_index, *getPreviousMicroDeformation( ), previousMicroConfigurations, previousInverseMicroConfigurations, true );

        start_index += ( ( *getNumConfigurations( ) ) - 1 ) * ( *getDimension( ) ) * ( *getDimension( ) );

        computeGradientMicroConfigurations( getPreviousStateVariables( ), start_index, *getConfigurations( ), microConfigurations,
                                            *getGradientMicroDeformation( ), gradientMicroConfigurations );

        computeGradientMicroConfigurations( getPreviousStateVariables( ), start_index, *getPreviousConfigurations( ), previousMicroConfigurations,
                                            *getPreviousGradientMicroDeformation( ), previousGradientMicroConfigurations );

        // Set the configurations

        setMicroConfigurations( microConfigurations );

        setInverseMicroConfigurations( inverseMicroConfigurations );

        setGradientMicroConfigurations( gradientMicroConfigurations );

        setPreviousMicroConfigurations( previousMicroConfigurations );

        setPreviousInverseMicroConfigurations( previousInverseMicroConfigurations );

        setPreviousGradientMicroConfigurations( previousGradientMicroConfigurations );

    }

    void hydraBaseMicromorphic::setMicroConfigurations( const floatMatrix &microConfigurations ){
        /*!
         * Set the micro-configurations
         * 
         * \param &microConfigurations: The list of micro-configurations in row-major matrix form
         */

        _microConfigurations.second = microConfigurations;

        _microConfigurations.first = true;

        addIterationData( &_microConfigurations );

    }

    void hydraBaseMicromorphic::setInverseMicroConfigurations( const floatMatrix &inverseMicroConfigurations ){
        /*!
         * Set the inverse micro-configurations
         * 
         * \param &inverseMicroConfigurations: The list of inverse micro-configurations in row-major matrix form
         */

        _inverseMicroConfigurations.second = inverseMicroConfigurations;

        _inverseMicroConfigurations.first = true;

        addIterationData( &_inverseMicroConfigurations );

    }

    void hydraBaseMicromorphic::setGradientMicroConfigurations( const floatMatrix &gradientMicroConfigurations ){
        /*!
         * Set the spatial gradients of the micro-configurations w.r.t. their reference configurations
         * 
         * \param &gradientMicroConfigurations: The list of the spatial gradients of the micro-configurations in row-major matrix form
         */

        _gradientMicroConfigurations.second = gradientMicroConfigurations;

        _gradientMicroConfigurations.first = true;

        addIterationData( &_gradientMicroConfigurations );

    }

    void hydraBaseMicromorphic::setPreviousMicroConfigurations( const floatMatrix &previousMicroConfigurations ){
        /*!
         * Set the previous micro-configurations
         * 
         * \param &previousMicroConfigurations: The list of previous micro-configurations in row-major matrix form
         */

        _previousMicroConfigurations.second = previousMicroConfigurations;

        _previousMicroConfigurations.first = true;

    }

    void hydraBaseMicromorphic::setPreviousInverseMicroConfigurations( const floatMatrix &previousInverseMicroConfigurations ){
        /*!
         * Set the previous inverse micro-configurations
         * 
         * \param &previousInverseMicroConfigurations: The list of previous inverse micro-configurations in row-major matrix form
         */

        _previousInverseMicroConfigurations.second = previousInverseMicroConfigurations;

        _previousInverseMicroConfigurations.first = true;

    }

    void hydraBaseMicromorphic::setPreviousGradientMicroConfigurations( const floatMatrix &previousGradientMicroConfigurations ){
        /*!
         * Set the spatial gradients of the previous micro-configurations w.r.t. their reference configurations
         * 
         * \param &previousGradientMicroConfigurations: The list of the spatial gradients of the previous micro-configurations in row-major matrix form
         */

        _previousGradientMicroConfigurations.second = previousGradientMicroConfigurations;

        _previousGradientMicroConfigurations.first = true;

    }

    floatVector hydraBaseMicromorphic::getSubMicroConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex ){
        /*!
         * Get a sub-micro configuration \f$\bf{\chi}^{sc}\f$ defined as
         *
         * \f$ \chi^{sc}_{iI} = \chi^{\text{lowerIndex}}_{i\hat{I}} \chi^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots \chi^{\text{upperIndex-1}}_{\bar{I}I} \f$
         * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
         * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
         *   Note, the configuration indicated by the index is NOT included in the sub-configuration
         */

        return getSubConfiguration( *getMicroConfigurations( ), lowerIndex, upperIndex );

    }

    floatVector hydraBaseMicromorphic::getPrecedingMicroConfiguration( const unsigned int &index ){
        /*!
         * Get the sub-micro configuration preceding but not including the index
         * 
         * \param &index: The index of the configuration immediately following the sub-micro configuration
         */

        return getSubMicroConfiguration( 0, index );

    }

    floatVector hydraBaseMicromorphic::getFollowingMicroConfiguration( const unsigned int &index ){
        /*!
         * Get the sub-micro configuration following but not including the index
         * 
         * \param &index: The index of the current configuration immediately before the sub-micro configuration
         */

        return getSubMicroConfiguration( index + 1, *getNumConfigurations( ) );

    }

    floatVector hydraBaseMicromorphic::getMicroConfiguration( const unsigned int &index ){
        /*!
         * Get the micro configuration indicated by the provided index
         * 
         * \param &index: The index of the current configuration to be extracted
         */

        return getSubMicroConfiguration( index, index + 1 );

    }

    floatVector hydraBaseMicromorphic::getPreviousSubMicroConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex ){
        /*!
         * Get a previous sub-micro configuration \f$\bf{\chi}^{sc}\f$ defined as
         *
         * \f$ \chi^{sc}_{iI} = \chi^{\text{lowerIndex}}_{i\hat{I}} \chi^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots \chi^{\text{upperIndex-1}}_{\bar{I}I} \f$
         * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
         * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
         *   Note, the configuration indicated by the index is NOT included in the sub-configuration
         */

        return getSubConfiguration( *getPreviousMicroConfigurations( ), lowerIndex, upperIndex );

    }

    floatVector hydraBaseMicromorphic::getPreviousPrecedingMicroConfiguration( const unsigned int &index ){
        /*!
         * Get the previous sub-micro configuration preceding but not including the index
         * 
         * \param &index: The index of the configuration immediately following the sub-micro configuration
         */

        return getPreviousSubMicroConfiguration( 0, index );

    }

    floatVector hydraBaseMicromorphic::getPreviousFollowingMicroConfiguration( const unsigned int &index ){
        /*!
         * Get the previous sub-micro configuration following but not including the index
         * 
         * \param &index: The index of the current configuration immediately before the sub-micro configuration
         */

        return getPreviousSubMicroConfiguration( index + 1, *getNumConfigurations( ) );

    }

    floatVector hydraBaseMicromorphic::getPreviousMicroConfiguration( const unsigned int &index ){
        /*!
         * Get the previous micro configuration indicated by the provided index
         * 
         * \param &index: The index of the current configuration to be extracted
         */

        return getPreviousSubMicroConfiguration( index, index + 1 );

    }

    floatMatrix hydraBaseMicromorphic::getSubMicroConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex ){
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

        return getSubConfigurationJacobian( *getMicroConfigurations( ), lowerIndex, upperIndex );

    }

    floatMatrix hydraBaseMicromorphic::getPrecedingMicroConfigurationJacobian( const unsigned int &index ){
        /*!
         * Get the jacobian of the sub-micro configuration preceding but not including the index with respect to the current configurations.
         * 
         * \param &index: The index of the configuration immediately following the sub-configuration
         */

        return getSubMicroConfigurationJacobian( 0, index );

    }

    floatMatrix hydraBaseMicromorphic::getFollowingMicroConfigurationJacobian( const unsigned int &index ){
        /*!
         * Get the jacobian of the sub-micro configuration following but not including the index with respect to the current configurations.
         * 
         * \param &index: The index of the current configuration immediately before the sub-configuration
         */

        return getSubMicroConfigurationJacobian( index + 1, *getNumConfigurations( ) );

    }

    floatMatrix hydraBaseMicromorphic::getPreviousSubMicroConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex ){
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

        return getSubConfigurationJacobian( *getPreviousMicroConfigurations( ), lowerIndex, upperIndex );

    }

    floatMatrix hydraBaseMicromorphic::getPreviousPrecedingMicroConfigurationJacobian( const unsigned int &index ){
        /*!
         * Get the jacobian of the previous sub-micro configuration preceding but not including the index with respect to the current configurations.
         * 
         * \param &index: The index of the configuration immediately following the sub-configuration
         */

        return getPreviousSubMicroConfigurationJacobian( 0, index );

    }

    floatMatrix hydraBaseMicromorphic::getPreviousFollowingMicroConfigurationJacobian( const unsigned int &index ){
        /*!
         * Get the jacobian of the previous sub-micro configuration following but not including the index with respect to the current configurations.
         * 
         * \param &index: The index of the current configuration immediately before the sub-configuration
         */

        return getPreviousSubMicroConfigurationJacobian( index + 1, *getNumConfigurations( ) );

    }

    void hydraBaseMicromorphic::setdChi1dChi( const floatMatrix &dChi1dChi ){
        /*!
         * Set the value of the derivative of the first micro-configuration w.r.t. the total micro configuration
         * 
         * \param &dChi1dChi: The value of the jacobian
         */

        _dChi1dChi.second = dChi1dChi;

        _dChi1dChi.first = true;

        addIterationData( &_dChi1dChi );

    }

    void hydraBaseMicromorphic::setdChi1dChin( const floatMatrix &dChi1dChin ){
        /*!
         * Set the value of the derivative of the first micro-configuration w.r.t. the remaining micro configurationn
         * 
         * \param &dChi1dChin: The value of the jacobian
         */

        _dChi1dChin.second = dChi1dChin;

        _dChi1dChin.first = true;

        addIterationData( &_dChi1dChin );

    }

    void hydraBaseMicromorphic::setPreviousdChi1dChi( const floatMatrix &previousdChi1dChi ){
        /*!
         * Set the value of the derivative of the previous first micro-configuration w.r.t. the total micro configuration
         * 
         * \param &previousdChi1dChi: The value of the jacobian
         */

        _previousdChi1dChi.second = previousdChi1dChi;

        _previousdChi1dChi.first = true;

    }

    void hydraBaseMicromorphic::setPreviousdChi1dChin( const floatMatrix &previousdChi1dChin ){
        /*!
         * Set the value of the derivative of the previous first micro-configuration w.r.t. the remaining micro configurationn
         * 
         * \param &previousdChi1dChin: The value of the jacobian
         */

        _previousdChi1dChin.second = previousdChi1dChin;

        _previousdChi1dChin.first = true;

    }

    void hydraBaseMicromorphic::setdGradChi1dFn( const floatMatrix &dGradChi1dFn ){
        /*!
         * Set the Jacobian of the spatial gradient of the first micro-configuration w.r.t. the configurations after the first.
         * 
         * \param &dGradChi1dFn: The value of the jacobian
         */

        _dGradChi1dFn.second = dGradChi1dFn;

        _dGradChi1dFn.first = true;

        addIterationData( &_dGradChi1dFn );

    }

    void hydraBaseMicromorphic::setdGradChi1dChi( const floatMatrix &dGradChi1dChi ){
        /*!
         * Set the Jacobian of the spatial gradient of the first micro-configuration w.r.t. the total micro-configuration
         * 
         * \param &dGradChi1dChi: The value of the jacobian
         */

        _dGradChi1dChi.second = dGradChi1dChi;

        _dGradChi1dChi.first = true;

        addIterationData( &_dGradChi1dChi );

    }

    void hydraBaseMicromorphic::setdGradChi1dChin( const floatMatrix &dGradChi1dChin ){
        /*!
         * Set the Jacobian of the spatial gradient of the first micro-configuration w.r.t. the sub micro-configurations
         * 
         * \param &dGradChi1dChin: The value of the jacobian
         */

        _dGradChi1dChin.second = dGradChi1dChin;

        _dGradChi1dChin.first = true;

        addIterationData( &_dGradChi1dChin );

    }

    void hydraBaseMicromorphic::setdGradChi1dGradChi( const floatMatrix &dGradChi1dGradChi ){
        /*!
         * Set the Jacobian of the spatial gradient of the first micro-configuration w.r.t. the total spatial gradient of the micro-configuration
         * 
         * \param &dGradChi1dGradChi: The value of the jacobian
         */

        _dGradChi1dGradChi.second = dGradChi1dGradChi;

        _dGradChi1dGradChi.first = true;

        addIterationData( &_dGradChi1dGradChi );

    }

    void hydraBaseMicromorphic::setdGradChi1dGradChin( const floatMatrix &dGradChi1dGradChin ){
        /*!
         * Set the Jacobian of the spatial gradient of the first micro-configuration w.r.t. the sub micro-configurations
         * 
         * \param &dGradChi1dGradChin: The value of the jacobian
         */

        _dGradChi1dGradChin.second = dGradChi1dGradChin;

        _dGradChi1dGradChin.first = true;

        addIterationData( &_dGradChi1dGradChin );

    }

    void hydraBaseMicromorphic::setPreviousdGradChi1dFn( const floatMatrix &previousdGradChi1dFn ){
        /*!
         * Set the Jacobian of the previous spatial gradient of the first micro-configuration w.r.t. the configurations after the first.
         * 
         * \param &previousdGradChi1dFn: The value of the jacobian
         */

        _previousdGradChi1dFn.second = previousdGradChi1dFn;

        _previousdGradChi1dFn.first = true;

    }

    void hydraBaseMicromorphic::setPreviousdGradChi1dChi( const floatMatrix &previousdGradChi1dChi ){
        /*!
         * Set the Jacobian of the previous spatial gradient of the first micro-configuration w.r.t. the total micro-configuration
         * 
         * \param &previousdGradChi1dChi: The value of the jacobian
         */

        _previousdGradChi1dChi.second = previousdGradChi1dChi;

        _previousdGradChi1dChi.first = true;

        addIterationData( &_previousdGradChi1dChi );

    }

    void hydraBaseMicromorphic::setPreviousdGradChi1dChin( const floatMatrix &previousdGradChi1dChin ){
        /*!
         * Set the Jacobian of the previous spatial gradient of the first micro-configuration w.r.t. the sub micro-configurations
         * 
         * \param &previousdGradChi1dChin: The value of the jacobian
         */

        _previousdGradChi1dChin.second = previousdGradChi1dChin;

        _previousdGradChi1dChin.first = true;

    }

    void hydraBaseMicromorphic::setPreviousdGradChi1dGradChi( const floatMatrix &previousdGradChi1dGradChi ){
        /*!
         * Set the Jacobian of the previous spatial gradient of the first micro-configuration w.r.t. the total spatial gradient of the micro-configuration
         * 
         * \param &previousdGradChi1dGradChi: The value of the jacobian
         */

        _previousdGradChi1dGradChi.second = previousdGradChi1dGradChi;

        _previousdGradChi1dGradChi.first = true;

    }

    void hydraBaseMicromorphic::setPreviousdGradChi1dGradChin( const floatMatrix &previousdGradChi1dGradChin ){
        /*!
         * Set the Jacobian of the previous spatial gradient of the first micro-configuration w.r.t. the sub micro-configurations
         * 
         * \param &previousdGradChi1dGradChin: The value of the jacobian
         */

        _previousdGradChi1dGradChin.second = previousdGradChi1dGradChin;

        _previousdGradChi1dGradChin.first = true;

    }

    void hydraBaseMicromorphic::setFirstMicroConfigurationJacobians( ){
        /*!
         * Set the Jacobians of the first micro configuration w.r.t. the total micro configuration and the remaining sub-micro configurations
         */

        floatMatrix dChi1dChi;

        floatMatrix dChi1dChin;

        calculateFirstConfigurationJacobians( *getMicroConfigurations( ), dChi1dChi, dChi1dChin );

        setdChi1dChi( dChi1dChi );

        setdChi1dChin( dChi1dChin );

    }

    void hydraBaseMicromorphic::setPreviousFirstMicroConfigurationJacobians( ){
        /*!
         * Set the Jacobians of the previous first micro configuration w.r.t. the total micro configuration and the remaining sub-micro configurations
         */

        floatMatrix previousdChi1dChi;

        floatMatrix previousdChi1dChin;

        calculateFirstConfigurationJacobians( *getPreviousMicroConfigurations( ), previousdChi1dChi, previousdChi1dChin );

        setPreviousdChi1dChi( previousdChi1dChi );

        setPreviousdChi1dChin( previousdChi1dChin );

    }

    void hydraBaseMicromorphic::setFirstGradientMicroConfigurationJacobians( ){
        /*!
         * Set the Jacobians of the gradient of the first micro configuration w.r.t. the total micro configuration and the remaining sub-micro configurations
         */

        floatMatrix dGradChi1dCn;

        floatMatrix dGradChi1dChi;

        floatMatrix dGradChi1dChin;

        floatMatrix dGradChi1dGradChi;

        floatMatrix dGradChi1dGradChin;

        calculateFirstConfigurationGradChiJacobian( *getConfigurations( ), *getMicroConfigurations( ),
                                                    *getGradientMicroDeformation( ), *getGradientMicroConfigurations( ),
                                                    *getdChi1dChi( ), *getdChi1dChin( ),
                                                    dGradChi1dCn, dGradChi1dChi, dGradChi1dChin, dGradChi1dGradChi, dGradChi1dGradChin );

        setdGradChi1dFn( dGradChi1dCn );

        setdGradChi1dChi( dGradChi1dChi );

        setdGradChi1dChin( dGradChi1dChin );

        setdGradChi1dGradChi( dGradChi1dGradChi );

        setdGradChi1dGradChin( dGradChi1dGradChin );

    }

    void hydraBaseMicromorphic::setPreviousFirstGradientMicroConfigurationJacobians( ){
        /*!
         * Set the Jacobians of the previous gradient of the first micro configuration w.r.t. the total micro configuration and the remaining sub-micro configurations
         */

        floatMatrix previousdGradChi1dCn;

        floatMatrix previousdGradChi1dChi;

        floatMatrix previousdGradChi1dChin;

        floatMatrix previousdGradChi1dGradChi;

        floatMatrix previousdGradChi1dGradChin;

        calculateFirstConfigurationGradChiJacobian( *getPreviousConfigurations( ), *getPreviousMicroConfigurations( ),
                                                    *getPreviousGradientMicroDeformation( ), *getPreviousGradientMicroConfigurations( ),
                                                    *getPreviousdChi1dChi( ), *getPreviousdChi1dChin( ),
                                                    previousdGradChi1dCn, previousdGradChi1dChi, previousdGradChi1dChin, previousdGradChi1dGradChi, previousdGradChi1dGradChin );

        setPreviousdGradChi1dFn( previousdGradChi1dCn );

        setPreviousdGradChi1dChi( previousdGradChi1dChi );

        setPreviousdGradChi1dChin( previousdGradChi1dChin );

        setPreviousdGradChi1dGradChi( previousdGradChi1dGradChi );

        setPreviousdGradChi1dGradChin( previousdGradChi1dGradChin );

    }

    const floatMatrix *hydraBaseMicromorphic::getdChi1dChi( ){
        /*!
         * Get the derivative of the first micro-configuration w.r.t. the total micro-configuration
         */

        if ( !_dChi1dChi.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setFirstMicroConfigurationJacobians( ) );

        }

        return &_dChi1dChi.second;

    }

    const floatMatrix *hydraBaseMicromorphic::getdChi1dChin( ){
        /*!
         * Get the derivative of the first micro-configuration w.r.t. the remaining micro-configurations
         */

        if ( !_dChi1dChin.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setFirstMicroConfigurationJacobians( ) );

        }

        return &_dChi1dChin.second;

    }

    const floatMatrix *hydraBaseMicromorphic::getPreviousdChi1dChi( ){
        /*!
         * Get the derivative of the previous first micro-configuration w.r.t. the total micro-configuration
         */

        if ( !_previousdChi1dChi.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousFirstMicroConfigurationJacobians( ) );

        }

        return &_previousdChi1dChi.second;

    }

    const floatMatrix *hydraBaseMicromorphic::getPreviousdChi1dChin( ){
        /*!
         * Get the derivative of the previous first micro-configuration w.r.t. the remaining micro-configurations
         */

        if ( !_previousdChi1dChin.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousFirstMicroConfigurationJacobians( ) );

        }

        return &_previousdChi1dChin.second;

    }

    const floatMatrix *hydraBaseMicromorphic::getdGradChi1dFn( ){
        /*!
         * Get the Jacobian of the spatial gradient of the first micro-configuration w.r.t. the sub configurations
         */

        if ( !_dGradChi1dFn.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setFirstGradientMicroConfigurationJacobians( ) );

        }

        return &_dGradChi1dFn.second;

    }

    const floatMatrix *hydraBaseMicromorphic::getdGradChi1dChi( ){
        /*!
         * Get the Jacobian of the spatial gradient of the first micro-configuration w.r.t. the micro-configuration
         */

        if ( !_dGradChi1dChi.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setFirstGradientMicroConfigurationJacobians( ) );

        }

        return &_dGradChi1dChi.second;

    }

    const floatMatrix *hydraBaseMicromorphic::getdGradChi1dChin( ){
        /*!
         * Get the Jacobian of the spatial gradient of the first micro-configuration w.r.t. the sub micro-configuration
         */

        if ( !_dGradChi1dChin.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setFirstGradientMicroConfigurationJacobians( ) );

        }

        return &_dGradChi1dChin.second;

    }

    const floatMatrix *hydraBaseMicromorphic::getdGradChi1dGradChi( ){
        /*!
         * Get the Jacobian of the spatial gradient of the first micro-configuration w.r.t. the spatial gradient of the micro-configuration
         */

        if ( !_dGradChi1dGradChi.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setFirstGradientMicroConfigurationJacobians( ) );

        }

        return &_dGradChi1dGradChi.second;

    }

    const floatMatrix *hydraBaseMicromorphic::getdGradChi1dGradChin( ){
        /*!
         * Get the Jacobian of the spatial gradient of the first micro-configuration w.r.t. the spatial gradient of the sub micro-configurations
         */

        if ( !_dGradChi1dGradChin.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setFirstGradientMicroConfigurationJacobians( ) );

        }

        return &_dGradChi1dGradChin.second;

    }

    const floatMatrix *hydraBaseMicromorphic::getPreviousdGradChi1dFn( ){
        /*!
         * Get the Jacobian of the previous spatial gradient of the first micro-configuration w.r.t. the sub configurations
         */

        if ( !_previousdGradChi1dFn.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousFirstGradientMicroConfigurationJacobians( ) );

        }

        return &_previousdGradChi1dFn.second;

    }

    const floatMatrix *hydraBaseMicromorphic::getPreviousdGradChi1dChi( ){
        /*!
         * Get the Jacobian of the previous spatial gradient of the first micro-configuration w.r.t. the micro-configuration
         */

        if ( !_previousdGradChi1dChi.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousFirstGradientMicroConfigurationJacobians( ) );

        }

        return &_previousdGradChi1dChi.second;

    }

    const floatMatrix *hydraBaseMicromorphic::getPreviousdGradChi1dChin( ){
        /*!
         * Get the Jacobian of the previous spatial gradient of the first micro-configuration w.r.t. the sub micro-configuration
         */

        if ( !_previousdGradChi1dChin.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousFirstGradientMicroConfigurationJacobians( ) );

        }

        return &_previousdGradChi1dChin.second;

    }

    const floatMatrix *hydraBaseMicromorphic::getPreviousdGradChi1dGradChi( ){
        /*!
         * Get the Jacobian of the previous spatial gradient of the first micro-configuration w.r.t. the spatial gradient of the micro-configuration
         */

        if ( !_previousdGradChi1dGradChi.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousFirstGradientMicroConfigurationJacobians( ) );

        }

        return &_previousdGradChi1dGradChi.second;

    }

    const floatMatrix *hydraBaseMicromorphic::getPreviousdGradChi1dGradChin( ){
        /*!
         * Get the Jacobian of the previous spatial gradient of the first micro-configuration w.r.t. the spatial gradient of the sub micro-configurations
         */

        if ( !_previousdGradChi1dGradChin.first ){

            TARDIGRADE_ERROR_TOOLS_CATCH( setPreviousFirstGradientMicroConfigurationJacobians( ) );

        }

        return &_previousdGradChi1dGradChin.second;

    }

    void hydraBaseMicromorphic::calculateFirstConfigurationGradChi( const floatMatrix &configurations, const floatMatrix &microConfigurations, const floatVector &gradientMicroConfiguration, floatMatrix &gradientMicroConfigurations ){
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

        const unsigned int *dim = getDimension( );

        // Compute the gradient in the reference configuration
        floatVector gradientChi1Reference = gradientMicroConfiguration; // Initialize to the total gradient in the reference configuration

        for ( unsigned int index = 1; index < *getNumConfigurations( ); index++ ){

            floatVector FFollow = getSubConfiguration( configurations, index + 1, *getNumConfigurations( ) );

            floatVector chiPrecede = getSubConfiguration( microConfigurations, 0, index );

            floatVector chiFollow  = getSubConfiguration( microConfigurations, index + 1, *getNumConfigurations( ) );

            // Add the contribution of the term
            for ( unsigned int i = 0; i < *dim; i++ ){

                for ( unsigned int I = 0; I < *dim; I++ ){

                    for ( unsigned int J = 0; J < *dim; J++ ){

                        for ( unsigned int j = 0; j < *dim; j++ ){

                            for ( unsigned int k = 0; k < *dim; k++ ){

                                for ( unsigned int l = 0; l < *dim; l++ ){

                                    gradientChi1Reference[ ( *dim ) * ( *dim ) * i + ( *dim ) * I + J ]
                                        -= chiPrecede[ ( *dim ) * i + j ] * chiFollow[ ( *dim ) * k + I ] * FFollow[ ( *dim ) * l + J ]
                                         * gradientMicroConfigurations[ index ][ ( *dim ) * ( *dim ) * j + ( *dim ) * k + l ];

                                }

                            }

                        }

                    }

                }

            }

        }

        // Map the gradient of the micro-configuration to the reference of the first configuration
        floatVector invChiFollow = tardigradeVectorTools::inverse( getSubConfiguration( microConfigurations, 1, *getNumConfigurations( ) ), *dim, *dim );

        floatVector invFFollow = tardigradeVectorTools::inverse( getSubConfiguration( configurations, 1, *getNumConfigurations( ) ), *dim, *dim );

        gradientMicroConfigurations[ 0 ] = floatVector( ( *dim ) * ( *dim ) * ( *dim ), 0 );

        for ( unsigned int i = 0; i < *dim; i++ ){

            for ( unsigned int I = 0; I < *dim; I++ ){

                for ( unsigned int J = 0; J < *dim; J++ ){

                    for ( unsigned int a = 0; a < *dim; a++ ){

                        for ( unsigned int b = 0; b < *dim; b++ ){

                            gradientMicroConfigurations[ 0 ][ ( *dim ) * ( *dim ) * i + ( *dim ) * I + J ] += gradientChi1Reference[ ( *dim ) * ( *dim ) * i + ( *dim ) * a + b ] * invChiFollow[ ( *dim ) * a + I ] * invFFollow[ ( *dim ) * b + J ];
   

                        }

                    }

                }

            }

        }

    }

    void hydraBaseMicromorphic::calculateFirstConfigurationGradChiJacobian( const floatMatrix &configurations, const floatMatrix &microConfigurations,
                                                                            const floatVector &gradientMicroConfiguration, const floatMatrix &gradientMicroConfigurations,
                                                                            const floatMatrix &dChi1dChi, const floatMatrix &dChi1dChin,
                                                                            floatMatrix &dGradChi1dCn,
                                                                            floatMatrix &dGradChi1dChi, floatMatrix &dGradChi1dChin,
                                                                            floatMatrix &dGradChi1dGradChi, floatMatrix &dGradChi1dGradChin ){
        /*!
         * Calculate the value of the jacobian of the gradient of the first micro-configuration given all of the configurations, the micro-configurations,
         * the spatial gradient of the micro deformation in the reference configuration, and the gradients of the micro-configurations
         * other than the first in their own reference configurations.
         * 
         * \param &configurations: The configuration matrix
         * \param &microConfigurations: The micro-configuration matrix
         * \param &gradientMicroConfiguration: The gradient of the micro-deformation in the reference configuration
         * \param &gradientMicroConfigurations: The gradient of the micro-deformations in their reference configurations
         * \param &dGradChi1dCn: The Jacobian of the gradient of the first micro-configuration w.r.t. the remaining configurations
         * \param &dGradChi1dChi: The Jacobian of the gradient of the first micro-configuration w.r.t. the total micro-configuration
         * \param &dGradChi1dChin: The Jacobian of the gradient of the first micro-configuration w.r.t. the remaining micro-configurations
         * \param &dGradChi1dGradChi: The Jacobian of the gradient of the first micro-configuration w.r.t. the gradient of the total micro-configuration
         * \param &dGradChi1dGradChin: The Jacobian of the gradient of the first micro-configuration w.r.t. the gradient of the remaining sub micro-configurations
         */

        const unsigned int *dim = getDimension( );

        floatVector eye( ( *dim ) * ( *dim ), 0 );
        tardigradeVectorTools::eye( eye );

        // Compute the gradient in the reference configuration
        floatVector gradientChi1Reference = gradientMicroConfiguration; // Initialize to the total gradient in the reference configuration

        floatMatrix dGradientChi1ReferencedCn( gradientMicroConfiguration.size( ), floatVector( ( ( *getNumConfigurations( ) ) - 1 ) * configurations[ 0 ].size( ), 0 ) );

        floatMatrix dGradientChi1ReferencedChi( gradientMicroConfiguration.size( ), floatVector( microConfigurations[ 0 ].size( ), 0 ) );

        floatMatrix dGradientChi1ReferencedChin( gradientMicroConfiguration.size( ), floatVector( ( ( *getNumConfigurations( ) ) - 1 ) * microConfigurations[ 0 ].size( ), 0 ) );

        floatMatrix dGradientChi1ReferencedGradChin( gradientMicroConfiguration.size( ), floatVector( ( ( *getNumConfigurations( ) ) - 1 ) * gradientMicroConfigurations[ 0 ].size( ), 0 ) );

        for ( unsigned int index = 1; index < *getNumConfigurations( ); index++ ){

            floatVector FFollow = getSubConfiguration( configurations, index + 1, *getNumConfigurations( ) );

            floatVector chiPrecede = getSubConfiguration( microConfigurations, 0, index );

            floatVector chiFollow  = getSubConfiguration( microConfigurations, index + 1, *getNumConfigurations( ) );

            // Set the Jacobians of the mapping terms
            floatMatrix dFFollowdCs = getSubConfigurationJacobian( configurations, index + 1, *getNumConfigurations( ) );

            floatMatrix dChiPrecededChis = getSubConfigurationJacobian( microConfigurations, 0, index );

            floatMatrix dChiPrecededChi( microConfigurations[ 0 ].size( ), floatVector( microConfigurations[ 0 ].size( ), 0 ) );

            floatMatrix dChiPrecededChin( microConfigurations[ 0 ].size( ), floatVector( ( ( *getNumConfigurations( ) ) - 1 ) * microConfigurations[ 0 ].size( ), 0 ) );

            for ( unsigned int i = 0; i < ( *dim ) * ( *dim ); i++ ){

                for ( unsigned int j = 0; j < ( *dim ) * ( *dim ); j++ ){

                    for ( unsigned int k = 0; k < ( *dim ) * ( *dim ); k++ ){

                        dChiPrecededChi[ i ][ j ] += dChiPrecededChis[ i ][ k ] * dChi1dChi[ k ][ j ];

                    }

                }
                for ( unsigned int j = 0; j < ( *dim ) * ( *dim ) * ( ( *getNumConfigurations( ) ) - 1 ); j++ ){

                    dChiPrecededChin[ i ][ j ] += dChiPrecededChis[ i ][ j + ( *dim ) * ( *dim ) ];

                    for ( unsigned int k = 0; k < ( *dim ) * ( *dim ); k++ ){

                        dChiPrecededChin[ i ][ j ] += dChiPrecededChis[ i ][ k ] * dChi1dChin[ k ][ j ];

                    }

                }

            }

            floatMatrix dChiFollowdChis = getSubConfigurationJacobian( microConfigurations, index + 1, *getNumConfigurations( ) );

            // Add the contribution of the term
            for ( unsigned int i = 0; i < *dim; i++ ){

                for ( unsigned int I = 0; I < *dim; I++ ){

                    for ( unsigned int J = 0; J < *dim; J++ ){

                        for ( unsigned int j = 0; j < *dim; j++ ){

                            for ( unsigned int k = 0; k < *dim; k++ ){

                                for ( unsigned int l = 0; l < *dim; l++ ){

                                    gradientChi1Reference[ ( *dim ) * ( *dim ) * i + ( *dim ) * I + J ]
                                        -= chiPrecede[ ( *dim ) * i + j ] * chiFollow[ ( *dim ) * k + I ] * FFollow[ ( *dim ) * l + J ]
                                         * gradientMicroConfigurations[ index ][ ( *dim ) * ( *dim ) * j + ( *dim ) * k + l ];

                                    dGradientChi1ReferencedGradChin[ ( *dim ) * ( *dim ) * i + ( *dim ) * I + J ][ gradientMicroConfiguration.size( ) * ( index - 1 ) + ( *dim ) * ( *dim ) * j + ( *dim ) * k + l ]
                                        -= chiPrecede[ ( *dim ) * i + j ] * chiFollow[ ( *dim ) * k + I ] * FFollow[ ( *dim ) * l + J ];

                                    for ( unsigned int A = 0; A < dFFollowdCs[ 0 ].size( ) - FFollow.size( ); A++ ){

                                        dGradientChi1ReferencedCn[ ( *dim ) * ( *dim ) * i + ( *dim ) * I + J ][ A ]
                                             -= chiPrecede[ ( *dim ) * i + j ] * chiFollow[ ( *dim ) * k + I ] * dFFollowdCs[ ( *dim ) * l + J ][ FFollow.size( ) + A ]
                                             * gradientMicroConfigurations[ index ][ ( *dim ) * ( *dim ) * j + ( *dim ) * k + l ];

                                    }

                                    for ( unsigned int A = 0; A < dChiPrecededChi[ 0 ].size( ); A++ ){

                                        dGradientChi1ReferencedChi[ ( *dim ) * ( *dim ) * i + ( *dim ) * I + J ][ A ]
                                            -= dChiPrecededChi[ ( *dim ) * i + j ][ A ] * chiFollow[ ( *dim ) * k + I ] * FFollow[ ( *dim ) * l + J ]
                                             * gradientMicroConfigurations[ index ][ ( *dim ) * ( *dim ) * j + ( *dim ) * k + l ];

                                    }

                                    for ( unsigned int A = 0; A < ( ( *getNumConfigurations( ) ) - 1 ) * ( *dim ) * ( *dim ); A++ ){

                                        dGradientChi1ReferencedChin[ ( *dim ) * ( *dim ) * i + ( *dim ) * I + J ][ A ]
                                            -= dChiPrecededChin[ ( *dim ) * i + j ][ A ] * chiFollow[ ( *dim ) * k + I ] * FFollow[ ( *dim ) * l + J ]
                                             * gradientMicroConfigurations[ index ][ ( *dim ) * ( *dim ) * j + ( *dim ) * k + l ]
                                             + chiPrecede[ ( *dim ) * i + j ] * dChiFollowdChis[ ( *dim ) * k + I ][ A + ( *dim ) * ( *dim ) ] * FFollow[ ( *dim ) * l + J ]
                                             * gradientMicroConfigurations[ index ][ ( *dim ) * ( *dim ) * j + ( *dim ) * k + l ];

                                    }

                                }

                            }

                        }

                    }

                }

            }

        }

        // Map the gradient of the micro-configuration to the reference of the first configuration
        floatVector ChiFollow = getSubConfiguration( microConfigurations, 1, *getNumConfigurations( ) );

        floatVector FFollow = getSubConfiguration( configurations, 1, *getNumConfigurations( ) );

        floatVector invChiFollow = tardigradeVectorTools::inverse( ChiFollow, *dim, *dim );

        floatVector invFFollow = tardigradeVectorTools::inverse( FFollow, *dim, *dim );

        floatMatrix dChiFollowdChis = getSubConfigurationJacobian( microConfigurations, 1, *getNumConfigurations( ) );

        floatMatrix dFFollowdFs = getSubConfigurationJacobian( configurations, 1, *getNumConfigurations( ) );

        floatMatrix dInvChiFollowdChiFollow( invChiFollow.size( ), floatVector( ChiFollow.size( ), 0 ) );

        floatMatrix dInvFFollowdFFollow( invFFollow.size( ), floatVector( FFollow.size( ), 0 ) );

        for ( unsigned int I = 0; I < ( *dim ); I++ ){

            for ( unsigned int i = 0; i < ( *dim ); i++ ){

                for ( unsigned int j = 0; j < ( *dim ); j++ ){

                    for ( unsigned int J = 0; J < ( *dim ); J++ ){

                        dInvChiFollowdChiFollow[ ( *dim ) * I + j ][ ( *dim ) * i + J ] -= invChiFollow[ ( *dim ) * I + i ] * invChiFollow[ ( *dim ) * J + j ];

                        dInvFFollowdFFollow[ ( *dim ) * I + j ][ ( *dim ) * i + J ] -= invFFollow[ ( *dim ) * I + i ] * invFFollow[ ( *dim ) * J + j ];

                    }

                }

            }

        }

        floatMatrix dInvChiFollowdChin( ( *dim ) * ( *dim ), floatVector( ( *dim ) * ( *dim ) * ( ( *getNumConfigurations( ) ) - 1 ), 0 ) );

        floatMatrix dInvFFollowdFn( ( *dim ) * ( *dim ), floatVector( ( *dim ) * ( *dim ) * ( ( *getNumConfigurations( ) ) - 1 ), 0 ) );

        for ( unsigned int i = 0; i < ( *dim ) * ( *dim ); i++ ){

            for ( unsigned int j = 0; j < ( ( *getNumConfigurations( ) ) - 1 ) * ( *dim ) * ( *dim ); j++ ){

                for ( unsigned int k = 0; k < ( *dim ) * ( *dim ); k++ ){

                    dInvChiFollowdChin[ i ][ j ] += dInvChiFollowdChiFollow[ i ][ k ] * dChiFollowdChis[ k ][ j + ( *dim ) * ( *dim ) ];

                    dInvFFollowdFn[ i ][ j ] += dInvFFollowdFFollow[ i ][ k ] * dFFollowdFs[ k ][ j + ( *dim ) * ( *dim ) ];

                }

            }

        }

        dGradChi1dChi = floatMatrix( ( *dim ) * ( *dim ) * ( *dim ), floatVector( ( *dim ) * ( *dim ), 0 ) );

        dGradChi1dChin = floatMatrix( ( *dim ) * ( *dim ) * ( *dim ), floatVector( ( ( *getNumConfigurations( ) ) - 1 ) * ( *dim ) * ( *dim ), 0 ) );

        dGradChi1dGradChi = floatMatrix( ( *dim ) * ( *dim ) * ( *dim ), floatVector( ( *dim ) * ( *dim ) * ( *dim ), 0 ) );

        dGradChi1dCn = floatMatrix( ( *dim ) * ( *dim ) * ( *dim ), floatVector( ( ( *getNumConfigurations( ) ) - 1 ) * ( *dim ) * ( *dim ), 0 ) );

        dGradChi1dGradChin = floatMatrix( ( *dim ) * ( *dim ) * ( *dim ), floatVector( ( ( *getNumConfigurations( ) ) - 1 ) * gradientMicroConfigurations[ 0 ].size( ), 0 ) );

        for ( unsigned int i = 0; i < *dim; i++ ){

            for ( unsigned int I = 0; I < *dim; I++ ){

                for ( unsigned int J = 0; J < *dim; J++ ){

                    for ( unsigned int a = 0; a < *dim; a++ ){

                        for ( unsigned int b = 0; b < *dim; b++ ){

                            for ( unsigned int k = 0; k < *dim; k++ ){

                                for ( unsigned int l = 0; l < *dim; l++ ){

                                    dGradChi1dChi[ ( *dim ) * ( *dim ) * i + ( *dim ) * I + J ][ ( *dim ) * a + b ]
                                        += dGradientChi1ReferencedChi[ ( *dim ) * ( *dim ) * i + ( *dim ) * k + l ][ ( *dim ) * a + b ] * invChiFollow[ ( *dim ) * k + I ] * invFFollow[ ( *dim ) * l + J ];

                                    for ( unsigned int index = 1; index < *getNumConfigurations( ); index++ ){

                                        dGradChi1dChin[ ( *dim ) * ( *dim ) * i + ( *dim ) * I + J ][ ( *dim ) * ( *dim ) * ( index - 1 ) + ( *dim ) * a + b ]
                                            += dGradientChi1ReferencedChin[ ( *dim ) * ( *dim ) * i + ( *dim ) * k + l ][ ( *dim ) * ( *dim ) * ( index - 1 ) + ( *dim ) * a + b ] * invChiFollow[ ( *dim ) * k + I ] * invFFollow[ ( *dim ) * l + J ]
                                             + gradientChi1Reference[ ( *dim ) * ( *dim ) * i + ( *dim ) * k + l ] * dInvChiFollowdChin[ ( *dim ) * k + I ][ ( *dim ) * ( *dim ) * ( index - 1 ) + ( *dim ) * a + b ] * invFFollow[ ( *dim ) * l + J ];

                                        dGradChi1dCn[ ( *dim ) * ( *dim ) * i + ( *dim ) * I + J ][ ( *dim ) * ( *dim ) * ( index - 1 ) + ( *dim ) * a + b ]
                                            += dGradientChi1ReferencedCn[ ( *dim ) * ( *dim ) * i + ( *dim ) * k + l ][ ( *dim ) * ( *dim ) * ( index - 1 ) + ( *dim ) * a + b ] * invChiFollow[ ( *dim ) * k + I ] * invFFollow[ ( *dim ) * l + J ]
                                             + gradientChi1Reference[ ( *dim ) * ( *dim ) * i + ( *dim ) * k + l ] * invChiFollow[ ( *dim ) * k + I ] * dInvFFollowdFn[ ( *dim ) * l + J ][ ( *dim ) * ( *dim ) * ( index - 1 ) + ( *dim ) * a + b ];

                                    }

                                }

                            }

                            for ( unsigned int c = 0; c < *dim; c++ ){

                                dGradChi1dGradChi[ ( *dim ) * ( *dim ) * i + ( *dim ) * I + J ][ ( *dim ) * ( *dim ) * a + ( *dim ) * b + c ]
                                    += eye[ ( *dim ) * i + a ] * invChiFollow[ ( *dim ) * b + I ] * invFFollow[ ( *dim ) * c + J ];

                                for ( unsigned int index = 1; index < *getNumConfigurations( ); index++ ){

                                    for ( unsigned int k = 0; k < *dim; k++ ){

                                        for ( unsigned int l = 0; l < *dim; l++ ){

                                            dGradChi1dGradChin[ ( *dim ) * ( *dim ) * i + ( *dim ) * I + J ][ gradientMicroConfiguration.size( ) * ( index - 1 ) + ( * dim ) * ( *dim ) * a + ( *dim ) * b + c ]
                                                += dGradientChi1ReferencedGradChin[ ( *dim ) * ( *dim ) * i + ( *dim ) * k + l ][ gradientMicroConfiguration.size( ) * ( index - 1 ) + ( *dim ) * ( *dim ) * a + ( *dim ) * b + c ] * invChiFollow[ ( *dim ) * k + I ] * invFFollow[ ( *dim ) * l + J ];

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

}
