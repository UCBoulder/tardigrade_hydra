/**
  ******************************************************************************
  * \file hydra.cpp
  ******************************************************************************
  * A C++ library for printing messages to stdout. Used as a stub repo example.
  ******************************************************************************
  */

#include<hydra.h>

namespace hydra{

    //Define hydra global constants in a place that Doxygen can pick up for documentation
    /** \brief Define the expected number of tensor spatial dimensions for the Abaqus interface. */
    const int spatialDimensions = 3;

    /** \brief Define required number of Abaqus material constants for the Abaqus interface. */
    const int nStateVariables = 2;

    /** \brief Define required number of Abaqus material constants for the Abaqus interface. */
    const int nMaterialParameters = 2;

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

        const unsigned int* dim = getDimension( );

        const unsigned int* nConfig = getNumConfigurations( );

        const unsigned int* nNLISV  = getNumNonLinearSolveStateVariables( );

        // Extract the previous configurations
        floatVector eye( ( *dim ) * ( *dim ) );
        vectorTools::eye( eye );

        if ( getPreviousStateVariables( )->size( ) < ( ( *nConfig ) * ( *dim ) * ( *dim ) + ( *nNLISV ) ) ){

            std::string message = "The number of state variables is less than required for the configurations and\n";
            message            += "non-linear state variables\n";
            message            += "  # previousStateVariables          : " + std::to_string( getPreviousStateVariables( )->size( ) ) + "\n";
            message            += "  # ( configurations - 1 ) * dim**2 : " + std::to_string( ( ( *nConfig ) - 1 ) * ( *dim ) * ( *dim ) ) + "\n";
            message            += "  # non-linear solve ISVs           : " + std::to_string( ( *nNLISV ) ) + "\n";
            message            += "  # minimum required ISVs           : " + std::to_string( ( *nConfig ) * ( *dim ) * ( *dim ) + ( *nNLISV ) );

            ERROR_TOOLS_CATCH( throw std::runtime_error( message ) );

        }

        _previousConfigurations.second = floatMatrix( *nConfig, floatVector( ( *dim ) * ( *dim ), 0 ) );

        _configurations.second = floatMatrix( *nConfig, floatVector( ( *dim ) * ( *dim ), 0 ) );

        // Initialize the first configuration with the total deformation gradient
        _configurations.second[ 0 ] = *getDeformationGradient( );

        _previousConfigurations.second[ 0 ] = *getPreviousDeformationGradient( );

        for ( unsigned int i = 1; i < *nConfig; i++ ){

            // Set the current configuration as being equal to the previous
            _configurations.second[ i ] = floatVector( getPreviousStateVariables( )->begin( ) + ( i - 1 ) * ( *dim ) * ( *dim ),
                                                       getPreviousStateVariables( )->begin( ) + i * ( *dim ) * ( *dim ) ) + eye;

            // Compute the inverse of the current configuration and store it
            _inverseConfigurations.second[ i ] = vectorTools::inverse( _configurations.second[ i ], ( *dim ), ( *dim ) );

            // Set the previous configuration
            _previousConfigurations.second[ i ] = _configurations.second[ i ];

            // Set the previous inverse configuration
            _previousInverseConfigurations.second[ i ] = _inverseConfigurations.second[ i ];

            // Add contribution of deformation gradient to the first configuration
            _configurations.second[ 0 ] = vectorTools::matrixMultiply( _configurations.second[ 0 ], _inverseConfigurations.second[ i ],
                                                                       ( *dim ), ( *dim ), ( *dim ), ( *dim ) );

            // Add the contribution of the deformation gradient to the previous configuration
            _previousConfigurations.second[ 0 ] = vectorTools::matrixMultiply( _previousConfigurations.second[ 0 ], _previousInverseConfigurations.second[ i ],
                                                                               ( *dim ), ( *dim ), ( *dim ), ( *dim ) );

        }

        _inverseConfigurations.second[ 0 ] = vectorTools::inverse( _configurations.second[ 0 ], ( *dim ), ( *dim ) );

        _previousInverseConfigurations.second[ 0 ] = _inverseConfigurations.second[ 0 ];

        // Extract the remaining state variables required for the non-linear solve
        _nonLinearSolveStateVariables.second = floatVector( getPreviousStateVariables( )->begin( ) + ( *nConfig ) * ( *dim ) * ( *dim ),
                                                            getPreviousStateVariables( )->begin( ) + ( *nConfig ) * ( *dim ) * ( *dim ) + *nNLISV );

        _previousNonLinearSolveStateVariables.second = _nonLinearSolveStateVariables.second;

        // Extract the additional state variables
        _additionalStateVariables.second = floatVector( getPreviousStateVariables( )->begin( ) + ( *nConfig ) * ( *dim ) * ( *dim ) + *nNLISV,
                                                        getPreviousStateVariables( )->end( ) );

        _previousAdditionalStateVariables.second = _additionalStateVariables.second;

        _configurations.first = true;

        _previousConfigurations.first = true;

        _inverseConfigurations.first = true;

        _previousInverseConfigurations.first = true;

        _nonLinearSolveStateVariables.first = true;

        _previousNonLinearSolveStateVariables.first = true;

        _additionalStateVariables.first = true;

        _previousAdditionalStateVariables.first = true;

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
        const std::string cmname( abaqusTools::FtoCString( 80, CMNAME ) );
        const std::vector< double > props( PROPS, PROPS + NPROPS );
        const std::vector< double > coords( COORDS, COORDS + spatialDimensions );
        const std::vector< int > jstep( JSTEP, JSTEP + 4 );
        //Fortran two-dimensional arrays require careful column to row major conversions to c++ types
        std::vector< std::vector< double > > ddsdde = abaqusTools::columnToRowMajor( DDSDDE, NTENS, NTENS );
        const std::vector< std::vector< double > > drot = abaqusTools::columnToRowMajor( DROT, spatialDimensions, spatialDimensions );
        const std::vector< std::vector< double > > dfgrd0 = abaqusTools::columnToRowMajor( DFGRD0, spatialDimensions, spatialDimensions );
        const std::vector< std::vector< double > > dfgrd1 = abaqusTools::columnToRowMajor( DFGRD1, spatialDimensions, spatialDimensions );

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
        abaqusTools::rowToColumnMajor( STRESS, stress, 1, NTENS );
        abaqusTools::rowToColumnMajor( DDSDDT, ddsddt, 1, NTENS );
        abaqusTools::rowToColumnMajor( DRPLDE, drplde, 1, NTENS );
        abaqusTools::rowToColumnMajor( STATEV, statev, 1, NSTATV );
        //Arrays require vector of vector to column major conversion
        abaqusTools::rowToColumnMajor( DDSDDE, ddsdde, NTENS, NTENS );

    }

}
