/**
  * \file test_tardigrade-hydra.cpp
  *
  * Tests for tardigrade-hydra
  */

#include<tardigrade-hydra.h>
#include<tardigrade-hydraLinearElasticity.h>
#include<constitutive_tools.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_tardigrade-hydra
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef errorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer)
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

namespace tardigradeHydra{

    namespace unit_test{

        class hydraBaseTester{

            public:

                static void checkTime( hydraBase &hydra ){
    
                    BOOST_CHECK( &hydra._time == hydra.getTime( ) );
    
                }
    
                static void checkDeltaTime( hydraBase &hydra ){
    
                    BOOST_CHECK( &hydra._deltaTime == hydra.getDeltaTime( ) );
    
                }
    
                static void checkTemperature( hydraBase &hydra ){
    
                    BOOST_CHECK( &hydra._temperature == hydra.getTemperature( ) );
    
                }
    
                static void checkPreviousTemperature( hydraBase &hydra ){
    
                    BOOST_CHECK( &hydra._previousTemperature == hydra.getPreviousTemperature( ) );
    
                }
    
                static void checkDeformationGradient( hydraBase &hydra ){
    
                    BOOST_CHECK( &hydra._deformationGradient == hydra.getDeformationGradient( ) );
    
                }
    
                static void checkPreviousDeformationGradient( hydraBase &hydra ){
    
                    BOOST_CHECK( &hydra._previousDeformationGradient == hydra.getPreviousDeformationGradient( ) );
    
                }
    
                static void checkPreviousStateVariables( hydraBase &hydra ){
    
                    BOOST_CHECK( &hydra._previousStateVariables == hydra.getPreviousStateVariables( ) );
    
                }
    
                static void checkParameters( hydraBase &hydra ){
    
                    BOOST_CHECK( &hydra._parameters == hydra.getParameters( ) );
    
                }

                static void checkNumConfigurations( hydraBase &hydra ){
    
                    BOOST_CHECK( &hydra._numConfigurations == hydra.getNumConfigurations( ) );
    
                }

                static void checkNumNonLinearSolveStateVariables( hydraBase &hydra ){
    
                    BOOST_CHECK( &hydra._numNonLinearSolveStateVariables == hydra.getNumNonLinearSolveStateVariables( ) );
    
                }

                static void checkDimension( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._dimension == hydra.getDimension( ) );

                }

                static void checkConfigurations( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._configurations.second == hydra.getConfigurations( ) );

                }

                static void checkPreviousConfigurations( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._previousConfigurations.second == hydra.getPreviousConfigurations( ) );

                }

                static void checkInverseConfigurations( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._inverseConfigurations.second == hydra.getInverseConfigurations( ) );

                }

                static void checkPreviousInverseConfigurations( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._previousInverseConfigurations.second == hydra.getPreviousInverseConfigurations( ) );

                }

                static void checkNonLinearSolveStateVariables( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._nonLinearSolveStateVariables.second == hydra.getNonLinearSolveStateVariables( ) );

                }

                static void checkPreviousNonLinearSolveStateVariables( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._previousNonLinearSolveStateVariables.second == hydra.getPreviousNonLinearSolveStateVariables( ) );

                }

                static void checkAdditionalStateVariables( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._additionalStateVariables.second == hydra.getAdditionalStateVariables( ) );

                }

                static void checkPreviousAdditionalStateVariables( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._previousAdditionalStateVariables.second == hydra.getPreviousAdditionalStateVariables( ) );

                }

                static void decomposeStateVariableVector( hydraBase &hydra ){

                    ERROR_TOOLS_CATCH( hydra.decomposeStateVariableVector( ) );

                }

                static void decomposeUnknownVector( hydraBase &hydra ){

                    ERROR_TOOLS_CATCH( hydra.decomposeUnknownVector( ) );

                }

                static void checkdF1dF( hydraBase &hydra ){

                    BOOST_CHECK( std::find( hydra._iterationData.begin( ), hydra._iterationData.end( ), &hydra._dF1dF ) != hydra._iterationData.end( ) );

                }

                static void checkdF1dFn( hydraBase &hydra ){

                    BOOST_CHECK( std::find( hydra._iterationData.begin( ), hydra._iterationData.end( ), &hydra._dF1dFn ) != hydra._iterationData.end( ) );

                }

                static void checkRelativeTolerance( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._tolr == hydra.getRelativeTolerance( ) );

                }

                static void checkAbsoluteTolerance( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._tola == hydra.getAbsoluteTolerance( ) );

                }

                static void checkLSAlpha( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._lsAlpha == hydra.getLSAlpha( ) );

                }

                static void formNonLinearProblem( hydraBase &hydra ){

                    BOOST_CHECK_NO_THROW( hydra.formNonLinearProblem( ) );

                }

                static void initializeUnknownVector( hydraBase &hydra ){

                    BOOST_CHECK_NO_THROW( hydra.initializeUnknownVector( ) );

                }

                static void set_residual( hydraBase &hydra, const floatVector &value ){

                    hydra._residual.second = value;
                    hydra._residual.first = true;

                    hydra.addIterationData( &hydra._residual );

                }

                static void set_unknownVector( hydraBase &hydra, const floatVector &value ){

                    hydra._X.second = value;
                    hydra._X.first = true;

                }

                static void set_flatJacobian( hydraBase &hydra, const floatVector &value ){

                    hydra._jacobian.second = value;
                    hydra._jacobian.first = true;

                    hydra.addIterationData( &hydra._jacobian );
                }

                static unsigned int get_iteration( hydraBase &hydra ){

                    return hydra._iteration;

                }

                static unsigned int get_LSIteration( hydraBase &hydra ){

                    return hydra._LSIteration;

                }

                static void solveNonLinearProblem( hydraBase &hydra ){

                    hydra.solveNonLinearProblem( );

                }

                static void resetIterationData( hydraBase &hydra ){

                    hydra.resetIterationData( );

                }

                static void updateUnknownVector( hydraBase &hydra, const floatVector &value ){

                    hydra.updateUnknownVector( value );

                }

        };

    }

}

BOOST_AUTO_TEST_CASE( testSayHello ){
    /*!
     * Test message printed to stdout in sayHello function
     */

    //Setup redirect variables for stdout
    std::stringbuf buffer;
    cout_redirect rd(&buffer);
    boost::test_tools::output_test_stream result;

    //Initialize test variables
    std::string message;
    std::string answer;
    errorOut error = NULL;

    cout_redirect guard( result.rdbuf() );

    //Check normal operation
    message = "World!";
    answer = "Hello World!\n";
    error = tardigradeHydra::sayHello(message);
    BOOST_CHECK( ! error );
    BOOST_CHECK( result.is_equal( answer ) );

    //Reset error code between tests
    error = NULL;

    //Check for "George" error
    message = "George";
    error = tardigradeHydra::sayHello(message);
    BOOST_CHECK( error );

}

BOOST_AUTO_TEST_CASE( testAbaqusInterface ){
    /*!
     * Test the tardigrade-hydra abaqus interface
     */

    double double_scalar = 0.0;
    int int_scalar = 0;

    //Create nominally correct variable holders that match expected Abaqus Fortran interface
    //TODO: fill out nominally correct variable shape and values
    //Strings
    char CMNAME[ ] = "tardigrade-hydra";
    //Scalar integers
    int NDI = 3;
    int NSHR = 3;
    int NTENS = 6;
    int NSTATV = 2;
    int NPROPS = 2;
    int NOEL = int_scalar;
    int NPT = int_scalar;
    int LAYER = int_scalar;
    int KSPT = int_scalar;
    int KINC = int_scalar;
    //Scalar doubles
    double SSE = double_scalar;
    double SPD = double_scalar;
    double SCD = double_scalar;
    double RPL = double_scalar;
    double DRPLDT = double_scalar;
    double DTIME = double_scalar;
    double TEMP = double_scalar;
    double DTEMP = double_scalar;
    double PNEWDT = double_scalar;
    double CELENT = double_scalar;
    //Fortan int column major arrays
    std::vector< int > jstep( 4 );
    int* JSTEP  = jstep.data( );
    //Fortan double column major arrays
    std::vector< double > stress( NTENS );
    double* STRESS = stress.data( );
    std::vector< double > statev( NSTATV );
    double* STATEV = statev.data( );
    std::vector< double > ddsdde( NTENS * NTENS );
    double* DDSDDE = ddsdde.data( );
    std::vector< double > ddsddt( NTENS );
    double* DDSDDT = ddsddt.data( );
    std::vector< double > drplde( NTENS );
    double* DRPLDE = drplde.data( );
    std::vector< double > strain( NTENS );
    double* STRAN  = strain.data( );
    std::vector< double > dstrain( NTENS );
    double* DSTRAN = dstrain.data( );
    std::vector< double > time( 2 );
    double* TIME   = time.data( );
    std::vector< double > predef( 1 );
    double* PREDEF = predef.data( );
    std::vector< double > dpred( 1 );
    double* DPRED  = dpred.data( );
    std::vector< double > props( NPROPS );
    double* PROPS  = props.data( );
    //TODO: figure out how to use tardigradeHydra::spatialDimensions here
    std::vector< double > coords( 3 );
    double* COORDS = coords.data( );
    //TODO: figure out how to use tardigradeHydra::spatialDimensions here
    std::vector< double > drot( 3 * 3);
    double* DROT   = drot.data( );
    //TODO: figure out how to use tardigradeHydra::spatialDimensions here
    std::vector< double > dfgrd0( 3 * 3);
    double* DFGRD0 = dfgrd0.data( );
    //TODO: figure out how to use tardigradeHydra::spatialDimensions here
    std::vector< double > dfgrd1( 3 * 3);
    double* DFGRD1 = dfgrd1.data( );

    //Sign of life test. Just run to see if any exceptions are thrown.
    tardigradeHydra::abaqusInterface(
        STRESS, STATEV, DDSDDE, SSE,    SPD,
        SCD,    RPL,    DDSDDT, DRPLDE, DRPLDT,
        STRAN,  DSTRAN, TIME,   DTIME,  TEMP,
        DTEMP,  PREDEF, DPRED,  CMNAME, NDI,
        NSHR,   NTENS,  NSTATV, PROPS,  NPROPS,
        COORDS, DROT,   PNEWDT, CELENT, DFGRD0,
        DFGRD1, NOEL,   NPT,    LAYER,  KSPT,
        JSTEP,  KINC );

    //Check for nStateVariables thrown exception
    std::vector< double > temp = { 1 };
    double* STATEV_incorrect = temp.data( );
    int NSTATV_incorrect = temp.size( );
    BOOST_CHECK_THROW(
        tardigradeHydra::abaqusInterface(
           STRESS, STATEV_incorrect, DDSDDE,           SSE,    SPD,
           SCD,    RPL,              DDSDDT,           DRPLDE, DRPLDT,
           STRAN,  DSTRAN,           TIME,             DTIME,  TEMP,
           DTEMP,  PREDEF,           DPRED,            CMNAME, NDI,
           NSHR,   NTENS,            NSTATV_incorrect, PROPS,  NPROPS,
           COORDS, DROT,             PNEWDT,           CELENT, DFGRD0,
           DFGRD1, NOEL,             NPT,              LAYER,  KSPT,
           JSTEP,  KINC ),
        std::exception );

    //Check for nMaterialParameters thrown exception
    temp = { 1 };
    double* PROPS_incorrect = temp.data( );
    int NPROPS_incorrect = temp.size( );

    BOOST_CHECK_THROW(
        tardigradeHydra::abaqusInterface(
           STRESS, STATEV, DDSDDE, SSE,             SPD,
           SCD,    RPL,    DDSDDT, DRPLDE,          DRPLDT,
           STRAN,  DSTRAN, TIME,   DTIME,           TEMP,
           DTEMP,  PREDEF, DPRED,  CMNAME,          NDI,
           NSHR,   NTENS,  NSTATV, PROPS_incorrect, NPROPS_incorrect,
           COORDS, DROT,   PNEWDT, CELENT,          DFGRD0,
           DFGRD1, NOEL,   NPT,    LAYER,           KSPT,
           JSTEP,  KINC ),
        std::exception );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getTime ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkTime( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getDeltaTime ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkDeltaTime( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getTemperature ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkTemperature( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousTemperature ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousTemperature( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getDeformationGradient ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkDeformationGradient( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousDeformationGradient ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousDeformationGradient( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousStateVariables ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getParameters ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkParameters( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getNumConfigurations ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkNumConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getNumNonLinearSolveStateVariables ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkNumNonLinearSolveStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getConfigurations ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousConfigurations ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getInverseConfigurations ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkInverseConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousInverseConfigurations ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousInverseConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getNonLinearSolveStateVariables ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkNonLinearSolveStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousNonLinearSolveStateVariables ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousNonLinearSolveStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getAdditionalStateVariables ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkAdditionalStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousAdditionalStateVariables ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousAdditionalStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getRelativeTolerance ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkRelativeTolerance( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getAbsoluteTolerance ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkAbsoluteTolerance( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getLSAlpha ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkLSAlpha( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_decomposeStateVariableVector ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatMatrix configurationsAnswer = {
                                           { 1.05936416, -0.30634264, -0.86204929,
                                             0.03274673,  0.17917379, -0.22403642,
                                             1.24895144, -0.21368066, -1.05360316 },
                                           { 1.53155137,  0.53182759,  0.63440096,
                                             0.84943179,  1.72445532,  0.61102351,
                                             0.72244338,  0.32295891,  1.36178866 },
                                           { 1.22826323,  0.29371405,  0.63097612,
                                             0.09210494,  1.43370117,  0.43086276,
                                             0.4936851 ,  0.42583029,  1.31226122 },
                                           { 1.42635131,  0.89338916,  0.94416002,
                                             0.50183668,  1.62395295,  0.1156184 ,
                                             0.31728548,  0.41482621,  1.86630916 }
                                       };

    floatMatrix previousConfigurationsAnswer = {
                                                   { -0.41803693, -0.10444357,  0.58891991,
                                                      0.33131016, -0.34276201, -0.04604603,
                                                      1.32949612, -0.42712653, -1.03631144 },
                                                   {  1.53155137,  0.53182759,  0.63440096,
                                                      0.84943179,  1.72445532,  0.61102351,
                                                      0.72244338,  0.32295891,  1.36178866 },
                                                   {  1.22826323,  0.29371405,  0.63097612,
                                                      0.09210494,  1.43370117,  0.43086276,
                                                      0.4936851 ,  0.42583029,  1.31226122 },
                                                   {  1.42635131,  0.89338916,  0.94416002,
                                                      0.50183668,  1.62395295,  0.1156184 ,
                                                      0.31728548,  0.41482621,  1.86630916 }
                                               };

    floatMatrix inverseConfigurationsAnswer( 4 );

    floatMatrix previousInverseConfigurationsAnswer( 4 );

    for ( unsigned int i = 0; i < 4; i++ ){

        inverseConfigurationsAnswer[ i ]         = vectorTools::inverse( configurationsAnswer[ i ], 3, 3 );
        previousInverseConfigurationsAnswer[ i ] = vectorTools::inverse( previousConfigurationsAnswer[ i ], 3, 3 );

    }

    floatVector nonLinearSolveStateVariablesAnswer = { 0.25045537, 0.48303426, 0.98555979,
                                                       0.51948512, 0.61289453 };

    floatVector previousNonLinearSolveStateVariablesAnswer = { 0.25045537, 0.48303426, 0.98555979,
                                                               0.51948512, 0.61289453 };

    floatVector additionalStateVariablesAnswer = { 0.12062867, 0.8263408 , 0.60306013,
                                                   0.54506801, 0.34276383, 0.30412079  };

    floatVector previousAdditionalStateVariablesAnswer = { 0.12062867, 0.8263408 , 0.60306013,
                                                           0.54506801, 0.34276383, 0.30412079  };

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::decomposeStateVariableVector( hydra );

    BOOST_CHECK( vectorTools::fuzzyEquals( configurationsAnswer, *hydra.getConfigurations( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( previousConfigurationsAnswer, *hydra.getPreviousConfigurations( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( inverseConfigurationsAnswer, *hydra.getInverseConfigurations( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( previousInverseConfigurationsAnswer, *hydra.getPreviousInverseConfigurations( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( nonLinearSolveStateVariablesAnswer, *hydra.getNonLinearSolveStateVariables( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( previousNonLinearSolveStateVariablesAnswer, *hydra.getPreviousNonLinearSolveStateVariables( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( additionalStateVariablesAnswer, *hydra.getAdditionalStateVariables( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( previousAdditionalStateVariablesAnswer, *hydra.getPreviousAdditionalStateVariables( ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getSubConfiguration ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getSubConfiguration( 0, 4 ), *hydra.getDeformationGradient( ) ) );

    floatVector answer1 = { 2.24332648, 1.48246714, 2.02801682,
                            1.50380989, 2.98203598, 2.08079721,
                            1.58939152, 1.2551092 , 2.38201794 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getSubConfiguration( 1, 3 ), answer1 ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPrecedingConfiguration ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getPrecedingConfiguration( 4 ), *hydra.getDeformationGradient( ) ) );

    floatVector answer1 = { 0.73947165, -0.24328161, -0.68904986,
                            0.04049558,  0.25403825, -0.1748363 ,
                            0.97015752, -0.04452644, -0.77301275 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getPrecedingConfiguration( 2 ), answer1 ) );

    floatVector answer2 = { 0.54568475, -0.42501821, -0.54244544,
                           -0.01317666,  0.30165847, -0.09442353,
                            0.80588282, -0.10806097, -0.42143322 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getPrecedingConfiguration( 3 ), answer2 ) );

    floatVector answer3 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getPrecedingConfiguration( 0 ), answer3 ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getFollowingConfiguration ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer1 = { 2.09953091, 1.83604029, 2.3712323 ,
                            0.98756433, 2.58928197, 1.05684715,
                            1.33422708, 1.67694162, 2.96443669 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getFollowingConfiguration( 1 ), answer1 ) );

    floatVector answer2 = { 1.42635131, 0.89338916, 0.94416002,
                            0.50183668, 1.62395295, 0.1156184 ,
                            0.31728548, 0.41482621, 1.86630916 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getFollowingConfiguration( 2 ), answer2 ) );

    floatVector answer3 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getFollowingConfiguration( 3 ), answer3 ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousSubConfiguration ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getPreviousSubConfiguration( 0, 4 ), *hydra.getPreviousDeformationGradient( ) ) );

    floatVector answer1 = { 2.24332648, 1.48246714, 2.02801682,
                            1.50380989, 2.98203598, 2.08079721,
                            1.58939152, 1.2551092 , 2.38201794 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getPreviousSubConfiguration( 1, 3 ), answer1 ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousPrecedingConfiguration ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getPreviousPrecedingConfiguration( 4 ), *hydra.getPreviousDeformationGradient( ) ) );

    floatVector answer1 = { -0.30350143, -0.21223491,  0.47296395,
                             0.18299993, -0.42974886, -0.06195712,
                             0.92470041, -0.36418391, -0.8287879 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getPreviousPrecedingConfiguration( 2 ), answer1 ) );

    floatVector answer2 = { -0.15883228, -0.1920217 ,  0.33770597,
                             0.15460279, -0.58876502, -0.15099813,
                             0.69307214, -0.60345639, -0.66103563 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getPreviousPrecedingConfiguration( 3 ), answer2 ) );

    floatVector answer3 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getPreviousPrecedingConfiguration( 0 ), answer3 ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousFollowingConfiguration ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer1 = { 2.09953091, 1.83604029, 2.3712323 ,
                            0.98756433, 2.58928197, 1.05684715,
                            1.33422708, 1.67694162, 2.96443669 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getPreviousFollowingConfiguration( 1 ), answer1 ) );

    floatVector answer2 = { 1.42635131, 0.89338916, 0.94416002,
                            0.50183668, 1.62395295, 0.1156184 ,
                            0.31728548, 0.41482621, 1.86630916 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getPreviousFollowingConfiguration( 2 ), answer2 ) );

    floatVector answer3 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getPreviousFollowingConfiguration( 3 ), answer3 ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getSubConfigurationGradient ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatMatrix configurations = *hydra.getConfigurations( );

    floatVector x = vectorTools::appendVectors( configurations );

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 0;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatMatrix delta( numConfigurations, floatVector( dimension * dimension, 0 ) );

        unsigned int config = i / 9;

        unsigned int index = i % 9;

        delta[ config ][ index ] = std::fabs( eps * configurations[ config ][ index ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ config ][ index ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradient, hydra.getSubConfigurationGradient( configurations, lower, upper ) ) );

    lower = 1;
    upper = 3;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatMatrix delta( numConfigurations, floatVector( dimension * dimension, 0 ) );

        unsigned int config = i / 9;

        unsigned int index = i % 9;

        delta[ config ][ index ] = std::fabs( eps * configurations[ config ][ index ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ config ][ index ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradient, hydra.getSubConfigurationGradient( configurations, lower, upper ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getSubConfigurationGradient2 ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatMatrix configurations = *hydra.getConfigurations( );

    floatVector x = vectorTools::appendVectors( configurations );

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 0;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatMatrix delta( numConfigurations, floatVector( dimension * dimension, 0 ) );

        unsigned int config = i / 9;

        unsigned int index = i % 9;

        delta[ config ][ index ] = std::fabs( eps * configurations[ config ][ index ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ config ][ index ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradient, hydra.getSubConfigurationGradient( lower, upper ) ) );

    lower = 1;
    upper = 3;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatMatrix delta( numConfigurations, floatVector( dimension * dimension, 0 ) );

        unsigned int config = i / 9;

        unsigned int index = i % 9;

        delta[ config ][ index ] = std::fabs( eps * configurations[ config ][ index ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ config ][ index ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradient, hydra.getSubConfigurationGradient( lower, upper ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPrecedingConfigurationGradient ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatMatrix configurations = *hydra.getConfigurations( );

    floatVector x = vectorTools::appendVectors( configurations );

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 0;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatMatrix delta( numConfigurations, floatVector( dimension * dimension, 0 ) );

        unsigned int config = i / 9;

        unsigned int index = i % 9;

        delta[ config ][ index ] = std::fabs( eps * configurations[ config ][ index ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ config ][ index ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradient, hydra.getPrecedingConfigurationGradient( upper ) ) );

    lower = 0;
    upper = 3;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatMatrix delta( numConfigurations, floatVector( dimension * dimension, 0 ) );

        unsigned int config = i / 9;

        unsigned int index = i % 9;

        delta[ config ][ index ] = std::fabs( eps * configurations[ config ][ index ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ config ][ index ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradient, hydra.getPrecedingConfigurationGradient( upper ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getFollowingConfigurationGradient ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatMatrix configurations = *hydra.getConfigurations( );

    floatVector x = vectorTools::appendVectors( configurations );

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 1;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatMatrix delta( numConfigurations, floatVector( dimension * dimension, 0 ) );

        unsigned int config = i / 9;

        unsigned int index = i % 9;

        delta[ config ][ index ] = std::fabs( eps * configurations[ config ][ index ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower + 1, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower + 1, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ config ][ index ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradient, hydra.getFollowingConfigurationGradient( lower ) ) );

    lower = 2;
    upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatMatrix delta( numConfigurations, floatVector( dimension * dimension, 0 ) );

        unsigned int config = i / 9;

        unsigned int index = i % 9;

        delta[ config ][ index ] = std::fabs( eps * configurations[ config ][ index ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower + 1, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower + 1, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ config ][ index ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradient, hydra.getFollowingConfigurationGradient( lower ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousSubConfigurationGradient ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatMatrix configurations = *hydra.getPreviousConfigurations( );

    floatVector x = vectorTools::appendVectors( configurations );

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 0;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatMatrix delta( numConfigurations, floatVector( dimension * dimension, 0 ) );

        unsigned int config = i / 9;

        unsigned int index = i % 9;

        delta[ config ][ index ] = std::fabs( eps * configurations[ config ][ index ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ config ][ index ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradient, hydra.getPreviousSubConfigurationGradient( lower, upper ) ) );

    lower = 1;
    upper = 3;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatMatrix delta( numConfigurations, floatVector( dimension * dimension, 0 ) );

        unsigned int config = i / 9;

        unsigned int index = i % 9;

        delta[ config ][ index ] = std::fabs( eps * configurations[ config ][ index ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ config ][ index ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradient, hydra.getPreviousSubConfigurationGradient( lower, upper ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousPrecedingConfigurationGradient ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatMatrix configurations = *hydra.getPreviousConfigurations( );

    floatVector x = vectorTools::appendVectors( configurations );

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 0;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatMatrix delta( numConfigurations, floatVector( dimension * dimension, 0 ) );

        unsigned int config = i / 9;

        unsigned int index = i % 9;

        delta[ config ][ index ] = std::fabs( eps * configurations[ config ][ index ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ config ][ index ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradient, hydra.getPreviousPrecedingConfigurationGradient( upper ) ) );

    lower = 0;
    upper = 3;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatMatrix delta( numConfigurations, floatVector( dimension * dimension, 0 ) );

        unsigned int config = i / 9;

        unsigned int index = i % 9;

        delta[ config ][ index ] = std::fabs( eps * configurations[ config ][ index ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ config ][ index ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradient, hydra.getPreviousPrecedingConfigurationGradient( upper ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousFollowingConfigurationGradient ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatMatrix configurations = *hydra.getPreviousConfigurations( );

    floatVector x = vectorTools::appendVectors( configurations );

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 1;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatMatrix delta( numConfigurations, floatVector( dimension * dimension, 0 ) );

        unsigned int config = i / 9;

        unsigned int index = i % 9;

        delta[ config ][ index ] = std::fabs( eps * configurations[ config ][ index ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower + 1, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower + 1, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ config ][ index ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradient, hydra.getPreviousFollowingConfigurationGradient( lower ) ) );

    lower = 2;
    upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatMatrix delta( numConfigurations, floatVector( dimension * dimension, 0 ) );

        unsigned int config = i / 9;

        unsigned int index = i % 9;

        delta[ config ][ index ] = std::fabs( eps * configurations[ config ][ index ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower + 1, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower + 1, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ config ][ index ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradient, hydra.getPreviousFollowingConfigurationGradient( lower ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraTest_setFirstConfigurationGradients ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatType eps = 1e-6;

    floatMatrix dF1dF_answer( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatMatrix dF1dFn_answer( deformationGradient.size( ), floatVector( deformationGradient.size( ) * ( numConfigurations - 1 ), 0 ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        tardigradeHydra::hydraBase hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                                  previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::hydraBase hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                                  previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        floatVector F1_p, F1_m;

        F1_p = ( *hydra_p.getConfigurations( ) )[ 0 ];

        F1_m = ( *hydra_m.getConfigurations( ) )[ 0 ];

        for ( unsigned int j = 0; j < F1_p.size( ); j++ ){

            dF1dF_answer[ j ][ i ] = ( F1_p[ j ] - F1_m[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dF1dF_answer, *hydra.getdF1dF( ) ) );

    tardigradeHydra::unit_test::hydraBaseTester::checkdF1dF( hydra );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ) * deformationGradient.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] += eps * std::fabs( previousStateVariables[ i ] ) + eps;

        tardigradeHydra::hydraBase hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                  previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::hydraBase hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                  previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        floatVector F1_p, F1_m;

        F1_p = ( *hydra_p.getConfigurations( ) )[ 0 ];

        F1_m = ( *hydra_m.getConfigurations( ) )[ 0 ];

        for ( unsigned int j = 0; j < F1_p.size( ); j++ ){

            dF1dFn_answer[ j ][ i ] = ( F1_p[ j ] - F1_m[ j ] ) / ( 2 * delta[ i ] );

        }
        
    }

    BOOST_CHECK( vectorTools::fuzzyEquals( dF1dFn_answer, *hydra.getdF1dFn( ) ) );

    tardigradeHydra::unit_test::hydraBaseTester::checkdF1dFn( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraTest_setPreviousFirstConfigurationGradients ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatType eps = 1e-6;

    floatMatrix previousdF1dF_answer( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatMatrix previousdF1dFn_answer( deformationGradient.size( ), floatVector( deformationGradient.size( ) * ( numConfigurations - 1 ), 0 ) );

    for ( unsigned int i = 0; i < previousDeformationGradient.size( ); i++ ){

        floatVector delta( previousDeformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        tardigradeHydra::hydraBase hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                                  previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::hydraBase hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                                  previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        floatVector F1_p, F1_m;

        F1_p = ( *hydra_p.getPreviousConfigurations( ) )[ 0 ];

        F1_m = ( *hydra_m.getPreviousConfigurations( ) )[ 0 ];

        for ( unsigned int j = 0; j < F1_p.size( ); j++ ){

            previousdF1dF_answer[ j ][ i ] = ( F1_p[ j ] - F1_m[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( previousdF1dF_answer, *hydra.getPreviousdF1dF( ) ) );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ) * deformationGradient.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] += eps * std::fabs( previousStateVariables[ i ] ) + eps;

        tardigradeHydra::hydraBase hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                  previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::hydraBase hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                  previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        floatVector F1_p, F1_m;

        F1_p = ( *hydra_p.getPreviousConfigurations( ) )[ 0 ];

        F1_m = ( *hydra_m.getPreviousConfigurations( ) )[ 0 ];

        for ( unsigned int j = 0; j < F1_p.size( ); j++ ){

            previousdF1dFn_answer[ j ][ i ] = ( F1_p[ j ] - F1_m[ j ] ) / ( 2 * delta[ i ] );

        }
        
    }

    BOOST_CHECK( vectorTools::fuzzyEquals( previousdF1dFn_answer, *hydra.getPreviousdF1dFn( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residualBase_residualBase ){

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    tardigradeHydra::residualBase residual( &hydra, numEquations );

    BOOST_CHECK( residual.hydra == &hydra );

    BOOST_CHECK( *residual.getNumEquations( ) == numEquations );

}

BOOST_AUTO_TEST_CASE( test_residualBase_checkDefaults ){

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    tardigradeHydra::residualBase residual( &hydra, numEquations );

    BOOST_CHECK_THROW( residual.setResidual( ), std::nested_exception );

    BOOST_CHECK_THROW( residual.setJacobian( ), std::nested_exception );

    BOOST_CHECK_THROW( residual.setdRdF( ), std::nested_exception );

    BOOST_CHECK_THROW( residual.setdRdT( ), std::nested_exception );

    BOOST_CHECK_NO_THROW( residual.setAdditionalDerivatives( ) );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setResidual ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            using tardigradeHydra::residualBase::setResidual;

            floatVector residual = { 1, 2, 3 };

            residualBaseMock( tardigradeHydra::hydraBase *hydra, unsigned int numEquations ) : residualBase( hydra, numEquations ){ }
    
            virtual void setResidual( ){
    
                setResidual( residual );
    
            }

    };

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    residualBaseMock residual( &hydra, numEquations );

    BOOST_CHECK( vectorTools::fuzzyEquals( *residual.getResidual( ), residual.residual ) );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setJacobian ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            using tardigradeHydra::residualBase::setJacobian;

            floatMatrix jacobian = { { 1, 2, 3 }, { 4, 5, 6 } };
    
            residualBaseMock( tardigradeHydra::hydraBase *hydra, unsigned int numEquations ) : residualBase( hydra, numEquations ){ }

            virtual void setJacobian( ){
    
                setJacobian( jacobian );
    
            }

    };

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    residualBaseMock residual( &hydra, numEquations );

    BOOST_CHECK( vectorTools::fuzzyEquals( *residual.getJacobian( ), residual.jacobian ) );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setdRdF ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            using tardigradeHydra::residualBase::setdRdF;

            floatMatrix dRdF = { { 1, 2, 3 }, { 4, 5, 6 } };

            residualBaseMock( tardigradeHydra::hydraBase *hydra, unsigned int numEquations ) : residualBase( hydra, numEquations ){ }
    
            virtual void setdRdF( ){
    
                setdRdF( dRdF );
    
            }

    };

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    residualBaseMock residual( &hydra, numEquations );

    BOOST_CHECK( vectorTools::fuzzyEquals( *residual.getdRdF( ), residual.dRdF ) );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setdRdT ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            using tardigradeHydra::residualBase::setdRdT;

            floatVector dRdT = { 4, 5, 6 };

            residualBaseMock( tardigradeHydra::hydraBase *hydra, unsigned int numEquations ) : residualBase( hydra, numEquations ){ }
    
            virtual void setdRdT( ){
    
                setdRdT( dRdT );
    
            }

    };

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    residualBaseMock residual( &hydra, numEquations );

    BOOST_CHECK( vectorTools::fuzzyEquals( *residual.getdRdT( ), residual.dRdT ) );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setAdditionalDerivatives ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            using tardigradeHydra::residualBase::setAdditionalDerivatives;

            floatMatrix additionalDerivatives = { { 4, 5, 6 } };

            residualBaseMock( tardigradeHydra::hydraBase *hydra, unsigned int numEquations ) : residualBase( hydra, numEquations ){ }
    
            virtual void setAdditionalDerivatives( ){
    
                setAdditionalDerivatives( additionalDerivatives );
    
            }

    };

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    residualBaseMock residual( &hydra, numEquations );

    BOOST_CHECK( vectorTools::fuzzyEquals( *residual.getAdditionalDerivatives( ), residual.additionalDerivatives ) );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setCauchyStress ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            using tardigradeHydra::residualBase::setCauchyStress;

            floatVector cauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            residualBaseMock( tardigradeHydra::hydraBase *hydra, unsigned int numEquations ) : residualBase( hydra, numEquations ){ }
    
            virtual void setCauchyStress( ){
    
                setCauchyStress( cauchyStress );
    
            }

    };

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    residualBaseMock residual( &hydra, numEquations );

    BOOST_CHECK( vectorTools::fuzzyEquals( *residual.getCauchyStress( ), residual.cauchyStress ) );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setCurrentAdditionalStateVariables ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            using tardigradeHydra::residualBase::setCurrentAdditionalStateVariables;

            floatVector currentAdditionalStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };

            residualBaseMock( tardigradeHydra::hydraBase *hydra, unsigned int numEquations ) : residualBase( hydra, numEquations ){ }
    
            virtual void setCurrentAdditionalStateVariables( ){
    
                setCurrentAdditionalStateVariables( currentAdditionalStateVariables );
    
            }

    };

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    residualBaseMock residual( &hydra, numEquations );

    BOOST_CHECK( vectorTools::fuzzyEquals( *residual.getCurrentAdditionalStateVariables( ), residual.currentAdditionalStateVariables ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setResidualClasses ){

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::residualBase r1;
        
            tardigradeHydra::residualBase r2;
        
            tardigradeHydra::residualBase r3;

            unsigned int s1 = 36;

            unsigned int s2 = 2;

            unsigned int s3 = 3;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                r1 = tardigradeHydra::residualBase( this, s1 );

                r2 = tardigradeHydra::residualBase( this, s2 );

                r3 = tardigradeHydra::residualBase( this, s3 );

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                residuals[ 0 ] = &r1;

                residuals[ 1 ] = &r2;

                residuals[ 2 ] = &r3;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_CHECK_NO_THROW( hydra.setResidualClasses( ) );

    BOOST_CHECK( *( *hydra.getResidualClasses( ) )[ 0 ]->getNumEquations( ) == hydra.s1 );

    BOOST_CHECK( *( *hydra.getResidualClasses( ) )[ 1 ]->getNumEquations( ) == hydra.s2 );

    BOOST_CHECK( *( *hydra.getResidualClasses( ) )[ 2 ]->getNumEquations( ) == hydra.s3 );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setResidualClasses2 ){

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::residualBase r1;
        
            tardigradeHydra::residualBase r2;
        
            tardigradeHydra::residualBase r3;

            unsigned int s1 = 2;

            unsigned int s2 = 4;

            unsigned int s3 = 3;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                r1 = tardigradeHydra::residualBase( this, s1 );

                r2 = tardigradeHydra::residualBase( this, s2 );

                r3 = tardigradeHydra::residualBase( this, s3 );

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                residuals[ 0 ] = &r1;

                residuals[ 1 ] = &r2;

                residuals[ 2 ] = &r3;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_CHECK_THROW( hydra.setResidualClasses( ), std::nested_exception );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_formNonLinearProblem ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            unsigned int numVariables = 41;

            using tardigradeHydra::residualBase::residualBase;

            using tardigradeHydra::residualBase::setResidual;

            using tardigradeHydra::residualBase::setJacobian;

            using tardigradeHydra::residualBase::setdRdF;

            using tardigradeHydra::residualBase::setdRdT;

            using tardigradeHydra::residualBase::setAdditionalDerivatives;

            virtual void setResidual( ){

                floatVector residual( *getNumEquations( ), 0 );

                for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                    residual[ i ] = i;

                }

                setResidual( residual );

            }

            virtual void setJacobian( ){

                floatMatrix jacobian( *getNumEquations( ), floatVector( numVariables, 0 ) );

                for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                    for ( unsigned int j = 0; j < numVariables; j++ ){

                        jacobian[ i ][ j ] = i + 0.1 * j;

                    }

                }

                setJacobian( jacobian );

            }

            virtual void setdRdF( ){

                floatMatrix dRdF( *getNumEquations( ), floatVector( 9, 0 ) );

                for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                    for ( unsigned int j = 0; j < 9; j++ ){

                        dRdF[ i ][ j ] = i - 0.1 * j;

                    }

                }

                setdRdF( dRdF );

            }

            virtual void setdRdT( ){

                floatVector dRdT( *getNumEquations( ), 0 );

                for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                    dRdT[ i ] = 0.3 * i;

                }

                setdRdT( dRdT );

            }

            virtual void setAdditionalDerivatives( ){

                floatMatrix additionalDerivatives( *getNumEquations( ), floatVector( 4,  0 ) );

                for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                    for ( unsigned int j = 0; j < 4; j++ ){

                        additionalDerivatives[ i ][ j ] = i - 0.7 * j;

                    }

                }

                setAdditionalDerivatives( additionalDerivatives );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            residualBaseMock r1;
        
            residualBaseMock r2;
        
            residualBaseMock r3;

            unsigned int s1 = 36;

            unsigned int s2 = 2;

            unsigned int s3 = 3;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                r1 = residualBaseMock( this, s1 );

                r2 = residualBaseMock( this, s2 );

                r3 = residualBaseMock( this, s3 );

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                residuals[ 0 ] = &r1;

                residuals[ 1 ] = &r2;

                residuals[ 2 ] = &r3;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    hydraBaseMock hydraGet( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    hydraGet.setResidualClasses( );

    floatVector residualAnswer = vectorTools::appendVectors( { *hydraGet.r1.getResidual( ),
                                                               *hydraGet.r2.getResidual( ),
                                                               *hydraGet.r3.getResidual( ) } );

    floatMatrix jacobianAnswer = vectorTools::appendVectors( { *hydraGet.r1.getJacobian( ),
                                                               *hydraGet.r2.getJacobian( ),
                                                               *hydraGet.r3.getJacobian( ) } );

    floatMatrix dRdFAnswer = vectorTools::appendVectors( { *hydraGet.r1.getdRdF( ),
                                                           *hydraGet.r2.getdRdF( ),
                                                           *hydraGet.r3.getdRdF( ) } );

    floatVector dRdTAnswer = vectorTools::appendVectors( { *hydraGet.r1.getdRdT( ),
                                                           *hydraGet.r2.getdRdT( ),
                                                           *hydraGet.r3.getdRdT( ) } );

    floatMatrix additionalDerivativesAnswer = vectorTools::appendVectors( { *hydraGet.r1.getAdditionalDerivatives( ),
                                                                            *hydraGet.r2.getAdditionalDerivatives( ),
                                                                            *hydraGet.r3.getAdditionalDerivatives( ) } );

    tardigradeHydra::unit_test::hydraBaseTester::formNonLinearProblem( hydra );

    BOOST_CHECK( vectorTools::fuzzyEquals( residualAnswer, *hydra.getResidual( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( jacobianAnswer, hydra.getJacobian( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dRdFAnswer, hydra.getdRdF( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( dRdTAnswer, *hydra.getdRdT( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( additionalDerivativesAnswer, hydra.getAdditionalDerivatives( ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_initializeUnknownVector ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            unsigned int numVariables = 41;

            using tardigradeHydra::residualBase::residualBase;

            using tardigradeHydra::residualBase::setResidual;

            using tardigradeHydra::residualBase::setJacobian;

            using tardigradeHydra::residualBase::setdRdF;

            using tardigradeHydra::residualBase::setdRdT;

            using tardigradeHydra::residualBase::setAdditionalDerivatives;

            virtual void setResidual( ){

                floatVector residual( *getNumEquations( ), 0 );

                for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                    residual[ i ] = i;

                }

                setResidual( residual );

            }

            virtual void setJacobian( ){

                floatMatrix jacobian( *getNumEquations( ), floatVector( numVariables, 0 ) );

                for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                    for ( unsigned int j = 0; j < numVariables; j++ ){

                        jacobian[ i ][ j ] = i + 0.1 * j;

                    }

                }

                setJacobian( jacobian );

            }

            virtual void setdRdF( ){

                floatMatrix dRdF( *getNumEquations( ), floatVector( 9, 0 ) );

                for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                    for ( unsigned int j = 0; j < 9; j++ ){

                        dRdF[ i ][ j ] = i - 0.1 * j;

                    }

                }

                setdRdF( dRdF );

            }

            virtual void setdRdT( ){

                floatVector dRdT( *getNumEquations( ), 0 );

                for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                    dRdT[ i ] = 0.3 * i;

                }

                setdRdT( dRdT );

            }

    };

    class residualBaseMockStress : public residualBaseMock{

        public:

            floatVector cauchyStress = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            using residualBaseMock::residualBaseMock;

            using residualBaseMock::setResidual;

            using residualBaseMock::setJacobian;

            using residualBaseMock::setdRdF;

            using residualBaseMock::setdRdT;

            using residualBaseMock::setAdditionalDerivatives;

            using tardigradeHydra::residualBase::setCauchyStress;

            virtual void setCauchyStress( ){

                setCauchyStress( cauchyStress );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            residualBaseMockStress r1;
        
            residualBaseMock r2;
        
            residualBaseMock r3;

            unsigned int s1 = 36;

            unsigned int s2 = 2;

            unsigned int s3 = 3;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                r1 = residualBaseMockStress( this, s1 );

                r2 = residualBaseMock( this, s2 );

                r3 = residualBaseMock( this, s3 );

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                residuals[ 0 ] = &r1;

                residuals[ 1 ] = &r2;

                residuals[ 2 ] = &r3;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    hydra.setResidualClasses( );

    tardigradeHydra::unit_test::hydraBaseTester::initializeUnknownVector( hydra );

    floatVector unknownVectorAnswer = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                        1.53155137, 0.53182759, 0.63440096, 0.84943179, 1.72445532, 0.61102351, 0.72244338, 0.32295891, 1.36178866,
                                        1.22826323, 0.29371405, 0.63097612, 0.09210494, 1.43370117, 0.43086276, 0.4936851 , 0.42583029, 1.31226122,
                                        1.42635131, 0.89338916, 0.94416002, 0.50183668, 1.62395295, 0.1156184 , 0.31728548, 0.41482621, 1.86630916,
                                        0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453 };
                                         

    BOOST_CHECK( vectorTools::fuzzyEquals( unknownVectorAnswer, *hydra.getUnknownVector( ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setTolerance ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector residual = { 1, 2, -3, 0 };

    floatVector unknownVector = { -2, 5, 10, 0.3 };

    floatVector toleranceAnswer = 1e-9 * (vectorTools::abs( residual ) + vectorTools::abs( unknownVector )) + 1e-9;

    tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra, unknownVector );

    BOOST_CHECK( vectorTools::fuzzyEquals( *hydra.getTolerance( ), toleranceAnswer ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_checkConvergence ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector residual = { 1, 2, -3, 0 };

    floatVector unknownVector = { -2, 5, 10, 0.3 };

    tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra, unknownVector );

    BOOST_CHECK( !hydra.checkConvergence( ) );

    residual = { 0, 2, -3, 0 };

    tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    BOOST_CHECK( !hydra.checkConvergence( ) );

    residual = { 0, 0, 0, 0 };

    tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    BOOST_CHECK( hydra.checkConvergence( ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getLSResidualNorm ){

    tardigradeHydra::hydraBase hydra;

    floatVector residual = { 1, 2, 3 };

    floatType lsResidualNormAnswer = vectorTools::l2norm( residual );
      tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    BOOST_CHECK( vectorTools::fuzzyEquals( *hydra.getLSResidualNorm( ), lsResidualNormAnswer ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_checkLSConvergence ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector residual = { 1, 2, -3, 0 };

    tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    BOOST_CHECK( !hydra.checkLSConvergence( ) );

    residual = { 0.1, 2, -3, 0 };

    tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    BOOST_CHECK( hydra.checkLSConvergence( ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_decomposeUnknownVector ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector unknownVector = { -0.15564097, -0.82143154, -0.27568448, -0.86079836,  0.07935105,
                                   0.14494499,  0.24999455,  0.43334209, -0.56655812, -0.30347709,
                                  -0.6423031 ,  0.31882307, -0.40934079, -0.68385209, -0.93550481,
                                   0.87222402, -0.44175239, -0.65098503, -0.70008955, -0.42273504,
                                   0.28778728,  0.04507432,  0.90587052, -0.23518004,  0.1069888 ,
                                  -0.63207718, -0.3026225 ,  0.21407967,  0.2462129 ,  0.62193457,
                                   0.88751802, -0.06915571,  0.81058704,  0.87534601,  0.58086916,
                                   0.68175081,  0.7558753 ,  0.90210319, -0.15050147, -0.83735103,
                                  -0.8684919 };

    floatVector cauchyStressAnswer = { -0.15564097, -0.82143154, -0.27568448, -0.86079836,  0.07935105,
                                        0.14494499,  0.24999455,  0.43334209, -0.56655812 };

    floatMatrix configurationsAnswer = {  { -2.51323148, -2.34583945,  1.64112676,  0.66198512,  1.02970651,  0.93612154, -1.30632269,  0.42502445,  2.20594616 },
                                          { -0.30347709, -0.6423031 ,  0.31882307, -0.40934079, -0.68385209, -0.93550481,  0.87222402, -0.44175239, -0.65098503 },
                                          { -0.70008955, -0.42273504,  0.28778728,  0.04507432,  0.90587052, -0.23518004,  0.1069888 , -0.63207718, -0.3026225  },
                                          {  0.21407967,  0.2462129 ,  0.62193457,  0.88751802, -0.06915571,  0.81058704,  0.87534601,  0.58086916,  0.68175081 } };

    floatVector isvAnswer = { 0.7558753 ,  0.90210319, -0.15050147, -0.83735103, -0.8684919 };

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra, unknownVector );

    tardigradeHydra::unit_test::hydraBaseTester::decomposeUnknownVector( hydra );

    BOOST_CHECK( vectorTools::fuzzyEquals( *hydra.getCauchyStress( ), cauchyStressAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *hydra.getConfigurations( ), configurationsAnswer ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( *hydra.getNonLinearSolveStateVariables( ), isvAnswer ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_solveNonLinearProblem ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector initialUnknownVector = { 1, 2, 3 };

            std::vector< unsigned int > numLSIterations = { 1, 2, 1, 3 };

            std::vector< std::vector< floatVector > > residual = {
                                                                     {
                                                                       { 1, 2, 3 } 
                                                                     },
                                                                     { 
                                                                       { 4, 5, 6 },
                                                                       { 7, 8, 9 }
                                                                     },
                                                                     { 
                                                                       { 10, 9, 8 }
                                                                     },
                                                                     { 
                                                                       { 7, 6, 5 },
                                                                       { 4, 3, 2 },
                                                                       { 1, 1, 1 }
                                                                     },
                                                                 };

            std::vector< std::vector< floatVector > > flatJacobian = {
                                                                         {
                                                                           {  0.99951474, -0.18854971, -0.59377647,  0.73798389,  0.50649978,
                                                                             -0.24789123, -0.32889876, -0.60029673, -0.14211995 } 
                                                                         },
                                                                         { 
                                                                           {  0.39531964, -0.92103617,  0.15529267, -0.88403283, -0.29325643,
                                                                              0.33220876, -0.86852642,  0.29966728,  0.83480685 },
                                                                           { -0.52407245,  0.35366894,  0.08551378,  0.70688925, -0.14848006,
                                                                              0.24923792,  0.67051005, -0.08228108,  0.65747464 }
                                                                         },
                                                                         { 
                                                                           { -0.81643355,  0.85279545, -0.18966435, -0.42616781,  0.53676336,
                                                                             -0.92714295, -0.624413  ,  0.09900722, -0.01636474 }
                                                                         },
                                                                         { 
                                                                           { -0.2448239 ,  0.74837591,  0.68140055, -0.14194502,  0.67937811,
                                                                             -0.28307126,  0.01624776,  0.02491045, -0.39841209 },
                                                                           { -0.65594132, -0.44378971, -0.43630527,  0.13856471,  0.95333808,
                                                                             -0.89484459, -0.26311552, -0.94030312,  0.69678361 },
                                                                           { -0.82400549, -0.05626004, -0.17170742,  0.34131108,  0.49817018,
                                                                             -0.78483747,  0.95727384,  0.23324588, -0.22667228 }
                                                                         },
                                                                     };

            std::vector< std::vector< floatVector > > expectedXVectors = {
                                                                             {
                                                                               { 17.15005309, -9.55454694, 35.53888217 } 
                                                                             },
                                                                             { 
                                                                               { 18.95259726, -5.57445709, 28.79822658 },
                                                                               { 18.05132518, -7.56450202, 32.16855437 }
                                                                             },
                                                                             { 
                                                                               { 31.01502518, -5.8216338 , 36.92596115 }
                                                                             },
                                                                             { 
                                                                               { 407.30865444,  77.76786468,  70.0478993  },
                                                                               { 219.16183981,  35.97311544,  53.48693023 },
                                                                               { 125.0884325 ,  15.07574082,  45.20644569 }
                                                                             },
                                                                         };

        private:

            using tardigradeHydra::hydraBase::getResidual;

            virtual void initializeUnknownVector( ){

                tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( *this, initialUnknownVector );

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual[ 0 ][ 0 ] );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, flatJacobian[ 0 ][ 0 ] );

            }

            virtual bool checkConvergence( ){

                unsigned int iteration = tardigradeHydra::unit_test::hydraBaseTester::get_iteration( *this );

                if ( iteration < residual.size( ) ){

                    return false;

                }

                return true;

            }

            virtual bool checkLSConvergence( ){

                unsigned int iteration = tardigradeHydra::unit_test::hydraBaseTester::get_iteration( *this );

                unsigned int LSIteration = tardigradeHydra::unit_test::hydraBaseTester::get_LSIteration( *this );

                getResidual( );

                if ( LSIteration < residual[ iteration ].size( ) - 1 ){

                    return false;

                }

                return true;

            }

            virtual void formNonLinearProblem( ){

                unsigned int iteration = tardigradeHydra::unit_test::hydraBaseTester::get_iteration( *this );

                unsigned int LSIteration = tardigradeHydra::unit_test::hydraBaseTester::get_LSIteration( *this );

                unsigned int iterationOffset = 0;

                unsigned int LSoffset = 0;

                if ( LSIteration < residual[ iteration ].size( ) - 1 ){

                    LSoffset += 1;

                }
                else if ( iteration < residual.size( ) - 1 ){

                    iterationOffset += 1;
                    LSIteration = 0;

                }

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual[ iteration + iterationOffset ][ LSIteration + LSoffset ] );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, flatJacobian[ iteration + iterationOffset ][ LSIteration + LSoffset ] );

            }

            virtual void updateUnknownVector( const floatVector &newUnknownVector ){

                tardigradeHydra::unit_test::hydraBaseTester::resetIterationData( *this );

                unsigned int iteration = tardigradeHydra::unit_test::hydraBaseTester::get_iteration( *this );

                unsigned int LSIteration = tardigradeHydra::unit_test::hydraBaseTester::get_LSIteration( *this );

                BOOST_CHECK( vectorTools::fuzzyEquals( expectedXVectors[ iteration ][ LSIteration ], newUnknownVector ) );

                tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( *this, newUnknownVector );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::solveNonLinearProblem( hydra );

}

BOOST_AUTO_TEST_CASE( test_residual_setdRdF ){

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::linearElasticity::residual elasticity;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            floatVector unknownVector = {   1,  1,  1,  1,  1,  1,  1,  1,  1 };

            virtual void setResidualClasses( ){

                elasticity = tardigradeHydra::linearElasticity::residual( this, elasticitySize, *getParameters( ) );

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                residuals[ 0 ] = &elasticity;

                setResidualClasses( residuals );

            }

        private:

            virtual void initializeUnknownVector( ){

                tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( *this, unknownVector );

            }


    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 1.05, 0, 0,
                                        0.00, 1, 0,
                                        0.00, 1, 1};

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { };

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    floatVector unknownVector = {   1,  1,  1,  1,  1,  1,  1,  1,  1 };

    floatVector PK2Answer = {  73.836  ,   0.     ,   0.     ,
                                0.     , 124.72425,  56.7    ,
                                0.     ,  56.7    ,  68.02425 };

    floatVector cauchyStressAnswer = {  77.5278,   0.    ,   0.    ,
                                         0.    , 118.785 , 172.785 ,
                                         0.    , 172.785 , 291.57   };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_CHECK_NO_THROW( hydra.evaluate( ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( PK2Answer, *hydra.elasticity.getPK2Stress( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( cauchyStressAnswer, *hydra.getCauchyStress( ) ) );

}
