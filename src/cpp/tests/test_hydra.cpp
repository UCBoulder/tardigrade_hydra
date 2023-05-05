/**
  * \file test_hydra.cpp
  *
  * Tests for hydra
  */

#include<hydra.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_hydra
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef errorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef hydra::floatType floatType; //!< Redefinition of the floating point type
typedef hydra::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef hydra::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

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

namespace hydra{

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

                    try{
                        hydra.decomposeStateVariableVector( );
                    }
                    catch(std::exception &e){
                        errorTools::printNestedExceptions( e );
                    }

                    ERROR_TOOLS_CATCH( hydra.decomposeStateVariableVector( ) );

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
    error = hydra::sayHello(message);
    BOOST_CHECK( ! error );
    BOOST_CHECK( result.is_equal( answer ) );

    //Reset error code between tests
    error = NULL;

    //Check for "George" error
    message = "George";
    error = hydra::sayHello(message);
    BOOST_CHECK( error );

}

BOOST_AUTO_TEST_CASE( testAbaqusInterface ){
    /*!
     * Test the hydra abaqus interface
     */

    double double_scalar = 0.0;
    int int_scalar = 0;

    //Create nominally correct variable holders that match expected Abaqus Fortran interface
    //TODO: fill out nominally correct variable shape and values
    //Strings
    char CMNAME[ ] = "hydra";
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
    //TODO: figure out how to use hydra::spatialDimensions here
    std::vector< double > coords( 3 );
    double* COORDS = coords.data( );
    //TODO: figure out how to use hydra::spatialDimensions here
    std::vector< double > drot( 3 * 3);
    double* DROT   = drot.data( );
    //TODO: figure out how to use hydra::spatialDimensions here
    std::vector< double > dfgrd0( 3 * 3);
    double* DFGRD0 = dfgrd0.data( );
    //TODO: figure out how to use hydra::spatialDimensions here
    std::vector< double > dfgrd1( 3 * 3);
    double* DFGRD1 = dfgrd1.data( );

    //Sign of life test. Just run to see if any exceptions are thrown.
    hydra::abaqusInterface(
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
        hydra::abaqusInterface(
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
        hydra::abaqusInterface(
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

BOOST_AUTO_TEST_CASE( test_hydra_getTime ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkTime( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getDeltaTime ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkDeltaTime( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getTemperature ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkTemperature( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getPreviousTemperature ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkPreviousTemperature( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getDeformationGradient ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkDeformationGradient( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getPreviousDeformationGradient ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkPreviousDeformationGradient( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getPreviousStateVariables ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkPreviousStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getParameters ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkParameters( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getNumConfigurations ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkNumConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getNumNonLinearSolveStateVariables ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkNumNonLinearSolveStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getConfigurations ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getPreviousConfigurations ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkPreviousConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getInverseConfigurations ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkInverseConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getPreviousInverseConfigurations ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkPreviousInverseConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getNonLinearSolveStateVariables ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkNonLinearSolveStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getPreviousNonLinearSolveStateVariables ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkPreviousNonLinearSolveStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getAdditionalStateVariables ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkAdditionalStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_getPreviousAdditionalStateVariables ){

    hydra::hydraBase hydra;

    hydra::unit_test::hydraBaseTester::checkPreviousAdditionalStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydra_decomposeStateVariableVector ){

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

    hydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    hydra::unit_test::hydraBaseTester::decomposeStateVariableVector( hydra );

    BOOST_CHECK( vectorTools::fuzzyEquals( configurationsAnswer, *hydra.getConfigurations( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( previousConfigurationsAnswer, *hydra.getPreviousConfigurations( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( inverseConfigurationsAnswer, *hydra.getInverseConfigurations( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( previousInverseConfigurationsAnswer, *hydra.getPreviousInverseConfigurations( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( nonLinearSolveStateVariablesAnswer, *hydra.getNonLinearSolveStateVariables( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( previousNonLinearSolveStateVariablesAnswer, *hydra.getPreviousNonLinearSolveStateVariables( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( additionalStateVariablesAnswer, *hydra.getAdditionalStateVariables( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( previousAdditionalStateVariablesAnswer, *hydra.getPreviousAdditionalStateVariables( ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydra_getSubConfiguration ){

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

    hydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getSubConfiguration( 0, 3 ), *hydra.getDeformationGradient( ) ) );

    floatVector answer1 = { 2.24332648, 1.48246714, 2.02801682,
                            1.50380989, 2.98203598, 2.08079721,
                            1.58939152, 1.2551092 , 2.38201794 };

    BOOST_CHECK( vectorTools::fuzzyEquals( hydra.getSubConfiguration( 1, 2 ), answer1 ) );

}
