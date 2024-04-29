/**
  * \file test_tardigrade_hydra.cpp
  *
  * Tests for tardigrade_hydra
  */

#include<tardigrade_hydra.h>
#include<tardigrade_hydraLinearElasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_tardigrade_hydra
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
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

                    BOOST_CHECK( hydra._dimension == hydra.getDimension( ) );

                }

                static void checkSOTDimension( hydraBase &hydra ){

                    BOOST_CHECK( hydra._dimension * hydra._dimension == hydra.getSOTDimension( ) );

                }

                static void checkTOTDimension( hydraBase &hydra ){

                    BOOST_CHECK( hydra._dimension * hydra._dimension * hydra._dimension == hydra.getTOTDimension( ) );

                }

                static void checkFOTDimension( hydraBase &hydra ){

                    BOOST_CHECK( hydra._dimension * hydra._dimension * hydra._dimension * hydra._dimension == hydra.getFOTDimension( ) );

                }

                static void checkConfigurations( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._configurations.second == hydra.get_configurations( ) );

                }

                static void checkPreviousConfigurations( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._previousConfigurations.second == hydra.get_previousConfigurations( ) );

                }

                static void checkInverseConfigurations( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._inverseConfigurations.second == hydra.get_inverseConfigurations( ) );

                }

                static void checkPreviousInverseConfigurations( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._previousInverseConfigurations.second == hydra.get_previousInverseConfigurations( ) );

                }

                static void checkNonLinearSolveStateVariables( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._nonLinearSolveStateVariables.second == hydra.get_nonLinearSolveStateVariables( ) );

                }

                static void checkPreviousNonLinearSolveStateVariables( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._previousNonLinearSolveStateVariables.second == hydra.get_previousNonLinearSolveStateVariables( ) );

                }

                static void checkAdditionalStateVariables( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._additionalStateVariables.second == hydra.get_additionalStateVariables( ) );

                }

                static void checkPreviousAdditionalStateVariables( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._previousAdditionalStateVariables.second == hydra.get_previousAdditionalStateVariables( ) );

                }

                static void decomposeStateVariableVector( hydraBase &hydra ){

                    TARDIGRADE_ERROR_TOOLS_CATCH( hydra.decomposeStateVariableVector( ) );

                }

                static void decomposeUnknownVector( hydraBase &hydra ){

                    TARDIGRADE_ERROR_TOOLS_CATCH( hydra.decomposeUnknownVector( ) );

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

                static void checkUsePreconditioner( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._use_preconditioner == hydra.getUsePreconditioner( ) );

                }

                static void checkPreconditionerType( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._preconditioner_type == hydra.getPreconditionerType( ) );

                }

                static void checkPreconditioner( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._preconditioner.second == hydra.getFlatPreconditioner( ) );

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

                static void set_flatPreconditioner( hydraBase &hydra, const floatVector &value ){

                    hydra._preconditioner.second = value;
                    hydra._preconditioner.first = true;

                    hydra.addIterationData( &hydra._preconditioner );
                }

                static void set_usePreconditioner( hydraBase &hydra, const bool &value ){

                    hydra._use_preconditioner = value;

                }

                static void set_preconditionerType( hydraBase &hydra, const unsigned int &value ){

                    hydra._preconditioner_type = value;

                }

                static void set_flatdRdF( hydraBase &hydra, const floatVector &value ){

                    hydra._dRdF.second = value;
                    hydra._dRdF.first = true;

                    hydra.addIterationData( &hydra._dRdF );
                }

                static void set_flatdRdT( hydraBase &hydra, const floatVector &value ){

                    hydra._dRdT.second = value;
                    hydra._dRdT.first = true;

                    hydra.addIterationData( &hydra._dRdT );
                }

                static void set_configuration_unknown_count( hydraBase &hydra, const unsigned int &value ){

                    hydra._configuration_unknown_count = value;

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

BOOST_AUTO_TEST_CASE( testAbaqusInterface ){
    /*!
     * Test the tardigrade_hydra abaqus interface
     */

    double double_scalar = 0.0;
    int int_scalar = 0;

    //Create nominally correct variable holders that match expected Abaqus Fortran interface
    //TODO: fill out nominally correct variable shape and values
    //Strings
    char CMNAME[ ] = "tardigrade_hydra";
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

BOOST_AUTO_TEST_CASE( test_hydraBase_getDimension ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkDimension( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getSOTDimension ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkSOTDimension( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getTOTDimension ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkTOTDimension( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getFOTDimension ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkFOTDimension( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_configurations ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_previousConfigurations ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_inverseConfigurations ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkInverseConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_previousInverseConfigurations ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousInverseConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_nonLinearSolveStateVariables ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkNonLinearSolveStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_previousNonLinearSolveStateVariables ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousNonLinearSolveStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_additionalStateVariables ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkAdditionalStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_previousAdditionalStateVariables ){

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

        inverseConfigurationsAnswer[ i ]         = tardigradeVectorTools::inverse( configurationsAnswer[ i ], 3, 3 );
        previousInverseConfigurationsAnswer[ i ] = tardigradeVectorTools::inverse( previousConfigurationsAnswer[ i ], 3, 3 );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( configurationsAnswer ), *hydra.get_configurations( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( previousConfigurationsAnswer ), *hydra.get_previousConfigurations( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( inverseConfigurationsAnswer ), *hydra.get_inverseConfigurations( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( previousInverseConfigurationsAnswer ), *hydra.get_previousInverseConfigurations( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( nonLinearSolveStateVariablesAnswer, *hydra.get_nonLinearSolveStateVariables( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousNonLinearSolveStateVariablesAnswer, *hydra.get_previousNonLinearSolveStateVariables( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( additionalStateVariablesAnswer, *hydra.get_additionalStateVariables( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousAdditionalStateVariablesAnswer, *hydra.get_previousAdditionalStateVariables( ) ) );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getSubConfiguration( 0, 4 ), *hydra.getDeformationGradient( ) ) );

    floatVector answer1 = { 2.24332648, 1.48246714, 2.02801682,
                            1.50380989, 2.98203598, 2.08079721,
                            1.58939152, 1.2551092 , 2.38201794 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getSubConfiguration( 1, 3 ), answer1 ) );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getPrecedingConfiguration( 4 ), *hydra.getDeformationGradient( ) ) );

    floatVector answer1 = { 0.73947165, -0.24328161, -0.68904986,
                            0.04049558,  0.25403825, -0.1748363 ,
                            0.97015752, -0.04452644, -0.77301275 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getPrecedingConfiguration( 2 ), answer1 ) );

    floatVector answer2 = { 0.54568475, -0.42501821, -0.54244544,
                           -0.01317666,  0.30165847, -0.09442353,
                            0.80588282, -0.10806097, -0.42143322 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getPrecedingConfiguration( 3 ), answer2 ) );

    floatVector answer3 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getPrecedingConfiguration( 0 ), answer3 ) );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getFollowingConfiguration( 1 ), answer1 ) );

    floatVector answer2 = { 1.42635131, 0.89338916, 0.94416002,
                            0.50183668, 1.62395295, 0.1156184 ,
                            0.31728548, 0.41482621, 1.86630916 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getFollowingConfiguration( 2 ), answer2 ) );

    floatVector answer3 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getFollowingConfiguration( 3 ), answer3 ) );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getPreviousSubConfiguration( 0, 4 ), *hydra.getPreviousDeformationGradient( ) ) );

    floatVector answer1 = { 2.24332648, 1.48246714, 2.02801682,
                            1.50380989, 2.98203598, 2.08079721,
                            1.58939152, 1.2551092 , 2.38201794 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getPreviousSubConfiguration( 1, 3 ), answer1 ) );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getPreviousPrecedingConfiguration( 4 ), *hydra.getPreviousDeformationGradient( ) ) );

    floatVector answer1 = { -0.30350143, -0.21223491,  0.47296395,
                             0.18299993, -0.42974886, -0.06195712,
                             0.92470041, -0.36418391, -0.8287879 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getPreviousPrecedingConfiguration( 2 ), answer1 ) );

    floatVector answer2 = { -0.15883228, -0.1920217 ,  0.33770597,
                             0.15460279, -0.58876502, -0.15099813,
                             0.69307214, -0.60345639, -0.66103563 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getPreviousPrecedingConfiguration( 3 ), answer2 ) );

    floatVector answer3 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getPreviousPrecedingConfiguration( 0 ), answer3 ) );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getPreviousFollowingConfiguration( 1 ), answer1 ) );

    floatVector answer2 = { 1.42635131, 0.89338916, 0.94416002,
                            0.50183668, 1.62395295, 0.1156184 ,
                            0.31728548, 0.41482621, 1.86630916 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getPreviousFollowingConfiguration( 2 ), answer2 ) );

    floatVector answer3 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getPreviousFollowingConfiguration( 3 ), answer3 ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getSubConfigurationJacobian ){

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

    floatVector configurations = *hydra.get_configurations( );

    floatVector x = configurations;

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 0;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension, 0 );

        delta[ i ] = std::fabs( eps * configurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( gradient ), hydra.getSubConfigurationJacobian( configurations, lower, upper ) ) );

    lower = 1;
    upper = 3;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension, 0 );

        delta[ i ] = std::fabs( eps * configurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( gradient ), hydra.getSubConfigurationJacobian( configurations, lower, upper ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getSubConfigurationJacobian2 ){

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

    floatVector configurations = *hydra.get_configurations( );

    floatVector x = configurations;

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 0;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension, 0 );

        delta[ i ] = std::fabs( eps * configurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( gradient ), hydra.getSubConfigurationJacobian( lower, upper ) ) );

    lower = 1;
    upper = 3;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension, 0 );

        delta[ i ] = std::fabs( eps * configurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( gradient ), hydra.getSubConfigurationJacobian( lower, upper ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPrecedingConfigurationJacobian ){

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

    floatVector configurations = *hydra.get_configurations( );

    floatVector x = configurations;

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 0;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension, 0 );

        delta[ i ] = std::fabs( eps * configurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( gradient ), hydra.getPrecedingConfigurationJacobian( upper ) ) );

    lower = 0;
    upper = 3;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension, 0 );

        delta[ i ] = std::fabs( eps * configurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( gradient ), hydra.getPrecedingConfigurationJacobian( upper ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getFollowingConfigurationJacobian ){

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

    floatVector configurations = *hydra.get_configurations( );

    floatVector x = configurations;

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 1;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension, 0 );

        delta[ i ] = std::fabs( eps * configurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower + 1, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower + 1, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( gradient ), hydra.getFollowingConfigurationJacobian( lower ) ) );

    lower = 2;
    upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension, 0 );

        delta[ i ] = std::fabs( eps * configurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower + 1, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower + 1, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( gradient ), hydra.getFollowingConfigurationJacobian( lower ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousSubConfigurationJacobian ){

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

    floatVector configurations = *hydra.get_previousConfigurations( );

    floatVector x = configurations;

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 0;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension, 0 );

        delta[ i ] = std::fabs( eps * configurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( gradient ), hydra.getPreviousSubConfigurationJacobian( lower, upper ) ) );

    lower = 1;
    upper = 3;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension, 0 );

        delta[ i ] = std::fabs( eps * configurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( gradient ), hydra.getPreviousSubConfigurationJacobian( lower, upper ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousPrecedingConfigurationJacobian ){

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

    floatVector configurations = *hydra.get_previousConfigurations( );

    floatVector x = configurations;

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 0;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension, 0 );

        delta[ i ] = std::fabs( eps * configurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( gradient ), hydra.getPreviousPrecedingConfigurationJacobian( upper ) ) );

    lower = 0;
    upper = 3;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension, 0 );

        delta[ i ] = std::fabs( eps * configurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( gradient ), hydra.getPreviousPrecedingConfigurationJacobian( upper ) ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousFollowingConfigurationJacobian ){

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

    floatVector configurations = *hydra.get_previousConfigurations( );

    floatVector x = configurations;

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 1;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension );

        delta[ i ] = std::fabs( eps * configurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower + 1, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower + 1, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( gradient ), hydra.getPreviousFollowingConfigurationJacobian( lower ) ) );

    lower = 2;
    upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension );

        delta[ i ] = std::fabs( eps * configurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( configurations + delta, lower + 1, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( configurations - delta, lower + 1, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( gradient ), hydra.getPreviousFollowingConfigurationJacobian( lower ) ) );

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

        F1_p = floatVector( hydra_p.get_configurations( )->begin( ),
                            hydra_p.get_configurations( )->begin( ) + 9 );

        F1_m = floatVector( hydra_m.get_configurations( )->begin( ),
                            hydra_m.get_configurations( )->begin( ) + 9 );

        for ( unsigned int j = 0; j < F1_p.size( ); j++ ){

            dF1dF_answer[ j ][ i ] = ( F1_p[ j ] - F1_m[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( dF1dF_answer ), *hydra.get_dF1dF( ) ) );

    tardigradeHydra::unit_test::hydraBaseTester::checkdF1dF( hydra );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ) * deformationGradient.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] += eps * std::fabs( previousStateVariables[ i ] ) + eps;

        tardigradeHydra::hydraBase hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                  previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::hydraBase hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                  previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        floatVector F1_p, F1_m;

        F1_p = floatVector( hydra_p.get_configurations( )->begin( ),
                            hydra_p.get_configurations( )->begin( ) + 9 );

        F1_m = floatVector( hydra_m.get_configurations( )->begin( ),
                            hydra_m.get_configurations( )->begin( ) + 9 );

        for ( unsigned int j = 0; j < F1_p.size( ); j++ ){

            dF1dFn_answer[ j ][ i ] = ( F1_p[ j ] - F1_m[ j ] ) / ( 2 * delta[ i ] );

        }
        
    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( dF1dFn_answer ), *hydra.get_dF1dFn( ) ) );

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

        F1_p = floatVector( hydra_p.get_previousConfigurations( )->begin( ),
                            hydra_p.get_previousConfigurations( )->begin( ) + 9 );

        F1_m = floatVector( hydra_m.get_previousConfigurations( )->begin( ),
                            hydra_m.get_previousConfigurations( )->begin( ) + 9 );

        for ( unsigned int j = 0; j < F1_p.size( ); j++ ){

            previousdF1dF_answer[ j ][ i ] = ( F1_p[ j ] - F1_m[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( previousdF1dF_answer ), *hydra.get_previousdF1dF( ) ) );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ) * deformationGradient.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] += eps * std::fabs( previousStateVariables[ i ] ) + eps;

        tardigradeHydra::hydraBase hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                  previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::hydraBase hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                  previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        floatVector F1_p, F1_m;

        F1_p = floatVector( hydra_p.get_previousConfigurations( )->begin( ),
                            hydra_p.get_previousConfigurations( )->begin( ) + 9 );

        F1_m = floatVector( hydra_m.get_previousConfigurations( )->begin( ),
                            hydra_m.get_previousConfigurations( )->begin( ) + 9 );

        for ( unsigned int j = 0; j < F1_p.size( ); j++ ){

            previousdF1dFn_answer[ j ][ i ] = ( F1_p[ j ] - F1_m[ j ] ) / ( 2 * delta[ i ] );

        }
        
    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( tardigradeVectorTools::appendVectors( previousdF1dFn_answer ), *hydra.get_previousdF1dFn( ) ) );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *residual.getResidual( ), residual.residual ) );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setJacobian ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            using tardigradeHydra::residualBase::setJacobian;

            floatVector jacobian = { 1, 2, 3, 4, 5, 6 };
    
            residualBaseMock( tardigradeHydra::hydraBase *hydra, unsigned int numEquations ) : residualBase( hydra, numEquations ){ }

            virtual void setJacobian( ){
    
                setJacobian( jacobian );
    
            }

    };

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    residualBaseMock residual( &hydra, numEquations );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *residual.getJacobian( ), residual.jacobian ) );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setdRdF ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            using tardigradeHydra::residualBase::setdRdF;

            floatVector dRdF = { 1, 2, 3, 4, 5, 6 };

            residualBaseMock( tardigradeHydra::hydraBase *hydra, unsigned int numEquations ) : residualBase( hydra, numEquations ){ }
    
            virtual void setdRdF( ){
    
                setdRdF( dRdF );
    
            }

    };

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    residualBaseMock residual( &hydra, numEquations );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *residual.getdRdF( ), residual.dRdF ) );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *residual.getdRdT( ), residual.dRdT ) );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setAdditionalDerivatives ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            using tardigradeHydra::residualBase::setAdditionalDerivatives;

            floatVector additionalDerivatives = { 4, 5, 6 };

            residualBaseMock( tardigradeHydra::hydraBase *hydra, unsigned int numEquations ) : residualBase( hydra, numEquations ){ }
    
            virtual void setAdditionalDerivatives( ){
    
                setAdditionalDerivatives( additionalDerivatives );
    
            }

    };

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    residualBaseMock residual( &hydra, numEquations );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *residual.getAdditionalDerivatives( ), residual.additionalDerivatives ) );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setStress ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            using tardigradeHydra::residualBase::setStress;

            floatVector cauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            residualBaseMock( tardigradeHydra::hydraBase *hydra, unsigned int numEquations ) : residualBase( hydra, numEquations ){ }
    
            virtual void setStress( ){
    
                setStress( cauchyStress );
    
            }

    };

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    residualBaseMock residual( &hydra, numEquations );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *residual.getStress( ), residual.cauchyStress ) );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setPreviousStress ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            using tardigradeHydra::residualBase::setPreviousStress;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            residualBaseMock( tardigradeHydra::hydraBase *hydra, unsigned int numEquations ) : residualBase( hydra, numEquations ){ }
    
            virtual void setPreviousStress( ){
    
                setPreviousStress( previousCauchyStress );
    
            }

    };

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    residualBaseMock residual( &hydra, numEquations );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *residual.getPreviousStress( ), residual.previousCauchyStress ) );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *residual.getCurrentAdditionalStateVariables( ), residual.currentAdditionalStateVariables ) );

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

                setJacobian( tardigradeVectorTools::appendVectors( jacobian ) );

            }

            virtual void setdRdF( ){

                floatMatrix dRdF( *getNumEquations( ), floatVector( 9, 0 ) );

                for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                    for ( unsigned int j = 0; j < 9; j++ ){

                        dRdF[ i ][ j ] = i - 0.1 * j;

                    }

                }

                setdRdF( tardigradeVectorTools::appendVectors( dRdF ) );

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

                setAdditionalDerivatives( tardigradeVectorTools::appendVectors( additionalDerivatives ) );

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

    floatVector residualAnswer = tardigradeVectorTools::appendVectors( { *hydraGet.r1.getResidual( ),
                                                               *hydraGet.r2.getResidual( ),
                                                               *hydraGet.r3.getResidual( ) } );

    floatMatrix jacobianAnswer = tardigradeVectorTools::appendVectors( { tardigradeVectorTools::inflate( *hydraGet.r1.getJacobian( ), 36, 41 ),
                                                                         tardigradeVectorTools::inflate( *hydraGet.r2.getJacobian( ),  2, 41 ),
                                                                         tardigradeVectorTools::inflate( *hydraGet.r3.getJacobian( ),  3, 41 ) } );

    floatMatrix dRdFAnswer = tardigradeVectorTools::appendVectors( { tardigradeVectorTools::inflate( *hydraGet.r1.getdRdF( ), 36, 9 ),
                                                                     tardigradeVectorTools::inflate( *hydraGet.r2.getdRdF( ),  2, 9 ),
                                                                     tardigradeVectorTools::inflate( *hydraGet.r3.getdRdF( ),  3, 9 ) } );

    floatVector dRdTAnswer = tardigradeVectorTools::appendVectors( { *hydraGet.r1.getdRdT( ),
                                                           *hydraGet.r2.getdRdT( ),
                                                           *hydraGet.r3.getdRdT( ) } );

    floatMatrix additionalDerivativesAnswer = tardigradeVectorTools::appendVectors( { tardigradeVectorTools::inflate( *hydraGet.r1.getAdditionalDerivatives( ), 36, 4 ),
                                                                                      tardigradeVectorTools::inflate( *hydraGet.r2.getAdditionalDerivatives( ),  2, 4 ),
                                                                                      tardigradeVectorTools::inflate( *hydraGet.r3.getAdditionalDerivatives( ),  3, 4 ) } );

    tardigradeHydra::unit_test::hydraBaseTester::formNonLinearProblem( hydra );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( residualAnswer, *hydra.getResidual( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( jacobianAnswer, hydra.getJacobian( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dRdFAnswer, hydra.getdRdF( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dRdTAnswer, *hydra.getdRdT( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( additionalDerivativesAnswer, hydra.getAdditionalDerivatives( ) ) );

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

                setJacobian( tardigradeVectorTools::appendVectors( jacobian ) );

            }

            virtual void setdRdF( ){

                floatMatrix dRdF( *getNumEquations( ), floatVector( 9, 0 ) );

                for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                    for ( unsigned int j = 0; j < 9; j++ ){

                        dRdF[ i ][ j ] = i - 0.1 * j;

                    }

                }

                setdRdF( tardigradeVectorTools::appendVectors( dRdF ) );

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

            using tardigradeHydra::residualBase::setStress;

            virtual void setStress( ){

                setStress( cauchyStress );

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
                                         

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( unknownVectorAnswer, *hydra.getUnknownVector( ) ) );

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

    floatVector toleranceAnswer = 1e-9 * (tardigradeVectorTools::abs( residual ) + tardigradeVectorTools::abs( unknownVector )) + 1e-9;

    tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra, unknownVector );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *hydra.getTolerance( ), toleranceAnswer ) );

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

    floatType lsResidualNormAnswer = tardigradeVectorTools::l2norm( residual );
      tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *hydra.getLSResidualNorm( ), lsResidualNormAnswer ) );

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

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousStress ){


    class residualMockGood : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            residualMockGood residual1;

            tardigradeHydra::residualBase residual2;

            tardigradeHydra::residualBase remainder;

            using tardigradeHydra::hydraBase::hydraBase;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                residual1 = residualMockGood( this, 9 );

                residual2 = tardigradeHydra::residualBase( this, 9 );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &residual1;

                residuals[ 1 ] = &residual2;

                residuals[ 2 ] = &remainder;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = { 1.1, 0.1, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0 };

    floatVector previousDeformationGradient = { 1.0, 0.0, 0.0,
                                                0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0 };

    floatVector previousStateVariables = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0.02, 0.03 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  4, 5, 6 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *hydra.getPreviousStress( ) ) );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *hydra.getStress( ), cauchyStressAnswer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *hydra.get_configurations( ), tardigradeVectorTools::appendVectors( configurationsAnswer ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *hydra.get_nonLinearSolveStateVariables( ), isvAnswer ) );

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

                BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( expectedXVectors[ iteration ][ LSIteration ], newUnknownVector ) );

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

    hydraBaseMock hydra_pre( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension,
                             9, 1e-9, 1e-9, 20, 5, 1e-4, true, 0 );

    tardigradeHydra::unit_test::hydraBaseTester::solveNonLinearProblem( hydra_pre );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( PK2Answer, *hydra.elasticity.get_PK2Stress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( cauchyStressAnswer, *hydra.getStress( ) ) );

}

BOOST_AUTO_TEST_CASE( test_getConfiguration ){
    /*!
     * Boost test of the get-configuration command
     */

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getConfiguration( 1 ), floatVector( hydra.get_configurations( )->begin( ) + 1 * 9,
                                                                                               hydra.get_configurations( )->begin( ) + 2 * 9 ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.getPreviousConfiguration( 3 ), floatVector( hydra.get_previousConfigurations( )->begin( ) + 3 * 9,
                                                                                                       hydra.get_previousConfigurations( )->begin( ) + 4 * 9 ) ) );

}

BOOST_AUTO_TEST_CASE( test_computeTangents ){
    /*!
     * Boost test of the get-configuration command
     */

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

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector residual = { 1, 2, 3, 4, 5, 6, 7, 8 };

            floatVector jacobian = {  1.020e+00, -2.100e-02, -2.700e-02,  5.000e-03,  2.200e-02,
                                     -8.000e-03,  4.800e-02,  1.800e-02, -2.000e-03,  9.890e-01,
                                     -1.600e-02,  2.300e-02, -6.000e-03, -4.400e-02, -1.000e-02,
                                      2.400e-02, -3.200e-02, -3.200e-02,  1.003e+00,  3.000e-03,
                                      1.300e-02,  3.500e-02,  2.200e-02,  1.100e-02,  2.200e-02,
                                     -1.800e-02, -1.400e-02,  9.730e-01, -2.100e-02,  1.300e-02,
                                     -4.100e-02, -7.000e-03, -7.000e-03, -1.000e-03, -7.000e-03,
                                     -1.900e-02,  9.930e-01,  3.900e-02,  4.400e-02,  0.000e+00,
                                      1.200e-02, -3.800e-02, -1.800e-02, -9.000e-03,  3.700e-02,
                                      9.750e-01, -2.000e-03,  4.900e-02,  2.000e-03,  1.100e-02,
                                     -3.800e-02,  3.300e-02,  1.000e-02,  5.000e-03,  9.840e-01,
                                     -2.000e-02, -8.000e-03,  1.800e-02,  3.800e-02,  1.000e-03,
                                      1.700e-02,  9.000e-03,  1.200e-02,  1.017e+00 };

            floatVector dRdF = { 0.842, 0.083, 0.764, 0.244, 0.194, 0.572, 0.096, 0.885, 0.627,
                                 0.723, 0.016, 0.594, 0.557, 0.159, 0.153, 0.696, 0.319, 0.692,
                                 0.554, 0.389, 0.925, 0.842, 0.357, 0.044, 0.305, 0.398, 0.705,
                                 0.995, 0.356, 0.763, 0.593, 0.692 };

            floatVector dRdT = { 0.151, 0.399, 0.241, 0.343, 0.513, 0.667, 0.106, 0.131 };

            virtual void formNonLinearProblem( ) override{

                tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( *this, residual );

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, jacobian );

                tardigradeHydra::unit_test::hydraBaseTester::set_configuration_unknown_count( *this, 4 );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatdRdF( *this, dRdF );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatdRdT( *this, dRdT );

            }


    };
    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector dXdF = { -0.82465846, -0.07170032, -0.6980051 , -0.20331326, -0.23335885,
                         -0.61372185, -0.10491576, -0.88596209, -0.61044411, -0.68740006,
                         -0.00136119, -0.58912653, -0.57592882, -0.20858956, -0.18456305,
                         -0.78968331, -0.29210434, -0.65513278, -0.52247678, -0.36671689,
                         -0.93810693, -0.84253114, -0.31654198, -0.05205265, -0.30854419,
                         -0.42011514, -0.71233187, -1.00584317, -0.31220525, -0.69069361,
                         -0.56654897, -0.62510316 };

    floatVector dXdT = { -0.14962987, -0.4307929 , -0.22433599, -0.36650931, -0.49577669,
                         -0.68301613, -0.09246248, -0.09819724 };

    const floatVector *flatdXdF = hydra.getFlatdXdF( );

    const floatVector *flatdXdT = hydra.getFlatdXdT( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dXdF, *flatdXdF ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dXdT, *flatdXdT ) );

    hydraBaseMock hydra_pre( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension,
                             9, 1e-9, 1e-9, 20, 5, 1e-4, true, 0 );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dXdF, *hydra_pre.getFlatdXdF( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dXdT, *hydra_pre.getFlatdXdT( ) ) );
}

BOOST_AUTO_TEST_CASE( test_formPreconditioner ){
    /*!
     * Boost test of the get-configuration command
     */

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

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector expected_preconditioner = { 1, 2, 3, 4, 5, 6, 7 };

            void setPreconditionerType( const unsigned int val ){

                tardigradeHydra::unit_test::hydraBaseTester::set_preconditionerType( *this, val );

            }

            virtual void formMaxRowPreconditioner( ) override{

                tardigradeHydra::unit_test::hydraBaseTester::set_flatPreconditioner( *this, expected_preconditioner );

            }

    };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( hydra.expected_preconditioner, *hydra.getFlatPreconditioner( ) ) );

    hydraBaseMock bad_hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    bad_hydra.setPreconditionerType( 7 );

    BOOST_CHECK_THROW( bad_hydra.getFlatPreconditioner( ), std::nested_exception );

}

BOOST_AUTO_TEST_CASE( test_formMaxRowPreconditioner ){
    /*!
     * Boost test of the get-configuration command
     */

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

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            floatVector jacobian = {  1.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                                      0.        ,  1.        ,  0.        ,  0.        ,  0.        ,
                                      0.        , 48.07641984,  1.        ,  0.        , -7.68935399,
                                      0.        ,  0.        , 18.48297386,  1.        ,  0.        ,
                                      0.        ,  0.        ,  0.        ,  0.        ,  1.        };

            using tardigradeHydra::hydraBase::hydraBase;

            void setPreconditionerType( const unsigned int val ){

                tardigradeHydra::unit_test::hydraBaseTester::set_preconditionerType( *this, val );

            }

            virtual void formMaxRowPreconditioner( ) override{

                tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( *this, floatVector( 5, 0 ) );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, jacobian );

                tardigradeHydra::hydraBase::formMaxRowPreconditioner( );

            }

    };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer = { 1.        , 1.        , 0.02080022, 0.05410385, 1.        };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *hydra.getFlatPreconditioner( ) ) );

}
