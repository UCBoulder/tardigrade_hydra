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

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

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

                static void checkAdditionalDOF( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._additionalDOF == hydra.getAdditionalDOF( ) );

                }

                static void checkPreviousAdditionalDOF( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._previousAdditionalDOF == hydra.getPreviousAdditionalDOF( ) );

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

                static void checkGradientSigma( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._gradientSigma == hydra.getGradientSigma( ) );

                }

                static void checkGradientBeta( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._gradientBeta == hydra.getGradientBeta( ) );

                }

                static void checkMaxGradientIterations( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._maxGradientIterations == hydra.getMaxGradientIterations( ) );

                }

                static void checkGradientRho( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._gradientRho == hydra.getGradientRho( ) );

                }

                static void checkGradientP( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._gradientP == hydra.getGradientP( ) );

                }

                static void checkLMMu( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._lm_mu == hydra.getLMMu( ) );

                }

                static void checkMuk( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._mu_k == hydra.getMuk( ) );

                }

                static void checkUseLevenbergMarquardt( hydraBase &hydra ){

                    BOOST_CHECK( &hydra._use_LM_step == hydra.getUseLevenbergMarquardt( ) );

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

                static void set_flatAdditionalDerivatives( hydraBase &hydra, const floatVector &value ){

                    hydra._additionalDerivatives.second = value;
                    hydra._additionalDerivatives.first = true;

                    hydra.addIterationData( &hydra._additionalDerivatives );
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

                static void set_flatdRdAdditionalDOF( hydraBase &hydra, const floatVector &value ){

                    hydra._dRdAdditionalDOF.second = value;
                    hydra._dRdAdditionalDOF.first = true;

                    hydra.addIterationData( &hydra._dRdAdditionalDOF );
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

                static unsigned int get_gradientIteration( hydraBase &hydra ){

                    return hydra._gradientIteration;

                }

                static void solveNewtonUpdate( hydraBase &hydra, floatVector &deltaX ){

                    hydra.solveNewtonUpdate( deltaX );

                }

                static void performPreconditionedSolve( hydraBase &hydra, floatVector &deltaX ){

                    hydra.performPreconditionedSolve( deltaX );

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

BOOST_AUTO_TEST_CASE( testAbaqusInterface, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

BOOST_AUTO_TEST_CASE( test_hydraBase_getTime, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkTime( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getDeltaTime, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkDeltaTime( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getTemperature, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkTemperature( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousTemperature, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousTemperature( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getDeformationGradient, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkDeformationGradient( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousDeformationGradient, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousDeformationGradient( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getAdditionalDOF, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkAdditionalDOF( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousAdditionalDOF, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousAdditionalDOF( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousStateVariables, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getParameters, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkParameters( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getNumConfigurations, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkNumConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getNumNonLinearSolveStateVariables, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkNumNonLinearSolveStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getDimension, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkDimension( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getSOTDimension, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkSOTDimension( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getTOTDimension, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkTOTDimension( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getFOTDimension, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkFOTDimension( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_configurations, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_previousConfigurations, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_inverseConfigurations, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkInverseConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_previousInverseConfigurations, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousInverseConfigurations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_nonLinearSolveStateVariables, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkNonLinearSolveStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_previousNonLinearSolveStateVariables, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousNonLinearSolveStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_additionalStateVariables, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkAdditionalStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_get_previousAdditionalStateVariables, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkPreviousAdditionalStateVariables( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getRelativeTolerance, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkRelativeTolerance( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getAbsoluteTolerance, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkAbsoluteTolerance( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getLSAlpha, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkLSAlpha( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getGradientSigma, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkGradientSigma( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getGradientBeta, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkGradientBeta( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getMaxGradientIterations, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkMaxGradientIterations( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getGradientRho, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkGradientRho( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getGradientP, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkGradientP( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getLMMu, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkLMMu( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getMuk, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkMuk( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getUseLevenbergMarquardt, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    tardigradeHydra::unit_test::hydraBaseTester::checkUseLevenbergMarquardt( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setGradientSigma, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    hydra.setGradientSigma( 123.4 );

    BOOST_TEST( 123.4 == *hydra.getGradientSigma( ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setGradientBeta, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    hydra.setGradientBeta( 123.4 );

    BOOST_TEST( 123.4 == *hydra.getGradientBeta( ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setMaxGradientIterations, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    hydra.setMaxGradientIterations( 123 );

    BOOST_TEST( 123 == *hydra.getMaxGradientIterations( ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setGradientRho, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    hydra.setGradientRho( 123.4 );

    BOOST_TEST( 123.4 == *hydra.getGradientRho( ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setGradientP, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    hydra.setGradientP( 123.4 );

    BOOST_TEST( 123.4 == *hydra.getGradientP( ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setLMMu, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    hydra.setLMMu( 123.4 );

    BOOST_TEST( 123.4 == *hydra.getLMMu( ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setMuk, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    hydra.setMuk( 123.4 );

    BOOST_TEST( 123.4 == *hydra.getMuk( ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setUseLevenbergMarquardt, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    hydra.setUseLevenbergMarquardt( true );

    BOOST_TEST( true == *hydra.getUseLevenbergMarquardt( ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_decomposeStateVariableVector, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector additionalDOF = { };

    floatVector previousAdditionalDOF = { };

    floatVector previousAdditionalStateVariablesAnswer = { 0.12062867, 0.8263408 , 0.60306013,
                                                           0.54506801, 0.34276383, 0.30412079  };

    floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                      additionalDOF, previousAdditionalDOF,
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::decomposeStateVariableVector( hydra );

    BOOST_TEST( tardigradeVectorTools::appendVectors( configurationsAnswer ) == *hydra.get_configurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousConfigurationsAnswer ) == *hydra.get_previousConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( inverseConfigurationsAnswer ) == *hydra.get_inverseConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousInverseConfigurationsAnswer ) == *hydra.get_previousInverseConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( nonLinearSolveStateVariablesAnswer == *hydra.get_nonLinearSolveStateVariables( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousNonLinearSolveStateVariablesAnswer == *hydra.get_previousNonLinearSolveStateVariables( ), CHECK_PER_ELEMENT );

    BOOST_TEST( additionalStateVariablesAnswer == *hydra.get_additionalStateVariables( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAdditionalStateVariablesAnswer == *hydra.get_previousAdditionalStateVariables( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getSubConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector additionalDOF = { };

    floatVector previousAdditionalDOF = { };

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
                                      additionalDOF, previousAdditionalDOF,
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_TEST( hydra.getSubConfiguration( 0, 4 ) == *hydra.getDeformationGradient( ), CHECK_PER_ELEMENT );

    floatVector answer1 = { 2.24332648, 1.48246714, 2.02801682,
                            1.50380989, 2.98203598, 2.08079721,
                            1.58939152, 1.2551092 , 2.38201794 };

    BOOST_TEST( hydra.getSubConfiguration( 1, 3 ) == answer1, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPrecedingConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector additionalDOF = { };

    floatVector previousAdditionalDOF = { };

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
                                      additionalDOF, previousAdditionalDOF,
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_TEST( hydra.getPrecedingConfiguration( 4 ) == *hydra.getDeformationGradient( ), CHECK_PER_ELEMENT );

    floatVector answer1 = { 0.73947165, -0.24328161, -0.68904986,
                            0.04049558,  0.25403825, -0.1748363 ,
                            0.97015752, -0.04452644, -0.77301275 };

    BOOST_TEST( hydra.getPrecedingConfiguration( 2 ) == answer1, CHECK_PER_ELEMENT );

    floatVector answer2 = { 0.54568475, -0.42501821, -0.54244544,
                           -0.01317666,  0.30165847, -0.09442353,
                            0.80588282, -0.10806097, -0.42143322 };

    BOOST_TEST( hydra.getPrecedingConfiguration( 3 ) == answer2, CHECK_PER_ELEMENT );

    floatVector answer3 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    BOOST_TEST( hydra.getPrecedingConfiguration( 0 ) == answer3, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getFollowingConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector additionalDOF = { };

    floatVector previousAdditionalDOF = { };

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
                                      additionalDOF, previousAdditionalDOF,
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer1 = { 2.09953091, 1.83604029, 2.3712323 ,
                            0.98756433, 2.58928197, 1.05684715,
                            1.33422708, 1.67694162, 2.96443669 };

    BOOST_TEST( hydra.getFollowingConfiguration( 1 ) == answer1, CHECK_PER_ELEMENT );

    floatVector answer2 = { 1.42635131, 0.89338916, 0.94416002,
                            0.50183668, 1.62395295, 0.1156184 ,
                            0.31728548, 0.41482621, 1.86630916 };

    BOOST_TEST( hydra.getFollowingConfiguration( 2 ) == answer2, CHECK_PER_ELEMENT );

    floatVector answer3 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    BOOST_TEST( hydra.getFollowingConfiguration( 3 ) == answer3, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousSubConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
    floatVector additionalDOF = { };

    floatVector previousAdditionalDOF = { };

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
                                      additionalDOF, previousAdditionalDOF,
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_TEST( hydra.getPreviousSubConfiguration( 0, 4 ) == *hydra.getPreviousDeformationGradient( ), CHECK_PER_ELEMENT );

    floatVector answer1 = { 2.24332648, 1.48246714, 2.02801682,
                            1.50380989, 2.98203598, 2.08079721,
                            1.58939152, 1.2551092 , 2.38201794 };

    BOOST_TEST( hydra.getPreviousSubConfiguration( 1, 3 ) == answer1, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousPrecedingConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
    floatVector additionalDOF = { };

    floatVector previousAdditionalDOF = { };

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
                                      additionalDOF, previousAdditionalDOF,
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_TEST( hydra.getPreviousPrecedingConfiguration( 4 ) == *hydra.getPreviousDeformationGradient( ), CHECK_PER_ELEMENT );

    floatVector answer1 = { -0.30350143, -0.21223491,  0.47296395,
                             0.18299993, -0.42974886, -0.06195712,
                             0.92470041, -0.36418391, -0.8287879 };

    BOOST_TEST( hydra.getPreviousPrecedingConfiguration( 2 ) == answer1, CHECK_PER_ELEMENT );

    floatVector answer2 = { -0.15883228, -0.1920217 ,  0.33770597,
                             0.15460279, -0.58876502, -0.15099813,
                             0.69307214, -0.60345639, -0.66103563 };

    BOOST_TEST( hydra.getPreviousPrecedingConfiguration( 3 ) == answer2, CHECK_PER_ELEMENT );

    floatVector answer3 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    BOOST_TEST( hydra.getPreviousPrecedingConfiguration( 0 ) == answer3, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousFollowingConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
    floatVector additionalDOF = { };

    floatVector previousAdditionalDOF = { };

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
                                      additionalDOF, previousAdditionalDOF,
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer1 = { 2.09953091, 1.83604029, 2.3712323 ,
                            0.98756433, 2.58928197, 1.05684715,
                            1.33422708, 1.67694162, 2.96443669 };

    BOOST_TEST( hydra.getPreviousFollowingConfiguration( 1 ) == answer1, CHECK_PER_ELEMENT );

    floatVector answer2 = { 1.42635131, 0.89338916, 0.94416002,
                            0.50183668, 1.62395295, 0.1156184 ,
                            0.31728548, 0.41482621, 1.86630916 };

    BOOST_TEST( hydra.getPreviousFollowingConfiguration( 2 ) == answer2, CHECK_PER_ELEMENT );

    floatVector answer3 = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

    BOOST_TEST( hydra.getPreviousFollowingConfiguration( 3 ) == answer3, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getSubConfigurationJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
    floatVector additionalDOF = { };

    floatVector previousAdditionalDOF = { };

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
                                      additionalDOF, previousAdditionalDOF,
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getSubConfigurationJacobian( configurations, lower, upper ), CHECK_PER_ELEMENT );

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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getSubConfigurationJacobian( configurations, lower, upper ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getSubConfigurationJacobian2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
    floatVector additionalDOF = { };

    floatVector previousAdditionalDOF = { };

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
                                      additionalDOF, previousAdditionalDOF,
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getSubConfigurationJacobian( lower, upper ), CHECK_PER_ELEMENT );

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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getSubConfigurationJacobian( lower, upper ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPrecedingConfigurationJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
    floatVector additionalDOF = { };

    floatVector previousAdditionalDOF = { };

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
                                      additionalDOF, previousAdditionalDOF,
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPrecedingConfigurationJacobian( upper ), CHECK_PER_ELEMENT );

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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPrecedingConfigurationJacobian( upper ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getFollowingConfigurationJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
    floatVector additionalDOF = { };

    floatVector previousAdditionalDOF = { };

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
                                      additionalDOF, previousAdditionalDOF,
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getFollowingConfigurationJacobian( lower ), CHECK_PER_ELEMENT );

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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getFollowingConfigurationJacobian( lower ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousSubConfigurationJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector additionalDOF = { };

    floatVector previousAdditionalDOF = { };

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
                                      additionalDOF, previousAdditionalDOF,
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPreviousSubConfigurationJacobian( lower, upper ), CHECK_PER_ELEMENT );

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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPreviousSubConfigurationJacobian( lower, upper ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousPrecedingConfigurationJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
    floatVector additionalDOF = { };

    floatVector previousAdditionalDOF = { };

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
                                      additionalDOF, previousAdditionalDOF,
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPreviousPrecedingConfigurationJacobian( upper ), CHECK_PER_ELEMENT );

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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPreviousPrecedingConfigurationJacobian( upper ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousFollowingConfigurationJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
    floatVector additionalDOF = { };

    floatVector previousAdditionalDOF = { };

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
                                      additionalDOF, previousAdditionalDOF,
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPreviousFollowingConfigurationJacobian( lower ), CHECK_PER_ELEMENT );

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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPreviousFollowingConfigurationJacobian( lower ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraTest_setFirstConfigurationGradients, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                      { }, { },
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatType eps = 1e-6;

    floatMatrix dF1dF_answer( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatMatrix dF1dFn_answer( deformationGradient.size( ), floatVector( deformationGradient.size( ) * ( numConfigurations - 1 ), 0 ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        tardigradeHydra::hydraBase hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                                            { }, { },
                                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::hydraBase hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                                            { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( dF1dF_answer ) == *hydra.get_dF1dF( ), CHECK_PER_ELEMENT );

    tardigradeHydra::unit_test::hydraBaseTester::checkdF1dF( hydra );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ) * deformationGradient.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] += eps * std::fabs( previousStateVariables[ i ] ) + eps;

        tardigradeHydra::hydraBase hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                            { }, { },
                                            previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::hydraBase hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                            { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( dF1dFn_answer ) == *hydra.get_dF1dFn( ), CHECK_PER_ELEMENT );

    tardigradeHydra::unit_test::hydraBaseTester::checkdF1dFn( hydra );

}

BOOST_AUTO_TEST_CASE( test_hydraTest_setPreviousFirstConfigurationGradients, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                      { }, { },
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatType eps = 1e-6;

    floatMatrix previousdF1dF_answer( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatMatrix previousdF1dFn_answer( deformationGradient.size( ), floatVector( deformationGradient.size( ) * ( numConfigurations - 1 ), 0 ) );

    for ( unsigned int i = 0; i < previousDeformationGradient.size( ); i++ ){

        floatVector delta( previousDeformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        tardigradeHydra::hydraBase hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                                            { }, { },
                                            previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::hydraBase hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                                            { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdF1dF_answer ) == *hydra.get_previousdF1dF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ) * deformationGradient.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] += eps * std::fabs( previousStateVariables[ i ] ) + eps;

        tardigradeHydra::hydraBase hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                            { }, { },
                                            previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::hydraBase hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                            { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdF1dFn_answer ) == *hydra.get_previousdF1dFn( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residualBase_residualBase, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    tardigradeHydra::residualBase residual( &hydra, numEquations );

    BOOST_CHECK( residual.hydra == &hydra );

    BOOST_CHECK( *residual.getNumEquations( ) == numEquations );

}

BOOST_AUTO_TEST_CASE( test_residualBase_checkDefaults, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    unsigned int numEquations = 3;

    tardigradeHydra::residualBase residual( &hydra, numEquations );

    BOOST_CHECK_THROW( residual.setResidual( ), std::nested_exception );

    BOOST_CHECK_THROW( residual.setJacobian( ), std::nested_exception );

    BOOST_CHECK_THROW( residual.setdRdF( ), std::nested_exception );

    BOOST_CHECK_THROW( residual.setdRdT( ), std::nested_exception );

    BOOST_CHECK_NO_THROW( residual.setAdditionalDerivatives( ) );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setResidual, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    BOOST_TEST( *residual.getResidual( ) == residual.residual, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    BOOST_TEST( *residual.getJacobian( ) == residual.jacobian, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setdRdF, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    BOOST_TEST( *residual.getdRdF( ) == residual.dRdF, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setdRdT, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    BOOST_TEST( *residual.getdRdT( ) == residual.dRdT, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setAdditionalDerivatives, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    BOOST_TEST( *residual.getAdditionalDerivatives( ) == residual.additionalDerivatives, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setStress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    BOOST_TEST( *residual.getStress( ) == residual.cauchyStress, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setPreviousStress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    BOOST_TEST( *residual.getPreviousStress( ) == residual.previousCauchyStress, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residualBase_setCurrentAdditionalStateVariables, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    BOOST_TEST( *residual.getCurrentAdditionalStateVariables( ) == residual.currentAdditionalStateVariables, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setResidualClasses, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_CHECK_NO_THROW( hydra.setResidualClasses( ) );

    BOOST_CHECK( *( *hydra.getResidualClasses( ) )[ 0 ]->getNumEquations( ) == hydra.s1 );

    BOOST_CHECK( *( *hydra.getResidualClasses( ) )[ 1 ]->getNumEquations( ) == hydra.s2 );

    BOOST_CHECK( *( *hydra.getResidualClasses( ) )[ 2 ]->getNumEquations( ) == hydra.s3 );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setResidualClasses2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_CHECK_THROW( hydra.setResidualClasses( ), std::nested_exception );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_formNonLinearProblem, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualBaseMock : public tardigradeHydra::residualBase{

        public:

            unsigned int numVariables = 41;

            using tardigradeHydra::residualBase::residualBase;

            using tardigradeHydra::residualBase::setResidual;

            using tardigradeHydra::residualBase::setJacobian;

            using tardigradeHydra::residualBase::setdRdF;

            using tardigradeHydra::residualBase::setdRdT;

            using tardigradeHydra::residualBase::setdRdAdditionalDOF;

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

            virtual void setdRdAdditionalDOF( ){

                floatVector dRdAdditionalDOF( ( *getNumEquations( ) ) * ( hydra->getAdditionalDOF( )->size( ) ), 0 );

                for ( unsigned int i = 0; i < *getNumEquations( ); i++ ){

                    for ( unsigned int j = 0; j < hydra->getAdditionalDOF( )->size( ); j++ ){

                        dRdAdditionalDOF[ hydra->getAdditionalDOF( )->size( ) * i + j ] = 0.6 * i - 0.1 * j;

                    }

                }

                setdRdAdditionalDOF( dRdAdditionalDOF );

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33 };

    floatVector previousAdditionalDOF = { 0.111, 0.222, 0.333 };

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
                         additionalDOF, additionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    hydraBaseMock hydraGet( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                            additionalDOF, additionalDOF,
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

    floatVector dRdAdditionalDOFAnswer = tardigradeVectorTools::appendVectors( { *hydraGet.r1.getdRdAdditionalDOF( ),
                                                                                 *hydraGet.r2.getdRdAdditionalDOF( ),
                                                                                 *hydraGet.r3.getdRdAdditionalDOF( ) } );

    floatMatrix additionalDerivativesAnswer = tardigradeVectorTools::appendVectors( { tardigradeVectorTools::inflate( *hydraGet.r1.getAdditionalDerivatives( ), 36, 4 ),
                                                                                      tardigradeVectorTools::inflate( *hydraGet.r2.getAdditionalDerivatives( ),  2, 4 ),
                                                                                      tardigradeVectorTools::inflate( *hydraGet.r3.getAdditionalDerivatives( ),  3, 4 ) } );

    tardigradeHydra::unit_test::hydraBaseTester::formNonLinearProblem( hydra );

    BOOST_TEST( residualAnswer == *hydra.getResidual( ), CHECK_PER_ELEMENT );

    BOOST_CHECK( jacobianAnswer.size( ) == hydra.getJacobian( ).size( ) );

    for ( unsigned int i = 0; i < jacobianAnswer.size( ); i++ ){

        BOOST_TEST( jacobianAnswer[ i ] == hydra.getJacobian( )[ i ], CHECK_PER_ELEMENT );

    }

    BOOST_CHECK( dRdFAnswer.size( ) == hydra.getdRdF( ).size( ) );

    for ( unsigned int i = 0; i < dRdFAnswer.size( ); i++ ){

        BOOST_TEST( dRdFAnswer[ i ] == hydra.getdRdF( )[ i ], CHECK_PER_ELEMENT );

    }

    BOOST_TEST( dRdTAnswer == *hydra.getdRdT( ), CHECK_PER_ELEMENT );


    BOOST_CHECK( additionalDerivativesAnswer.size( ) == hydra.getAdditionalDerivatives( ).size( ) );

    for ( unsigned int i = 0; i < additionalDerivativesAnswer.size( ); i++ ){

        BOOST_TEST( additionalDerivativesAnswer[ i ] == hydra.getAdditionalDerivatives( )[ i ], CHECK_PER_ELEMENT );

    }

}

BOOST_AUTO_TEST_CASE( test_hydraBase_initializeUnknownVector, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    hydra.setResidualClasses( );

    tardigradeHydra::unit_test::hydraBaseTester::initializeUnknownVector( hydra );

    floatVector unknownVectorAnswer = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                        1.53155137, 0.53182759, 0.63440096, 0.84943179, 1.72445532, 0.61102351, 0.72244338, 0.32295891, 1.36178866,
                                        1.22826323, 0.29371405, 0.63097612, 0.09210494, 1.43370117, 0.43086276, 0.4936851 , 0.42583029, 1.31226122,
                                        1.42635131, 0.89338916, 0.94416002, 0.50183668, 1.62395295, 0.1156184 , 0.31728548, 0.41482621, 1.86630916,
                                        0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453 };
                                         

    BOOST_TEST( unknownVectorAnswer == *hydra.getUnknownVector( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_initializeUnknownVector_2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    class residualBaseMockSuggest : public residualBaseMock{

        public:

            unsigned int config_suggest = 1;

            using residualBaseMock::residualBaseMock;

            virtual void suggestInitialIterateValues( std::vector< unsigned int > &indices, std::vector< floatType > &values ) override{

                indices = std::vector< unsigned int >( 5, 0 );

                values = std::vector< floatType >( 5, 0 );

                for ( unsigned int i = 0; i < indices.size( ); i++ ){

                    indices[ i ] = 9 * config_suggest + i;

                    values[ i ] = i;

                }

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
        
            residualBaseMockSuggest r3;

            unsigned int s1 = 36;

            unsigned int s2 = 2;

            unsigned int s3 = 3;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                r1 = residualBaseMockStress( this, s1 );

                r2 = residualBaseMock( this, s2 );

                r3 = residualBaseMockSuggest( this, s3 );

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    hydra.setResidualClasses( );

    tardigradeHydra::unit_test::hydraBaseTester::initializeUnknownVector( hydra );

    floatVector unknownVectorAnswer = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                        0,          1,          2,          3,          4,          0.61102351, 0.72244338, 0.32295891, 1.36178866,
                                        1.22826323, 0.29371405, 0.63097612, 0.09210494, 1.43370117, 0.43086276, 0.4936851 , 0.42583029, 1.31226122,
                                        1.42635131, 0.89338916, 0.94416002, 0.50183668, 1.62395295, 0.1156184 , 0.31728548, 0.41482621, 1.86630916,
                                        0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453 };
                                         
    BOOST_TEST( unknownVectorAnswer == *hydra.getUnknownVector( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setTolerance, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector residual = { 1, 2, -3, 0 };

    floatVector unknownVector = { -2, 5, 10, 0.3 };

    floatVector toleranceAnswer = 1e-9 * (tardigradeVectorTools::abs( residual ) + tardigradeVectorTools::abs( unknownVector )) + 1e-9;

    tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra, unknownVector );

    BOOST_TEST( *hydra.getTolerance( ) == toleranceAnswer, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_checkConvergence, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
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

BOOST_AUTO_TEST_CASE( test_hydraBase_getLSResidualNorm, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::hydraBase hydra;

    floatVector residual = { 1, 2, 3 };

    floatType lsResidualNormAnswer = tardigradeVectorTools::l2norm( residual );
      tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    BOOST_TEST( *hydra.getLSResidualNorm( ) == lsResidualNormAnswer );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_checkLSConvergence, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector residual = { 1, 2, -3, 0 };

    tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    BOOST_CHECK( !hydra.checkLSConvergence( ) );

    residual = { 0.1, 2, -3, 0 };

    tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    BOOST_CHECK( hydra.checkLSConvergence( ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getPreviousStress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){


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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    BOOST_TEST( answer == *hydra.getPreviousStress( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_decomposeUnknownVector, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
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

    BOOST_TEST( *hydra.getStress( ) == cauchyStressAnswer, CHECK_PER_ELEMENT );

    BOOST_TEST( *hydra.get_configurations( ) == tardigradeVectorTools::appendVectors( configurationsAnswer ), CHECK_PER_ELEMENT );

    BOOST_TEST( *hydra.get_nonLinearSolveStateVariables( ) == isvAnswer, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_solveNewtonUpdate, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector flatJacobian = { 1, 2, 3, 4, 5, 6, 7, 8, 2 };

            floatVector residual = { 1, 2, 3 };

            virtual const unsigned int getNumUnknowns( ) override{ return residual.size( ); }

        protected:

            virtual void initializeUnknownVector( ){

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, flatJacobian );

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer = { 1./3, -2./3, 0 };

    floatVector result( 3, 0 );

    tardigradeHydra::unit_test::hydraBaseTester::initializeUnknownVector( hydra );
    tardigradeHydra::unit_test::hydraBaseTester::solveNewtonUpdate( hydra, result );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    result = { 0, 0, 0 };

    hydraBaseMock hydra_pre( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension,
                             9, 1e-9, 1e-9, 20, 5, 1e-4, true, 0 );

    tardigradeHydra::unit_test::hydraBaseTester::initializeUnknownVector( hydra_pre );
    tardigradeHydra::unit_test::hydraBaseTester::solveNewtonUpdate( hydra_pre, result );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_performPreconditionedSolve, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector flatJacobian = { 1, 2, 3, 4, 5, 6, 7, 8, 2 };

            floatVector residual = { 1, 2, 3 };

            virtual const unsigned int getNumUnknowns( ) override{ return residual.size( ); }

        protected:

            virtual void initializeUnknownVector( ){

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, flatJacobian );

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer = { 1./3, -2./3, 0 };

    floatVector result( 3, 0 );

    tardigradeHydra::unit_test::hydraBaseTester::initializeUnknownVector( hydra );
    tardigradeHydra::unit_test::hydraBaseTester::performPreconditionedSolve( hydra, result );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_solveNonLinearProblem, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

            std::vector< floatType > expectedBaseResidualNorms = { 14, 77, 245, 110 };

            std::vector< floatVector > expectedBasedResidualNormdXs = {
                                                                         {   2.97757248,  -1.95288068,  -3.03183756 },
                                                                         { -16.10008822,  -6.7048463 ,  14.58211116 },
                                                                         { -33.99029958,  28.301765  , -20.74369594 },
                                                                         {  -4.96839724,  18.87890456,   2.15863168 }
                                                                      };

            std::vector< floatType > expectedMuk = { 7e-8, 7e-8, 7e-8, 7e-8 };

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
                                                                               { 33.25014131, -2.84970064, 20.95677101 },
                                                                               { 31.64013249, -3.52018527, 22.41498213 }
                                                                             },
                                                                             { 
                                                                               { 44.6038325 , -1.77731706, 27.17238891 }
                                                                             },
                                                                             { 
                                                                               { 420.89746176,  81.81218142,  60.29432705 },
                                                                               { 232.75064713,  40.01743218,  43.73335798 },
                                                                               { 138.67723982,  19.12005756,  35.45287345 }
                                                                             },
                                                                         };

        virtual const unsigned int getNumUnknowns( ) override{ return initialUnknownVector.size( ); }

        unsigned int in_gradient_convergence = 0;

        private:

            using tardigradeHydra::hydraBase::getResidual;

            virtual void initializeUnknownVector( ){

                tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( *this, initialUnknownVector );

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual[ 0 ][ 0 ] );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, flatJacobian[ 0 ][ 0 ] );

            }

            virtual bool checkConvergence( ) override{

                unsigned int iteration = tardigradeHydra::unit_test::hydraBaseTester::get_iteration( *this );

                if ( iteration < residual.size( ) ){

                    return false;

                }

                return true;

            }

            virtual bool checkLSConvergence( ) override{

                unsigned int iteration = tardigradeHydra::unit_test::hydraBaseTester::get_iteration( *this );

                unsigned int LSIteration = tardigradeHydra::unit_test::hydraBaseTester::get_LSIteration( *this );

                getResidual( );

                if ( LSIteration < residual[ iteration ].size( ) - 1 ){

                    return false;

                }

                return true;

            }

            virtual bool checkDescentDirection( const floatVector &dx ) override{

                unsigned int iteration = tardigradeHydra::unit_test::hydraBaseTester::get_iteration( *this );

                if ( iteration == 1 ){

                    return false;

                }

                return true;

            }

            virtual void performGradientStep( const floatVector &X0 ) override{

                in_gradient_convergence = 1;

                tardigradeHydra::hydraBase::performGradientStep( X0 );

            }

            virtual bool checkGradientConvergence( const floatVector &X0 ) override{

                unsigned int iteration = tardigradeHydra::unit_test::hydraBaseTester::get_iteration( *this );

                unsigned int gradientIteration = tardigradeHydra::unit_test::hydraBaseTester::get_gradientIteration( *this );

                in_gradient_convergence = 1;

                getResidual( );

                if ( gradientIteration < residual[ iteration ].size( ) - 1 ){

                    return false;

                }

                in_gradient_convergence = 0;

                return true;

            }

            virtual void formNonLinearProblem( ) override{

                unsigned int iteration = tardigradeHydra::unit_test::hydraBaseTester::get_iteration( *this );

                unsigned int LSIteration = tardigradeHydra::unit_test::hydraBaseTester::get_LSIteration( *this );

                unsigned int gradIteration = tardigradeHydra::unit_test::hydraBaseTester::get_gradientIteration( *this );

                unsigned int iterationOffset = 0;

                unsigned int LSoffset = 0;

                if ( ( LSIteration < residual[ iteration ].size( ) - 1 ) && ( gradIteration < residual[ iteration ].size( ) - 2 ) ){

                    LSoffset += 1;

                }
                else if ( iteration < residual.size( ) - 1 ){

                    iterationOffset += 1;
                    LSIteration = 0;
                    gradIteration = 0;

                }

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual[ iteration + iterationOffset ][ LSIteration + LSoffset + gradIteration ] );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, flatJacobian[ iteration + iterationOffset ][ LSIteration + LSoffset + gradIteration ] );

            }

            virtual void updateUnknownVector( const floatVector &newUnknownVector ) override{

                tardigradeHydra::unit_test::hydraBaseTester::resetIterationData( *this );

                unsigned int iteration = tardigradeHydra::unit_test::hydraBaseTester::get_iteration( *this );

                unsigned int LSIteration = tardigradeHydra::unit_test::hydraBaseTester::get_LSIteration( *this );

                unsigned int gradIteration = tardigradeHydra::unit_test::hydraBaseTester::get_gradientIteration( *this );

                unsigned int subIteration = LSIteration + gradIteration + in_gradient_convergence;

                BOOST_TEST( expectedXVectors[ iteration ][ subIteration ] == newUnknownVector, CHECK_PER_ELEMENT );

                BOOST_TEST( expectedBaseResidualNorms[ iteration ] == *get_baseResidualNorm( ) );

                BOOST_TEST( expectedBasedResidualNormdXs[ iteration ] == *get_basedResidualNormdX( ), CHECK_PER_ELEMENT );

                BOOST_TEST( expectedMuk[ iteration ] == *getMuk( ) );

                tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( *this, newUnknownVector );

                in_gradient_convergence = 0;

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::solveNonLinearProblem( hydra );

    BOOST_TEST( hydra.getNumNewton( ) == 2 );

    BOOST_TEST( hydra.getNumLS( ) == 1 );

    BOOST_TEST( hydra.getNumGrad( ) == 1 );

    hydraBaseMock hydra_pre( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension,
                             9, 1e-9, 1e-9, 20, 5, 1e-4, true, 0 );

    tardigradeHydra::unit_test::hydraBaseTester::solveNonLinearProblem( hydra_pre );

    BOOST_TEST( hydra_pre.getNumNewton( ) == 2 );

    BOOST_TEST( hydra_pre.getNumLS( ) == 1 );

    BOOST_TEST( hydra_pre.getNumGrad( ) == 1 );

}

BOOST_AUTO_TEST_CASE( test_residual_setdRdF, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_CHECK_NO_THROW( hydra.evaluate( ) );

    BOOST_TEST( PK2Answer == *hydra.elasticity.get_PK2Stress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( cauchyStressAnswer == *hydra.getStress( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_getConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
                                      { }, { },
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_TEST( hydra.getConfiguration( 1 ) == floatVector( hydra.get_configurations( )->begin( ) + 1 * 9,
                                                            hydra.get_configurations( )->begin( ) + 2 * 9 ), CHECK_PER_ELEMENT );

    BOOST_TEST( hydra.getPreviousConfiguration( 3 ) == floatVector( hydra.get_previousConfigurations( )->begin( ) + 3 * 9,
                                                                    hydra.get_previousConfigurations( )->begin( ) + 4 * 9 ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_computeTangents, * boost::unit_test::tolerance( 2e-6 ) ){
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

            virtual const unsigned int getNumUnknowns( ) override{ return residual.size( ); }

    };
    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector dXdF = { -0.82465846, -0.07170032, -0.6980051 , -0.20331326, -0.23335885,
                         -0.61372185, -0.10491576, -0.88596209, -0.61044411, -0.68740006,
                         -0.0013611877, -0.58912653, -0.57592882, -0.20858956, -0.18456305,
                         -0.78968331, -0.29210434, -0.65513278, -0.52247678, -0.36671689,
                         -0.93810693, -0.84253114, -0.31654198, -0.05205265, -0.30854419,
                         -0.42011514, -0.71233187, -1.00584317, -0.31220525, -0.69069361,
                         -0.56654897, -0.62510316 };

    floatVector dXdT = { -0.14962987, -0.4307929 , -0.22433599, -0.36650931, -0.49577669,
                         -0.68301613, -0.09246248, -0.09819724 };

    const floatVector *flatdXdF = hydra.getFlatdXdF( );

    const floatVector *flatdXdT = hydra.getFlatdXdT( );

    BOOST_TEST( dXdF == *flatdXdF, CHECK_PER_ELEMENT );

    BOOST_TEST( dXdT == *flatdXdT, CHECK_PER_ELEMENT );

    hydraBaseMock hydra_pre( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension,
                             9, 1e-9, 1e-9, 20, 5, 1e-4, true, 0 );

    BOOST_TEST( dXdF == *hydra_pre.getFlatdXdF( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dXdT == *hydra_pre.getFlatdXdT( ), CHECK_PER_ELEMENT );
}

BOOST_AUTO_TEST_CASE( test_computeFlatdXdAdditionalDOF, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

            floatVector dRdAdditionalDOF = {  0.39293837, -0.42772133, -0.54629709,  0.10262954,  0.43893794,
                                             -0.15378708,  0.9615284 ,  0.36965948, -0.0381362 , -0.21576496,
                                             -0.31364397,  0.45809941, -0.12285551, -0.88064421, -0.20391149,
                                              0.47599081, -0.63501654, -0.64909649,  0.06310275,  0.06365517,
                                              0.26880192,  0.69886359,  0.44891065,  0.22204702,  0.44488677,
                                             -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225,
                                             -0.81579012, -0.13259765, -0.13827447, -0.0126298 , -0.14833942,
                                             -0.37547755, -0.14729739,  0.78677833,  0.88832004,  0.00367335 };

            virtual void formNonLinearProblem( ) override{

                tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( *this, residual );

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, jacobian );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatdRdAdditionalDOF( *this, dRdAdditionalDOF );

            }

            virtual const unsigned int getNumUnknowns( ) override{ return residual.size( ); }

    };
    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { 1, 2, 3, 4, 5 }, { 1, 2, 3, 4, 5 },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector dXdAdditionalDOF = { -0.4084987 ,  0.39168929,  0.5515543 , -0.05139529, -0.41998415,
                                      0.18443685, -0.98977328, -0.34171132,  0.09720518,  0.21059288,
                                      0.27548634, -0.48194765,  0.11825706,  0.87472897,  0.20725542,
                                     -0.44595879,  0.60609488,  0.63118768, -0.06644165, -0.04874605,
                                     -0.33283885, -0.70895038, -0.46424604, -0.2398367 , -0.44833713,
                                      0.37340906,  0.25434963,  0.60171606,  0.49794981, -0.2340264 ,
                                      0.86191456,  0.11588005,  0.11268818,  0.02929378,  0.16441709,
                                      0.34495751,  0.19107884, -0.76717361, -0.90895772, -0.01071376 };

    const floatVector *flatdXdAdditionalDOF = hydra.getFlatdXdAdditionalDOF( );

    BOOST_TEST( dXdAdditionalDOF == *flatdXdAdditionalDOF, CHECK_PER_ELEMENT );

    hydraBaseMock hydra_pre( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             { 1, 2, 3, 4, 5 }, { 1, 2, 3, 4, 5 },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension,
                             9, 1e-9, 1e-9, 20, 5, 1e-4, true, 0 );

    BOOST_TEST( dXdAdditionalDOF == *hydra_pre.getFlatdXdAdditionalDOF( ), CHECK_PER_ELEMENT );

}


BOOST_AUTO_TEST_CASE( test_formPreconditioner, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    BOOST_TEST( hydra.expected_preconditioner == *hydra.getFlatPreconditioner( ), CHECK_PER_ELEMENT );

    hydraBaseMock bad_hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    bad_hydra.setPreconditionerType( 7 );

    BOOST_CHECK_THROW( bad_hydra.getFlatPreconditioner( ), std::nested_exception );

}

BOOST_AUTO_TEST_CASE( test_formMaxRowPreconditioner, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

            virtual const unsigned int getNumUnknowns( ) override{ return 5; }

    };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer = { 1.        , 1.        , 0.02080022, 0.05410385, 1.        };

    BOOST_TEST( answer == *hydra.getFlatPreconditioner( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setResidualNorm, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the norm of the residual
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

            floatVector jacobian = { -0.15378708,  0.9615284 ,  0.36965948, -0.0381362 , -0.21576496,
                                     -0.31364397,  0.45809941, -0.12285551, -0.88064421, -0.20391149,
                                      0.47599081, -0.63501654, -0.64909649,  0.06310275,  0.06365517,
                                      0.26880192,  0.69886359,  0.44891065,  0.22204702,  0.44488677,
                                     -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225 };

            using tardigradeHydra::hydraBase::hydraBase;

            virtual void formNonLinearProblem( ) override{

                floatVector residual( 5, 0 );

                const floatVector *X = getUnknownVector( );

                for ( unsigned int i = 0; i < 5; i++ ){
                    for ( unsigned int j = 0; j < 5; j++ ){
                        residual[ i ] += jacobian[5*i+j] * ( *X )[j];
                    }
                }

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, jacobian );

            }

            virtual const unsigned int getNumUnknowns( ) override{ return 5; }

    };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector unknownVector = { 0.39293837, -0.42772133, -0.54629709,  0.10262954,  0.43893794 };

    floatType answer = 1.6716509825117496;

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra, unknownVector );

    BOOST_TEST( answer == *hydra.get_residualNorm( ) );

    floatVector dResidualNormdX( 5, 0 );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 5; i++ ){

        floatType delta = eps * std::fabs( unknownVector[ i ] ) + eps;

        floatVector xp = unknownVector;

        floatVector xm = unknownVector;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydrap, xp );

        tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydram, xm );

        dResidualNormdX[ i ] = ( *hydrap.get_residualNorm( ) - *hydram.get_residualNorm( ) ) / ( 2 * delta );

    }

    BOOST_TEST( dResidualNormdX == *hydra.get_dResidualNormdX( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_checkGradientConvergence, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test checking the gradient convergence
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

            floatVector X0 = { -0.81579012, -0.13259765, -0.13827447, -0.0126298 , -0.14833942 };

            floatType baseResidualNorm = 0.4459139462561169;

            floatVector basedResidualNormdX = { -0.86442794, -0.34410741, -0.58249594, -0.97271835, -0.32478706 };

            floatVector jacobian = { -0.15378708,  0.9615284 ,  0.36965948, -0.0381362 , -0.21576496,
                                     -0.31364397,  0.45809941, -0.12285551, -0.88064421, -0.20391149,
                                      0.47599081, -0.63501654, -0.64909649,  0.06310275,  0.06365517,
                                      0.26880192,  0.69886359,  0.44891065,  0.22204702,  0.44488677,
                                     -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225 };

            using tardigradeHydra::hydraBase::hydraBase;

            virtual void mockInitialize( ){

                set_baseResidualNorm( baseResidualNorm );

                set_basedResidualNormdX( basedResidualNormdX );

            }

            virtual void formNonLinearProblem( ) override{

                floatVector residual( 5, 0 );

                const floatVector *X = getUnknownVector( );

                for ( unsigned int i = 0; i < 5; i++ ){
                    for ( unsigned int j = 0; j < 5; j++ ){
                        residual[ i ] += jacobian[5*i+j] * ( *X )[j];
                    }
                }

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, jacobian );

            }

            virtual const unsigned int getNumUnknowns( ) override{ return 5; }

    };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector unknownVector = { 0.39293837, -0.42772133, -0.54629709,  0.10262954,  0.43893794 };

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra, unknownVector );

    hydra.mockInitialize( );

    BOOST_TEST( !hydra.checkGradientConvergence( hydra.X0 ) );

    hydraBaseMock hydra2( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra2, 0.5 * hydra.X0 );

    hydra2.mockInitialize( );

    BOOST_TEST( hydra2.checkGradientConvergence( hydra2.X0 ) );

}

BOOST_AUTO_TEST_CASE( test_performGradientStep, * boost::unit_test::tolerance( 1e-5 ) ){
    /*!
     * Test checking the gradient convergence
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

            floatVector X0 = { -0.81579012, -0.13259765, -0.13827447, -0.0126298 , -0.14833942 };

            floatType baseResidualNorm = 0.2408779076031648;

            floatVector basedResidualNormdX = { -1.17899799,  0.07843952, -0.01708813, -0.01779959, -0.06410942 };

            floatVector A = { -0.15378708,  0.9615284 ,  0.36965948, -0.0381362 , -0.21576496,
                              -0.31364397,  0.45809941, -0.12285551, -0.88064421, -0.20391149,
                               0.47599081, -0.63501654, -0.64909649,  0.06310275,  0.06365517,
                               0.26880192,  0.69886359,  0.44891065,  0.22204702,  0.44488677,
                              -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225 };

            using tardigradeHydra::hydraBase::hydraBase;

            virtual void mockInitialize( ){

                set_baseResidualNorm( baseResidualNorm );

                set_basedResidualNormdX( basedResidualNormdX );

            }

            virtual void runGradientStep( ){

                performGradientStep( X0 );

            }

            virtual void formNonLinearProblem( ) override{

                floatVector residual( 5, 0 );

                floatVector jacobian( 25, 0 );

                const floatVector *X = getUnknownVector( );

                for ( unsigned int i = 0; i < 5; i++ ){
                    for ( unsigned int j = 0; j < 5; j++ ){
                        residual[ i ] += A[ 5 * i + j ] * ( ( *X )[ j ] * ( *X )[ j ] );
                        jacobian[ 5 * i + j ] += 2 * A[ 5 * i + j ] * ( *X )[ j ];
                    }
                }

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, jacobian );

            }

            virtual void decomposeUnknownVector( ) override{ return; }
            virtual const unsigned int getNumUnknowns( ) override{ return 5; }

    };

    floatVector answer = { 0.36320787, -0.21103717, -0.12118634,  0.00516978, -0.08423 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector unknownVector = { 0.39293837, -0.42772133, -0.54629709,  0.10262954,  0.43893794 };

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra, unknownVector );

    hydra.mockInitialize( );

    hydra.runGradientStep( );

    BOOST_TEST( answer == *hydra.getUnknownVector( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_checkDescentDirection, * boost::unit_test::tolerance( 1e-5 ) ){
    /*!
     * Test checking the descent direction
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

            floatType residualNorm = 0.2408779076031648;

            floatVector dResidualNormdX = { -1.17899799,  0.07843952, -0.01708813, -0.01779959, -0.06410942 };

            using tardigradeHydra::hydraBase::hydraBase;

            virtual void mockInitialize( ){

                set_residualNorm( residualNorm );

                set_dResidualNormdX( dResidualNormdX );

            }

            virtual bool runCheckDescentDirection( const floatVector &dx ){

                return checkDescentDirection( dx );

            }

            virtual const unsigned int getNumUnknowns( ) override{ return 5; }

    };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector dx = hydra.dResidualNormdX;

    hydra.mockInitialize( );

    BOOST_TEST( !hydra.runCheckDescentDirection( dx ) );

    dx = -hydra.dResidualNormdX;

    BOOST_TEST( hydra.runCheckDescentDirection( dx ) );

}

BOOST_AUTO_TEST_CASE( test_getNonlinearLMTerms, * boost::unit_test::tolerance( 1e-5 ) ){
    /*!
     * Test checking the gradient convergence
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

            floatVector residual = { 1, 2, 3, 4, 5 };

            floatVector jacobian = { -0.15378708,  0.9615284 ,  0.36965948, -0.0381362 , -0.21576496,
                                     -0.31364397,  0.45809941, -0.12285551, -0.88064421, -0.20391149,
                                      0.47599081, -0.63501654, -0.64909649,  0.06310275,  0.06365517,
                                      0.26880192,  0.69886359,  0.44891065,  0.22204702,  0.44488677,
                                     -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225 };

            floatType mu_k = 1.34;

            using tardigradeHydra::hydraBase::hydraBase;

            virtual void formNonLinearProblem( ) override{

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, jacobian );

                setMuk( mu_k );

            }

            virtual void decomposeUnknownVector( ) override{ return; }
            virtual const unsigned int getNumUnknowns( ) override{ return 5; }

    };

    floatVector answerRHS = { -0.04830576,  1.38601851, -2.74506611, -2.78478784,  2.6566859 };

    floatVector answerLHS = {  1.88621891, -0.30808058, -0.01417759,  0.51788095,  0.15427055,
                              -0.30808058,  3.44245776,  1.17530079, -0.21093811, -0.10279234,
                              -0.01417759,  1.17530079,  2.40995212,  0.377036  , -0.03867597,
                               0.51788095, -0.21093811,  0.377036  ,  2.34049101,  0.18253039,
                               0.15427055, -0.10279234, -0.03867597,  0.18253039,  1.69872961 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector unknownVector = { 0.39293837, -0.42772133, -0.54629709,  0.10262954,  0.43893794 };

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra, unknownVector );

    hydra.setUseLevenbergMarquardt( false );

    BOOST_TEST( hydra.residual == *hydra.getNonlinearRHS( ), CHECK_PER_ELEMENT );

    BOOST_TEST( hydra.jacobian == *hydra.getFlatNonlinearLHS( ), CHECK_PER_ELEMENT );

    hydra.setUseLevenbergMarquardt( true );

    BOOST_TEST( answerRHS == *hydra.getNonlinearRHS( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answerLHS == *hydra.getFlatNonlinearLHS( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setBaseQuantities, * boost::unit_test::tolerance( 1e-5 ) ){
    /*!
     * Test checking the gradient convergence
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

            floatVector residual = { 1, 2, 3, 4, 5 };

            floatVector jacobian = { -0.15378708,  0.9615284 ,  0.36965948, -0.0381362 , -0.21576496,
                                     -0.31364397,  0.45809941, -0.12285551, -0.88064421, -0.20391149,
                                      0.47599081, -0.63501654, -0.64909649,  0.06310275,  0.06365517,
                                      0.26880192,  0.69886359,  0.44891065,  0.22204702,  0.44488677,
                                     -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225 };

            floatType mu_k = 1.34;

            using tardigradeHydra::hydraBase::hydraBase;

            floatType rnorm = 10.3;

            floatVector dRNormdX = { 1, 2, 3 };

            void runSetBaseQuantities( ){ setBaseQuantities( ); }

            const floatType *getBaseResidualNorm( ){ return get_baseResidualNorm( ); }

            const floatVector *getBasedResidualNormdX( ){ return get_basedResidualNormdX( ); }

            virtual void setResidualNorm( ) override{

                set_residualNorm( rnorm );

            }

            virtual void setdResidualNormdX( ) override{

                set_dResidualNormdX( dRNormdX );;

            }

            virtual void formNonLinearProblem( ) override{

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, jacobian );

                setMuk( mu_k );

            }

            virtual void decomposeUnknownVector( ) override{ return; }
            virtual const unsigned int getNumUnknowns( ) override{ return 5; }

    };

    floatType answer1 = 0.5 * 1e-8 * 10.3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    hydra.runSetBaseQuantities( );

    BOOST_TEST( answer1 == *hydra.getMuk( ) );

    BOOST_TEST( hydra.rnorm == *hydra.getBaseResidualNorm( ) );

    BOOST_TEST( hydra.dRNormdX == *hydra.getBasedResidualNormdX( ) );

    hydra.rnorm = 1e-9;

    hydra.setResidualNorm( );

    hydra.runSetBaseQuantities( );

    BOOST_TEST( hydra.rnorm == *hydra.getMuk( ) );

}
