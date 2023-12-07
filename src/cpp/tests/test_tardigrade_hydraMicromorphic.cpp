/**
  * \file test_tardigrade_hydraMicromorphic.cpp
  *
  * Tests for tardigrade_hydraMicromorphic
  */

#include<tardigrade_hydraMicromorphic.h>

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

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_constructor ){

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    floatVector previousDeformationGradient = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

    floatVector microDeformation = { 19, 20, 21, 22, 23, 24, 25, 26, 27 };

    floatVector previousMicroDeformation = { 28, 29, 30, 31, 32, 33, 34, 35, 36 };

    floatVector gradientMicroDeformation = { 36, 37, 38, 39, 40, 41, 42, 43, 44,
                                             45, 46, 47, 48, 49, 50, 51, 52, 53,
                                             54, 55, 56, 57, 58, 59, 60, 61, 62 };

    floatVector previousGradientMicroDeformation = { 63, 64, 65, 66, 67, 68, 69, 70, 71,
                                                     72, 73, 74, 75, 76, 77, 78, 79, 80,
                                                     81, 82, 83, 84, 85, 86, 87, 88, 89 };

    floatVector previousStateVariables = { 0.01,  0.02,  0.03,  0.04,  0.05,  0.06,  0.07,  0.08,  0.09,
                                          -0.01, -0.02, -0.03, -0.04, -0.05, -0.06, -0.07, -0.08, -0.09,
                                           1.00,  2.00,  3.00,  4.00,  5.00,  6.00,  7.00,  8.00,  9.00,
                                          10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00, 17.00, 18.00,
                                          19.00, 20.00, 21.00, 22.00, 23.00, 24.00, 25.00, 26.00, 27.00 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 6;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 27;

    floatType tolr = 1e-2;

    floatType tola = 1e-3;

    unsigned int maxIterations = 24;

    unsigned int maxLSIterations = 45;

    floatType lsAlpha = 2.3;

    class hydraBaseMicromorphicMock : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            unsigned int ndecomp = 0;
    
            unsigned int nsrc = 0;

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

    };

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( time, *hydra.getTime( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( deltaTime, *hydra.getDeltaTime( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( temperature, *hydra.getTemperature( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousTemperature, *hydra.getPreviousTemperature( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( deformationGradient, *hydra.getDeformationGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousDeformationGradient, *hydra.getPreviousDeformationGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( microDeformation, *hydra.getMicroDeformation( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousMicroDeformation, *hydra.getPreviousMicroDeformation( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradientMicroDeformation, *hydra.getGradientMicroDeformation( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousGradientMicroDeformation, *hydra.getPreviousGradientMicroDeformation( ) ) );

}
