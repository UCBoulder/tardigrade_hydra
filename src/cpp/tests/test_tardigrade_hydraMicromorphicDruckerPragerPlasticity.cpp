/**
  * \file test_tardigrade_hydraMicromorphicDruckerPragerPlasticity.cpp
  *
  * Tests for tardigrade_hydraMicromorphicDruckerPragerPlasticity
  */

#include<tardigrade_hydraMicromorphicDruckerPragerPlasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_micromorphic_tools.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_tardigrade_hydraMicromorphicDruckerPragerPlasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::variableType variableType; //!< Redefinition of the variable type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::variableVector variableVector; //!< Redefinition of the vector of variable types
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::variableMatrix variableMatrix; //!< Redefinition of the matrix of variable types

typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::parameterType parameterType; //!< Redefinition of the parameter type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::parameterVector parameterVector; //!< Redefinition of the vector of parameters

typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::constantType constantType; //!< Redefinition of the constant type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::constantVector constantVector; //!< Redefinition of the vector of constants
typedef tardigradeHydra::micromorphicDruckerPragerPlasticity::constantMatrix constantMatrix; //!< Redefinition of the matrix of constants

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

                static void updateUnknownVector( tardigradeHydra::hydraBase &hydra, const floatVector &value ){

                    BOOST_CHECK_NO_THROW( hydra.updateUnknownVector( value ) );

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

BOOST_AUTO_TEST_CASE( test_setDrivingStresses ){
    /*!
     * Test of the driving stresses
     */

    // Form a hydra class and residual class
    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    variableVector deformationGradient = {  1.00157757,  0.00159138,  0.00672005,
                                            0.01747159,  1.01122277,  0.00555118,
                                            0.01112217, -0.00885205,  0.99308943 };

    variableVector microDeformation = { 0.99460588, -0.0078411 ,  0.01145249,
                                       -0.00307139,  0.97798389, -0.00509779,
                                        0.01189977, -0.01587541,  0.98377259 };

    variableVector gradientMicroDeformation = { -1.35868385e-02, -1.03142977e-02,  6.54880619e-03, -2.03947530e-02,
                                                -3.31494137e-03, -3.45686183e-03, -3.15745117e-04, -3.70848549e-03,
                                                -9.38693885e-03, -3.68243465e-03,  1.96694582e-02,  2.22080009e-02,
                                                 9.18337942e-05,  6.19764759e-03, -1.92190802e-02, -9.13572591e-03,
                                                -4.25868940e-03,  1.83154579e-02, -1.24772317e-02, -8.48286787e-04,
                                                 2.42779893e-02,  9.74255963e-04,  5.64472629e-03, -1.89685667e-02,
                                                 1.63170400e-02,  5.15300642e-03,  2.25340032e-03 };

    floatVector previousDeformationGradient = { 9.94656270e-01,  4.82152400e-02,  3.31984800e-02,  2.81918700e-02,
         1.02086536e+00, -1.77592100e-02, -2.24798000e-03, -1.28410000e-04,
         9.77165250e-01 };

    floatVector previousMicroDeformation = { 0.96917405, -0.01777599,  0.00870406, -0.02163002,  0.9998683 ,
        -0.01669352,  0.03355217,  0.04427456,  1.01778466 };

    floatVector previousGradientMicroDeformation = { 0.05043761,  0.02160516, -0.0565408 ,  0.01218304, -0.05851034,
         0.00485749, -0.00962607, -0.03455912,  0.04490067,  0.01552915,
        -0.02878364,  0.00595866,  0.04750406, -0.02377005, -0.05041534,
        -0.02922214,  0.06280788,  0.02850865, -0.00226005,  0.0146049 ,
         0.01560184,  0.03224767,  0.05822091, -0.05294424, -0.03518206,
         0.01831308,  0.03774438 };

    floatVector previousStateVariables = {-0.02495446, -0.00169657,  0.04855598,  0.00194851,  0.01128945,
       -0.03793713,  0.03263408,  0.01030601,  0.0045068 , -0.01572362,
       -0.01958792, -0.00829778,  0.01813008,  0.03754568,  0.00104223,
        0.01693138,  0.00859366,  0.01249035,  0.01746891,  0.03423424,
       -0.0416805 ,  0.02636828, -0.02563336, -0.0305777 ,  0.0072457 ,
       -0.04042875,  0.03853268,  0.0127249 ,  0.02234164, -0.04838708,
        0.00944319,  0.00567852, -0.03410404, -0.03469295,  0.01955295,
       -0.01812336,  0.01919703,  0.00543832, -0.01110494,  0.04251325,
        0.034167  , -0.01426024, -0.04564085, -0.01952319, -0.01018143,
        0.1, 0.2, 0.3, 0.4, 0.5 };

    floatVector parameters = { 2, 0.1, 0.2, 5, 0.3, 0.4, 0.5, 0.6, 0.7, 11, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2, 1.9, 2.0 };

    floatVector unknownVector = {  4.31202550e-01, -1.78960429e-01, -6.17986089e-01,  9.34988614e-01,
                                   3.01500733e-01,  7.30919703e-01, -9.49515284e-01, -4.66188370e-01,
                                   4.14220065e-03, -8.65102730e-01,  9.86066522e-01, -5.27075208e-01,
                                  -2.51415635e-01, -5.71976170e-01, -7.89108268e-01, -5.35040429e-01,
                                  -3.98779729e-01,  2.68884536e-01, -4.37530437e-01, -2.75446478e-01,
                                  -9.88114313e-01, -2.68561748e-01,  6.77719634e-02, -6.75968326e-01,
                                   1.94866217e-01, -4.13695063e-01,  2.64100990e-01, -9.47606789e-01,
                                   7.75186921e-01, -9.67762739e-01, -7.46083938e-01,  5.54324923e-01,
                                  -9.08209536e-01,  4.21997387e-01,  9.42092281e-01,  7.43365866e-01,
                                   4.20323303e-01,  9.17019486e-01, -1.40373324e-01,  7.45757829e-01,
                                  -2.88084664e-01,  8.59527306e-01, -7.02444688e-01,  8.80058030e-01,
                                   6.65432395e-01,  1.01730274e+00, -1.88038495e-02,  4.82434492e-03,
                                  -2.41803760e-02,  1.01105922e+00, -2.46131243e-02, -2.07588861e-02,
                                  -1.37250795e-02,  1.01875623e+00,  9.93178816e-01,  1.99799676e-03,
                                   3.40516069e-03, -1.37268320e-02,  1.00360734e+00,  8.04758975e-03,
                                  -1.00877303e-02, -4.06865705e-03,  9.97654446e-01,  2.16175331e-02,
                                   4.37468737e-03,  2.24126186e-02,  2.80173769e-03,  2.80710425e-05,
                                  -2.48233895e-02, -9.55547806e-04,  2.13727499e-02, -1.50817155e-02,
                                  -2.23954433e-02, -4.66105533e-03, -6.38017597e-03,  1.78576529e-02,
                                  -2.36694442e-02,  2.10074615e-02,  9.04514995e-03,  2.02112997e-02,
                                   5.37645354e-03,  1.55976656e-02, -8.22280632e-03, -7.52168860e-03,
                                  -5.50628848e-03,  1.27398541e-02, -6.53544128e-03, -1.28890097e-02,
                                   2.18834178e-02,  2.04005542e-02, 0.01, 0.02, 0.03, 0.04, 0.05 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 5;

    std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5 };

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

    floatType tolr = 1e-2;

    floatType tola = 1e-3;

    unsigned int maxIterations = 24;

    unsigned int maxLSIterations = 45;

    floatType lsAlpha = 2.3;

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticity::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::residual;

    };

    class hydraBaseMicromorphicMock : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

    };

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 45, 1, stateVariableIndices, parameters );

}
