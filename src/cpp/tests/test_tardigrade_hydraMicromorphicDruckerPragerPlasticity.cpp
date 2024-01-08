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

BOOST_AUTO_TEST_CASE( test_extractParameters ){
    /*!
     * Test of the extraction of the parameters
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
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

    floatVector parameters = { 2, 0.53895133, 0.37172145,
                               2, 0.37773052, 0.92739145,
                               2, 0.53186824, 0.75454313,
                               3, 0.95338442, 0.74042148, 0.09916127,
                               3, 0.38093104, 0.49241325, 0.46187452,
                               5, 0.82121039, 0.90566759, 0.50466975, 0.04830311, 0.85951495,
                               3, 0.01166325, 0.05331896, 0.28081774,
                               3, 0.32982199, 0.60161431, 0.33157768,
                               5, 0.58881096, 0.11473813, 0.58001078, 0.83382529, 0.3260217 };

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
                                   2.18834178e-02,  2.04005542e-02,
                                   0.01, 0.02, 0.03, 0.04, 0.05,
                                   0.06, 0.07, 0.08, 0.09, 0.10 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 10;

    std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5 };

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

    floatType tolr = 1e-2;

    floatType tola = 1e-3;

    unsigned int maxIterations = 24;

    unsigned int maxLSIterations = 45;

    floatType lsAlpha = 2.3;

    class stressMock : public tardigradeHydra::residualBaseMicromorphic{

        public:

            using tardigradeHydra::residualBaseMicromorphic::residualBaseMicromorphic;

    };

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticity::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::residual;

    };

    class hydraBaseMicromorphicMock : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

            stressMock elasticity;

            residualMock plasticity;

            std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4 };

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              3, 0.95338442, 0.74042148, 0.09916127,
                                              3, 0.38093104, 0.49241325, 0.46187452,
                                              5, 0.82121039, 0.90566759, 0.50466975, 0.04830311, 0.85951495,
                                              3, 0.01166325, 0.05331896, 0.28081774,
                                              3, 0.32982199, 0.60161431, 0.33157768,
                                              5, 0.58881096, 0.11473813, 0.58001078, 0.83382529, 0.3260217 };

        private:

            using tardigradeHydra::hydraBaseMicromorphic::setResidualClasses;

            virtual void setResidualClasses( ) override{

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = stressMock( this, 45 );

                plasticity = residualMock( this, 55, 1, stateVariableIndices, plasticParameters );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                setResidualClasses( residuals );

            }

    };

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 55, 1, stateVariableIndices, parameters );

    floatVector answer1 = { 0.53895133, 0.37172145 };

    floatVector answer2 = { 0.37773052, 0.92739145 };

    floatVector answer3 = { 0.53186824, 0.75454313 };

    floatVector answer4 = { 0.95338442, 0.74042148, 0.09916127 };

    floatVector answer5 = { 0.38093104, 0.49241325, 0.46187452 };

    floatVector answer6 = { 0.82121039, 0.90566759, 0.50466975, 0.04830311, 0.85951495 };

    floatVector answer7 = { 0.01166325, 0.05331896, 0.28081774 };

    floatVector answer8 = { 0.32982199, 0.60161431, 0.33157768 };

    floatVector answer9 = { 0.58881096, 0.11473813, 0.58001078, 0.83382529, 0.3260217 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer1, *R.get_macroHardeningParameters( )         ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *R.get_microHardeningParameters( )         ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer3, *R.get_microGradientHardeningParameters( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer4, *R.get_macroFlowParameters( )              ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer5, *R.get_microFlowParameters( )              ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer6, *R.get_microGradientFlowParameters( )      ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer7, *R.get_macroYieldParameters( )             ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer8, *R.get_microYieldParameters( )             ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer9, *R.get_microGradientYieldParameters( )     ) );

}

BOOST_AUTO_TEST_CASE( test_set_state_variables ){
    /*!
     * Test of the extraction of the state variables
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
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

    floatVector parameters = { 2, 0.53895133, 0.37172145,
                               2, 0.37773052, 0.92739145,
                               2, 0.53186824, 0.75454313,
                               3, 0.95338442, 0.74042148, 0.09916127,
                               3, 0.38093104, 0.49241325, 0.46187452,
                               5, 0.82121039, 0.90566759, 0.50466975, 0.04830311, 0.85951495,
                               3, 0.01166325, 0.05331896, 0.28081774,
                               3, 0.32982199, 0.60161431, 0.33157768,
                               5, 0.58881096, 0.11473813, 0.58001078, 0.83382529, 0.3260217 };

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
                                   2.18834178e-02,  2.04005542e-02,
                                   0.01, 0.02, 0.03, 0.04, 0.05,
                                   0.06, 0.07, 0.08, 0.09, 0.10 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 10;

    std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

    floatType tolr = 1e-2;

    floatType tola = 1e-3;

    unsigned int maxIterations = 24;

    unsigned int maxLSIterations = 45;

    floatType lsAlpha = 2.3;

    class stressMock : public tardigradeHydra::residualBaseMicromorphic{

        public:

            using tardigradeHydra::residualBaseMicromorphic::residualBaseMicromorphic;

    };

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticity::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::residual;

    };

    class hydraBaseMicromorphicMock : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

            stressMock elasticity;

            residualMock plasticity;

            std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              3, 0.95338442, 0.74042148, 0.09916127,
                                              3, 0.38093104, 0.49241325, 0.46187452,
                                              5, 0.82121039, 0.90566759, 0.50466975, 0.04830311, 0.85951495,
                                              3, 0.01166325, 0.05331896, 0.28081774,
                                              3, 0.32982199, 0.60161431, 0.33157768,
                                              5, 0.58881096, 0.11473813, 0.58001078, 0.83382529, 0.3260217 };

        private:

            using tardigradeHydra::hydraBaseMicromorphic::setResidualClasses;

            virtual void setResidualClasses( ) override{

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = stressMock( this, 45 );

                plasticity = residualMock( this, 55, 1, stateVariableIndices, plasticParameters );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                setResidualClasses( residuals );

            }

    };

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 55, 1, stateVariableIndices, parameters );

    floatVector stateVariablesAnswer = { 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10 };

    floatVector previousStateVariablesAnswer = { 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00 };

    floatVector plasticMultipliersAnswer = { 0.01, 0.02, 0.03, 0.04, 0.05 };

    floatVector previousPlasticMultipliersAnswer = { 0.10, 0.20, 0.30, 0.40, 0.50 };

    floatVector plasticStrainLikeISVsAnswer = { 0.06, 0.07, 0.08, 0.09, 0.10 };

    floatVector previousPlasticStrainLikeISVsAnswer = { 0.60, 0.70, 0.80, 0.90, 1.00 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( stateVariablesAnswer, *R.get_plasticStateVariables( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousStateVariablesAnswer, *R.get_previousPlasticStateVariables( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( plasticMultipliersAnswer, *R.get_plasticMultipliers( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousPlasticMultipliersAnswer, *R.get_previousPlasticMultipliers( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( plasticStrainLikeISVsAnswer, *R.get_plasticStrainLikeISVs( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousPlasticStrainLikeISVsAnswer, *R.get_previousPlasticStrainLikeISVs( ) ) );

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
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

    floatVector parameters = { 2, 0.53895133, 0.37172145,
                               2, 0.37773052, 0.92739145,
                               2, 0.53186824, 0.75454313,
                               3, 0.95338442, 0.74042148, 0.09916127,
                               3, 0.38093104, 0.49241325, 0.46187452,
                               5, 0.82121039, 0.90566759, 0.50466975, 0.04830311, 0.85951495,
                               3, 0.01166325, 0.05331896, 0.28081774,
                               3, 0.32982199, 0.60161431, 0.33157768,
                               5, 0.58881096, 0.11473813, 0.58001078, 0.83382529, 0.3260217 };

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
                                   2.18834178e-02,  2.04005542e-02,
                                   0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 10;

    std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

    floatType tolr = 1e-2;

    floatType tola = 1e-3;

    unsigned int maxIterations = 24;

    unsigned int maxLSIterations = 45;

    floatType lsAlpha = 2.3;

    class stressMock : public tardigradeHydra::residualBaseMicromorphic{

        public:

            using tardigradeHydra::residualBaseMicromorphic::residualBaseMicromorphic;

            floatVector previousPK2Stress = {  1,  2,  3,  4,  5,  6,  7,  8,  9 };

            floatVector previousSigma     = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

            floatVector previousM         = { 19, 20, 21, 22, 23, 24, 25, 26, 27,
                                              28, 29, 30, 31, 32, 33, 34, 35, 36,
                                              37, 38, 39, 40, 41, 42, 43, 44, 45 };

        protected:

            using tardigradeHydra::residualBase::setPreviousStress;

            virtual void setPreviousStress( ) override{

                setPreviousStress( tardigradeVectorTools::appendVectors( { previousPK2Stress, previousSigma, previousM } ) );

            }

    };

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticity::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::residual;

    };

    class hydraBaseMicromorphicMock : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

            stressMock elasticity;

            residualMock plasticity;

            std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              3, 0.95338442, 0.74042148, 0.09916127,
                                              3, 0.38093104, 0.49241325, 0.46187452,
                                              5, 0.82121039, 0.90566759, 0.50466975, 0.04830311, 0.85951495,
                                              3, 0.01166325, 0.05331896, 0.28081774,
                                              3, 0.32982199, 0.60161431, 0.33157768,
                                              5, 0.58881096, 0.11473813, 0.58001078, 0.83382529, 0.3260217 };

            floatVector _local_deltaPK2Stress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            floatVector _local_deltaSigma     = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            floatVector _local_deltaM         = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                  0, 0, 0, 0, 0, 0, 0, 0, 0,
                                                  0, 0, 0, 0, 0, 0, 0, 0, 0 };

        private:

            using tardigradeHydra::hydraBaseMicromorphic::setResidualClasses;

            virtual void setResidualClasses( ) override{

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = stressMock( this, 45 );

                elasticity.previousPK2Stress += _local_deltaPK2Stress;

                elasticity.previousSigma     += _local_deltaSigma;

                elasticity.previousM         += _local_deltaM;

                plasticity = residualMock( this, 55, 1, stateVariableIndices, plasticParameters );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                setResidualClasses( residuals );

            }

    };

    // compute the expected current answer

    floatVector PK2Stress( unknownVector.begin( ),
                           unknownVector.begin( ) + dimension * dimension );

    floatVector referenceMicroStress( unknownVector.begin( ) + dimension * dimension,
                                      unknownVector.begin( ) + 2 * dimension * dimension );

    floatVector referenceHigherOrderStress( unknownVector.begin( ) + 2 * dimension * dimension,
                                            unknownVector.begin( ) + configuration_unknown_count );

    floatVector F( unknownVector.begin( ) + configuration_unknown_count,
                   unknownVector.begin( ) + configuration_unknown_count + dimension * dimension );

    floatVector chi( unknownVector.begin( ) + configuration_unknown_count + dimension * dimension,
                     unknownVector.begin( ) + configuration_unknown_count + 2 * dimension * dimension );

    floatVector answerMacroStress;

    floatVector answerMicroStress;

    floatVector answerHigherOrderStress;

    tardigradeMicromorphicTools::pushForwardPK2Stress(                          PK2Stress, F,      answerMacroStress       );

    tardigradeMicromorphicTools::pushForwardReferenceMicroStress(    referenceMicroStress, F,      answerMicroStress       );

    tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, F, chi, answerHigherOrderStress );

    // compute the expected previous answer
    floatVector eye( dimension * dimension );

    tardigradeVectorTools::eye( eye );

    floatVector previousPK2Stress                  = {  1,  2,  3,  4,  5,  6,  7,  8,  9 };

    floatVector previousReferenceMicroStress       = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

    floatVector previousReferenceHigherOrderStress = { 19, 20, 21, 22, 23, 24, 25, 26, 27,
                                                       28, 29, 30, 31, 32, 33, 34, 35, 36,
                                                       37, 38, 39, 40, 41, 42, 43, 44, 45 };

    floatVector previousF( previousStateVariables.begin( ),
                           previousStateVariables.begin( ) + dimension * dimension );

    previousF += eye;

    floatVector previousChi( previousStateVariables.begin( ) + dimension * dimension,
                             previousStateVariables.begin( ) + 2 * dimension * dimension );

    previousChi += eye;

    floatVector answerPreviousMacroStress;

    floatVector answerPreviousMicroStress;

    floatVector answerPreviousHigherOrderStress;

    tardigradeMicromorphicTools::pushForwardPK2Stress(                          previousPK2Stress, previousF,              answerPreviousMacroStress       );

    tardigradeMicromorphicTools::pushForwardReferenceMicroStress(    previousReferenceMicroStress, previousF,              answerPreviousMicroStress       );

    tardigradeMicromorphicTools::pushForwardHigherOrderStress( previousReferenceHigherOrderStress, previousF, previousChi, answerPreviousHigherOrderStress );

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 55, 1, stateVariableIndices, parameters );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMacroStress,               *R.get_macroDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMicroStress,               *R.get_symmetricMicroDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerHigherOrderStress,         *R.get_higherOrderDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousMacroStress,       *R.get_previousMacroDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousMicroStress,       *R.get_previousSymmetricMicroDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousHigherOrderStress, *R.get_previousHigherOrderDrivingStress( ) ) );

    // Check the Jacobians

    floatType eps = 1e-6;

    floatMatrix dMacrodX(     9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dMicrodX(     9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dHigherdX(   27, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dMacrodF(     9, floatVector( 9, 0 ) );

    floatMatrix dMicrodF(     9, floatVector( 9, 0 ) );

    floatMatrix dHigherdF(   27, floatVector( 9, 0 ) );

    floatMatrix dMacrodChi(   9, floatVector( 9, 0 ) );

    floatMatrix dMicrodChi(   9, floatVector( 9, 0 ) );

    floatMatrix dHigherdChi( 27, floatVector( 9, 0 ) );

    floatMatrix previousdMacrodX(     9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix previousdMicrodX(     9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix previousdHigherdX(   27, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix previousdMacrodF(     9, floatVector( 9, 0 ) );

    floatMatrix previousdMicrodF(     9, floatVector( 9, 0 ) );

    floatMatrix previousdHigherdF(   27, floatVector( 9, 0 ) );

    floatMatrix previousdMacrodChi(   9, floatVector( 9, 0 ) );

    floatMatrix previousdMicrodChi(   9, floatVector( 9, 0 ) );

    floatMatrix previousdHigherdChi( 27, floatVector( 9, 0 ) );

    // Check the Jacobians of the deformation gradient
    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );
    
        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm( &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_macroDrivingStress( );

        floatVector vm = *Rm.get_macroDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMacrodF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_symmetricMicroDrivingStress( );

        vm = *Rm.get_symmetricMicroDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicrodF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_higherOrderDrivingStress( );

        vm = *Rm.get_higherOrderDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dHigherdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals(  dMacrodF, floatMatrix(  9, floatVector( 9, 0 ) ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals(  dMicrodF, floatMatrix(  9, floatVector( 9, 0 ) ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dHigherdF, floatMatrix( 27, floatVector( 9, 0 ) ) ) );

    // Check the Jacobians of the micro deformation
    for ( unsigned int i = 0; i < microDeformation.size( ); i++ ){

        floatVector delta( microDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation + delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation - delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );
    
        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm( &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_macroDrivingStress( );

        floatVector vm = *Rm.get_macroDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMacrodChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_symmetricMicroDrivingStress( );

        vm = *Rm.get_symmetricMicroDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicrodChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_higherOrderDrivingStress( );

        vm = *Rm.get_higherOrderDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dHigherdChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals(  dMacrodChi, floatMatrix(  9, floatVector( 9, 0 ) ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals(  dMicrodChi, floatMatrix(  9, floatVector( 9, 0 ) ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dHigherdChi, floatMatrix( 27, floatVector( 9, 0 ) ) ) );

    // Check the Jacobians of the unknown vector
    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );
    
        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        residualMock Rp( &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm( &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_macroDrivingStress( );

        floatVector vm = *Rm.get_macroDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMacrodX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_symmetricMicroDrivingStress( );

        vm = *Rm.get_symmetricMicroDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicrodX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_higherOrderDrivingStress( );

        vm = *Rm.get_higherOrderDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dHigherdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    // Assemble the results
    floatMatrix result_dMacrodX(   9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix result_dMicrodX(   9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix result_dHigherdX( 27, floatVector( unknownVector.size( ), 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        unsigned int col = 0;

        for ( auto v  = ( *R.get_dMacroDrivingStressdMacroStress( ) )[ i ].begin( );
                   v != ( *R.get_dMacroDrivingStressdMacroStress( ) )[ i ].end( );
                   v++ ){

            result_dMacrodX[ i ][ col ] = *v;

            col++;

        }

        col = configuration_unknown_count;

        for ( auto v  = ( *R.get_dMacroDrivingStressdFn( ) )[ i ].begin( );
                   v != ( *R.get_dMacroDrivingStressdFn( ) )[ i ].end( );
                   v++ ){

            result_dMacrodX[ i ][ col ] = *v;

            col++;

        }

    }

    for ( unsigned int i = 0; i < 9; i++ ){

        unsigned int col = dimension * dimension;

        for ( auto v  = ( *R.get_dSymmetricMicroDrivingStressdMicroStress( ) )[ i ].begin( );
                   v != ( *R.get_dSymmetricMicroDrivingStressdMicroStress( ) )[ i ].end( );
                   v++ ){

            result_dMicrodX[ i ][ col ] = *v;

            col++;

        }

        col = configuration_unknown_count;

        for ( auto v  = ( *R.get_dSymmetricMicroDrivingStressdFn( ) )[ i ].begin( );
                   v != ( *R.get_dSymmetricMicroDrivingStressdFn( ) )[ i ].end( );
                   v++ ){

            result_dMicrodX[ i ][ col ] = *v;

            col++;

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        unsigned int col = 2 * dimension * dimension;

        for ( auto v  = ( *R.get_dHigherOrderDrivingStressdHigherOrderStress( ) )[ i ].begin( );
                   v != ( *R.get_dHigherOrderDrivingStressdHigherOrderStress( ) )[ i ].end( );
                   v++ ){

            result_dHigherdX[ i ][ col ] = *v;

            col++;

        }

        col = configuration_unknown_count;

        for ( auto v  = ( *R.get_dHigherOrderDrivingStressdFn( ) )[ i ].begin( );
                   v != ( *R.get_dHigherOrderDrivingStressdFn( ) )[ i ].end( );
                   v++ ){

            result_dHigherdX[ i ][ col ] = *v;

            col++;

        }

        for ( auto v  = ( *R.get_dHigherOrderDrivingStressdChin( ) )[ i ].begin( );
                   v != ( *R.get_dHigherOrderDrivingStressdChin( ) )[ i ].end( );
                   v++ ){

            result_dHigherdX[ i ][ col ] = *v;

            col++;

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals(  dMacrodX, result_dMacrodX  ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals(  dMicrodX, result_dMicrodX  ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dHigherdX, result_dHigherdX ) );




    // Check the Jacobians of the deformation gradient
    floatVector dvec = tardigradeVectorTools::appendVectors( { hydra.elasticity.previousPK2Stress, hydra.elasticity.previousSigma, hydra.elasticity.previousM, previousStateVariables } );

    for ( unsigned int i = 0; i < previousDeformationGradient.size( ); i++ ){

        floatVector delta( previousDeformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );
    
        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm( &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousMacroDrivingStress( );

        floatVector vm = *Rm.get_previousMacroDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMacrodF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousSymmetricMicroDrivingStress( );

        vm = *Rm.get_previousSymmetricMicroDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicrodF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousHigherOrderDrivingStress( );

        vm = *Rm.get_previousHigherOrderDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdHigherdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals(  previousdMacrodF, floatMatrix(  9, floatVector( 9, 0 ) ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals(  previousdMicrodF, floatMatrix(  9, floatVector( 9, 0 ) ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdHigherdF, floatMatrix( 27, floatVector( 9, 0 ) ) ) );

    // Check the Jacobians of the micro deformation
    for ( unsigned int i = 0; i < previousMicroDeformation.size( ); i++ ){

        floatVector delta( previousMicroDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation + delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation - delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );
    
        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm( &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousMacroDrivingStress( );

        floatVector vm = *Rm.get_previousMacroDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMacrodChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousSymmetricMicroDrivingStress( );

        vm = *Rm.get_previousSymmetricMicroDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicrodChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousHigherOrderDrivingStress( );

        vm = *Rm.get_previousHigherOrderDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdHigherdChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals(  previousdMacrodChi, floatMatrix(  9, floatVector( 9, 0 ) ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals(  previousdMicrodChi, floatMatrix(  9, floatVector( 9, 0 ) ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdHigherdChi, floatMatrix( 27, floatVector( 9, 0 ) ) ) );

    // Check the Jacobians of the unknown vector
    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatVector delta( dvec.size( ), 0 );

        delta[ i ] = eps * std::fabs( dvec[ i ] ) + eps;

        floatVector dPK2(   delta.begin( ),      delta.begin( ) + 9  );
        floatVector dSIGMA( delta.begin( ) + 9,  delta.begin( ) + 18 );
        floatVector dM(     delta.begin( ) + 18, delta.begin( ) + 45 );
        floatVector dXi(    delta.begin( ) + 45, delta.end( )     );

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables + dXi, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables - dXi, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );
    
        hydrap._local_deltaPK2Stress = dPK2;

        hydrap._local_deltaSigma     = dSIGMA;

        hydrap._local_deltaM         = dM;

        hydram._local_deltaPK2Stress = -dPK2;

        hydram._local_deltaSigma     = -dSIGMA;

        hydram._local_deltaM         = -dM;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm( &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousMacroDrivingStress( );

        floatVector vm = *Rm.get_previousMacroDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMacrodX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousSymmetricMicroDrivingStress( );

        vm = *Rm.get_previousSymmetricMicroDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicrodX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousHigherOrderDrivingStress( );

        vm = *Rm.get_previousHigherOrderDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdHigherdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    // Assemble the results
    floatMatrix result_previousdMacrodX(   9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix result_previousdMicrodX(   9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix result_previousdHigherdX( 27, floatVector( unknownVector.size( ), 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        unsigned int col = 0;

        for ( auto v  = ( *R.get_previousdMacroDrivingStressdMacroStress( ) )[ i ].begin( );
                   v != ( *R.get_previousdMacroDrivingStressdMacroStress( ) )[ i ].end( );
                   v++ ){

            result_previousdMacrodX[ i ][ col ] = *v;

            col++;

        }

        col = configuration_unknown_count;

        for ( auto v  = ( *R.get_previousdMacroDrivingStressdFn( ) )[ i ].begin( );
                   v != ( *R.get_previousdMacroDrivingStressdFn( ) )[ i ].end( );
                   v++ ){

            result_previousdMacrodX[ i ][ col ] = *v;

            col++;

        }

    }

    for ( unsigned int i = 0; i < 9; i++ ){

        unsigned int col = dimension * dimension;

        for ( auto v  = ( *R.get_previousdSymmetricMicroDrivingStressdMicroStress( ) )[ i ].begin( );
                   v != ( *R.get_previousdSymmetricMicroDrivingStressdMicroStress( ) )[ i ].end( );
                   v++ ){

            result_previousdMicrodX[ i ][ col ] = *v;

            col++;

        }

        col = configuration_unknown_count;

        for ( auto v  = ( *R.get_previousdSymmetricMicroDrivingStressdFn( ) )[ i ].begin( );
                   v != ( *R.get_previousdSymmetricMicroDrivingStressdFn( ) )[ i ].end( );
                   v++ ){

            result_previousdMicrodX[ i ][ col ] = *v;

            col++;

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        unsigned int col = 2 * dimension * dimension;

        for ( auto v  = ( *R.get_previousdHigherOrderDrivingStressdHigherOrderStress( ) )[ i ].begin( );
                   v != ( *R.get_previousdHigherOrderDrivingStressdHigherOrderStress( ) )[ i ].end( );
                   v++ ){

            result_previousdHigherdX[ i ][ col ] = *v;

            col++;

        }

        col = configuration_unknown_count;

        for ( auto v  = ( *R.get_previousdHigherOrderDrivingStressdFn( ) )[ i ].begin( );
                   v != ( *R.get_previousdHigherOrderDrivingStressdFn( ) )[ i ].end( );
                   v++ ){

            result_previousdHigherdX[ i ][ col ] = *v;

            col++;

        }

        for ( auto v  = ( *R.get_previousdHigherOrderDrivingStressdChin( ) )[ i ].begin( );
                   v != ( *R.get_previousdHigherOrderDrivingStressdChin( ) )[ i ].end( );
                   v++ ){

            result_previousdHigherdX[ i ][ col ] = *v;

            col++;

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals(  previousdMacrodX, result_previousdMacrodX  ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals(  previousdMicrodX, result_previousdMicrodX  ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdHigherdX, result_previousdHigherdX ) );

}

BOOST_AUTO_TEST_CASE( test_setCohesion ){
    /*!
     * Test setting the cohesion values
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
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

    floatVector parameters = { 2, 0.53895133, 0.37172145,
                               2, 0.37773052, 0.92739145,
                               2, 0.53186824, 0.75454313,
                               2, 0.95338442, 0.74042148,
                               2, 0.38093104, 0.49241325,
                               2, 0.82121039, 0.90566759,
                               2, 0.01166325, 0.05331896,
                               2, 0.32982199, 0.60161431,
                               2, 0.58881096, 0.11473813 };

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
                                   2.18834178e-02,  2.04005542e-02,
                                   0.01, 0.02, 0.03, 0.04, 0.05,
                                   0.06, 0.07, 0.08, 0.09, 0.10 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 10;

    std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

    floatType tolr = 1e-2;

    floatType tola = 1e-3;

    unsigned int maxIterations = 24;

    unsigned int maxLSIterations = 45;

    floatType lsAlpha = 2.3;

    class stressMock : public tardigradeHydra::residualBaseMicromorphic{

        public:

            using tardigradeHydra::residualBaseMicromorphic::residualBaseMicromorphic;

    };

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticity::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::residual;

    };

    class hydraBaseMicromorphicMock : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

            stressMock elasticity;

            residualMock plasticity;

            std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4 };

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              2, 0.95338442, 0.74042148, 
                                              2, 0.38093104, 0.49241325, 
                                              2, 0.82121039, 0.90566759,
                                              2, 0.01166325, 0.05331896, 
                                              2, 0.32982199, 0.60161431, 
                                              2, 0.58881096, 0.11473813 };

        private:

            using tardigradeHydra::hydraBaseMicromorphic::setResidualClasses;

            virtual void setResidualClasses( ) override{

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = stressMock( this, 45 );

                plasticity = residualMock( this, 55, 1, stateVariableIndices, plasticParameters );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                setResidualClasses( residuals );

            }

    };

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 55, 1, stateVariableIndices, parameters );

    floatType answer1 = 0.53895133 + 0.37172145 * 0.06;

    floatType answer2 = 0.37773052 + 0.92739145 * 0.07;

    floatVector temp = { 0.08, 0.09, 0.10 };

    floatVector answer3 = 0.53186824 + 0.75454313 * temp;

    floatType answer4 = 0.53895133 + 0.37172145 * 0.6;

    floatType answer5 = 0.37773052 + 0.92739145 * 0.7;

    temp = { 0.8, 0.9, 1.00 };

    floatVector answer6 = 0.53186824 + 0.75454313 * temp;

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer1, *R.get_macroCohesion( )                 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *R.get_microCohesion( )                 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer3, *R.get_microGradientCohesion( )         ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer4, *R.get_previousMacroCohesion( )         ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer5, *R.get_previousMicroCohesion( )         ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer6, *R.get_previousMicroGradientCohesion( ) ) );

    // Test the Jacobian
    floatType eps = 1e-6;

    floatVector dMacroCohesiondX( unknownVector.size( ), 0 );

    floatVector dMicroCohesiondX( unknownVector.size( ), 0 );

    floatMatrix dMicroGradientCohesiondX( 3, floatVector( unknownVector.size( ), 0 ) );

    floatVector previousdMacroCohesiondX( previousStateVariables.size( ), 0 );

    floatVector previousdMicroCohesiondX( previousStateVariables.size( ), 0 );

    floatMatrix previousdMicroGradientCohesiondX( 3, floatVector( previousStateVariables.size( ), 0 ) );

    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );
    
        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );
    
        residualMock Rp( &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm( &hydram, 55, 1, stateVariableIndices, parameters );

        floatType sp = *Rp.get_macroCohesion( );

        floatType sm = *Rm.get_macroCohesion( );

        for ( unsigned int j = 0; j < 1; j++ ){

            dMacroCohesiondX[ i ] = ( sp - sm ) / ( 2 * delta[ i ] );

        }

        sp = *Rp.get_microCohesion( );

        sm = *Rm.get_microCohesion( );

        for ( unsigned int j = 0; j < 1; j++ ){

            dMicroCohesiondX[ i ] = ( sp - sm ) / ( 2 * delta[ i ] );

        }

        floatVector vp = *Rp.get_microGradientCohesion( );

        floatVector vm = *Rm.get_microGradientCohesion( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicroGradientCohesiondX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatVector dMacroCohesiondXAssembled( unknownVector.size( ), 0 );

    floatVector dMicroCohesiondXAssembled( unknownVector.size( ), 0 );

    floatMatrix dMicroGradientCohesiondXAssembled( 3, floatVector( unknownVector.size( ), 0 ) );

    for ( unsigned int i = 0; i < R.get_dMacroCohesiondStateVariables( )->size( ); i++ ){

        dMacroCohesiondXAssembled[ i + numConfigurations * configuration_unknown_count ] = ( *R.get_dMacroCohesiondStateVariables( ) )[ i ];

        dMicroCohesiondXAssembled[ i + numConfigurations * configuration_unknown_count ] = ( *R.get_dMicroCohesiondStateVariables( ) )[ i ];

    }

    for ( unsigned int i = 0; i < R.get_dMicroGradientCohesiondStateVariables( )->size( ); i++ ){

        for ( unsigned int j = 0; j < ( *R.get_dMicroGradientCohesiondStateVariables( ) )[ i ].size( ); j++ ){

            dMicroGradientCohesiondXAssembled[ i ][ j + numConfigurations * configuration_unknown_count ] = ( *R.get_dMicroGradientCohesiondStateVariables( ) )[ i ][ j ];

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMacroCohesiondXAssembled, dMacroCohesiondX ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroCohesiondXAssembled, dMicroCohesiondX ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroGradientCohesiondXAssembled, dMicroGradientCohesiondX ) );

    for ( unsigned int i = 0; i < previousStateVariables.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables - delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );
    
        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );
    
        residualMock Rp( &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm( &hydram, 55, 1, stateVariableIndices, parameters );

        floatType sp = *Rp.get_previousMacroCohesion( );

        floatType sm = *Rm.get_previousMacroCohesion( );

        for ( unsigned int j = 0; j < 1; j++ ){

            previousdMacroCohesiondX[ i ] = ( sp - sm ) / ( 2 * delta[ i ] );

        }

        sp = *Rp.get_previousMicroCohesion( );

        sm = *Rm.get_previousMicroCohesion( );

        for ( unsigned int j = 0; j < 1; j++ ){

            previousdMicroCohesiondX[ i ] = ( sp - sm ) / ( 2 * delta[ i ] );

        }

        floatVector vp = *Rp.get_previousMicroGradientCohesion( );

        floatVector vm = *Rm.get_previousMicroGradientCohesion( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroGradientCohesiondX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatVector previousdMacroCohesiondXAssembled( previousStateVariables.size( ), 0 );

    floatVector previousdMicroCohesiondXAssembled( previousStateVariables.size( ), 0 );

    floatMatrix previousdMicroGradientCohesiondXAssembled( 3, floatVector( previousStateVariables.size( ), 0 ) );

    for ( unsigned int i = 0; i < R.get_dMacroCohesiondStateVariables( )->size( ); i++ ){

        previousdMacroCohesiondXAssembled[ i + ( numConfigurations - 1 ) * configuration_unknown_count ] = ( *R.get_dMacroCohesiondStateVariables( ) )[ i ];

        previousdMicroCohesiondXAssembled[ i + ( numConfigurations - 1 ) * configuration_unknown_count ] = ( *R.get_dMicroCohesiondStateVariables( ) )[ i ];

    }

    for ( unsigned int i = 0; i < R.get_dMicroGradientCohesiondStateVariables( )->size( ); i++ ){

        for ( unsigned int j = 0; j < ( *R.get_dMicroGradientCohesiondStateVariables( ) )[ i ].size( ); j++ ){

            previousdMicroGradientCohesiondXAssembled[ i ][ j + ( numConfigurations - 1 ) * configuration_unknown_count ] = ( *R.get_dMicroGradientCohesiondStateVariables( ) )[ i ][ j ];

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMacroCohesiondXAssembled, previousdMacroCohesiondX ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroCohesiondXAssembled, previousdMicroCohesiondX ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroGradientCohesiondXAssembled, previousdMicroGradientCohesiondX ) );

}

BOOST_AUTO_TEST_CASE( testComputeDruckerPragerInternalParameters ){
    /*!
     * Test the computation of the internal parameters for the Drucker-Prager plasticity.
     *
     */

    parameterType frictionAngle = 0.25;
    parameterType beta = .423;

    parameterType Aanswer = 1.528893501990677;
    parameterType Banswer = 0.39039060414065774;

    parameterType Aresult, Bresult;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeDruckerPragerInternalParameters( frictionAngle, beta, Aresult, Bresult );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( Aresult, Aanswer ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( Bresult, Banswer ) );

}

BOOST_AUTO_TEST_CASE( testComputeSecondOrderDruckerPragerYieldEquation ){
    /*!
     * Test the computation of the second order stress Drucker-Prager yield equation.
     *
     */

    variableVector S = { 0.77189588, -0.84417528,  0.95929231,
                        -0.50465708,  0.50576944,  0.05335127,
                         0.81510751,  0.76814059, -0.82146208 };

    variableVector precedingF = { 1.03482346, 0.01430697, 0.01134257,
                                  0.02756574, 1.03597345, 0.02115532,
                                  0.04903821, 0.03424149, 1.0240466 };

    variableType cohesion = 0.51;
    parameterType frictionAngle = 0.46;
    parameterType beta = 0.34;

    constantType answer = 1.5860256857059118;
    variableType result;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, precedingF,
                                                                                                        frictionAngle, beta, result );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    //Test the Jacobian
    variableType resultJ, dFdCohesion;
    variableVector dFdS, dFdPrecedingF;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, precedingF,
                                                                                frictionAngle, beta, resultJ,
                                                                                dFdS, dFdCohesion, dFdPrecedingF );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    //Test the Jacobian
    variableType resultJ2, dFdcohesionJ2;
    variableVector dFdSJ2, dFdPrecedingFJ2;
    variableMatrix d2FdStress2J2;
    variableMatrix d2FdStressdElasticRCGJ2;
    variableMatrix d2FdS2J2, d2FdSdPrecedingFJ2;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, precedingF,
                                                                                frictionAngle, beta, resultJ2,
                                                                                dFdSJ2, dFdcohesionJ2, dFdPrecedingFJ2, d2FdS2J2,
                                                                                d2FdSdPrecedingFJ2);

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ2, answer ) );

    //Test dFdStress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < S.size(); i++ ){
        constantVector delta( S.size(), 0 );
        delta[i] = eps * fabs( S[i] ) + eps;

        floatType rp, rm;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S + delta, cohesion, precedingF,
                                                                                                            frictionAngle, beta, rp );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S - delta, cohesion, precedingF,
                                                                                                            frictionAngle, beta, rm );

        constantType gradCol = ( rp - rm ) / ( 2 * delta[i] );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol, dFdS[i] ) );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol, dFdSJ2[i] ) );

    }

    //Test dFdPrecedingF
    for ( unsigned int i = 0; i < precedingF.size(); i++ ){
        constantVector delta( precedingF.size(), 0 );
        delta[i] = eps * fabs( precedingF[i] ) + eps;

        floatType rp, rm;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, precedingF + delta,
                                                                                    frictionAngle, beta, rp );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, precedingF - delta,
                                                                                    frictionAngle, beta, rm );

        constantType gradCol = ( rp - rm ) / ( 2 * delta[i] );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol, dFdPrecedingF[i] ) );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol, dFdPrecedingFJ2[i] ) );

    }

    //Test dFdcohesion
    for ( unsigned int i = 0; i < 1; i++ ){
        constantType delta = eps * fabs( cohesion ) + eps;

        floatType rp, rm;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion + delta, precedingF,
                                                                                                            frictionAngle, beta, rp );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion - delta, precedingF,
                                                                                                                    frictionAngle, beta, rm );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( ( rp - rm ) / ( 2 * delta ), dFdCohesion ) );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( ( rp - rm ) / ( 2 * delta ), dFdcohesionJ2 ) );

    }

    //Test d2FdStress2
    for ( unsigned int i = 0; i < S.size(); i++ ){
        constantVector delta( S.size(), 0 );
        delta[i] = eps * fabs( S[i] ) + eps;

        floatType rp, rm;

        variableVector dFdSp, dFdSm, dFdPrecedingFp, dFdPrecedingFm;
        variableType dFdcohesionp, dFdcohesionm;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S + delta, cohesion, precedingF,
                                                                                    frictionAngle, beta, rp,
                                                                                    dFdSp, dFdcohesionp, dFdPrecedingFp );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S - delta, cohesion, precedingF,
                                                                                    frictionAngle, beta, rm,
                                                                                    dFdSm, dFdcohesionm, dFdPrecedingFm );

        constantVector gradCol = ( dFdSp - dFdSm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], d2FdS2J2[j][i] ) );
        }
    }

    //Test d2FdStressdPrecedingF
    for ( unsigned int i = 0; i < precedingF.size(); i++ ){
        constantVector delta( precedingF.size(), 0 );
        delta[i] = eps * fabs( precedingF[i] ) + eps;

        floatType rp, rm;
        variableVector dFdSp, dFdSm, dFdPrecedingFp, dFdPrecedingFm;
        variableType dFdcohesionp, dFdcohesionm;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, precedingF + delta,
                                                                                    frictionAngle, beta, rp,
                                                                                    dFdSp, dFdcohesionp, dFdPrecedingFp );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, precedingF - delta,
                                                                                    frictionAngle, beta, rm,
                                                                                    dFdSm, dFdcohesionm, dFdPrecedingFm );

        constantVector gradCol = ( dFdSp - dFdSm ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], d2FdSdPrecedingFJ2[j][i] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( testComputeHigherOrderDruckerPragerYieldEquation ){
    /*!
     * Test the computation of the higher order stress Drucker-Prager yield equation.
     *
     */

    variableVector M = { 0.80732114,  0.79202055,  0.17990022,  0.97454675,  0.703207  ,
                        -0.58236697,  0.53324571, -0.93438873, -0.40650796,  0.14071918,
                         0.66933708, -0.67854069, -0.30317772, -0.93821882,  0.97270622,
                         0.00295302, -0.12441126,  0.30539971, -0.0580227 ,  0.89696105,
                         0.17567709, -0.9592962 ,  0.63535407,  0.95437804, -0.64531877,
                         0.69978907,  0.81327586 };

    variableVector precedingF = { 1.03482346, 0.01430697, 0.01134257,
                                  0.02756574, 1.03597345, 0.02115532,
                                  0.04903821, 0.03424149, 1.0240466 };

    variableVector cohesion = { 0.51, 0.71, .14 };
    parameterType frictionAngle = 0.46;
    parameterType beta = 0.34;

    constantVector answer = { 1.08862863, 1.39180044, 1.82194477 };
    variableVector result;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, precedingF,
                                                                                                        frictionAngle, beta, result );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( result, answer ) );

    //Test the Jacobians

    variableVector resultJ;
    variableMatrix dFdStress, dFdc, dFdPrecedingF;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, precedingF,
                                                                                                        frictionAngle, beta, resultJ,
                                                                                                        dFdStress, dFdc, dFdPrecedingF );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ, answer ) );

    variableVector resultJ2;
    variableMatrix dFdStressJ2, dFdcJ2, dFdPrecedingFJ2;
    variableMatrix d2FdStress2J2, d2FdStressdPrecedingFJ2;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, precedingF,
                                                                                                        frictionAngle, beta, resultJ2,
                                                                                                        dFdStressJ2, dFdcJ2, dFdPrecedingFJ2,
                                                                                                        d2FdStress2J2, d2FdStressdPrecedingFJ2 );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( resultJ2, answer ) );

    //Test derivatives w.r.t stress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < M.size(); i++ ){
        constantVector delta( M.size(), 0 );
        delta[i] = eps * fabs( M[i] ) + eps;

        variableVector resultP, resultM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M + delta, cohesion, precedingF,
                                                                                                            frictionAngle, beta, resultP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M - delta, cohesion, precedingF,
                                                                                                            frictionAngle, beta, resultM );

        constantVector gradCol = ( resultP - resultM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dFdStress[j][i] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dFdStressJ2[j][i] ) );
        }

        variableMatrix dFdStressP, dFdcP, dFdPrecedingFP;
        variableMatrix dFdStressM, dFdcM, dFdPrecedingFM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M + delta, cohesion, precedingF,
                                                                                                            frictionAngle, beta, resultP,
                                                                                                            dFdStressP, dFdcP, dFdPrecedingFP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M - delta, cohesion, precedingF,
                                                                                                            frictionAngle, beta, resultM,
                                                                                                            dFdStressM, dFdcM, dFdPrecedingFM );

        constantMatrix gradMat = ( dFdStressP - dFdStressM ) / ( 2 * delta[i] );

        unsigned int n, o, p;

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        n = ( int )( i / 9 );
                        o = ( int )( (i - 9 * n ) / 3 );
                        p = ( i - 9 * n - 3 * o ) % 3;
                        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradMat[ j ][ 9 * k + 3 * l + m ],
                                                        d2FdStress2J2[ j ][ 243 * k + 81 * l + 27 * m + 9 * n + 3 * o + p ] ) );
                    }
                }
            }
        }
    }

    //Test derivatives w.r.t. the cohesion
    for ( unsigned int i = 0; i < cohesion.size(); i++ ){
        constantVector delta( cohesion.size(), 0 );
        delta[i] = eps * fabs( cohesion[i] ) + eps;

        variableVector resultP, resultM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion + delta, precedingF,
                                                                                                            frictionAngle, beta, resultP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion - delta, precedingF,
                                                                                                            frictionAngle, beta, resultM );

        constantVector gradCol = ( resultP - resultM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dFdc[j][i] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dFdcJ2[j][i] ) );
        }
    }

    //Test derivatives w.r.t. the preceding deformation gradient
    for ( unsigned int i = 0; i < precedingF.size(); i++ ){
        constantVector delta( precedingF.size(), 0 );
        delta[i] = eps * fabs( precedingF[i] ) + eps;

        variableVector resultP, resultM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, precedingF + delta,
                                                                                                            frictionAngle, beta, resultP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, precedingF - delta,
                                                                                                            frictionAngle, beta, resultM );

        constantVector gradCol = ( resultP - resultM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dFdPrecedingF[j][i] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dFdPrecedingF[j][i] ) );
        }

        variableMatrix dFdStressP, dFdcP, dFdPrecedingFP;
        variableMatrix dFdStressM, dFdcM, dFdPrecedingFM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, precedingF + delta,
                                                                                                            frictionAngle, beta, resultP,
                                                                                                            dFdStressP, dFdcP, dFdPrecedingFP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, precedingF - delta,
                                                                                                            frictionAngle, beta, resultM,
                                                                                                            dFdStressM, dFdcM, dFdPrecedingFM );

        constantMatrix gradMat = ( dFdStressP - dFdStressM ) / ( 2 * delta[i] );

        unsigned int n, o;

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        n = ( int )( i / 3 );
                        o = ( i % 3 );
                        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradMat[ j ][ 9 * k + 3 * l + m ],
                                                        d2FdStressdPrecedingFJ2[ j ][ 81 * k + 27 * l + 9 * m + 3 * n + o ] ) );
                    }
                }
            }
        }
    }

}

BOOST_AUTO_TEST_CASE( test_setFlowDerivatives ){
    /*!
     * Test setting the cohesion values
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
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

    floatVector parameters = { 2, 0.53895133, 0.37172145,
                               2, 0.37773052, 0.92739145,
                               2, 0.53186824, 0.75454313,
                               2, 0.95338442, 0.74042148,
                               2, 0.38093104, 0.49241325,
                               2, 0.82121039, 0.90566759,
                               2, 0.01166325, 0.05331896,
                               2, 0.32982199, 0.60161431,
                               2, 0.58881096, 0.11473813 };

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
                                   2.18834178e-02,  2.04005542e-02,
                                   0.01, 0.02, 0.03, 0.04, 0.05,
                                   0.06, 0.07, 0.08, 0.09, 0.10 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 10;

    std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

    floatType tolr = 1e-2;

    floatType tola = 1e-3;

    unsigned int maxIterations = 24;

    unsigned int maxLSIterations = 45;

    floatType lsAlpha = 2.3;

    class stressMock : public tardigradeHydra::residualBaseMicromorphic{

        public:

            using tardigradeHydra::residualBaseMicromorphic::residualBaseMicromorphic;

    };

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticity::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::residual;

            floatType macroCohesion = 1.23;

            floatType microCohesion = 2.34;

            floatVector microGradientCohesion = { 3.4, 4.5, 5.6 };

            floatVector macroDrivingStress = {  1,  2,  3,  4,  5,  6,  7,  8,  9 };

            floatVector microDrivingStress = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

            floatVector microGradientDrivingStress = { 19, 20, 21, 22, 23, 24, 25, 26, 27,
                                                       28, 29, 30, 31, 32, 33, 34, 35, 36,
                                                       37, 38, 39, 40, 41, 42, 43, 44, 45 };

            floatType previousMacroCohesion = -1.23;

            floatType previousMicroCohesion = -2.34;

            floatVector previousMicroGradientCohesion = { -3.4, -4.5, -5.6 };

            floatVector previousMacroDrivingStress = {  -1,  -2,  -3,  -4,  -5,  -6,  -7,  -8,  -9 };

            floatVector previousMicroDrivingStress = { -10, -11, -12, -13, -14, -15, -16, -17, -18 };

            floatVector previousMicroGradientDrivingStress = { -19, -20, -21, -22, -23, -24, -25, -26, -27,
                                                               -28, -29, -30, -31, -32, -33, -34, -35, -36,
                                                               -37, -38, -39, -40, -41, -42, -43, -44, -45 };

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              2, 0.95338442, 0.74042148,
                                              2, 0.38093104, 0.49241325,
                                              2, 0.82121039, 0.90566759,
                                              2, 0.01166325, 0.05331896,
                                              2, 0.32982199, 0.60161431,
                                              2, 0.58881096, 0.11473813 };

            virtual void setCohesions( const bool isPrevious ) override{

                if ( isPrevious ){

                    set_previousMacroCohesion( previousMacroCohesion );

                    set_previousMicroCohesion( previousMicroCohesion );

                    set_previousMicroGradientCohesion( previousMicroGradientCohesion );

                }
                else{

                    set_macroCohesion( macroCohesion );

                    set_microCohesion( microCohesion );

                    set_microGradientCohesion( microGradientCohesion );

                }

            }

            virtual void setDrivingStresses( const bool isPrevious ) override{

                if ( isPrevious ){

                    set_previousMacroDrivingStress( previousMacroDrivingStress );

                    set_previousSymmetricMicroDrivingStress( previousMicroDrivingStress );

                    set_previousHigherOrderDrivingStress( previousMicroGradientDrivingStress );

                }
                else{

                    set_macroDrivingStress( macroDrivingStress );

                    set_symmetricMicroDrivingStress( microDrivingStress );

                    set_higherOrderDrivingStress( microGradientDrivingStress );

                }

            }

    };

    class hydraBaseMicromorphicMock : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

            stressMock elasticity;

            residualMock plasticity;

            std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              2, 0.95338442, 0.74042148, 
                                              2, 0.38093104, 0.49241325, 
                                              2, 0.82121039, 0.90566759,
                                              2, 0.01166325, 0.05331896, 
                                              2, 0.32982199, 0.60161431, 
                                              2, 0.58881096, 0.11473813 };

        private:

            using tardigradeHydra::hydraBaseMicromorphic::setResidualClasses;

            virtual void setResidualClasses( ) override{

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = stressMock( this, 45 );

                plasticity = residualMock( this, 55, 1, stateVariableIndices, plasticParameters );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                setResidualClasses( residuals );

            }

    };

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 55, 1, stateVariableIndices, parameters );

    residualMock RJ( &hydra, 55, 1, stateVariableIndices, parameters );

    floatType tempYield, dMacroFlowdCohesion, dMicroFlowdCohesion;

    floatVector tempVectorYield, dMacroFlowdDrivingStress, dMicroFlowdDrivingStress, dMacroFlowdPrecedingF, dMicroFlowdPrecedingF;

    floatMatrix dMicroGradientFlowdDrivingStress, dMicroGradientFlowdCohesion, dMicroGradientFlowdPrecedingF;

    floatType previousdMacroFlowdCohesion, previousdMicroFlowdCohesion;

    floatVector previousdMacroFlowdDrivingStress, previousdMicroFlowdDrivingStress, previousdMacroFlowdPrecedingF, previousdMicroFlowdPrecedingF;

    floatMatrix previousdMicroGradientFlowdDrivingStress, previousdMicroGradientFlowdCohesion, previousdMicroGradientFlowdPrecedingF;

    floatVector Fp = floatVector( unknownVector.begin( ) + configuration_unknown_count,
                                  unknownVector.begin( ) + configuration_unknown_count + 9 );

    floatVector eye( 9 );
    tardigradeVectorTools::eye( eye );

    floatVector previousFp = floatVector( previousStateVariables.begin( ),
                                          previousStateVariables.begin( ) + 9 ) + eye;

    floatVector precedingDeformationGradient = tardigradeVectorTools::matrixMultiply( deformationGradient, tardigradeVectorTools::inverse( Fp, 3, 3 ), 3, 3, 3, 3 );

    floatVector previousPrecedingDeformationGradient = tardigradeVectorTools::matrixMultiply( previousDeformationGradient, tardigradeVectorTools::inverse( previousFp, 3, 3 ), 3, 3, 3, 3 );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( R.macroDrivingStress, R.macroCohesion, precedingDeformationGradient,
                                                                                                        R.plasticParameters[ 10 ], R.plasticParameters[ 11 ],
                                                                                                        tempYield, dMacroFlowdDrivingStress, dMacroFlowdCohesion, dMacroFlowdPrecedingF );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( R.microDrivingStress, R.microCohesion, precedingDeformationGradient,
                                                                                                        R.plasticParameters[ 13 ], R.plasticParameters[ 14 ],
                                                                                                        tempYield, dMicroFlowdDrivingStress, dMicroFlowdCohesion, dMicroFlowdPrecedingF );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( R.microGradientDrivingStress, R.microGradientCohesion, precedingDeformationGradient,
                                                                                                        R.plasticParameters[ 16 ], R.plasticParameters[ 17 ],
                                                                                                        tempVectorYield, dMicroGradientFlowdDrivingStress, dMicroGradientFlowdCohesion, dMicroGradientFlowdPrecedingF );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( R.previousMacroDrivingStress, R.previousMacroCohesion, previousPrecedingDeformationGradient,
                                                                                                        R.plasticParameters[ 10 ], R.plasticParameters[ 11 ],
                                                                                                        tempYield, previousdMacroFlowdDrivingStress, previousdMacroFlowdCohesion, previousdMacroFlowdPrecedingF );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( R.previousMicroDrivingStress, R.previousMicroCohesion, previousPrecedingDeformationGradient,
                                                                                                        R.plasticParameters[ 13 ], R.plasticParameters[ 14 ],
                                                                                                        tempYield, previousdMicroFlowdDrivingStress, previousdMicroFlowdCohesion, previousdMicroFlowdPrecedingF );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( R.previousMicroGradientDrivingStress, R.previousMicroGradientCohesion, previousPrecedingDeformationGradient,
                                                                                                        R.plasticParameters[ 16 ], R.plasticParameters[ 17 ],
                                                                                                        tempVectorYield, previousdMicroGradientFlowdDrivingStress, previousdMicroGradientFlowdCohesion, previousdMicroGradientFlowdPrecedingF );

    RJ.get_d2MacroFlowdDrivingStressdStress( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMacroFlowdCohesion,                      *R.get_dMacroFlowdc( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroFlowdCohesion,                      *R.get_dMicroFlowdc( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroGradientFlowdCohesion,              *R.get_dMicroGradientFlowdc( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMacroFlowdDrivingStress,                 *R.get_dMacroFlowdDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroFlowdDrivingStress,                 *R.get_dMicroFlowdDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroGradientFlowdDrivingStress,         *R.get_dMicroGradientFlowdDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMacroFlowdCohesion,              *R.get_previousdMacroFlowdc( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroFlowdCohesion,              *R.get_previousdMicroFlowdc( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroGradientFlowdCohesion,      *R.get_previousdMicroGradientFlowdc( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMacroFlowdDrivingStress,         *R.get_previousdMacroFlowdDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroFlowdDrivingStress,         *R.get_previousdMicroFlowdDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroGradientFlowdDrivingStress, *R.get_previousdMicroGradientFlowdDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMacroFlowdCohesion,                      *RJ.get_dMacroFlowdc( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroFlowdCohesion,                      *RJ.get_dMicroFlowdc( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroGradientFlowdCohesion,              *RJ.get_dMicroGradientFlowdc( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMacroFlowdDrivingStress,                 *RJ.get_dMacroFlowdDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroFlowdDrivingStress,                 *RJ.get_dMicroFlowdDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroGradientFlowdDrivingStress,         *RJ.get_dMicroGradientFlowdDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMacroFlowdCohesion,              *RJ.get_previousdMacroFlowdc( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroFlowdCohesion,              *RJ.get_previousdMicroFlowdc( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroGradientFlowdCohesion,      *RJ.get_previousdMicroGradientFlowdc( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMacroFlowdDrivingStress,         *RJ.get_previousdMacroFlowdDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroFlowdDrivingStress,         *RJ.get_previousdMicroFlowdDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroGradientFlowdDrivingStress, *RJ.get_previousdMicroGradientFlowdDrivingStress( ) ) );
}
