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

            floatMatrix dMacroDrivingStressdMacroStress = initialize( 9, 9 );

            floatMatrix dMacroDrivingStressdF           = initialize( 9, 9 );

            floatMatrix dMacroDrivingStressdFn          = initialize( 9, 9 );

            floatMatrix dMicroDrivingStressdMicroStress = initialize( 9, 9 );

            floatMatrix dMicroDrivingStressdF           = initialize( 9, 9 );

            floatMatrix dMicroDrivingStressdFn          = initialize( 9, 9 );

            floatMatrix dMicroGradientDrivingStressdHigherOrderStress = initialize( 27, 27 );

            floatMatrix dMicroGradientDrivingStressdF                 = initialize( 27,  9 );

            floatMatrix dMicroGradientDrivingStressdFn                = initialize( 27,  9 );

            floatMatrix dMicroGradientDrivingStressdChi               = initialize( 27,  9 );

            floatMatrix dMicroGradientDrivingStressdChin              = initialize( 27,  9 );

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

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

            }

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

            virtual void setDrivingStressesJacobians( const bool isPrevious ) override{

                setDrivingStresses( isPrevious );

                if ( isPrevious ){

                    set_previousdMacroDrivingStressdMacroStress( dMacroDrivingStressdMacroStress );

                    set_previousdMacroDrivingStressdF( dMacroDrivingStressdF );

                    set_previousdMacroDrivingStressdFn( dMacroDrivingStressdFn );

                    set_previousdSymmetricMicroDrivingStressdMicroStress( dMicroDrivingStressdMicroStress );

                    set_previousdSymmetricMicroDrivingStressdF( dMicroDrivingStressdF );

                    set_previousdSymmetricMicroDrivingStressdFn( dMicroDrivingStressdFn );

                    set_previousdHigherOrderDrivingStressdHigherOrderStress( dMicroGradientDrivingStressdHigherOrderStress );

                    set_previousdHigherOrderDrivingStressdF( dMicroGradientDrivingStressdF );

                    set_previousdHigherOrderDrivingStressdFn( dMicroGradientDrivingStressdFn );

                    set_previousdHigherOrderDrivingStressdChi( dMicroGradientDrivingStressdChi );

                    set_previousdHigherOrderDrivingStressdChin( dMicroGradientDrivingStressdChin );

                }
                else{

                    set_dMacroDrivingStressdMacroStress( dMacroDrivingStressdMacroStress );

                    set_dMacroDrivingStressdF( dMacroDrivingStressdF );

                    set_dMacroDrivingStressdFn( dMacroDrivingStressdFn );

                    set_dSymmetricMicroDrivingStressdMicroStress( dMicroDrivingStressdMicroStress );

                    set_dSymmetricMicroDrivingStressdF( dMicroDrivingStressdF );

                    set_dSymmetricMicroDrivingStressdFn( dMicroDrivingStressdFn );

                    set_dHigherOrderDrivingStressdHigherOrderStress( dMicroGradientDrivingStressdHigherOrderStress );

                    set_dHigherOrderDrivingStressdF( dMicroGradientDrivingStressdF );

                    set_dHigherOrderDrivingStressdFn( dMicroGradientDrivingStressdFn );

                    set_dHigherOrderDrivingStressdChi( dMicroGradientDrivingStressdChi );

                    set_dHigherOrderDrivingStressdChin( dMicroGradientDrivingStressdChin );
                
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
    RJ.get_previousd2MacroFlowdDrivingStressdStress( );

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

BOOST_AUTO_TEST_CASE( test_setFlowDerivatives2 ){
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

            floatVector previousPK2Stress = {  1,  2,  3,  4,  5,  6,  7,  8,  9 };

            floatVector previousSIGMA     = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

            floatVector previousM         = { 19, 20, 21, 22, 23, 24, 25, 26, 27,
                                              28, 29, 30, 31, 32, 33, 34, 35, 36,
                                              37, 38, 39, 40, 41, 42, 43, 44, 45 };

        protected:

            using tardigradeHydra::residualBase::setPreviousStress;

            virtual void setPreviousStress( ) override{

                setPreviousStress( tardigradeVectorTools::appendVectors( { previousPK2Stress, previousSIGMA, previousM } ) );

            }

    };

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticity::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::residual;

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

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

            floatVector _local_deltaPK2   = initialize(  9 );

            floatVector _local_deltaSIGMA = initialize(  9 );

            floatVector _local_deltaM     = initialize( 27 );

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

            }

        private:

            using tardigradeHydra::hydraBaseMicromorphic::setResidualClasses;

            virtual void setResidualClasses( ) override{

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = stressMock( this, 45 );

                elasticity.previousPK2Stress += _local_deltaPK2;

                elasticity.previousSIGMA     += _local_deltaSIGMA;

                elasticity.previousM         += _local_deltaM;

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

    // Jacobian tests

    floatMatrix d2MacroFlowdDrivingStressdX( 9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix d2MicroFlowdDrivingStressdX( 9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix d2MicroGradientFlowdDrivingStressdX( 3 * 27, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix d2MacroFlowdDrivingStressdF( 9, floatVector( 9, 0 ) );

    floatMatrix d2MicroFlowdDrivingStressdF( 9, floatVector( 9, 0 ) );

    floatMatrix d2MicroGradientFlowdDrivingStressdF( 3 * 27, floatVector( 9, 0 ) );

    floatMatrix d2MacroFlowdDrivingStressdChi( 9, floatVector( 9, 0 ) );

    floatMatrix d2MicroFlowdDrivingStressdChi( 9, floatVector( 9, 0 ) );

    floatMatrix d2MicroGradientFlowdDrivingStressdChi( 3 * 27, floatVector( 9, 0 ) );

    floatMatrix previousd2MacroFlowdDrivingStressdMacroStress( 9, floatVector( 9, 0 ) );

    floatMatrix previousd2MicroFlowdDrivingStressdMicroStress( 9, floatVector( 9, 0 ) );

    floatMatrix previousd2MicroGradientFlowdDrivingStressdHigherOrderStress( 3 * 27, floatVector( 27, 0 ) );

    floatMatrix previousd2MacroFlowdDrivingStressdX( 9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix previousd2MicroFlowdDrivingStressdX( 9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix previousd2MicroGradientFlowdDrivingStressdX( 3 * 27, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix previousd2MacroFlowdDrivingStressdF( 9, floatVector( 9, 0 ) );

    floatMatrix previousd2MicroFlowdDrivingStressdF( 9, floatVector( 9, 0 ) );

    floatMatrix previousd2MicroGradientFlowdDrivingStressdF( 3 * 27, floatVector( 9, 0 ) );

    floatType eps = 1e-6;

    // derivatives w.r.t. current stresses, Fn, Chin, gradchin, and ISVs
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

        floatVector vp = *Rp.get_dMacroFlowdDrivingStress( );

        floatVector vm = *Rm.get_dMacroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            d2MacroFlowdDrivingStressdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_dMicroFlowdDrivingStress( );

        vm = *Rm.get_dMicroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            d2MicroFlowdDrivingStressdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = tardigradeVectorTools::appendVectors( *Rp.get_dMicroGradientFlowdDrivingStress( ) );

        vm = tardigradeVectorTools::appendVectors( *Rm.get_dMicroGradientFlowdDrivingStress( ) );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            d2MicroGradientFlowdDrivingStressdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_d2MacroFlowdDrivingStressdX(              9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix assembled_d2MicroFlowdDrivingStressdX(              9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix assembled_d2MicroGradientFlowdDrivingStressdX( 3 * 27, floatVector( unknownVector.size( ), 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_d2MacroFlowdDrivingStressdX[ i ][ j ] = ( *R.get_d2MacroFlowdDrivingStressdStress( ) )[ i ][ j ];

        }

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_d2MacroFlowdDrivingStressdX[ i ][ j + configuration_unknown_count ] = ( *R.get_d2MacroFlowdDrivingStressdFn( ) )[ i ][ j ];

        }

    }

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_d2MicroFlowdDrivingStressdX[ i ][ j + 9 ] = ( *R.get_d2MicroFlowdDrivingStressdStress( ) )[ i ][ j ];

        }

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_d2MicroFlowdDrivingStressdX[ i ][ j + configuration_unknown_count ] = ( *R.get_d2MicroFlowdDrivingStressdFn( ) )[ i ][ j ];

        }

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 27; j++ ){

            for ( unsigned int k = 0; k < 27; k++ ){

                assembled_d2MicroGradientFlowdDrivingStressdX[ 27 * i + j ][ k + 18 ] = ( *R.get_d2MicroGradientFlowdDrivingStressdStress( ) )[ i ][ 27 * j + k ];

            }

            for ( unsigned int k = 0; k < 9; k++ ){

                assembled_d2MicroGradientFlowdDrivingStressdX[ 27 * i + j ][ k + configuration_unknown_count     ] = ( *R.get_d2MicroGradientFlowdDrivingStressdFn( ) )[ i ][ 9 * j + k ];

                assembled_d2MicroGradientFlowdDrivingStressdX[ 27 * i + j ][ k + configuration_unknown_count + 9 ] = ( *R.get_d2MicroGradientFlowdDrivingStressdChin( ) )[ i ][ 9 * j + k ];

            }

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assembled_d2MacroFlowdDrivingStressdX, d2MacroFlowdDrivingStressdX ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assembled_d2MicroFlowdDrivingStressdX, d2MicroFlowdDrivingStressdX ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assembled_d2MicroGradientFlowdDrivingStressdX, d2MicroGradientFlowdDrivingStressdX ) );

    // derivatives w.r.t. current F
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

        floatVector vp = *Rp.get_dMacroFlowdDrivingStress( );

        floatVector vm = *Rm.get_dMacroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            d2MacroFlowdDrivingStressdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_dMicroFlowdDrivingStress( );

        vm = *Rm.get_dMicroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            d2MicroFlowdDrivingStressdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = tardigradeVectorTools::appendVectors( *Rp.get_dMicroGradientFlowdDrivingStress( ) );

        vm = tardigradeVectorTools::appendVectors( *Rm.get_dMicroGradientFlowdDrivingStress( ) );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            d2MicroGradientFlowdDrivingStressdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_d2MicroGradientFlowdDrivingStressdF( 3 * 27, floatVector( 9, 0 ) );

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 27; j++ ){

            for ( unsigned int k = 0; k < 9; k++ ){

                assembled_d2MicroGradientFlowdDrivingStressdF[ 27 * i + j ][ k ] = ( *R.get_d2MicroGradientFlowdDrivingStressdF( ) )[ i ][ 9 * j + k ];

            }

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( d2MacroFlowdDrivingStressdF, *R.get_d2MacroFlowdDrivingStressdF( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( d2MicroFlowdDrivingStressdF, *R.get_d2MicroFlowdDrivingStressdF( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( d2MicroGradientFlowdDrivingStressdF, assembled_d2MicroGradientFlowdDrivingStressdF ) );

    // derivatives w.r.t. current chi
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

        floatVector vp = *Rp.get_dMacroFlowdDrivingStress( );

        floatVector vm = *Rm.get_dMacroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            d2MacroFlowdDrivingStressdChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_dMicroFlowdDrivingStress( );

        vm = *Rm.get_dMicroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            d2MicroFlowdDrivingStressdChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = tardigradeVectorTools::appendVectors( *Rp.get_dMicroGradientFlowdDrivingStress( ) );

        vm = tardigradeVectorTools::appendVectors( *Rm.get_dMicroGradientFlowdDrivingStress( ) );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            d2MicroGradientFlowdDrivingStressdChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_d2MicroGradientFlowdDrivingStressdChi( 3 * 27, floatVector( 9, 0 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( d2MacroFlowdDrivingStressdChi, floatMatrix( 9, floatVector( 9, 0 ) ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( d2MicroFlowdDrivingStressdChi, floatMatrix( 9, floatVector( 9, 0 ) ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( d2MicroGradientFlowdDrivingStressdChi, assembled_d2MicroGradientFlowdDrivingStressdChi ) );

    // derivatives w.r.t. previous macro stress
    for ( unsigned int i = 0; i < hydra.elasticity.previousPK2Stress.size( ); i++ ){

        floatVector delta( hydra.elasticity.previousPK2Stress.size( ), 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousPK2Stress[ i ] ) + eps;

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

        hydrap._local_deltaPK2 =  delta;

        hydram._local_deltaPK2 = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm( &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousdMacroFlowdDrivingStress( );

        floatVector vm = *Rm.get_previousdMacroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousd2MacroFlowdDrivingStressdMacroStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousdMicroFlowdDrivingStress( );

        vm = *Rm.get_previousdMicroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        vp = tardigradeVectorTools::appendVectors( *Rp.get_previousdMicroGradientFlowdDrivingStress( ) );

        vm = tardigradeVectorTools::appendVectors( *Rm.get_previousdMicroGradientFlowdDrivingStress( ) );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousd2MacroFlowdDrivingStressdMacroStress, *R.get_previousd2MacroFlowdDrivingStressdStress( ) ) );

    // derivatives w.r.t. previous micro stress
    for ( unsigned int i = 0; i < hydra.elasticity.previousSIGMA.size( ); i++ ){

        floatVector delta( hydra.elasticity.previousSIGMA.size( ), 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousSIGMA[ i ] ) + eps;

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

        hydrap._local_deltaSIGMA =  delta;

        hydram._local_deltaSIGMA = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm( &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousdMacroFlowdDrivingStress( );

        floatVector vm = *Rm.get_previousdMacroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        vp = *Rp.get_previousdMicroFlowdDrivingStress( );

        vm = *Rm.get_previousdMicroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousd2MicroFlowdDrivingStressdMicroStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = tardigradeVectorTools::appendVectors( *Rp.get_previousdMicroGradientFlowdDrivingStress( ) );

        vm = tardigradeVectorTools::appendVectors( *Rm.get_previousdMicroGradientFlowdDrivingStress( ) );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousd2MicroFlowdDrivingStressdMicroStress, *R.get_previousd2MicroFlowdDrivingStressdStress( ) ) );

    // derivatives w.r.t. previous micro gradient stress
    for ( unsigned int i = 0; i < hydra.elasticity.previousM.size( ); i++ ){

        floatVector delta( hydra.elasticity.previousM.size( ), 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousM[ i ] ) + eps;

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

        hydrap._local_deltaM =  delta;

        hydram._local_deltaM = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm( &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousdMacroFlowdDrivingStress( );

        floatVector vm = *Rm.get_previousdMacroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        vp = *Rp.get_previousdMicroFlowdDrivingStress( );

        vm = *Rm.get_previousdMicroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        vp = tardigradeVectorTools::appendVectors( *Rp.get_previousdMicroGradientFlowdDrivingStress( ) );

        vm = tardigradeVectorTools::appendVectors( *Rm.get_previousdMicroGradientFlowdDrivingStress( ) );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousd2MicroGradientFlowdDrivingStressdHigherOrderStress[ j ][ i ] =  ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_previousd2MicroGradientFlowdDrivingStressdStress( 3 * 27, floatVector( 27, 0 ) );

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 27; j++ ){

            for ( unsigned int k = 0; k < 27; k++ ){

                assembled_previousd2MicroGradientFlowdDrivingStressdStress[ 27 * i + j ][ k ] += ( *R.get_previousd2MicroGradientFlowdDrivingStressdStress( ) )[ i ][ 27 * j + k ];

            }

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousd2MicroGradientFlowdDrivingStressdHigherOrderStress, assembled_previousd2MicroGradientFlowdDrivingStressdStress ) );

    // derivatives w.r.t. previous, Fn, Chin, gradchin, and ISVs
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

        floatVector vp = *Rp.get_previousdMacroFlowdDrivingStress( );

        floatVector vm = *Rm.get_previousdMacroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousd2MacroFlowdDrivingStressdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousdMicroFlowdDrivingStress( );

        vm = *Rm.get_previousdMicroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousd2MicroFlowdDrivingStressdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = tardigradeVectorTools::appendVectors( *Rp.get_previousdMicroGradientFlowdDrivingStress( ) );

        vm = tardigradeVectorTools::appendVectors( *Rm.get_previousdMicroGradientFlowdDrivingStress( ) );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousd2MicroGradientFlowdDrivingStressdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_previousd2MacroFlowdDrivingStressdX(              9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix assembled_previousd2MicroFlowdDrivingStressdX(              9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix assembled_previousd2MicroGradientFlowdDrivingStressdX( 3 * 27, floatVector( previousStateVariables.size( ), 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_previousd2MacroFlowdDrivingStressdX[ i ][ j ] = ( *R.get_previousd2MacroFlowdDrivingStressdFn( ) )[ i ][ j ];

        }

    }

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_previousd2MicroFlowdDrivingStressdX[ i ][ j ] = ( *R.get_previousd2MicroFlowdDrivingStressdFn( ) )[ i ][ j ];

        }

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 27; j++ ){

            for ( unsigned int k = 0; k < 9; k++ ){

                assembled_previousd2MicroGradientFlowdDrivingStressdX[ 27 * i + j ][ k     ] = ( *R.get_previousd2MicroGradientFlowdDrivingStressdFn( ) )[ i ][ 9 * j + k ];

                assembled_previousd2MicroGradientFlowdDrivingStressdX[ 27 * i + j ][ k + 9 ] = ( *R.get_previousd2MicroGradientFlowdDrivingStressdChin( ) )[ i ][ 9 * j + k ];

            }

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assembled_previousd2MacroFlowdDrivingStressdX, previousd2MacroFlowdDrivingStressdX ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assembled_previousd2MicroFlowdDrivingStressdX, previousd2MicroFlowdDrivingStressdX ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assembled_previousd2MicroGradientFlowdDrivingStressdX, previousd2MicroGradientFlowdDrivingStressdX ) );

    // derivatives w.r.t. previous F
    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

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

        floatVector vp = *Rp.get_previousdMacroFlowdDrivingStress( );

        floatVector vm = *Rm.get_previousdMacroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousd2MacroFlowdDrivingStressdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousdMicroFlowdDrivingStress( );

        vm = *Rm.get_previousdMicroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousd2MicroFlowdDrivingStressdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = tardigradeVectorTools::appendVectors( *Rp.get_previousdMicroGradientFlowdDrivingStress( ) );

        vm = tardigradeVectorTools::appendVectors( *Rm.get_previousdMicroGradientFlowdDrivingStress( ) );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousd2MicroGradientFlowdDrivingStressdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_previousd2MicroGradientFlowdDrivingStressdF( 3 * 27, floatVector( 9, 0 ) );

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 27; j++ ){

            for ( unsigned int k = 0; k < 9; k++ ){

                assembled_previousd2MicroGradientFlowdDrivingStressdF[ 27 * i + j ][ k ] = ( *R.get_previousd2MicroGradientFlowdDrivingStressdF( ) )[ i ][ 9 * j + k ];

            }

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousd2MacroFlowdDrivingStressdF, *R.get_previousd2MacroFlowdDrivingStressdF( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousd2MicroFlowdDrivingStressdF, *R.get_previousd2MicroFlowdDrivingStressdF( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousd2MicroGradientFlowdDrivingStressdF, assembled_previousd2MicroGradientFlowdDrivingStressdF ) );

}

BOOST_AUTO_TEST_CASE( test_setPlasticStrainLikeISVEvolutionRates ){
    /*!
     * Test setting the strain-like ISV evolution rates
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

            floatType dMacroFlowdc = 1.23;

            floatType dMicroFlowdc = 2.34;

            floatMatrix dMicroGradientFlowdc = { { 0.34, 0.45, 0.56 },
                                                 { 0.67, 0.78, 0.89 },
                                                 { 0.91, 1.01, 1.12 } };

            floatType previousdMacroFlowdc = -1.03;

            floatType previousdMicroFlowdc = -2.04;

            floatMatrix previousdMicroGradientFlowdc = { { -3.4,  -4.5,  -5.6 },
                                                         { -6.7,  -7.8,  -8.9 },
                                                         { -9.1, -10.1, -11.2 } };

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              2, 0.95338442, 0.74042148,
                                              2, 0.38093104, 0.49241325,
                                              2, 0.82121039, 0.90566759,
                                              2, 0.01166325, 0.05331896,
                                              2, 0.32982199, 0.60161431,
                                              2, 0.58881096, 0.11473813 };

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

            }

            virtual void setFlowPotentialGradients( const bool isPrevious ) override{


                if ( isPrevious ){

                    set_previousdMacroFlowdc( previousdMacroFlowdc );

                    set_previousdMicroFlowdc( previousdMicroFlowdc );

                    set_previousdMicroGradientFlowdc( previousdMicroGradientFlowdc );

                }
                else{

                    set_dMacroFlowdc( dMacroFlowdc );

                    set_dMicroFlowdc( dMicroFlowdc );

                    set_dMicroGradientFlowdc( dMicroGradientFlowdc );

                }

            }

            virtual void setFlowPotentialGradientsJacobians( const bool isPrevious ) override{


                if ( isPrevious ){

                    set_previousdMacroFlowdc( previousdMacroFlowdc );

                    set_previousdMicroFlowdc( previousdMicroFlowdc );

                    set_previousdMicroGradientFlowdc( previousdMicroGradientFlowdc );

                }
                else{

                    set_dMacroFlowdc( dMacroFlowdc );

                    set_dMicroFlowdc( dMicroFlowdc );

                    set_dMicroGradientFlowdc( dMicroGradientFlowdc );

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

    floatVector evolutionRates = { -1.23 * 0.01,
                                   -2.34 * 0.02,
                                   -( 0.34 * 0.03 + 0.45 * 0.04 + 0.56 * 0.05 ),
                                   -( 0.67 * 0.03 + 0.78 * 0.04 + 0.89 * 0.05 ),
                                   -( 0.91 * 0.03 + 1.01 * 0.04 + 1.12 * 0.05 ) };

    floatVector previousEvolutionRates = { 1.03 * 0.1,
                                           2.04 * 0.2,
                                           ( 3.4 * 0.3 +  4.5 * 0.4 +  5.6 * 0.5 ),
                                           ( 6.7 * 0.3 +  7.8 * 0.4 +  8.9 * 0.5 ),
                                           ( 9.1 * 0.3 + 10.1 * 0.4 + 11.2 * 0.5 ) };

    RJ.get_dPlasticStrainLikeISVEvolutionRatesdStateVariables( );
    RJ.get_previousdPlasticStrainLikeISVEvolutionRatesdStateVariables( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( evolutionRates, *R.get_plasticStrainLikeISVEvolutionRates( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousEvolutionRates, *R.get_previousPlasticStrainLikeISVEvolutionRates( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( evolutionRates, *RJ.get_plasticStrainLikeISVEvolutionRates( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousEvolutionRates, *RJ.get_previousPlasticStrainLikeISVEvolutionRates( ) ) );

}

BOOST_AUTO_TEST_CASE( test_setPlasticStrainLikeISVEvolutionRates2 ){
    /*!
     * Test setting the strain-like ISV evolution rates
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

            floatVector previousPK2Stress = {  1,  2,  3,  4,  5,  6,  7,  8,  9 };

            floatVector previousSIGMA     = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

            floatVector previousM         = { 19, 20, 21, 22, 23, 24, 25, 26, 27,
                                              28, 29, 30, 31, 32, 33, 34, 35, 36,
                                              37, 38, 39, 40, 41, 42, 43, 44, 45 };

        protected:

            using tardigradeHydra::residualBase::setPreviousStress;

            virtual void setPreviousStress( ) override{

                setPreviousStress( tardigradeVectorTools::appendVectors( { previousPK2Stress, previousSIGMA, previousM } ) );

            }

    };

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticity::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::residual;

            floatType dMacroFlowdc = 1.23;

            floatType dMicroFlowdc = 2.34;

            floatMatrix dMicroGradientFlowdc = { { 0.34, 0.45, 0.56 },
                                                 { 0.67, 0.78, 0.89 },
                                                 { 0.91, 1.01, 1.12 } };

            floatType previousdMacroFlowdc = -1.03;

            floatType previousdMicroFlowdc = -2.04;

            floatMatrix previousdMicroGradientFlowdc = { { -3.4,  -4.5,  -5.6 },
                                                         { -6.7,  -7.8,  -8.9 },
                                                         { -9.1, -10.1, -11.2 } };

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              2, 0.95338442, 0.74042148,
                                              2, 0.38093104, 0.49241325,
                                              2, 0.82121039, 0.90566759,
                                              2, 0.01166325, 0.05331896,
                                              2, 0.32982199, 0.60161431,
                                              2, 0.58881096, 0.11473813 };

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

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

    // Set the Jacobians

    floatMatrix dEvolutionRatesdX( 5, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix previousdEvolutionRatesdX( 5, floatVector( previousStateVariables.size( ), 0 ) );

    floatType eps = 1e-6;

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

        floatVector vp = *Rp.get_plasticStrainLikeISVEvolutionRates( );

        floatVector vm = *Rm.get_plasticStrainLikeISVEvolutionRates( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dEvolutionRatesdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_dEvolutionRatesdX( 5, floatVector( unknownVector.size( ), 0 ) );

    for ( unsigned int i = 0; i < 5; i++ ){

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dEvolutionRatesdX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dPlasticStrainLikeISVEvolutionRatesdStateVariables( ) )[ i ][ j ];

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dEvolutionRatesdX, assembled_dEvolutionRatesdX ) );

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

        floatVector vp = *Rp.get_previousPlasticStrainLikeISVEvolutionRates( );

        floatVector vm = *Rm.get_previousPlasticStrainLikeISVEvolutionRates( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdEvolutionRatesdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_previousdEvolutionRatesdX( 5, floatVector( previousStateVariables.size( ), 0 ) );

    for ( unsigned int i = 0; i < 5; i++ ){

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_previousdEvolutionRatesdX[ i ][ j + configuration_unknown_count ] = ( *R.get_previousdPlasticStrainLikeISVEvolutionRatesdStateVariables( ) )[ i ][ j ];

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdEvolutionRatesdX, assembled_previousdEvolutionRatesdX ) );

}

BOOST_AUTO_TEST_CASE( test_setUpdatedPlasticStrainLikeISVs ){
    /*!
     * Test setting the strain-like ISV evolution rates
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

            floatVector evolutionRates = { 1, 2, 3, 4, 5 };

            floatVector previousEvolutionRates = { 0.1, 0.2, 0.3, 0.4, 0.5 };

            floatMatrix dPlasticStrainLikeISVEvolutionRatesdStateVariables = initialize( 5, 10 );

            floatMatrix previousdPlasticStrainLikeISVEvolutionRatesdStateVariables = initialize( 5, 10 );

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              2, 0.95338442, 0.74042148,
                                              2, 0.38093104, 0.49241325,
                                              2, 0.82121039, 0.90566759,
                                              2, 0.01166325, 0.05331896,
                                              2, 0.32982199, 0.60161431,
                                              2, 0.58881096, 0.11473813 };

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

            }

            virtual void setPlasticStrainLikeISVEvolutionRates( const bool isPrevious ) override{

                if ( isPrevious ){

                    set_previousPlasticStrainLikeISVEvolutionRates( previousEvolutionRates );

                }
                else{

                    set_plasticStrainLikeISVEvolutionRates( evolutionRates );

                }

            }

            virtual void setPlasticStrainLikeISVEvolutionRatesJacobians( const bool isPrevious ) override{

                if ( isPrevious ){

                    set_previousPlasticStrainLikeISVEvolutionRates( previousEvolutionRates );

                    set_previousdPlasticStrainLikeISVEvolutionRatesdStateVariables( previousdPlasticStrainLikeISVEvolutionRatesdStateVariables );

                }
                else{

                    set_plasticStrainLikeISVEvolutionRates( evolutionRates );

                    set_dPlasticStrainLikeISVEvolutionRatesdStateVariables( dPlasticStrainLikeISVEvolutionRatesdStateVariables );

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

    RJ.get_dUpdatedPlasticStrainLikeISVsdStateVariables( );

    floatVector updatedPlasticISVs = 0.5 * ( *R.get_plasticStrainLikeISVEvolutionRates( ) + *R.get_previousPlasticStrainLikeISVEvolutionRates( ) ) * deltaTime
                                   + floatVector( previousStateVariables.end( ) - 5, previousStateVariables.end( ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( updatedPlasticISVs, *R.get_updatedPlasticStrainLikeISVs( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( updatedPlasticISVs, *RJ.get_updatedPlasticStrainLikeISVs( ) ) );

}

BOOST_AUTO_TEST_CASE( test_setUpdatedPlasticStrainLikeISVs2 ){
    /*!
     * Test setting the strain-like ISV evolution rates
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

            floatVector previousPK2Stress = {  1,  2,  3,  4,  5,  6,  7,  8,  9 };

            floatVector previousSIGMA     = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

            floatVector previousM         = { 19, 20, 21, 22, 23, 24, 25, 26, 27,
                                              28, 29, 30, 31, 32, 33, 34, 35, 36,
                                              37, 38, 39, 40, 41, 42, 43, 44, 45 };

        protected:

            using tardigradeHydra::residualBase::setPreviousStress;

            virtual void setPreviousStress( ) override{

                setPreviousStress( tardigradeVectorTools::appendVectors( { previousPK2Stress, previousSIGMA, previousM } ) );

            }

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

            floatMatrix dMacroDrivingStressdMacroStress = initialize( 9, 9 );

            floatMatrix dMacroDrivingStressdF           = initialize( 9, 9 );

            floatMatrix dMacroDrivingStressdFn          = initialize( 9, 9 );

            floatMatrix dMicroDrivingStressdMicroStress = initialize( 9, 9 );

            floatMatrix dMicroDrivingStressdF           = initialize( 9, 9 );

            floatMatrix dMicroDrivingStressdFn          = initialize( 9, 9 );

            floatMatrix dMicroGradientDrivingStressdHigherOrderStress = initialize( 27, 27 );

            floatMatrix dMicroGradientDrivingStressdF                 = initialize( 27,  9 );

            floatMatrix dMicroGradientDrivingStressdFn                = initialize( 27,  9 );

            floatMatrix dMicroGradientDrivingStressdChi               = initialize( 27,  9 );

            floatMatrix dMicroGradientDrivingStressdChin              = initialize( 27,  9 );

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

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

            }

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

            virtual void setDrivingStressesJacobians( const bool isPrevious ) override{

                setDrivingStresses( isPrevious );

                if ( isPrevious ){

                    set_previousdMacroDrivingStressdMacroStress( dMacroDrivingStressdMacroStress );

                    set_previousdMacroDrivingStressdF( dMacroDrivingStressdF );

                    set_previousdMacroDrivingStressdFn( dMacroDrivingStressdFn );

                    set_previousdSymmetricMicroDrivingStressdMicroStress( dMicroDrivingStressdMicroStress );

                    set_previousdSymmetricMicroDrivingStressdF( dMicroDrivingStressdF );

                    set_previousdSymmetricMicroDrivingStressdFn( dMicroDrivingStressdFn );

                    set_previousdHigherOrderDrivingStressdHigherOrderStress( dMicroGradientDrivingStressdHigherOrderStress );

                    set_previousdHigherOrderDrivingStressdF( dMicroGradientDrivingStressdF );

                    set_previousdHigherOrderDrivingStressdFn( dMicroGradientDrivingStressdFn );

                    set_previousdHigherOrderDrivingStressdChi( dMicroGradientDrivingStressdChi );

                    set_previousdHigherOrderDrivingStressdChin( dMicroGradientDrivingStressdChin );

                }
                else{

                    set_dMacroDrivingStressdMacroStress( dMacroDrivingStressdMacroStress );

                    set_dMacroDrivingStressdF( dMacroDrivingStressdF );

                    set_dMacroDrivingStressdFn( dMacroDrivingStressdFn );

                    set_dSymmetricMicroDrivingStressdMicroStress( dMicroDrivingStressdMicroStress );

                    set_dSymmetricMicroDrivingStressdF( dMicroDrivingStressdF );

                    set_dSymmetricMicroDrivingStressdFn( dMicroDrivingStressdFn );

                    set_dHigherOrderDrivingStressdHigherOrderStress( dMicroGradientDrivingStressdHigherOrderStress );

                    set_dHigherOrderDrivingStressdF( dMicroGradientDrivingStressdF );

                    set_dHigherOrderDrivingStressdFn( dMicroGradientDrivingStressdFn );

                    set_dHigherOrderDrivingStressdChi( dMicroGradientDrivingStressdChi );

                    set_dHigherOrderDrivingStressdChin( dMicroGradientDrivingStressdChin );
                
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

            floatVector _local_deltaPK2   = initialize(  9 );

            floatVector _local_deltaSIGMA = initialize(  9 );

            floatVector _local_deltaM     = initialize( 27 );

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

            }

        private:

            using tardigradeHydra::hydraBaseMicromorphic::setResidualClasses;

            virtual void setResidualClasses( ) override{

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = stressMock( this, 45 );

                elasticity.previousPK2Stress += _local_deltaPK2;

                elasticity.previousSIGMA     += _local_deltaSIGMA;

                elasticity.previousM         += _local_deltaM;

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

    // Test of the Jacobians
    
    floatMatrix dUpdatedPlasticStrainLikeISVsdX( 5, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dUpdatedPlasticStrainLikeISVsdPreviousISVs( 5, floatVector( previousStateVariables.size( ), 0 ) );

    floatType eps = 1e-6;

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

        floatVector vp = *Rp.get_updatedPlasticStrainLikeISVs( );

        floatVector vm = *Rm.get_updatedPlasticStrainLikeISVs( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticStrainLikeISVsdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_dUpdatedPlasticStrainLikeISVsdX( 5, floatVector( unknownVector.size( ), 0 ) );

    for ( unsigned int i = 0; i < 5; i++ ){

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dUpdatedPlasticStrainLikeISVsdX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dUpdatedPlasticStrainLikeISVsdStateVariables( ) )[ i ][ j ];

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dUpdatedPlasticStrainLikeISVsdX, assembled_dUpdatedPlasticStrainLikeISVsdX ) );

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

        floatVector vp = *Rp.get_updatedPlasticStrainLikeISVs( );

        floatVector vm = *Rm.get_updatedPlasticStrainLikeISVs( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticStrainLikeISVsdPreviousISVs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_dUpdatedPlasticStrainLikeISVsdPreviousISVs( 5, floatVector( previousStateVariables.size( ), 0 ) );

    for ( unsigned int i = 0; i < 5; i++ ){

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dUpdatedPlasticStrainLikeISVsdPreviousISVs[ i ][ j + configuration_unknown_count ] = ( *R.get_dUpdatedPlasticStrainLikeISVsdPreviousStateVariables( ) )[ i ][ j ];

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dUpdatedPlasticStrainLikeISVsdPreviousISVs, assembled_dUpdatedPlasticStrainLikeISVsdPreviousISVs ) );

}

BOOST_AUTO_TEST_CASE( test_setYield ){
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

            floatVector dMacroCohesiondStateVariables         = initialize( 10 );

            floatVector dMicroCohesiondStateVariables         = initialize( 10 );

            floatMatrix dMicroGradientCohesiondStateVariables = initialize( 3, 10 );

            floatMatrix dMacroDrivingStressdMacroStress = initialize( 9, 9 );

            floatMatrix dMacroDrivingStressdF           = initialize( 9, 9 );

            floatMatrix dMacroDrivingStressdFn          = initialize( 9, 9 );

            floatMatrix dMicroDrivingStressdMicroStress = initialize( 9, 9 );

            floatMatrix dMicroDrivingStressdF           = initialize( 9, 9 );

            floatMatrix dMicroDrivingStressdFn          = initialize( 9, 9 );

            floatMatrix dMicroGradientDrivingStressdHigherOrderStress = initialize( 27, 27 );

            floatMatrix dMicroGradientDrivingStressdF                 = initialize( 27,  9 );

            floatMatrix dMicroGradientDrivingStressdFn                = initialize( 27,  9 );

            floatMatrix dMicroGradientDrivingStressdChi               = initialize( 27,  9 );

            floatMatrix dMicroGradientDrivingStressdChin              = initialize( 27,  9 );

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

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

            }

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

            virtual void setCohesionsJacobians( const bool isPrevious ) override{

                setCohesions( isPrevious );

                if ( isPrevious ){

                    set_previousdMacroCohesiondStateVariables( dMacroCohesiondStateVariables );

                    set_previousdMicroCohesiondStateVariables( dMacroCohesiondStateVariables );

                    set_previousdMicroGradientCohesiondStateVariables( dMicroGradientCohesiondStateVariables );

                }
                else{

                    set_dMacroCohesiondStateVariables( dMacroCohesiondStateVariables );

                    set_dMicroCohesiondStateVariables( dMacroCohesiondStateVariables );

                    set_dMicroGradientCohesiondStateVariables( dMicroGradientCohesiondStateVariables );

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

            virtual void setDrivingStressesJacobians( const bool isPrevious ) override{

                setDrivingStresses( isPrevious );

                if ( isPrevious ){

                    set_previousdMacroDrivingStressdMacroStress( dMacroDrivingStressdMacroStress );

                    set_previousdMacroDrivingStressdF( dMacroDrivingStressdF );

                    set_previousdMacroDrivingStressdFn( dMacroDrivingStressdFn );

                    set_previousdSymmetricMicroDrivingStressdMicroStress( dMicroDrivingStressdMicroStress );

                    set_previousdSymmetricMicroDrivingStressdF( dMicroDrivingStressdF );

                    set_previousdSymmetricMicroDrivingStressdFn( dMicroDrivingStressdFn );

                    set_previousdHigherOrderDrivingStressdHigherOrderStress( dMicroGradientDrivingStressdHigherOrderStress );

                    set_previousdHigherOrderDrivingStressdF( dMicroGradientDrivingStressdF );

                    set_previousdHigherOrderDrivingStressdFn( dMicroGradientDrivingStressdFn );

                    set_previousdHigherOrderDrivingStressdChi( dMicroGradientDrivingStressdChi );

                    set_previousdHigherOrderDrivingStressdChin( dMicroGradientDrivingStressdChin );

                }
                else{

                    set_dMacroDrivingStressdMacroStress( dMacroDrivingStressdMacroStress );

                    set_dMacroDrivingStressdF( dMacroDrivingStressdF );

                    set_dMacroDrivingStressdFn( dMacroDrivingStressdFn );

                    set_dSymmetricMicroDrivingStressdMicroStress( dMicroDrivingStressdMicroStress );

                    set_dSymmetricMicroDrivingStressdF( dMicroDrivingStressdF );

                    set_dSymmetricMicroDrivingStressdFn( dMicroDrivingStressdFn );

                    set_dHigherOrderDrivingStressdHigherOrderStress( dMicroGradientDrivingStressdHigherOrderStress );

                    set_dHigherOrderDrivingStressdF( dMicroGradientDrivingStressdF );

                    set_dHigherOrderDrivingStressdFn( dMicroGradientDrivingStressdFn );

                    set_dHigherOrderDrivingStressdChi( dMicroGradientDrivingStressdChi );

                    set_dHigherOrderDrivingStressdChin( dMicroGradientDrivingStressdChin );
                
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

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    residualMock RJ( &hydra, 55, 1, stateVariableIndices, parameters );

    floatType macroYield, microYield;

    floatVector microGradientYield;

    floatType previousMacroYield, previousMicroYield;

    floatVector previousMicroGradientYield;

    floatVector Fp = floatVector( unknownVector.begin( ) + configuration_unknown_count,
                                  unknownVector.begin( ) + configuration_unknown_count + 9 );

    floatVector eye( 9 );
    tardigradeVectorTools::eye( eye );

    floatVector previousFp = floatVector( previousStateVariables.begin( ),
                                          previousStateVariables.begin( ) + 9 ) + eye;

    floatVector precedingDeformationGradient = tardigradeVectorTools::matrixMultiply( deformationGradient, tardigradeVectorTools::inverse( Fp, 3, 3 ), 3, 3, 3, 3 );

    floatVector previousPrecedingDeformationGradient = tardigradeVectorTools::matrixMultiply( previousDeformationGradient, tardigradeVectorTools::inverse( previousFp, 3, 3 ), 3, 3, 3, 3 );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( R.macroDrivingStress, R.macroCohesion, precedingDeformationGradient,
                                                                                                        R.plasticParameters[ 19 ], R.plasticParameters[ 20 ], macroYield );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( R.microDrivingStress, R.microCohesion, precedingDeformationGradient,
                                                                                                        R.plasticParameters[ 22 ], R.plasticParameters[ 23 ], microYield );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( R.microGradientDrivingStress, R.microGradientCohesion, precedingDeformationGradient,
                                                                                                        R.plasticParameters[ 25 ], R.plasticParameters[ 26 ], microGradientYield );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( R.previousMacroDrivingStress, R.previousMacroCohesion, previousPrecedingDeformationGradient,
                                                                                                        R.plasticParameters[ 19 ], R.plasticParameters[ 20 ], previousMacroYield );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( R.previousMicroDrivingStress, R.previousMicroCohesion, previousPrecedingDeformationGradient,
                                                                                                        R.plasticParameters[ 22 ], R.plasticParameters[ 23 ], previousMicroYield );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( R.previousMicroGradientDrivingStress, R.previousMicroGradientCohesion, previousPrecedingDeformationGradient,
                                                                                                        R.plasticParameters[ 25 ], R.plasticParameters[ 26 ], previousMicroGradientYield );

    RJ.get_dMacroYielddStress( );
    RJ.get_previousdMacroYielddStress( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( macroYield,                 *R.get_macroYield( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( microYield,                 *R.get_microYield( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( microGradientYield,         *R.get_microGradientYield( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousMacroYield,         *R.get_previousMacroYield( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousMicroYield,         *R.get_previousMicroYield( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousMicroGradientYield, *R.get_previousMicroGradientYield( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( macroYield,                 *RJ.get_macroYield( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( microYield,                 *RJ.get_microYield( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( microGradientYield,         *RJ.get_microGradientYield( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousMacroYield,         *RJ.get_previousMacroYield( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousMicroYield,         *RJ.get_previousMicroYield( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousMicroGradientYield, *RJ.get_previousMicroGradientYield( ) ) );

}

BOOST_AUTO_TEST_CASE( test_setYield2 ){
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

            floatVector previousPK2       = {  1,  2,  3,  4,  5,  6,  7,  8,  9 };

            floatVector previousSIGMA     = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

            floatVector previousM         = { 19, 20, 21, 22, 23, 24, 25, 26, 27,
                                              28, 29, 30, 31, 32, 33, 34, 35, 36,
                                              37, 38, 39, 40, 41, 42, 43, 44, 45 };

        protected:

            using tardigradeHydra::residualBase::setPreviousStress;

            virtual void setPreviousStress( ) override{

                setPreviousStress( tardigradeVectorTools::appendVectors( { previousPK2, previousSIGMA, previousM } ) );

            }

    };

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticity::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::residual;

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              2, 0.95338442, 0.74042148,
                                              2, 0.38093104, 0.49241325,
                                              2, 0.82121039, 0.90566759,
                                              2, 0.01166325, 0.05331896,
                                              2, 0.32982199, 0.60161431,
                                              2, 0.58881096, 0.11473813 };

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

            }

    };

    class hydraBaseMicromorphicMock : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

            stressMock elasticity;

            residualMock plasticity;

            std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector _local_deltaPK2    = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            floatVector _local_deltaSIGMA  = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            floatVector _local_deltaM      = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0 };

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

                elasticity.previousPK2   += _local_deltaPK2;

                elasticity.previousSIGMA += _local_deltaSIGMA;

                elasticity.previousM     += _local_deltaM;

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

    // Set the Jacobians

    floatVector dMacroYielddX( unknownVector.size( ), 0 );

    floatVector dMicroYielddX( unknownVector.size( ), 0 );

    floatMatrix dMicroGradientYielddX( 3, floatVector( unknownVector.size( ), 0 ) );

    floatVector dMacroYielddF( deformationGradient.size( ), 0 );

    floatVector dMicroYielddF( deformationGradient.size( ), 0 );

    floatMatrix dMicroGradientYielddF( 3, floatVector( deformationGradient.size( ), 0 ) );

    floatMatrix dMicroGradientYielddChi( 3, floatVector( deformationGradient.size( ), 0 ) );

    floatVector previousdMacroYielddStress( 45, 0 );

    floatVector previousdMicroYielddStress( 45, 0 );

    floatMatrix previousdMicroGradientYielddStress( 3, floatVector( 45, 0 ) );

    floatVector previousdMacroYielddStateVariables( previousStateVariables.size( ), 0 );

    floatVector previousdMicroYielddStateVariables( previousStateVariables.size( ), 0 );

    floatMatrix previousdMicroGradientYielddStateVariables( 3, floatVector( previousStateVariables.size( ), 0 ) );

    floatVector previousdMacroYielddF( deformationGradient.size( ), 0 );

    floatVector previousdMicroYielddF( deformationGradient.size( ), 0 );

    floatMatrix previousdMicroGradientYielddF( 3, floatVector( deformationGradient.size( ), 0 ) );

    floatMatrix previousdMicroGradientYielddChi( 3, floatVector( deformationGradient.size( ), 0 ) );

    floatType eps = 1e-6;

    // Jacobians w.r.t. stress, Fn, Chin, and state variables
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

        floatType sp = *Rp.get_macroYield( );

        floatType sm = *Rm.get_macroYield( );

        dMacroYielddX[ i ] = ( sp - sm ) / ( 2 * delta[ i ] );

        sp = *Rp.get_microYield( );

        sm = *Rm.get_microYield( );

        dMicroYielddX[ i ] = ( sp - sm ) / ( 2 * delta[ i ] );

        floatVector vp = *Rp.get_microGradientYield( );

        floatVector vm = *Rm.get_microGradientYield( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicroGradientYielddX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatVector assemble_dMacroYielddX( unknownVector.size( ), 0 );

    floatVector assemble_dMicroYielddX( unknownVector.size( ), 0 );

    floatMatrix assemble_dMicroGradientYielddX( 3, floatVector( unknownVector.size( ), 0 ) );

    // Stress derivatives
    for ( unsigned int i = 0; i < 9; i++ ){

        assemble_dMacroYielddX[ i ] = ( *R.get_dMacroYielddStress( ) )[ i ];

    }

    for ( unsigned int i = 0; i < 9; i++ ){

        assemble_dMicroYielddX[ i + 9 ] = ( *R.get_dMicroYielddStress( ) )[ i ];

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 27; j++ ){

            assemble_dMicroGradientYielddX[ i ][ j + 18 ] = ( *R.get_dMicroGradientYielddStress( ) )[ i ][ j ];

        }

    }

    // Fn derivatives
    for ( unsigned int i = 0; i < 9; i++ ){

        assemble_dMacroYielddX[ i + configuration_unknown_count ] = ( *R.get_dMacroYielddFn( ) )[ i ];

    }

    for ( unsigned int i = 0; i < 9; i++ ){

        assemble_dMicroYielddX[ i + configuration_unknown_count ] = ( *R.get_dMicroYielddFn( ) )[ i ];

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assemble_dMicroGradientYielddX[ i ][ j + configuration_unknown_count ] = ( *R.get_dMicroGradientYielddFn( ) )[ i ][ j ];

        }

    }

    // Chin derivatives
    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assemble_dMicroGradientYielddX[ i ][ j + configuration_unknown_count + 9 ] = ( *R.get_dMicroGradientYielddChin( ) )[ i ][ j ];

        }

    }

    // State variable derivatives
    for ( unsigned int i = 0; i < 10; i++ ){

        assemble_dMacroYielddX[ i + 2 * configuration_unknown_count ] = ( *R.get_dMacroYielddStateVariables( ) )[ i ];

    }

    for ( unsigned int i = 0; i < 10; i++ ){

        assemble_dMicroYielddX[ i + 2 * configuration_unknown_count ] = ( *R.get_dMicroYielddStateVariables( ) )[ i ];

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 10; j++ ){

            assemble_dMicroGradientYielddX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dMicroGradientYielddStateVariables( ) )[ i ][ j ];

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assemble_dMacroYielddX, dMacroYielddX ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assemble_dMicroYielddX, dMicroYielddX ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assemble_dMicroGradientYielddX, dMicroGradientYielddX ) );

    // Jacobians w.r.t. F
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

        floatType sp = *Rp.get_macroYield( );

        floatType sm = *Rm.get_macroYield( );

        dMacroYielddF[ i ] = ( sp - sm ) / ( 2 * delta[ i ] );

        sp = *Rp.get_microYield( );

        sm = *Rm.get_microYield( );

        dMicroYielddF[ i ] = ( sp - sm ) / ( 2 * delta[ i ] );

        floatVector vp = *Rp.get_microGradientYield( );

        floatVector vm = *Rm.get_microGradientYield( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicroGradientYielddF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMacroYielddF, *R.get_dMacroYielddF( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroYielddF, *R.get_dMicroYielddF( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroGradientYielddF, *R.get_dMicroGradientYielddF( ) ) );

    // Jacobians w.r.t. Chi
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

        floatType sp = *Rp.get_macroYield( );

        floatType sm = *Rm.get_macroYield( );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( sp - sm ) / ( 2 * delta[ i ] ) ) );

        sp = *Rp.get_microYield( );

        sm = *Rm.get_microYield( );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( sp - sm ) / ( 2 * delta[ i ] ) ) );

        floatVector vp = *Rp.get_microGradientYield( );

        floatVector vm = *Rm.get_microGradientYield( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicroGradientYielddChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroGradientYielddChi, *R.get_dMicroGradientYielddChi( ) ) );

    // Jacobians w.r.t. the stress
    floatVector stress = tardigradeVectorTools::appendVectors( { hydra.elasticity.previousPK2, hydra.elasticity.previousSIGMA, hydra.elasticity.previousM } );

    for ( unsigned int i = 0; i < stress.size( ); i++ ){

        floatVector delta( stress.size( ), 0 );

        delta[ i ] = eps * std::fabs( stress[ i ] ) + eps;

        floatVector dpk2(   delta.begin( ) +  0, delta.begin( ) +  9 );

        floatVector dsigma( delta.begin( ) +  9, delta.begin( ) + 18 );

        floatVector dm(     delta.begin( ) + 18, delta.begin( ) + 45 );

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

        hydrap._local_deltaPK2   = dpk2;

        hydrap._local_deltaSIGMA = dsigma;

        hydrap._local_deltaM     = dm;

        hydram._local_deltaPK2   = -dpk2;

        hydram._local_deltaSIGMA = -dsigma;

        hydram._local_deltaM     = -dm;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm( &hydram, 55, 1, stateVariableIndices, parameters );

        floatType sp = *Rp.get_previousMacroYield( );

        floatType sm = *Rm.get_previousMacroYield( );

        previousdMacroYielddStress[ i ] = ( sp - sm ) / ( 2 * delta[ i ] );

        sp = *Rp.get_previousMicroYield( );

        sm = *Rm.get_previousMicroYield( );

        previousdMicroYielddStress[ i ] = ( sp - sm ) / ( 2 * delta[ i ] );

        floatVector vp = *Rp.get_previousMicroGradientYield( );

        floatVector vm = *Rm.get_previousMicroGradientYield( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroGradientYielddStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatVector assemble_previousdMacroYielddStress( 45, 0 );

    floatVector assemble_previousdMicroYielddStress( 45, 0 );

    floatMatrix assemble_previousdMicroGradientYielddStress( 3, floatVector( 45, 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        assemble_previousdMacroYielddStress[ i ] = ( *R.get_previousdMacroYielddStress( ) )[ i ];

    }

    for ( unsigned int i = 0; i < 9; i++ ){

        assemble_previousdMicroYielddStress[ i + 9 ] = ( *R.get_previousdMicroYielddStress( ) )[ i ];

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 27; j++ ){

            assemble_previousdMicroGradientYielddStress[ i ][ j + 18 ] = ( *R.get_previousdMicroGradientYielddStress( ) )[ i ][ j ];

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assemble_previousdMacroYielddStress, previousdMacroYielddStress ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assemble_previousdMicroYielddStress, previousdMicroYielddStress ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assemble_previousdMicroGradientYielddStress, previousdMicroGradientYielddStress ) );

    // Jacobians w.r.t. the previous Fn, Chin, and state variables
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

        floatType sp = *Rp.get_previousMacroYield( );

        floatType sm = *Rm.get_previousMacroYield( );

        previousdMacroYielddStateVariables[ i ] = ( sp - sm ) / ( 2 * delta[ i ] );

        sp = *Rp.get_previousMicroYield( );

        sm = *Rm.get_previousMicroYield( );

        previousdMicroYielddStateVariables[ i ] = ( sp - sm ) / ( 2 * delta[ i ] );

        floatVector vp = *Rp.get_previousMicroGradientYield( );

        floatVector vm = *Rm.get_previousMicroGradientYield( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroGradientYielddStateVariables[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatVector assemble_previousdMacroYielddStateVariables( previousStateVariables.size( ), 0 );

    floatVector assemble_previousdMicroYielddStateVariables( previousStateVariables.size( ), 0 );

    floatMatrix assemble_previousdMicroGradientYielddStateVariables( 3, floatVector( previousStateVariables.size( ), 0 ) );

    // Fn derivatives
    for ( unsigned int i = 0; i < 9; i++ ){

        assemble_previousdMacroYielddStateVariables[ i ] = ( *R.get_previousdMacroYielddFn( ) )[ i ];

    }

    for ( unsigned int i = 0; i < 9; i++ ){

        assemble_previousdMicroYielddStateVariables[ i ] = ( *R.get_previousdMicroYielddFn( ) )[ i ];

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assemble_previousdMicroGradientYielddStateVariables[ i ][ j ] = ( *R.get_previousdMicroGradientYielddFn( ) )[ i ][ j ];

        }

    }

    // Chin derivatives
    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assemble_previousdMicroGradientYielddStateVariables[ i ][ j + 9 ] = ( *R.get_previousdMicroGradientYielddChin( ) )[ i ][ j ];

        }

    }

    // State variable derivatives
    for ( unsigned int i = 0; i < 10; i++ ){

        assemble_previousdMacroYielddStateVariables[ i + configuration_unknown_count ] = ( *R.get_previousdMacroYielddStateVariables( ) )[ i ];

    }

    for ( unsigned int i = 0; i < 10; i++ ){

        assemble_previousdMicroYielddStateVariables[ i + configuration_unknown_count ] = ( *R.get_previousdMicroYielddStateVariables( ) )[ i ];

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 10; j++ ){

            assemble_previousdMicroGradientYielddStateVariables[ i ][ j + configuration_unknown_count ] = ( *R.get_previousdMicroGradientYielddStateVariables( ) )[ i ][ j ];

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assemble_previousdMacroYielddStateVariables, previousdMacroYielddStateVariables ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assemble_previousdMicroYielddStateVariables, previousdMicroYielddStateVariables ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( assemble_previousdMicroGradientYielddStateVariables, previousdMicroGradientYielddStateVariables ) );

    // Jacobians w.r.t. previous F
    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

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

        floatType sp = *Rp.get_previousMacroYield( );

        floatType sm = *Rm.get_previousMacroYield( );

        previousdMacroYielddF[ i ] = ( sp - sm ) / ( 2 * delta[ i ] );

        sp = *Rp.get_previousMicroYield( );

        sm = *Rm.get_previousMicroYield( );

        previousdMicroYielddF[ i ] = ( sp - sm ) / ( 2 * delta[ i ] );

        floatVector vp = *Rp.get_previousMicroGradientYield( );

        floatVector vm = *Rm.get_previousMicroGradientYield( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroGradientYielddF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMacroYielddF, *R.get_previousdMacroYielddF( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroYielddF, *R.get_previousdMicroYielddF( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroGradientYielddF, *R.get_previousdMicroGradientYielddF( ) ) );

    // Jacobians w.r.t. Chi
    for ( unsigned int i = 0; i < microDeformation.size( ); i++ ){

        floatVector delta( microDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

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

        floatType sp = *Rp.get_previousMacroYield( );

        floatType sm = *Rm.get_previousMacroYield( );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( sp - sm ) / ( 2 * delta[ i ] ) ) );

        sp = *Rp.get_microYield( );

        sm = *Rm.get_microYield( );

        BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( sp - sm ) / ( 2 * delta[ i ] ) ) );

        floatVector vp = *Rp.get_previousMicroGradientYield( );

        floatVector vm = *Rm.get_previousMicroGradientYield( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroGradientYielddChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroGradientYielddChi, *R.get_previousdMicroGradientYielddChi( ) ) );

}

BOOST_AUTO_TEST_CASE( test_setPrecedingDeformationGradient ){
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

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              2, 0.95338442, 0.74042148,
                                              2, 0.38093104, 0.49241325,
                                              2, 0.82121039, 0.90566759,
                                              2, 0.01166325, 0.05331896,
                                              2, 0.32982199, 0.60161431,
                                              2, 0.58881096, 0.11473813 };

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

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

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    residualMock RJ( &hydra, 55, 1, stateVariableIndices, parameters );

    floatVector Fp( unknownVector.begin( ) + configuration_unknown_count, unknownVector.begin( ) + configuration_unknown_count + 9 );

    floatVector eye( deformationGradient.size( ) );
    tardigradeVectorTools::eye( eye );

    floatVector previousFp( previousStateVariables.begin( ), previousStateVariables.begin( ) + 9 );
    previousFp += eye;

    floatVector answerFe         = tardigradeVectorTools::matrixMultiply( deformationGradient,         tardigradeVectorTools::inverse( Fp, 3, 3 ),         3, 3, 3, 3 );

    floatVector answerPreviousFe = tardigradeVectorTools::matrixMultiply( previousDeformationGradient, tardigradeVectorTools::inverse( previousFp, 3, 3 ), 3, 3, 3, 3 );

    RJ.get_dPrecedingDeformationGradientdF( );
    RJ.get_previousdPrecedingDeformationGradientdF( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerFe,         *R.get_precedingDeformationGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousFe, *R.get_previousPrecedingDeformationGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerFe,         *RJ.get_precedingDeformationGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousFe, *RJ.get_previousPrecedingDeformationGradient( ) ) );

    // Test the jacobians

    floatMatrix dPrecedingFdF(  9, floatVector( 9, 0 ) );

    floatMatrix dPrecedingFdFn( 9, floatVector( 9, 0 ) );

    floatMatrix previousdPrecedingFdF(  9, floatVector( 9, 0 ) );

    floatMatrix previousdPrecedingFdFn( 9, floatVector( 9, 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_precedingDeformationGradient( );

        floatVector vm = *Rm.get_precedingDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dPrecedingFdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPrecedingFdF, *R.get_dPrecedingDeformationGradientdF( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + configuration_unknown_count ] = eps * std::fabs( unknownVector[ i + configuration_unknown_count ] ) + eps;

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_precedingDeformationGradient( );

        floatVector vm = *Rm.get_precedingDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dPrecedingFdFn[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + configuration_unknown_count ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPrecedingFdFn, *R.get_dPrecedingDeformationGradientdFn( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousPrecedingDeformationGradient( );

        floatVector vm = *Rm.get_previousPrecedingDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdPrecedingFdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdPrecedingFdF, *R.get_previousdPrecedingDeformationGradientdF( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousPrecedingDeformationGradient( );

        floatVector vm = *Rm.get_previousPrecedingDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdPrecedingFdFn[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdPrecedingFdFn, *R.get_previousdPrecedingDeformationGradientdFn( ) ) );

}

BOOST_AUTO_TEST_CASE( test_setPrecedingMicroDeformation ){
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

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              2, 0.95338442, 0.74042148,
                                              2, 0.38093104, 0.49241325,
                                              2, 0.82121039, 0.90566759,
                                              2, 0.01166325, 0.05331896,
                                              2, 0.32982199, 0.60161431,
                                              2, 0.58881096, 0.11473813 };

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

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

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    residualMock RJ( &hydra, 55, 1, stateVariableIndices, parameters );

    floatVector chip( unknownVector.begin( ) + configuration_unknown_count + 9, unknownVector.begin( ) + configuration_unknown_count + 18 );

    floatVector eye( deformationGradient.size( ) );
    tardigradeVectorTools::eye( eye );

    floatVector previousChip( previousStateVariables.begin( ) + 9, previousStateVariables.begin( ) + 18 );
    previousChip += eye;

    floatVector answerChie         = tardigradeVectorTools::matrixMultiply( microDeformation,         tardigradeVectorTools::inverse( chip, 3, 3 ),         3, 3, 3, 3 );

    floatVector answerPreviousChie = tardigradeVectorTools::matrixMultiply( previousMicroDeformation, tardigradeVectorTools::inverse( previousChip, 3, 3 ), 3, 3, 3, 3 );

    RJ.get_dPrecedingMicroDeformationdChi( );
    RJ.get_previousdPrecedingMicroDeformationdChi( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerChie,         *R.get_precedingMicroDeformation( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousChie, *R.get_previousPrecedingMicroDeformation( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerChie,         *RJ.get_precedingMicroDeformation( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousChie, *RJ.get_previousPrecedingMicroDeformation( ) ) );

    // Test the jacobians

    floatMatrix dPrecedingChidChi(  9, floatVector( 9, 0 ) );

    floatMatrix dPrecedingChidChin( 9, floatVector( 9, 0 ) );

    floatMatrix previousdPrecedingChidChi(  9, floatVector( 9, 0 ) );

    floatMatrix previousdPrecedingChidChin( 9, floatVector( 9, 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_precedingMicroDeformation( );

        floatVector vm = *Rm.get_precedingMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dPrecedingChidChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPrecedingChidChi, *R.get_dPrecedingMicroDeformationdChi( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + configuration_unknown_count + 9 ] = eps * std::fabs( unknownVector[ i + configuration_unknown_count + 9 ] ) + eps;

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_precedingMicroDeformation( );

        floatVector vm = *Rm.get_precedingMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dPrecedingChidChin[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + configuration_unknown_count + 9 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPrecedingChidChin, *R.get_dPrecedingMicroDeformationdChin( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousPrecedingMicroDeformation( );

        floatVector vm = *Rm.get_previousPrecedingMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdPrecedingChidChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdPrecedingChidChi, *R.get_previousdPrecedingMicroDeformationdChi( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( previousStateVariables[ i + 9 ] ) + eps;

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousPrecedingMicroDeformation( );

        floatVector vm = *Rm.get_previousPrecedingMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdPrecedingChidChin[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdPrecedingChidChin, *R.get_previousdPrecedingMicroDeformationdChin( ) ) );

}

BOOST_AUTO_TEST_CASE( test_setPrecedingGradientMicroDeformation ){
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

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              2, 0.95338442, 0.74042148,
                                              2, 0.38093104, 0.49241325,
                                              2, 0.82121039, 0.90566759,
                                              2, 0.01166325, 0.05331896,
                                              2, 0.32982199, 0.60161431,
                                              2, 0.58881096, 0.11473813 };

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

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

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    residualMock RJ( &hydra, 55, 1, stateVariableIndices, parameters );

    floatVector answerGradChie         = ( *hydra.get_gradientMicroConfigurations( ) )[ 0 ];

    floatVector answerPreviousGradChie = ( *hydra.get_previousGradientMicroConfigurations( ) )[ 0 ];

    RJ.get_dPrecedingMicroDeformationdChi( );
    RJ.get_previousdPrecedingMicroDeformationdChi( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerGradChie,         *R.get_precedingGradientMicroDeformation( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousGradChie, *R.get_previousPrecedingGradientMicroDeformation( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerGradChie,         *RJ.get_precedingGradientMicroDeformation( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousGradChie, *RJ.get_previousPrecedingGradientMicroDeformation( ) ) );

    // Test the jacobians

    floatMatrix dPrecedingGradChidFn(       27, floatVector(  9, 0 ) );

    floatMatrix dPrecedingGradChidChi(      27, floatVector(  9, 0 ) );

    floatMatrix dPrecedingGradChidChin(     27, floatVector(  9, 0 ) );

    floatMatrix dPrecedingGradChidGradChi(  27, floatVector( 27, 0 ) );

    floatMatrix dPrecedingGradChidGradChin( 27, floatVector( 27, 0 ) );

    floatMatrix previousdPrecedingGradChidFn(       27, floatVector(  9, 0 ) );

    floatMatrix previousdPrecedingGradChidChi(      27, floatVector(  9, 0 ) );

    floatMatrix previousdPrecedingGradChidChin(     27, floatVector(  9, 0 ) );

    floatMatrix previousdPrecedingGradChidGradChi(  27, floatVector( 27, 0 ) );

    floatMatrix previousdPrecedingGradChidGradChin( 27, floatVector( 27, 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + configuration_unknown_count ] = eps * std::fabs( unknownVector[ i + configuration_unknown_count ] ) + eps;

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_precedingGradientMicroDeformation( );

        floatVector vm = *Rm.get_precedingGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dPrecedingGradChidFn[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + configuration_unknown_count ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPrecedingGradChidFn, *R.get_dPrecedingGradientMicroDeformationdFn( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_precedingGradientMicroDeformation( );

        floatVector vm = *Rm.get_precedingGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dPrecedingGradChidChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPrecedingGradChidChi, *R.get_dPrecedingGradientMicroDeformationdChi( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + configuration_unknown_count + 9 ] = eps * std::fabs( unknownVector[ i + configuration_unknown_count + 9 ] ) + eps;

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_precedingGradientMicroDeformation( );

        floatVector vm = *Rm.get_precedingGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dPrecedingGradChidChin[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + configuration_unknown_count + 9 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPrecedingGradChidChin, *R.get_dPrecedingGradientMicroDeformationdChin( ) ) );

    for ( unsigned int i = 0; i < 27; i++ ){

        floatVector delta( 27, 0 );

        delta[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation + delta, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation - delta, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_precedingGradientMicroDeformation( );

        floatVector vm = *Rm.get_precedingGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dPrecedingGradChidGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPrecedingGradChidGradChi, *R.get_dPrecedingGradientMicroDeformationdGradChi( ) ) );

    for ( unsigned int i = 0; i < 27; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + configuration_unknown_count + 18 ] = eps * std::fabs( unknownVector[ i + configuration_unknown_count + 18 ] ) + eps;

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_precedingGradientMicroDeformation( );

        floatVector vm = *Rm.get_precedingGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dPrecedingGradChidGradChin[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + configuration_unknown_count + 18 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPrecedingGradChidGradChin, *R.get_dPrecedingGradientMicroDeformationdGradChin( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousPrecedingGradientMicroDeformation( );

        floatVector vm = *Rm.get_previousPrecedingGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdPrecedingGradChidFn[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdPrecedingGradChidFn, *R.get_previousdPrecedingGradientMicroDeformationdFn( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousPrecedingGradientMicroDeformation( );

        floatVector vm = *Rm.get_previousPrecedingGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdPrecedingGradChidChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdPrecedingGradChidChi, *R.get_previousdPrecedingGradientMicroDeformationdChi( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( previousStateVariables[ i + 9 ] ) + eps;

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousPrecedingGradientMicroDeformation( );

        floatVector vm = *Rm.get_previousPrecedingGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdPrecedingGradChidChin[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdPrecedingGradChidChin, *R.get_previousdPrecedingGradientMicroDeformationdChin( ) ) );

    for ( unsigned int i = 0; i < 27; i++ ){

        floatVector delta( 27, 0 );

        delta[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation + delta,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation - delta,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousPrecedingGradientMicroDeformation( );

        floatVector vm = *Rm.get_previousPrecedingGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdPrecedingGradChidGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdPrecedingGradChidGradChi, *R.get_previousdPrecedingGradientMicroDeformationdGradChi( ) ) );

    for ( unsigned int i = 0; i < 27; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + 18 ] = eps * std::fabs( previousStateVariables[ i + 18 ] ) + eps;

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousPrecedingGradientMicroDeformation( );

        floatVector vm = *Rm.get_previousPrecedingGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdPrecedingGradChidGradChin[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + 18 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdPrecedingGradChidGradChin, *R.get_previousdPrecedingGradientMicroDeformationdGradChin( ) ) );

}

BOOST_AUTO_TEST_CASE( testComputePlasticMacroVelocityGradient ){
    /*!
     * Test the computation of the plastic macro velocity gradient.
     *
     */

    variableType macroGamma = 0.08166694603978908;
    variableType microGamma = 0.8652174130049269;

    variableVector Ce = { 15.81870565,  0.8392615 ,  2.09805203,
                           0.8392615 ,  2.68322729, -6.62003948,
                           2.09805203, -6.62003948, 24.82920808 };

    variableVector inverseCe = tardigradeVectorTools::inverse( Ce, 3, 3 );

    variableVector macroFlowDirection = { 0.78884638, 0.19639211, 0.15523073,
                                          0.47307595, 0.28241451, 0.66404732,
                                          0.1634089 , 0.92452471, 0.77390742 };

    variableVector microFlowDirection = { 0.86300151, 0.95736394, 0.61329255,
                                          0.27780339, 0.26054793, 0.33313753,
                                          0.34289169, 0.57971261, 0.51536929 };

    variableVector answerMacroLp = { -0.05489573, -0.01980382, -0.06060589,
                                      1.1610081 ,  0.4002548 ,  0.86866858,
                                      0.33607202,  0.12218348,  0.25723268 };

    variableVector resultMacroLp;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseCe,
                                                                                               macroFlowDirection, microFlowDirection,
                                                                                               resultMacroLp );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMacroLp, resultMacroLp ) );

    //Tests of the Jacobians
    variableVector resultMacroLpJ;
    variableVector dMacroLdMacroGammaJ, dMacroLdMicroGammaJ;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseCe,
                                                                                               macroFlowDirection, microFlowDirection,
                                                                                               resultMacroLpJ, dMacroLdMacroGammaJ,
                                                                                               dMacroLdMicroGammaJ );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMacroLp, resultMacroLpJ ) );

    variableVector resultMacroLpJ2;
    variableVector dMacroLdMacroGammaJ2, dMacroLdMicroGammaJ2;
    variableMatrix dMacroLdElasticRCG, dMacroLdMacroFlowDirection, dMacroLdMicroFlowDirection;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseCe,
                                                                       macroFlowDirection, microFlowDirection,
                                                                       resultMacroLpJ2, dMacroLdMacroGammaJ2,
                                                                       dMacroLdMicroGammaJ2, dMacroLdElasticRCG,
                                                                       dMacroLdMacroFlowDirection,
                                                                       dMacroLdMicroFlowDirection );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMacroLp, resultMacroLpJ2 ) );

    //Tests Jacobians w.r.t. macroGamma
    constantType eps = 1e-6;
    constantType scalarDelta = eps * fabs( macroGamma) + eps;

    variableVector resultMacroLpP, resultMacroLpM;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma + scalarDelta, microGamma,
                                                                       inverseCe, macroFlowDirection, microFlowDirection,
                                                                       resultMacroLpP );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma - scalarDelta, microGamma,
                                                                       inverseCe, macroFlowDirection, microFlowDirection,
                                                                       resultMacroLpM );

    variableVector gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * scalarDelta );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol, dMacroLdMacroGammaJ ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol, dMacroLdMacroGammaJ2 ) );

    //Test Jacobians w.r.t. microGamma
    scalarDelta = eps * fabs( microGamma) + eps;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma + scalarDelta,
                                                                       inverseCe, macroFlowDirection, microFlowDirection,
                                                                       resultMacroLpP );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma - scalarDelta,
                                                                       inverseCe, macroFlowDirection, microFlowDirection,
                                                                       resultMacroLpM );

    gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * scalarDelta );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol, dMacroLdMicroGammaJ ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol, dMacroLdMicroGammaJ2 ) );

    //Test Jacobians w.r.t. the right Cauchy-Green deformation tensor
    for ( unsigned int i = 0; i < Ce.size(); i++ ){

        constantVector delta( Ce.size(), 0 );

        delta[i] = eps * fabs( Ce[i] ) + eps;

        variableVector inverseCeTemp = tardigradeVectorTools::inverse( Ce + delta, 3, 3 );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                           inverseCeTemp, macroFlowDirection,
                                                                           microFlowDirection, resultMacroLpP );

        inverseCeTemp = tardigradeVectorTools::inverse( Ce - delta, 3, 3 );
    
        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                           inverseCeTemp, macroFlowDirection,
                                                                           microFlowDirection, resultMacroLpM );
    
        gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dMacroLdElasticRCG[ j ][ i ] ) );

        }

    }

    //Test Jacobians w.r.t. the macro flow direction
    for ( unsigned int i = 0; i < macroFlowDirection.size(); i++ ){

        constantVector delta( macroFlowDirection.size(), 0 );

        delta[i] = eps * fabs( macroFlowDirection[i] ) + eps;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                           inverseCe, macroFlowDirection + delta,
                                                                           microFlowDirection, resultMacroLpP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                           inverseCe, macroFlowDirection - delta,
                                                                           microFlowDirection, resultMacroLpM );

        gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dMacroLdMacroFlowDirection[ j ][ i ] ) );

        }

    }

    //Test Jacobians w.r.t. the micro flow direction
    for ( unsigned int i = 0; i < microFlowDirection.size(); i++ ){

        constantVector delta( microFlowDirection.size(), 0 );

        delta[i] = eps * fabs( microFlowDirection[i] ) + eps;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                           inverseCe, macroFlowDirection,
                                                                           microFlowDirection + delta, resultMacroLpP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma,
                                                                           inverseCe, macroFlowDirection,
                                                                           microFlowDirection - delta, resultMacroLpM );

        gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dMacroLdMicroFlowDirection[ j ][ i ] ) );

        }

    }

}

BOOST_AUTO_TEST_CASE( testComputePlasticMicroVelocityGradient ){
    /*!
     * Test the computation of the plastic micro velocity gradient.
     *
     */

    variableType microGamma = 0.8652174130049269;

    variableVector Ce = { 17.06524293,   5.38540351, -20.25732698,
                           5.38540351,   7.29152741,  -2.58538828,
                         -20.25732698,  -2.58538828,  33.09421588 };

    variableVector Psie = { 8.80605252,   2.3131131 , -19.28700431,
                            3.95982925,  -1.71986455,  -5.62211322,
                           -0.28252834,  11.47370888,   5.77705312 };

    variableVector invPsie = tardigradeVectorTools::inverse( Psie, 3, 3 );

    variableVector microFlowDirection = { 0.86300151, 0.95736394, 0.61329255,
                                          0.27780339, 0.26054793, 0.33313753,
                                          0.34289169, 0.57971261, 0.51536929 };

// NOTE: Replaced this with what I think is the right value
//    variableVector answerMicroLp = { -1.14307645,  1.41916637,  2.9980591 ,
//                                      0.03471568, -0.06540477, -0.10626039,
//                                     -0.42684692,  0.51754693,  1.10823355 };

    variableVector answerMicroLp = { 0.84711259,  0.41190088, -0.89488435,
                                     0.00935583,  0.00110041,  0.00441756,
                                     0.31536975,  0.15755706, -0.32577516 };

    variableVector resultMicroLp;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                                               microFlowDirection, resultMicroLp );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMicroLp, resultMicroLp ) );

    //Test the Jacobians
    variableVector resultMicroLpJ;
    variableVector dMicroLpdMicroGamma;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                                               microFlowDirection, resultMicroLpJ,
                                                                                               dMicroLpdMicroGamma );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMicroLp, resultMicroLpJ ) );

    variableVector resultMicroLpJ2;
    variableVector dMicroLpdMicroGamma2;
    variableMatrix dMicroLpdMicroRCG, dMicroLpdPsie, dMicroLpdMicroFlowDirection;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                                               microFlowDirection, resultMicroLpJ2,
                                                                                               dMicroLpdMicroGamma2,
                                                                                               dMicroLpdMicroRCG, dMicroLpdPsie,
                                                                                               dMicroLpdMicroFlowDirection );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMicroLp, resultMicroLpJ ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMicroLp, resultMicroLpJ2 ) );

    constantType eps = 1e-6;
    constantType scalarDelta = eps * fabs( microGamma ) + eps;

    variableVector resultMicroLpP, resultMicroLpM;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma + scalarDelta, Ce, Psie, invPsie,
                                                                                               microFlowDirection, resultMicroLpP );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma - scalarDelta, Ce, Psie, invPsie,
                                                                                               microFlowDirection, resultMicroLpM );

    variableVector gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * scalarDelta );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol, dMicroLpdMicroGamma ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol, dMicroLpdMicroGamma2 ) );

    //Test Jacobian w.r.t. the elastic micro right Cauchy-Green deformation tensor
    for ( unsigned int i = 0; i < Ce.size(); i++ ){

        constantVector delta( Ce.size(), 0 );

        delta[i] = eps * fabs( Ce[i] ) + eps;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce + delta, Psie, invPsie,
                                                                                                   microFlowDirection, resultMicroLpP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce - delta, Psie, invPsie,
                                                                                                   microFlowDirection, resultMicroLpM );

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dMicroLpdMicroRCG[ j ][ i ] ) );

        }

    }

    //Test Jacobian w.r.t. the elastic micro deformation measure Psi
    for ( unsigned int i = 0; i < Psie.size(); i++ ){

        constantVector delta( Psie.size(), 0 );

        delta[i] = eps * fabs( Psie[i] ) + eps;

        variableVector PsieTemp = Psie + delta;

        variableVector invPsieTemp = tardigradeVectorTools::inverse( PsieTemp, 3, 3 );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, PsieTemp, invPsieTemp,
                                                                                                   microFlowDirection, resultMicroLpP );

        PsieTemp = Psie - delta;

        invPsieTemp = tardigradeVectorTools::inverse( PsieTemp, 3, 3 );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, PsieTemp, invPsieTemp,
                                                                                                   microFlowDirection, resultMicroLpM );

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dMicroLpdPsie[ j ][ i ] ) );

        }

    }

    for ( unsigned int i = 0; i < microFlowDirection.size(); i++ ){

        constantVector delta( microFlowDirection.size(), 0 );

        delta[i] = eps * fabs( microFlowDirection[i] ) + eps;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                                                   microFlowDirection + delta, resultMicroLpP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                                                   microFlowDirection - delta, resultMicroLpM );

        gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[j], dMicroLpdMicroFlowDirection[ j ][ i ] ) );

        }

    }

}

BOOST_AUTO_TEST_CASE( testComputePlasticMicroGradientVelocityGradient ){
    /*!
     * Test the computation of the plastic micro gradient velocity gradient.
     *
     */

    variableVector microGradientGamma = { 0.17245016, 0.92420274, 0.28114459 };

    variableVector Psie = { 8.80605252,   2.3131131 , -19.28700431,
                            3.95982925,  -1.71986455,  -5.62211322,
                           -0.28252834,  11.47370888,   5.77705312 };

    variableVector invPsie = tardigradeVectorTools::inverse( Psie, 3, 3 );

    variableVector elasticGamma = { 4.80644001,  0.07861661,  1.64980155,  1.23072776, -0.5292324 ,
                                   -3.06821432,  6.87892124,  3.03299854,  0.04146446, -1.08196034,
                                    1.02647393, -2.6741583 , -0.07702067,  1.53487528, -1.46782133,
                                   -2.79814493, -3.08707902,  0.29650483,  7.95112472, -0.0823429 ,
                                    9.86433536,  0.55102384, -3.97123001,  1.26600849, 14.19808301,
                                    8.33368016,  0.57102355 };

    variableVector microGradientFlowDirection = { 0.38320117, 0.00147635, 0.22526135, 0.24857347, 0.44904944,
                                                  0.39175461, 0.94088825, 0.04088633, 0.95042374, 0.44676197,
                                                  0.33100061, 0.79806506, 0.05883935, 0.20924962, 0.83681153,
                                                  0.12116776, 0.39737069, 0.07417313, 0.5859491 , 0.28899583,
                                                  0.91967175, 0.413024  , 0.97723212, 0.81694258, 0.92037483,
                                                  0.84857389, 0.74623422, 0.65442987, 0.15706966, 0.03580793,
                                                  0.98977654, 0.6414159 , 0.03345668, 0.73436727, 0.25417675,
                                                  0.594925  , 0.4871345 , 0.27395216, 0.23644903, 0.42902409,
                                                  0.24760169, 0.16207352, 0.68475097, 0.86214768, 0.9734798 ,
                                                  0.86141159, 0.98250926, 0.25056881, 0.8315578 , 0.95970017,
                                                  0.62180382, 0.52207192, 0.66811873, 0.06083854, 0.59855098,
                                                  0.41784728, 0.41193658, 0.3161969 , 0.75697096, 0.2172361 ,
                                                  0.5170385 , 0.52482239, 0.55849978, 0.60039656, 0.38358062,
                                                  0.66214191, 0.22829067, 0.10781315, 0.40531347, 0.25340843,
                                                  0.89016033, 0.85477638, 0.43630125, 0.35699992, 0.3784267 ,
                                                  0.12262464, 0.38383612, 0.12695384, 0.74207569, 0.58531619,
                                                  0.08294492 };

    variableVector microLp = { -85.67983387, -16.91839826, 127.3318347 ,
                                 0.65035144,   0.1459189 ,  -0.71988301,
                               -36.05794838,  -7.86041652,  52.33737079 };

    variableVector answerMicroGradLp = { -83.50670143,  -11.14831106,  -39.09529065,  -16.58901227,
                                          -8.71869628,   19.76266969,   64.74989223,    3.55092183,
                                          59.30524632,  126.06867513,   26.68572242,   83.03661657,
                                          25.93560443,    6.13134257,   16.20263835, -185.220318  ,
                                         -38.17442863, -123.4797927 ,  -56.62958098,   -7.64862742,
                                         -16.00279814,  -11.26027965,   -4.19458173,    8.34519958,
                                          58.1455205 ,    5.42232797,   23.44933638 };

    variableVector answerSkewTerm = { -84.06988022,  -11.39160049,  -39.5285613 ,  -17.07306624,
                                       -8.94361366,   19.17011533,   64.59515902,    3.2768447 ,
                                       58.94068694,  126.07596462,   26.65609256,   83.02552879,
                                       25.99081849,    6.07289457,   16.21498441, -185.26655691,
                                      -38.22947657, -123.43856309,  -56.84233212,   -7.7271727 ,
                                      -16.14907393,  -11.46103628,   -4.28260582,    8.13100026,
                                       58.07906172,    5.31870592,   23.31357689 };

    variableVector resultMicroGradLp;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                                       microLp, resultMicroGradLp );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLp ) );

    variableVector resultMicroGradLp1;
    variableVector resultSkewTerm;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                                       microLp, resultMicroGradLp1,
                                                                                                       resultSkewTerm );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLp1 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerSkewTerm, resultSkewTerm ) );

    //Test the Jacobians
    variableVector resultMicroGradLpJ;
    variableMatrix dPlasticMicroGradientLdMicroGradientGamma;
    variableMatrix dPlasticMicroGradientLdPlasticMicroL;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                                       microLp, resultMicroGradLpJ,
                                                                                                       dPlasticMicroGradientLdMicroGradientGamma,
                                                                                                       dPlasticMicroGradientLdPlasticMicroL );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLpJ ) );

    variableVector resultMicroGradLpJ2, resultSkewTerm2;
    variableMatrix dPlasticMicroGradientLdMicroGradientGamma2;
    variableMatrix dPlasticMicroGradientLdPlasticMicroL2;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                                       microLp, resultMicroGradLpJ2, resultSkewTerm2,
                                                                                                       dPlasticMicroGradientLdMicroGradientGamma2,
                                                                                                       dPlasticMicroGradientLdPlasticMicroL2 );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLpJ2 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerSkewTerm, resultSkewTerm2 ) );

    variableVector resultMicroGradLpJ3;
    variableMatrix dPlasticMicroGradientLdMicroGradientGamma3;
    variableMatrix dPlasticMicroGradientLdPlasticMicroL3;
    variableMatrix dPlasticMicroGradientLdElasticPsi;
    variableMatrix dPlasticMicroGradientLdElasticGamma;
    variableMatrix dPlasticMicroGradientLdMicroGradientFlowDirection;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                                       microLp, resultMicroGradLpJ3,
                                                                                                       dPlasticMicroGradientLdMicroGradientGamma3,
                                                                                                       dPlasticMicroGradientLdPlasticMicroL3,
                                                                                                       dPlasticMicroGradientLdElasticPsi,
                                                                                                       dPlasticMicroGradientLdElasticGamma,
                                                                                                       dPlasticMicroGradientLdMicroGradientFlowDirection );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMicroGradLp, resultMicroGradLpJ3 ) );

    //Test computation of Jacobians w.r.t. microGradientGamma
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < microGradientGamma.size(); i++ ){
        constantVector delta( microGradientGamma.size(), 0 );
        delta[i] = eps * fabs( microGradientGamma[i] ) + eps;

        variableVector resultMicroGradLpP, resultMicroGradLpM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma + delta, Psie, invPsie,
                                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                                           microLp, resultMicroGradLpP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma - delta, Psie, invPsie,
                                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                                           microLp, resultMicroGradLpM );

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroGradientGamma[ j ][ i ] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroGradientGamma2[ j ][ i ] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroGradientGamma3[ j ][ i ] ) );
        }
    }

    //Test computation of Jacobians w.r.t. the plastic micro velocity gradient
    for ( unsigned int i = 0; i < microLp.size(); i++ ){
        constantVector delta( microLp.size(), 0 );
        delta[i] = eps * fabs( microLp[i] ) + eps;

        variableVector resultMicroGradLpP, resultMicroGradLpM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                                           microLp + delta, resultMicroGradLpP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                                           microLp - delta, resultMicroGradLpM );

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdPlasticMicroL[ j ][ i ] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdPlasticMicroL2[ j ][ i ] ) );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdPlasticMicroL3[ j ][ i ] ) );
        }
    }

    //Test computation of Jacobian w.r.t. the micro deformation measure Psi
    for ( unsigned int i = 0; i < Psie.size(); i++ ){
        constantVector delta( Psie.size(), 0 );
        delta[i] = eps * fabs( Psie[i] ) + eps;

        variableVector resultMicroGradLpP, resultMicroGradLpM;

        variableVector PsieTemp = Psie + delta;
        variableVector invPsieTemp = tardigradeVectorTools::inverse( PsieTemp, 3, 3 );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, PsieTemp, invPsieTemp,
                                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                                           microLp, resultMicroGradLpP );

        PsieTemp = Psie - delta;
        invPsieTemp = tardigradeVectorTools::inverse( PsieTemp, 3, 3 );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, PsieTemp, invPsieTemp,
                                                                                                           elasticGamma, microGradientFlowDirection,
                                                                                                           microLp, resultMicroGradLpM );

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdElasticPsi[ j ][ i ] ) );
        }
    }

    //Test computation of Jacobian w.r.t. the elastic higher order deformation metric Gamma
    for ( unsigned int i = 0; i < elasticGamma.size(); i++ ){
        constantVector delta( elasticGamma.size(), 0 );
        delta[i] = eps * fabs( elasticGamma[i] ) + eps;

        variableVector resultMicroGradLpP, resultMicroGradLpM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                           elasticGamma + delta,
                                                                                                           microGradientFlowDirection,
                                                                                                           microLp, resultMicroGradLpP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                           elasticGamma - delta,
                                                                                                           microGradientFlowDirection,
                                                                                                           microLp, resultMicroGradLpM );

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdElasticGamma[ j ][ i ] ) );
        }
    }

    //Test computation of Jacobian w.r.t. the elastic higher order deformation metric Gamma
    for ( unsigned int i = 0; i < microGradientFlowDirection.size(); i++ ){
        constantVector delta( microGradientFlowDirection.size(), 0 );
        delta[ i ] = eps * fabs( microGradientFlowDirection[ i ] ) + eps;

        variableVector resultMicroGradLpP, resultMicroGradLpM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                           elasticGamma,
                                                                                                           microGradientFlowDirection + delta,
                                                                                                           microLp, resultMicroGradLpP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                           elasticGamma,
                                                                                                           microGradientFlowDirection - delta,
                                                                                                           microLp, resultMicroGradLpM );

        variableVector gradCol = ( resultMicroGradLpP - resultMicroGradLpM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradCol[ j ], dPlasticMicroGradientLdMicroGradientFlowDirection[ j ][ i ] ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( test_setPlasticVelocityGradients ){
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

            floatVector precedingDeformationGradient = { 1, .2, .3, .4, 1.1, 0.3, 0.1, 0.2, 1.2 };

            floatVector previousPrecedingDeformationGradient = { 1, -0.2, -0.3, -0.4, 0.5, -0.6, -0.7, -0.8, 0.1 };

            floatMatrix dPrecedingDeformationGradientdF = initialize( 9, 9 );

            floatMatrix dPrecedingDeformationGradientdFn = initialize( 9, 9 );

            floatMatrix previousdPrecedingDeformationGradientdF = initialize( 9, 9 );

            floatMatrix previousdPrecedingDeformationGradientdFn = initialize( 9, 9 );

            floatVector precedingMicroDeformation = { 0.8, -0.1, 0.2, 0.3, 1.1, 0.4, -0.2, 0.5, 1.2 };

            floatVector previousPrecedingMicroDeformation = { 0.7, 0.11, -0.2, 0.2, 0.9, 0.2, 0.2, -0.1, 1.1 };

            floatVector precedingGradientMicroDeformation = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7 };

            floatVector previousPrecedingGradientMicroDeformation = { 0.2, 0.2, 0.7, 0.4, 1.5, 0.6, 0.7, 0.2, 0.9, 1.4, 1.1, 1.2, 1.8, 1.4, 1.5, 1.6, 1.7, 1.8, 5.9, 2.0, 2.1, 2.2, 2.3, 1.4, 2.5, 2.6, 2.7 };

            floatMatrix dPrecedingMicroDeformationdChi = initialize( 9, 9 );

            floatMatrix dPrecedingMicroDeformationdChin = initialize( 9, 9 );

            floatMatrix previousdPrecedingMicroDeformationdChi = initialize( 9, 9 );

            floatMatrix previousdPrecedingMicroDeformationdChin = initialize( 9, 9 );

            floatMatrix dPrecedingGradientMicroDeformationdFn = initialize( 27, 9 );

            floatMatrix dPrecedingGradientMicroDeformationdChi = initialize( 27, 9 );

            floatMatrix dPrecedingGradientMicroDeformationdChin = initialize( 27, 9 );

            floatMatrix dPrecedingGradientMicroDeformationdGradChi = initialize( 27, 27 );

            floatMatrix dPrecedingGradientMicroDeformationdGradChin = initialize( 27, 27 );

            floatMatrix previousdPrecedingGradientMicroDeformationdFn = initialize( 27, 9 );

            floatMatrix previousdPrecedingGradientMicroDeformationdChi = initialize( 27, 9 );

            floatMatrix previousdPrecedingGradientMicroDeformationdChin = initialize( 27, 9 );

            floatMatrix previousdPrecedingGradientMicroDeformationdGradChi = initialize( 27, 27 );

            floatMatrix previousdPrecedingGradientMicroDeformationdGradChin = initialize( 27, 27 );

            floatVector plasticMultipliers = { 1.1, 1.2, 1.3, 1.4, 1.5 };

            floatVector previousPlasticMultipliers = { 1.01, 1.02, 1.03, 1.04, 1.05 };

            floatVector dMacroFlowdDrivingStress = { 0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.09 };

            floatVector dMicroFlowdDrivingStress = { 0.21, 0.32, 0.43, 0.54, 0.15, 0.26, 0.37, 0.18, 0.29 };

            floatMatrix dMicroGradientFlowdDrivingStress = { {  .1,  .2,  .6,  .4,  .5,  .6,  .7,  .8,  .9, 1.0, 1.1, 1.2, 2.3, 1.4, 1.5, 1.6, 1.2, 1.8, 1.9, 2.0, 2.1, 1.2, 2.3, 2.4, 2.5, 2.6, 1.7 },
                                                             { 2.8, 5.9, 2.0, 1.1, 3.2, 4.3, 3.4, 3.6, 2.6, 3.1, 3.6, 3.9, 4.0, 1.1, 4.2, 4.3, 4.4, 5.5, 4.6, 4.7, 8.8, 4.4, 5.0, 2.1, 5.6, 5.3, 5.4 },
                                                             { 1.5, 5.9, 5.7, 5.8, 2.9, 6.0, 3.1, 6.2, 1.3, 5.4, 6.5, 3.6, 6.7, 6.8, 8.9, 7.0, 7.1, 7.2, 1.3, 7.4, 7.5, 7.6, 3.7, 7.8, 7.9, 3.2, 8.1 } };

            floatVector previousdMacroFlowdDrivingStress = { 0.11, -0.22, 0.33, 0.44, -0.55, 0.66, 0.77, -0.88, 0.99 };

            floatVector previousdMicroFlowdDrivingStress = { 0.21, 0.32, 0.43, 0.54, -0.15, -0.26, 0.37, -0.18, 0.29 };

            floatMatrix previousdMicroGradientFlowdDrivingStress = { {  .1,  .2,  .3,  .4,  .5,  .6,  .7,  .8,  .9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7 },
                                                                     { 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4 },
                                                                     { 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1 } };

            floatMatrix d2MacroFlowdDrivingStressdStress = initialize( 9, 9 );

            floatMatrix previousd2MacroFlowdDrivingStressdStress = initialize( 9, 9 );

            floatMatrix d2MacroFlowdDrivingStressdF = initialize( 9, 9 );

            floatMatrix previousd2MacroFlowdDrivingStressdF = initialize( 9, 9 );

            floatMatrix d2MacroFlowdDrivingStressdFn = initialize( 9, 9 );

            floatMatrix previousd2MacroFlowdDrivingStressdFn = initialize( 9, 9 );

            floatMatrix d2MicroFlowdDrivingStressdStress = initialize( 9, 9 );

            floatMatrix previousd2MicroFlowdDrivingStressdStress = initialize( 9, 9 );

            floatMatrix d2MicroFlowdDrivingStressdF = initialize( 9, 9 );

            floatMatrix previousd2MicroFlowdDrivingStressdF = initialize( 9, 9 );

            floatMatrix d2MicroFlowdDrivingStressdFn = initialize( 9, 9 );

            floatMatrix previousd2MicroFlowdDrivingStressdFn = initialize( 9, 9 );

            floatMatrix d2MicroGradientFlowdDrivingStressdStress   = initialize( 3, 27 * 27 );

            floatMatrix d2MicroGradientFlowdDrivingStressdF        = initialize( 3, 27 *  9 );

            floatMatrix d2MicroGradientFlowdDrivingStressdFn       = initialize( 3, 27 *  9 );

            floatMatrix d2MicroGradientFlowdDrivingStressdChi      = initialize( 3, 27 *  9 );

            floatMatrix d2MicroGradientFlowdDrivingStressdChin     = initialize( 3, 27 *  9 );

            floatMatrix previousd2MicroGradientFlowdDrivingStressdStress   = initialize( 3, 27 * 27 );

            floatMatrix previousd2MicroGradientFlowdDrivingStressdF        = initialize( 3, 27 *  9 );

            floatMatrix previousd2MicroGradientFlowdDrivingStressdFn       = initialize( 3, 27 *  9 );

            floatMatrix previousd2MicroGradientFlowdDrivingStressdChi      = initialize( 3, 27 *  9 );

            floatMatrix previousd2MicroGradientFlowdDrivingStressdChin     = initialize( 3, 27 *  9 );

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              2, 0.95338442, 0.74042148,
                                              2, 0.38093104, 0.49241325,
                                              2, 0.82121039, 0.90566759,
                                              2, 0.01166325, 0.05331896,
                                              2, 0.32982199, 0.60161431,
                                              2, 0.58881096, 0.11473813 };

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

            }

            virtual void setPlasticMultipliers( const bool isPrevious ) override{

                if ( isPrevious ){

                    set_previousPlasticMultipliers( previousPlasticMultipliers );

                }
                else{

                    set_plasticMultipliers( plasticMultipliers );

                }

            }

            virtual void setPrecedingDeformationGradient( const bool isPrevious ) override{

                if ( isPrevious ){

                    set_previousPrecedingDeformationGradient( previousPrecedingDeformationGradient );

                }
                else{

                    set_precedingDeformationGradient( precedingDeformationGradient );

                }

            }

            virtual void setPrecedingDeformationGradientJacobians( const bool isPrevious ) override{

                setPrecedingDeformationGradient( isPrevious );

                if ( isPrevious ){

                    set_previousdPrecedingDeformationGradientdF( previousdPrecedingDeformationGradientdF );

                    set_previousdPrecedingDeformationGradientdFn( previousdPrecedingDeformationGradientdFn );

                }
                else{

                    set_dPrecedingDeformationGradientdF( dPrecedingDeformationGradientdF );

                    set_dPrecedingDeformationGradientdFn( dPrecedingDeformationGradientdFn );

                }

            }

            virtual void setPrecedingMicroDeformation( const bool isPrevious ) override{

                if ( isPrevious ){

                    set_previousPrecedingMicroDeformation( previousPrecedingMicroDeformation );

                }
                else{

                    set_precedingMicroDeformation( precedingMicroDeformation );

                }

            }

            virtual void setPrecedingMicroDeformationJacobians( const bool isPrevious ) override{

                setPrecedingMicroDeformation( isPrevious );

                if ( isPrevious ){

                    set_previousdPrecedingMicroDeformationdChi( previousdPrecedingMicroDeformationdChi );

                    set_previousdPrecedingMicroDeformationdChin( previousdPrecedingMicroDeformationdChin );

                }
                else{

                    set_dPrecedingMicroDeformationdChi( dPrecedingMicroDeformationdChi );

                    set_dPrecedingMicroDeformationdChin( dPrecedingMicroDeformationdChin );

                }

            }

            virtual void setPrecedingGradientMicroDeformation( const bool isPrevious ) override{

                if ( isPrevious ){

                    set_previousPrecedingGradientMicroDeformation( previousPrecedingGradientMicroDeformation );

                }
                else{

                    set_precedingGradientMicroDeformation( precedingGradientMicroDeformation );

                }

            }

            virtual void setPrecedingGradientMicroDeformationJacobians( const bool isPrevious ) override{

                setPrecedingGradientMicroDeformation( isPrevious );

                if ( isPrevious ){

                    set_previousdPrecedingGradientMicroDeformationdFn(       previousdPrecedingGradientMicroDeformationdFn );

                    set_previousdPrecedingGradientMicroDeformationdChi(      previousdPrecedingGradientMicroDeformationdChi );

                    set_previousdPrecedingGradientMicroDeformationdChin(     previousdPrecedingGradientMicroDeformationdChin );

                    set_previousdPrecedingGradientMicroDeformationdGradChi(  previousdPrecedingGradientMicroDeformationdGradChi );

                    set_previousdPrecedingGradientMicroDeformationdGradChin( previousdPrecedingGradientMicroDeformationdGradChin );

                }
                else{

                    set_dPrecedingGradientMicroDeformationdFn(       dPrecedingGradientMicroDeformationdFn );

                    set_dPrecedingGradientMicroDeformationdChi(      dPrecedingGradientMicroDeformationdChi );

                    set_dPrecedingGradientMicroDeformationdChin(     dPrecedingGradientMicroDeformationdChin );

                    set_dPrecedingGradientMicroDeformationdGradChi(  dPrecedingGradientMicroDeformationdGradChi );

                    set_dPrecedingGradientMicroDeformationdGradChin( dPrecedingGradientMicroDeformationdGradChin );

                }

            }

            virtual void setFlowPotentialGradients( const bool isPrevious ) override{

                if ( isPrevious ){

                    set_previousdMacroFlowdDrivingStress( previousdMacroFlowdDrivingStress );

                    set_previousdMicroFlowdDrivingStress( previousdMicroFlowdDrivingStress );

                    set_previousdMicroGradientFlowdDrivingStress( previousdMicroGradientFlowdDrivingStress );

                }
                else{

                    set_dMacroFlowdDrivingStress( dMacroFlowdDrivingStress );

                    set_dMicroFlowdDrivingStress( dMicroFlowdDrivingStress );

                    set_dMicroGradientFlowdDrivingStress( dMicroGradientFlowdDrivingStress );

                }

            }

            virtual void setFlowPotentialGradientsJacobians( const bool isPrevious ) override{

                setFlowPotentialGradients( isPrevious );

                if ( isPrevious ){

                    set_previousd2MacroFlowdDrivingStressdStress( previousd2MacroFlowdDrivingStressdStress );

                    set_previousd2MacroFlowdDrivingStressdF( previousd2MacroFlowdDrivingStressdF );

                    set_previousd2MacroFlowdDrivingStressdFn( previousd2MacroFlowdDrivingStressdFn );

                    set_previousd2MicroFlowdDrivingStressdStress( previousd2MicroFlowdDrivingStressdStress );

                    set_previousd2MicroFlowdDrivingStressdF( previousd2MicroFlowdDrivingStressdF );

                    set_previousd2MicroFlowdDrivingStressdFn( previousd2MicroFlowdDrivingStressdFn );

                    set_previousd2MicroGradientFlowdDrivingStressdStress(   previousd2MicroGradientFlowdDrivingStressdStress );

                    set_previousd2MicroGradientFlowdDrivingStressdF(        previousd2MicroGradientFlowdDrivingStressdF );

                    set_previousd2MicroGradientFlowdDrivingStressdFn(       previousd2MicroGradientFlowdDrivingStressdFn );

                    set_previousd2MicroGradientFlowdDrivingStressdChi(      previousd2MicroGradientFlowdDrivingStressdChi );

                    set_previousd2MicroGradientFlowdDrivingStressdChin(     previousd2MicroGradientFlowdDrivingStressdChin );

                }
                else{

                    set_d2MacroFlowdDrivingStressdStress( d2MacroFlowdDrivingStressdStress );

                    set_d2MacroFlowdDrivingStressdF( d2MacroFlowdDrivingStressdF );

                    set_d2MacroFlowdDrivingStressdFn( d2MacroFlowdDrivingStressdFn );

                    set_d2MicroFlowdDrivingStressdStress( d2MicroFlowdDrivingStressdStress );

                    set_d2MicroFlowdDrivingStressdF( d2MicroFlowdDrivingStressdF );

                    set_d2MicroFlowdDrivingStressdFn( d2MicroFlowdDrivingStressdFn );

                    set_d2MicroGradientFlowdDrivingStressdStress(   d2MicroGradientFlowdDrivingStressdStress );

                    set_d2MicroGradientFlowdDrivingStressdF(        d2MicroGradientFlowdDrivingStressdF );

                    set_d2MicroGradientFlowdDrivingStressdFn(       d2MicroGradientFlowdDrivingStressdFn );

                    set_d2MicroGradientFlowdDrivingStressdChi(      d2MicroGradientFlowdDrivingStressdChi );

                    set_d2MicroGradientFlowdDrivingStressdChin(     d2MicroGradientFlowdDrivingStressdChin );

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

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    residualMock RJ( &hydra, 55, 1, stateVariableIndices, parameters );

    floatVector answerMacroL = { -0.05321195,  0.78092663,  0.88794544,
                                  0.29484406,  0.02862809,  0.58429782,
                                  0.44566833,  0.36929872, -0.24728282 };

    floatVector answerPreviousMacroL = { 0.2990453 ,  0.71640881,  0.95746865,
                                         0.87698387, -0.67728052, -0.33595266,
                                         2.37263201,  0.60392039,  2.8536219 };

    floatVector answerMicroL = { -6.21501491e-04,  7.93199942e-01,  4.83983583e-01,
                                  1.35508754e-01, -3.65731890e-01, -1.70009412e-01,
                                  1.53596752e-01,  2.97504168e-01,  2.66968338e-01 };

    floatVector answerPreviousMicroL = { -0.92471373,  0.16655198, -0.51178446,
                                          0.53976276,  0.79789782,  0.70459609,
                                          0.12024267,  0.33339548,  0.51643547 };

    floatVector answerGradientMicroL = {  3.8507487 ,  9.28081242,  4.46027602, 16.50809086, 13.53689506,
                                         18.95857114, 11.62882252,  7.02811906, 24.30162626,  4.99834674,
                                          7.21236049,  7.80149465, -4.85671956, -0.43475823,  0.63767785,
                                          8.87730438,  9.98386271, -2.38221783,  3.06628593,  6.65983251,
                                          9.29346899, 10.46800609,  9.98772098,  7.34645701, -1.88625885,
                                          6.42755711, 12.98624717 };

    floatVector answerPreviousGradientMicroL = { -21.46773019, -25.17130606, -29.16693711, -21.86601317,
                                                 -28.93581213, -31.13912857, -22.02294257, -27.50418022,
                                                 -32.17041273,  -6.1639954 , -10.27462586, -12.96246109,
                                                  -8.64986898, -10.8323775 , -14.64827127,  -8.87965286,
                                                 -10.86375779, -13.52206104,  -8.80618738, -17.72393425,
                                                 -21.79970553, -17.01989502, -20.89282785, -25.07806712,
                                                 -14.26432742, -20.14073019, -24.09169252 };

    RJ.get_dPlasticMacroVelocityGradientdMacroStress( );

    RJ.get_previousdPlasticMacroVelocityGradientdMacroStress( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMacroL,                 *R.get_plasticMacroVelocityGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousMacroL,         *R.get_previousPlasticMacroVelocityGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMicroL,                 *R.get_plasticMicroVelocityGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousMicroL,         *R.get_previousPlasticMicroVelocityGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerGradientMicroL,         *R.get_plasticGradientMicroVelocityGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousGradientMicroL, *R.get_previousPlasticGradientMicroVelocityGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMacroL,                 *RJ.get_plasticMacroVelocityGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousMacroL,         *RJ.get_previousPlasticMacroVelocityGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerMicroL,                 *RJ.get_plasticMicroVelocityGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousMicroL,         *RJ.get_previousPlasticMicroVelocityGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerGradientMicroL,         *RJ.get_plasticGradientMicroVelocityGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answerPreviousGradientMicroL, *RJ.get_previousPlasticGradientMicroVelocityGradient( ) ) );

}

BOOST_AUTO_TEST_CASE( test_setPlasticVelocityGradients2 ){
    /*!
     * Test setting the cohesion values
     */

    std::cout << "starting gradient check\n";

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

            floatVector previousPK2       = {  1,  2,  3,  4,  5,  6,  7,  8,  9 };

            floatVector previousSIGMA     = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

            floatVector previousM         = { 19, 20, 21, 22, 23, 24, 25, 26, 27,
                                              28, 29, 30, 31, 32, 33, 34, 35, 36,
                                              37, 38, 39, 40, 41, 42, 43, 44, 45 };

        protected:

            using tardigradeHydra::residualBase::setPreviousStress;

            virtual void setPreviousStress( ) override{

                setPreviousStress( tardigradeVectorTools::appendVectors( { previousPK2, previousSIGMA, previousM } ) );

            }

    };

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticity::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::residual;

            floatVector plasticParameters = { 2, 0.53895133, 0.37172145,
                                              2, 0.37773052, 0.92739145,
                                              2, 0.53186824, 0.75454313,
                                              2, 0.95338442, 0.74042148,
                                              2, 0.38093104, 0.49241325,
                                              2, 0.82121039, 0.90566759,
                                              2, 0.01166325, 0.05331896,
                                              2, 0.32982199, 0.60161431,
                                              2, 0.58881096, 0.11473813 };

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

            }

    };

    class hydraBaseMicromorphicMock : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

            stressMock elasticity;

            residualMock plasticity;

            std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector _local_deltaPK2    = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            floatVector _local_deltaSIGMA  = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            floatVector _local_deltaM      = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0,
                                               0, 0, 0, 0, 0, 0, 0, 0, 0 };

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

                elasticity.previousPK2   += _local_deltaPK2;

                elasticity.previousSIGMA += _local_deltaSIGMA;

                elasticity.previousM     += _local_deltaM;

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

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    // Test the Jacobians

    floatMatrix dMacroLdX(                9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dMicroLdX(                9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dMicroGradientLdX(       27, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dMacroLdF(                9, floatVector(  9, 0 ) );

    floatMatrix dMicroLdF(                9, floatVector(  9, 0 ) );

    floatMatrix dMicroGradientLdF(       27, floatVector(  9, 0 ) );

    floatMatrix dMicroLdChi(              9, floatVector(  9, 0 ) );

    floatMatrix dMicroGradientLdChi(     27, floatVector(  9, 0 ) );

    floatMatrix dMicroGradientLdGradChi( 27, floatVector( 27, 0 ) );

    floatMatrix previousdMacroLdStress(            9, floatVector( configuration_unknown_count, 0 ) );

    floatMatrix previousdMicroLdStress(            9, floatVector( configuration_unknown_count, 0 ) );

    floatMatrix previousdMicroGradientLdStress(   27, floatVector( configuration_unknown_count, 0 ) );

    floatMatrix previousdMacroLdISVs(              9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix previousdMicroLdISVs(              9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix previousdMicroGradientLdISVs(     27, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix previousdMacroLdF(                 9, floatVector(  9, 0 ) );

    floatMatrix previousdMicroLdF(                 9, floatVector(  9, 0 ) );

    floatMatrix previousdMicroGradientLdF(        27, floatVector(  9, 0 ) );

    floatMatrix previousdMicroLdChi(               9, floatVector(  9, 0 ) );

    floatMatrix previousdMicroGradientLdChi(      27, floatVector(  9, 0 ) );

    floatMatrix previousdMicroGradientLdGradChi(  27, floatVector( 27, 0 ) );

    floatType eps = 1e-6;

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_plasticMacroVelocityGradient( );

        floatVector vm = *Rm.get_plasticMacroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMacroLdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_plasticMicroVelocityGradient( );

        vm = *Rm.get_plasticMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicroLdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_plasticGradientMicroVelocityGradient( );

        vm = *Rm.get_plasticGradientMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicroGradientLdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_dMacroLdX( 9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix assembled_dMicroLdX( 9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix assembled_dMicroGradientLdX( 27, floatVector( unknownVector.size( ), 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dMacroLdX[ i ][ j ]     = ( *R.get_dPlasticMacroVelocityGradientdMacroStress( ) )[ i ][ j ];

            assembled_dMacroLdX[ i ][ j + 9 ] = ( *R.get_dPlasticMacroVelocityGradientdMicroStress( ) )[ i ][ j ];

            assembled_dMicroLdX[ i ][ j + 9 ] = ( *R.get_dPlasticMicroVelocityGradientdMicroStress( ) )[ i ][ j ];

            assembled_dMacroLdX[ i ][ j + configuration_unknown_count ] = ( *R.get_dPlasticMacroVelocityGradientdFn( ) )[ i ][ j ];

            assembled_dMicroLdX[ i ][ j + configuration_unknown_count ] = ( *R.get_dPlasticMicroVelocityGradientdFn( ) )[ i ][ j ];

            assembled_dMicroLdX[ i ][ j + configuration_unknown_count + 9 ] = ( *R.get_dPlasticMicroVelocityGradientdChin( ) )[ i ][ j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dMacroLdX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dPlasticMacroVelocityGradientdStateVariables( ) )[ i ][ j ];

            assembled_dMicroLdX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dPlasticMicroVelocityGradientdStateVariables( ) )[ i ][ j ];

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dMicroGradientLdX[ i ][ j + 9 ] = ( *R.get_dPlasticGradientMicroVelocityGradientdMicroStress( ) )[ i ][ j ];

            assembled_dMicroGradientLdX[ i ][ j + configuration_unknown_count ] = ( *R.get_dPlasticGradientMicroVelocityGradientdFn( ) )[ i ][ j ];

            assembled_dMicroGradientLdX[ i ][ j + configuration_unknown_count + 9 ] = ( *R.get_dPlasticGradientMicroVelocityGradientdChin( ) )[ i ][ j ];

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            assembled_dMicroGradientLdX[ i ][ j + 18 ] = ( *R.get_dPlasticGradientMicroVelocityGradientdHigherOrderStress( ) )[ i ][ j ];

            assembled_dMicroGradientLdX[ i ][ j + configuration_unknown_count + 18 ] = ( *R.get_dPlasticGradientMicroVelocityGradientdGradChin( ) )[ i ][ j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dMicroGradientLdX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dPlasticGradientMicroVelocityGradientdStateVariables( ) )[ i ][ j ];

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMacroLdX, assembled_dMacroLdX ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroLdX, assembled_dMicroLdX ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroGradientLdX, assembled_dMicroGradientLdX ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_plasticMacroVelocityGradient( );

        floatVector vm = *Rm.get_plasticMacroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMacroLdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_plasticMicroVelocityGradient( );

        vm = *Rm.get_plasticMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicroLdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_plasticGradientMicroVelocityGradient( );

        vm = *Rm.get_plasticGradientMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicroGradientLdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMacroLdF, *R.get_dPlasticMacroVelocityGradientdF( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroLdF, *R.get_dPlasticMicroVelocityGradientdF( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroGradientLdF, *R.get_dPlasticGradientMicroVelocityGradientdF( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_plasticMacroVelocityGradient( );

        floatVector vm = *Rm.get_plasticMacroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        vp = *Rp.get_plasticMicroVelocityGradient( );

        vm = *Rm.get_plasticMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicroLdChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_plasticGradientMicroVelocityGradient( );

        vm = *Rm.get_plasticGradientMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicroGradientLdChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroLdChi, *R.get_dPlasticMicroVelocityGradientdChi( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroGradientLdChi, *R.get_dPlasticGradientMicroVelocityGradientdChi( ) ) );

    for ( unsigned int i = 0; i < 27; i++ ){

        floatVector delta( 27, 0 );

        delta[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation + delta, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation - delta, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_plasticMacroVelocityGradient( );

        floatVector vm = *Rm.get_plasticMacroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        vp = *Rp.get_plasticMicroVelocityGradient( );

        vm = *Rm.get_plasticMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        vp = *Rp.get_plasticGradientMicroVelocityGradient( );

        vm = *Rm.get_plasticGradientMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicroGradientLdGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroGradientLdGradChi, *R.get_dPlasticGradientMicroVelocityGradientdGradChi( ) ) );

    // Check previous Jacobians
    for ( unsigned int i = 0; i < configuration_unknown_count; i++ ){

        floatVector delta( configuration_unknown_count, 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        floatVector dpk2(   delta.begin( ) +  0, delta.begin( ) +  9 );
        floatVector dsigma( delta.begin( ) +  9, delta.begin( ) + 18 );
        floatVector dm(     delta.begin( ) + 18, delta.end( )        );

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydrap._local_deltaPK2   = dpk2;

        hydrap._local_deltaSIGMA = dsigma;

        hydrap._local_deltaM     = dm;

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydram._local_deltaPK2   = -dpk2;

        hydram._local_deltaSIGMA = -dsigma;

        hydram._local_deltaM     = -dm;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousPlasticMacroVelocityGradient( );

        floatVector vm = *Rm.get_previousPlasticMacroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMacroLdStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousPlasticMicroVelocityGradient( );

        vm = *Rm.get_previousPlasticMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroLdStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousPlasticGradientMicroVelocityGradient( );

        vm = *Rm.get_previousPlasticGradientMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroGradientLdStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_previousdMacroLdStress( 9, floatVector( configuration_unknown_count, 0 ) );

    floatMatrix assembled_previousdMicroLdStress( 9, floatVector( configuration_unknown_count, 0 ) );

    floatMatrix assembled_previousdMicroGradientLdStress( 27, floatVector( configuration_unknown_count, 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_previousdMacroLdStress[ i ][ j ]     = ( *R.get_previousdPlasticMacroVelocityGradientdMacroStress( ) )[ i ][ j ];

            assembled_previousdMacroLdStress[ i ][ j + 9 ] = ( *R.get_previousdPlasticMacroVelocityGradientdMicroStress( ) )[ i ][ j ];

            assembled_previousdMicroLdStress[ i ][ j + 9 ] = ( *R.get_previousdPlasticMicroVelocityGradientdMicroStress( ) )[ i ][ j ];

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_previousdMicroGradientLdStress[ i ][ j + 9 ] = ( *R.get_previousdPlasticGradientMicroVelocityGradientdMicroStress( ) )[ i ][ j ];

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            assembled_previousdMicroGradientLdStress[ i ][ j + 18 ] = ( *R.get_previousdPlasticGradientMicroVelocityGradientdHigherOrderStress( ) )[ i ][ j ];

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMacroLdStress, assembled_previousdMacroLdStress ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroLdStress, assembled_previousdMicroLdStress ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroGradientLdStress, assembled_previousdMicroGradientLdStress ) );

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousPlasticMacroVelocityGradient( );

        floatVector vm = *Rm.get_previousPlasticMacroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMacroLdISVs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousPlasticMicroVelocityGradient( );

        vm = *Rm.get_previousPlasticMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroLdISVs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousPlasticGradientMicroVelocityGradient( );

        vm = *Rm.get_previousPlasticGradientMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroGradientLdISVs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_previousdMacroLdISVs( 9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix assembled_previousdMicroLdISVs( 9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix assembled_previousdMicroGradientLdISVs( 27, floatVector( previousStateVariables.size( ), 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_previousdMacroLdISVs[ i ][ j ] = ( *R.get_previousdPlasticMacroVelocityGradientdFn( ) )[ i ][ j ];

            assembled_previousdMicroLdISVs[ i ][ j ] = ( *R.get_previousdPlasticMicroVelocityGradientdFn( ) )[ i ][ j ];

            assembled_previousdMicroLdISVs[ i ][ j + 9 ] = ( *R.get_previousdPlasticMicroVelocityGradientdChin( ) )[ i ][ j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_previousdMacroLdISVs[ i ][ j + configuration_unknown_count ] = ( *R.get_previousdPlasticMacroVelocityGradientdStateVariables( ) )[ i ][ j ];

            assembled_previousdMicroLdISVs[ i ][ j + configuration_unknown_count ] = ( *R.get_previousdPlasticMicroVelocityGradientdStateVariables( ) )[ i ][ j ];

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_previousdMicroGradientLdISVs[ i ][ j ] = ( *R.get_previousdPlasticGradientMicroVelocityGradientdFn( ) )[ i ][ j ];

            assembled_previousdMicroGradientLdISVs[ i ][ j + 9 ] = ( *R.get_previousdPlasticGradientMicroVelocityGradientdChin( ) )[ i ][ j ];

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            assembled_previousdMicroGradientLdISVs[ i ][ j + 18 ] = ( *R.get_previousdPlasticGradientMicroVelocityGradientdGradChin( ) )[ i ][ j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_previousdMicroGradientLdISVs[ i ][ j + configuration_unknown_count ] = ( *R.get_previousdPlasticGradientMicroVelocityGradientdStateVariables( ) )[ i ][ j ];

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMacroLdISVs, assembled_previousdMacroLdISVs ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroLdISVs, assembled_previousdMicroLdISVs ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dMicroGradientLdX, assembled_dMicroGradientLdX ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousPlasticMacroVelocityGradient( );

        floatVector vm = *Rm.get_previousPlasticMacroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMacroLdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousPlasticMicroVelocityGradient( );

        vm = *Rm.get_previousPlasticMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroLdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousPlasticGradientMicroVelocityGradient( );

        vm = *Rm.get_previousPlasticGradientMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroGradientLdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMacroLdF, *R.get_previousdPlasticMacroVelocityGradientdF( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroLdF, *R.get_previousdPlasticMicroVelocityGradientdF( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroGradientLdF, *R.get_previousdPlasticGradientMicroVelocityGradientdF( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

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

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousPlasticMacroVelocityGradient( );

        floatVector vm = *Rm.get_previousPlasticMacroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        vp = *Rp.get_previousPlasticMicroVelocityGradient( );

        vm = *Rm.get_previousPlasticMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroLdChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousPlasticGradientMicroVelocityGradient( );

        vm = *Rm.get_previousPlasticGradientMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroGradientLdChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroLdChi, *R.get_previousdPlasticMicroVelocityGradientdChi( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroGradientLdChi, *R.get_previousdPlasticGradientMicroVelocityGradientdChi( ) ) );

    for ( unsigned int i = 0; i < 27; i++ ){

        floatVector delta( 27, 0 );

        delta[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation + delta,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation - delta,
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_previousPlasticMacroVelocityGradient( );

        floatVector vm = *Rm.get_previousPlasticMacroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        vp = *Rp.get_previousPlasticMicroVelocityGradient( );

        vm = *Rm.get_previousPlasticMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0., ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) ) );

        }

        vp = *Rp.get_previousPlasticGradientMicroVelocityGradient( );

        vm = *Rm.get_previousPlasticGradientMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroGradientLdGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousdMicroGradientLdGradChi, *R.get_previousdPlasticGradientMicroVelocityGradientdGradChi( ) ) );

}
