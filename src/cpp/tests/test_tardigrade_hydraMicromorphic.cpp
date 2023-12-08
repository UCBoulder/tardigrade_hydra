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

    floatVector deformationGradient = { 1, 2, 2, 4, 5, 6, 7, 8, 9 };

    floatVector previousDeformationGradient = { 10, 11, 12, 13, 14, 15, 16, 16, 18 };

    floatVector microDeformation = { 18, 20, 21, 22, 23, 24, 25, 26, 27 };

    floatVector previousMicroDeformation = { 28, 29, 30, 31, 32, 33, 34, 35, 33 };

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
                                          19.00, 20.00, 21.00, 22.00, 23.00, 24.00, 25.00, 26.00, 27.00,
                                          -1.00, -2.00, -3.00, -4.00, -5.00, -6.00 };

    floatMatrix configurationsAnswer = { { 0.80151542, 1.75875283, 1.71599025, 3.43145793, 4.30238634, 5.17331475, 6.12262672, 6.91517157, 7.70771643 },
                                         { 1.01, 0.02, 0.03, 0.04, 1.05, 0.06, 0.07, 0.08, 1.09 } };

    floatMatrix previousConfigurationsAnswer = { { 8.81379551,  9.5279568 , 10.2421181 , 11.50496429, 12.14074203, 12.77651977, 14.23044766, 13.79655112, 15.36265459 },
                                                 { 1.01, 0.02, 0.03, 0.04, 1.05, 0.06, 0.07, 0.08, 1.09 } };

    floatMatrix inverseConfigurationsAnswer = { { -1.        , -0.64666667,  0.65666667,  2.        , -1.65666667,         0.66666667, -1.        ,  2.        , -0.99       },
                                                {  0.99259711, -0.01689601, -0.02638913, -0.03431458,  0.95697614,        -0.05173315, -0.06122627, -0.06915172,  0.92292284 } };

    floatMatrix previousInverseConfigurationsAnswer = { { -1.96      ,  0.97      ,  0.5       , -0.97      ,  1.98      ,        -1.        ,  2.68666667, -2.67666667,  0.5        },
                                                        {  0.99259711, -0.01689601, -0.02638913, -0.03431458,  0.95697614,        -0.05173315, -0.06122627, -0.06915172,  0.92292284 } };

    floatMatrix microConfigurationsAnswer = { { 20.92702193, 23.6257958 , 25.32456968, 25.37137468, 27.1869842 , 29.00259373, 28.8021693 , 30.72388588, 32.64560245 },
                                              { 0.99, -0.02, -0.03, -0.04,  0.95, -0.06, -0.07, -0.08,  0.91 } };

    floatMatrix inverseMicroConfigurationsAnswer = { { -1.        ,   1.92      ,  -0.93      ,   2.        , -13.07      ,  10.06      ,  -1.        ,  10.60666667,  -8.61666667 },
                                                     {  1.01355812, 0.02428672, 0.03501533, 0.04786607, 1.05965574, 0.07144541, 0.08217402, 0.09502476, 1.1078755 } };

    floatMatrix previousMicroConfigurationsAnswer = { { 32.23296392, 34.26078755, 36.28861118, 35.66375855, 37.79768922, 39.9316199 , 38.84803112, 41.04951662, 40.25100212 },
                                                      { 0.99, -0.02, -0.03, -0.04,  0.95, -0.06, -0.07, -0.08,  0.91 } };

    floatMatrix previousInverseMicroConfigurationsAnswer = { { -11.1       ,  10.42333333,  -0.33333333,  10.91      , -10.58666667,   0.66666667,  -0.41333333,   0.73666667,  -0.33333333 },
                                                             { 1.01355812, 0.02428672, 0.03501533, 0.04786607, 1.05965574, 0.07144541, 0.08217402, 0.09502476, 1.1078755 } };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 6;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( configurationsAnswer, *hydra.getConfigurations( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousConfigurationsAnswer, *hydra.getPreviousConfigurations( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( inverseConfigurationsAnswer, *hydra.getInverseConfigurations( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousInverseConfigurationsAnswer, *hydra.getPreviousInverseConfigurations( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( microConfigurationsAnswer, *hydra.getMicroConfigurations( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousMicroConfigurationsAnswer, *hydra.getPreviousMicroConfigurations( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( inverseMicroConfigurationsAnswer, *hydra.getInverseMicroConfigurations( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousInverseMicroConfigurationsAnswer, *hydra.getPreviousInverseMicroConfigurations( ) ) );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getSubMicroConfiguration ){

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 1.19646919, -0.21386067, -0.27314855,
                                        0.05131477,  1.21946897, -0.07689354,
                                        0.4807642 ,  0.18482974,  0.9809319  };

    floatVector previousDeformationGradient = { 0.89211752, -0.15682198,  0.22904971,
                                               -0.06142776,  0.5596779 , -0.10195574,
                                                0.23799541, -0.31750827,  0.67545176 };

    floatVector microDeformation = { 1.03155137,  0.03182759,  0.13440096,
                                     0.34943179,  1.22445532,  0.11102351,
                                     0.22244338, -0.17704109,  0.86178866 };

    floatVector previousMicroDeformation = { 0.72826323, -0.20628595,  0.13097612,
                                            -0.40789506,  0.93370117, -0.06913724,
                                            -0.0063149 , -0.07416971,  0.81226122 };

    floatVector gradientMicroDeformation = {-0.07364869,  0.39338916,  0.44416002,  0.00183668,  0.12395295,
                                            -0.3843816 , -0.18271452, -0.08517379,  0.36630916, -0.24954463,
                                            -0.01696574,  0.48555979,  0.01948512,  0.11289453, -0.37937133,
                                             0.3263408 ,  0.10306013,  0.04506801, -0.15723617, -0.19587921,
                                            -0.08297779,  0.18130077,  0.37545684,  0.01042234,  0.16931378,
                                             0.08593655,  0.1249035 };

    floatVector previousGradientMicroDeformation = { 0.17468905,  0.34234244, -0.41680501,  0.26368284, -0.25633363,
                                                    -0.30577704,  0.07245696, -0.40428748,  0.38532683,  0.12724897,
                                                     0.22341636, -0.48387079,  0.09443188,  0.05678519, -0.34104036,
                                                    -0.34692948,  0.19552953, -0.18123357,  0.1919703 ,  0.05438325,
                                                    -0.11104943,  0.42513249,  0.34167   , -0.14260243, -0.45640854,
                                                    -0.19523193, -0.10181432};

    floatVector previousStateVariables = { 0.20495883,  0.49535848, -0.14408513,  0.26254781,  0.09317692,
                                           0.1917018 , -0.34887255, -0.10112371, -0.2591441 , -0.15654399,
                                           0.01312815,  0.16662455, -0.39409151, -0.36910505, -0.17801939,
                                           0.16156434,  0.34650623,  0.05325734,  0.35445249, -0.11516219,
                                          -0.1832121 , -0.14573532, -0.32891817,  0.32911263, -0.16132915,
                                           0.05237008,  0.07855147,  0.02153306, -0.49731194,  0.48834542,
                                           0.40534158, -0.29236414, -0.20751059,  0.02001015,  0.40191137,
                                           0.48363088, -0.24245794,  0.06435904,  0.30696868, -0.10562995,
                                           0.23107304, -0.33893099,  0.10069857,  0.36586446,  0.48352161,
                                          -0.42063421, -0.07165273, -0.29545714, -0.04936351,  0.04776357,
                                          -0.40667329, -0.20313922,  0.42758424,  0.06900373, -0.042588  ,
                                           0.25352599,  0.24186215, -0.45142097,  0.2086974 ,  0.33924335,
                                          -0.33406212,  0.28099794, -0.21346338, -0.19353025,  0.16526147,
                                          -0.38860783,  0.16487245,  0.38785679,  0.19631127, -0.05967212,
                                          -0.06178562,  0.2650961 ,  0.065642  , -0.41509584,  0.08267109,
                                           0.3148437 , -0.16293362,  0.42757658,  0.250717  ,  0.07406383,
                                           0.25164399, -0.42085104,  0.35938908,  0.32150411,  0.40987166,
                                          -0.3713688 , -0.41821991, -0.36158443, -0.10062129, -0.07569314,
                                           0.06221838, -0.37775645, -0.2986005 ,  0.31164435, -0.03201243,
                                           0.30793821, -0.49257362,  0.05159273,  0.43193215,  0.08217546,
                                          -0.29390427,  0.21775756, -0.12101415,  0.16838395, -0.47068028,
                                           0.13590036, -0.46780207,  0.24478066, -0.027087  , -0.37824564,
                                           0.04263593, -0.43322556,  0.15336487,  0.49608633,  0.26939734,
                                           0.07377411, -0.39736474,  0.19983407,  0.16116787, -0.45090287,
                                           0.2922993 ,  0.01871659, -0.07413231,  0.28818717, -0.08843078,
                                          -0.01897372, -0.31837116, -0.1786811 ,  0.345533  , -0.31309625,
                                          -0.08270894,  0.48903451, -0.26340019,  0.41683233,  0.41839747,
                                          -0.40870366, -0.03634728,  0.00221634, -0.18633105, -0.45266046,
                                           0.24168564,  0.09552964,  0.23824991 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

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

    floatVector answer = { 2.63306967, -0.16982965, -0.15109712,
                           1.26964411,  0.59017992,  0.35717977,
                           0.72297358, -0.4912966 ,  0.28183864 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, hydra.getSubMicroConfiguration( 0, 2 ) ) );

    floatVector answer2 = { 0.87556098, -0.36781419,  1.20660407,
                            0.21141974,  0.82131808, -0.42325899,
                            0.12210416,  1.0388779 ,  2.07092075 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, hydra.getSubMicroConfiguration( 1, 3 ) ) );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getPreceedingMicroConfiguration ){

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 1.19646919, -0.21386067, -0.27314855,
                                        0.05131477,  1.21946897, -0.07689354,
                                        0.4807642 ,  0.18482974,  0.9809319  };

    floatVector previousDeformationGradient = { 0.89211752, -0.15682198,  0.22904971,
                                               -0.06142776,  0.5596779 , -0.10195574,
                                                0.23799541, -0.31750827,  0.67545176 };

    floatVector microDeformation = { 1.03155137,  0.03182759,  0.13440096,
                                     0.34943179,  1.22445532,  0.11102351,
                                     0.22244338, -0.17704109,  0.86178866 };

    floatVector previousMicroDeformation = { 0.72826323, -0.20628595,  0.13097612,
                                            -0.40789506,  0.93370117, -0.06913724,
                                            -0.0063149 , -0.07416971,  0.81226122 };

    floatVector gradientMicroDeformation = {-0.07364869,  0.39338916,  0.44416002,  0.00183668,  0.12395295,
                                            -0.3843816 , -0.18271452, -0.08517379,  0.36630916, -0.24954463,
                                            -0.01696574,  0.48555979,  0.01948512,  0.11289453, -0.37937133,
                                             0.3263408 ,  0.10306013,  0.04506801, -0.15723617, -0.19587921,
                                            -0.08297779,  0.18130077,  0.37545684,  0.01042234,  0.16931378,
                                             0.08593655,  0.1249035 };

    floatVector previousGradientMicroDeformation = { 0.17468905,  0.34234244, -0.41680501,  0.26368284, -0.25633363,
                                                    -0.30577704,  0.07245696, -0.40428748,  0.38532683,  0.12724897,
                                                     0.22341636, -0.48387079,  0.09443188,  0.05678519, -0.34104036,
                                                    -0.34692948,  0.19552953, -0.18123357,  0.1919703 ,  0.05438325,
                                                    -0.11104943,  0.42513249,  0.34167   , -0.14260243, -0.45640854,
                                                    -0.19523193, -0.10181432};

    floatVector previousStateVariables = { 0.20495883,  0.49535848, -0.14408513,  0.26254781,  0.09317692,
                                           0.1917018 , -0.34887255, -0.10112371, -0.2591441 , -0.15654399,
                                           0.01312815,  0.16662455, -0.39409151, -0.36910505, -0.17801939,
                                           0.16156434,  0.34650623,  0.05325734,  0.35445249, -0.11516219,
                                          -0.1832121 , -0.14573532, -0.32891817,  0.32911263, -0.16132915,
                                           0.05237008,  0.07855147,  0.02153306, -0.49731194,  0.48834542,
                                           0.40534158, -0.29236414, -0.20751059,  0.02001015,  0.40191137,
                                           0.48363088, -0.24245794,  0.06435904,  0.30696868, -0.10562995,
                                           0.23107304, -0.33893099,  0.10069857,  0.36586446,  0.48352161,
                                          -0.42063421, -0.07165273, -0.29545714, -0.04936351,  0.04776357,
                                          -0.40667329, -0.20313922,  0.42758424,  0.06900373, -0.042588  ,
                                           0.25352599,  0.24186215, -0.45142097,  0.2086974 ,  0.33924335,
                                          -0.33406212,  0.28099794, -0.21346338, -0.19353025,  0.16526147,
                                          -0.38860783,  0.16487245,  0.38785679,  0.19631127, -0.05967212,
                                          -0.06178562,  0.2650961 ,  0.065642  , -0.41509584,  0.08267109,
                                           0.3148437 , -0.16293362,  0.42757658,  0.250717  ,  0.07406383,
                                           0.25164399, -0.42085104,  0.35938908,  0.32150411,  0.40987166,
                                          -0.3713688 , -0.41821991, -0.36158443, -0.10062129, -0.07569314,
                                           0.06221838, -0.37775645, -0.2986005 ,  0.31164435, -0.03201243,
                                           0.30793821, -0.49257362,  0.05159273,  0.43193215,  0.08217546,
                                          -0.29390427,  0.21775756, -0.12101415,  0.16838395, -0.47068028,
                                           0.13590036, -0.46780207,  0.24478066, -0.027087  , -0.37824564,
                                           0.04263593, -0.43322556,  0.15336487,  0.49608633,  0.26939734,
                                           0.07377411, -0.39736474,  0.19983407,  0.16116787, -0.45090287,
                                           0.2922993 ,  0.01871659, -0.07413231,  0.28818717, -0.08843078,
                                          -0.01897372, -0.31837116, -0.1786811 ,  0.345533  , -0.31309625,
                                          -0.08270894,  0.48903451, -0.26340019,  0.41683233,  0.41839747,
                                          -0.40870366, -0.03634728,  0.00221634, -0.18633105, -0.45266046,
                                           0.24168564,  0.09552964,  0.23824991 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

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

    floatVector answer = { 2.63306967, -0.16982965, -0.15109712,
                           1.26964411,  0.59017992,  0.35717977,
                           0.72297358, -0.4912966 ,  0.28183864 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, hydra.getPrecedingMicroConfiguration( 2 ) ) );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getFollowingMicroConfiguration ){

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 1.19646919, -0.21386067, -0.27314855,
                                        0.05131477,  1.21946897, -0.07689354,
                                        0.4807642 ,  0.18482974,  0.9809319  };

    floatVector previousDeformationGradient = { 0.89211752, -0.15682198,  0.22904971,
                                               -0.06142776,  0.5596779 , -0.10195574,
                                                0.23799541, -0.31750827,  0.67545176 };

    floatVector microDeformation = { 1.03155137,  0.03182759,  0.13440096,
                                     0.34943179,  1.22445532,  0.11102351,
                                     0.22244338, -0.17704109,  0.86178866 };

    floatVector previousMicroDeformation = { 0.72826323, -0.20628595,  0.13097612,
                                            -0.40789506,  0.93370117, -0.06913724,
                                            -0.0063149 , -0.07416971,  0.81226122 };

    floatVector gradientMicroDeformation = {-0.07364869,  0.39338916,  0.44416002,  0.00183668,  0.12395295,
                                            -0.3843816 , -0.18271452, -0.08517379,  0.36630916, -0.24954463,
                                            -0.01696574,  0.48555979,  0.01948512,  0.11289453, -0.37937133,
                                             0.3263408 ,  0.10306013,  0.04506801, -0.15723617, -0.19587921,
                                            -0.08297779,  0.18130077,  0.37545684,  0.01042234,  0.16931378,
                                             0.08593655,  0.1249035 };

    floatVector previousGradientMicroDeformation = { 0.17468905,  0.34234244, -0.41680501,  0.26368284, -0.25633363,
                                                    -0.30577704,  0.07245696, -0.40428748,  0.38532683,  0.12724897,
                                                     0.22341636, -0.48387079,  0.09443188,  0.05678519, -0.34104036,
                                                    -0.34692948,  0.19552953, -0.18123357,  0.1919703 ,  0.05438325,
                                                    -0.11104943,  0.42513249,  0.34167   , -0.14260243, -0.45640854,
                                                    -0.19523193, -0.10181432};

    floatVector previousStateVariables = { 0.20495883,  0.49535848, -0.14408513,  0.26254781,  0.09317692,
                                           0.1917018 , -0.34887255, -0.10112371, -0.2591441 , -0.15654399,
                                           0.01312815,  0.16662455, -0.39409151, -0.36910505, -0.17801939,
                                           0.16156434,  0.34650623,  0.05325734,  0.35445249, -0.11516219,
                                          -0.1832121 , -0.14573532, -0.32891817,  0.32911263, -0.16132915,
                                           0.05237008,  0.07855147,  0.02153306, -0.49731194,  0.48834542,
                                           0.40534158, -0.29236414, -0.20751059,  0.02001015,  0.40191137,
                                           0.48363088, -0.24245794,  0.06435904,  0.30696868, -0.10562995,
                                           0.23107304, -0.33893099,  0.10069857,  0.36586446,  0.48352161,
                                          -0.42063421, -0.07165273, -0.29545714, -0.04936351,  0.04776357,
                                          -0.40667329, -0.20313922,  0.42758424,  0.06900373, -0.042588  ,
                                           0.25352599,  0.24186215, -0.45142097,  0.2086974 ,  0.33924335,
                                          -0.33406212,  0.28099794, -0.21346338, -0.19353025,  0.16526147,
                                          -0.38860783,  0.16487245,  0.38785679,  0.19631127, -0.05967212,
                                          -0.06178562,  0.2650961 ,  0.065642  , -0.41509584,  0.08267109,
                                           0.3148437 , -0.16293362,  0.42757658,  0.250717  ,  0.07406383,
                                           0.25164399, -0.42085104,  0.35938908,  0.32150411,  0.40987166,
                                          -0.3713688 , -0.41821991, -0.36158443, -0.10062129, -0.07569314,
                                           0.06221838, -0.37775645, -0.2986005 ,  0.31164435, -0.03201243,
                                           0.30793821, -0.49257362,  0.05159273,  0.43193215,  0.08217546,
                                          -0.29390427,  0.21775756, -0.12101415,  0.16838395, -0.47068028,
                                           0.13590036, -0.46780207,  0.24478066, -0.027087  , -0.37824564,
                                           0.04263593, -0.43322556,  0.15336487,  0.49608633,  0.26939734,
                                           0.07377411, -0.39736474,  0.19983407,  0.16116787, -0.45090287,
                                           0.2922993 ,  0.01871659, -0.07413231,  0.28818717, -0.08843078,
                                          -0.01897372, -0.31837116, -0.1786811 ,  0.345533  , -0.31309625,
                                          -0.08270894,  0.48903451, -0.26340019,  0.41683233,  0.41839747,
                                          -0.40870366, -0.03634728,  0.00221634, -0.18633105, -0.45266046,
                                           0.24168564,  0.09552964,  0.23824991 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

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

    floatVector answer = { 0.37335959,  0.14440808,  0.07815635,
                          -0.05311829,  1.15252061, -0.83175389,
                          -0.26108048,  1.01045459,  1.40735072 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, hydra.getFollowingMicroConfiguration( 1 ) ) );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getMicroConfiguration ){

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 1.19646919, -0.21386067, -0.27314855,
                                        0.05131477,  1.21946897, -0.07689354,
                                        0.4807642 ,  0.18482974,  0.9809319  };

    floatVector previousDeformationGradient = { 0.89211752, -0.15682198,  0.22904971,
                                               -0.06142776,  0.5596779 , -0.10195574,
                                                0.23799541, -0.31750827,  0.67545176 };

    floatVector microDeformation = { 1.03155137,  0.03182759,  0.13440096,
                                     0.34943179,  1.22445532,  0.11102351,
                                     0.22244338, -0.17704109,  0.86178866 };

    floatVector previousMicroDeformation = { 0.72826323, -0.20628595,  0.13097612,
                                            -0.40789506,  0.93370117, -0.06913724,
                                            -0.0063149 , -0.07416971,  0.81226122 };

    floatVector gradientMicroDeformation = {-0.07364869,  0.39338916,  0.44416002,  0.00183668,  0.12395295,
                                            -0.3843816 , -0.18271452, -0.08517379,  0.36630916, -0.24954463,
                                            -0.01696574,  0.48555979,  0.01948512,  0.11289453, -0.37937133,
                                             0.3263408 ,  0.10306013,  0.04506801, -0.15723617, -0.19587921,
                                            -0.08297779,  0.18130077,  0.37545684,  0.01042234,  0.16931378,
                                             0.08593655,  0.1249035 };

    floatVector previousGradientMicroDeformation = { 0.17468905,  0.34234244, -0.41680501,  0.26368284, -0.25633363,
                                                    -0.30577704,  0.07245696, -0.40428748,  0.38532683,  0.12724897,
                                                     0.22341636, -0.48387079,  0.09443188,  0.05678519, -0.34104036,
                                                    -0.34692948,  0.19552953, -0.18123357,  0.1919703 ,  0.05438325,
                                                    -0.11104943,  0.42513249,  0.34167   , -0.14260243, -0.45640854,
                                                    -0.19523193, -0.10181432};

    floatVector previousStateVariables = { 0.20495883,  0.49535848, -0.14408513,  0.26254781,  0.09317692,
                                           0.1917018 , -0.34887255, -0.10112371, -0.2591441 , -0.15654399,
                                           0.01312815,  0.16662455, -0.39409151, -0.36910505, -0.17801939,
                                           0.16156434,  0.34650623,  0.05325734,  0.35445249, -0.11516219,
                                          -0.1832121 , -0.14573532, -0.32891817,  0.32911263, -0.16132915,
                                           0.05237008,  0.07855147,  0.02153306, -0.49731194,  0.48834542,
                                           0.40534158, -0.29236414, -0.20751059,  0.02001015,  0.40191137,
                                           0.48363088, -0.24245794,  0.06435904,  0.30696868, -0.10562995,
                                           0.23107304, -0.33893099,  0.10069857,  0.36586446,  0.48352161,
                                          -0.42063421, -0.07165273, -0.29545714, -0.04936351,  0.04776357,
                                          -0.40667329, -0.20313922,  0.42758424,  0.06900373, -0.042588  ,
                                           0.25352599,  0.24186215, -0.45142097,  0.2086974 ,  0.33924335,
                                          -0.33406212,  0.28099794, -0.21346338, -0.19353025,  0.16526147,
                                          -0.38860783,  0.16487245,  0.38785679,  0.19631127, -0.05967212,
                                          -0.06178562,  0.2650961 ,  0.065642  , -0.41509584,  0.08267109,
                                           0.3148437 , -0.16293362,  0.42757658,  0.250717  ,  0.07406383,
                                           0.25164399, -0.42085104,  0.35938908,  0.32150411,  0.40987166,
                                          -0.3713688 , -0.41821991, -0.36158443, -0.10062129, -0.07569314,
                                           0.06221838, -0.37775645, -0.2986005 ,  0.31164435, -0.03201243,
                                           0.30793821, -0.49257362,  0.05159273,  0.43193215,  0.08217546,
                                          -0.29390427,  0.21775756, -0.12101415,  0.16838395, -0.47068028,
                                           0.13590036, -0.46780207,  0.24478066, -0.027087  , -0.37824564,
                                           0.04263593, -0.43322556,  0.15336487,  0.49608633,  0.26939734,
                                           0.07377411, -0.39736474,  0.19983407,  0.16116787, -0.45090287,
                                           0.2922993 ,  0.01871659, -0.07413231,  0.28818717, -0.08843078,
                                          -0.01897372, -0.31837116, -0.1786811 ,  0.345533  , -0.31309625,
                                          -0.08270894,  0.48903451, -0.26340019,  0.41683233,  0.41839747,
                                          -0.40870366, -0.03634728,  0.00221634, -0.18633105, -0.45266046,
                                           0.24168564,  0.09552964,  0.23824991 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

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

    floatVector answer = { 2.00093808,  1.48052569, -0.5533865 ,
                           0.74087198,  1.25662963,  0.17264556,
                           0.75663016, -0.11948556, -0.07579564 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, hydra.getMicroConfiguration( 0 ) ) );

    floatVector answer2 = { 0.75754206,  0.06435904,  0.30696868,
                           -0.10562995,  1.23107304, -0.33893099,
                            0.10069857,  0.36586446,  1.48352161 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, hydra.getMicroConfiguration( 2 ) ) );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getPreviousSubMicroConfiguration ){

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 1.19646919, -0.21386067, -0.27314855,
                                        0.05131477,  1.21946897, -0.07689354,
                                        0.4807642 ,  0.18482974,  0.9809319  };

    floatVector previousDeformationGradient = { 0.89211752, -0.15682198,  0.22904971,
                                               -0.06142776,  0.5596779 , -0.10195574,
                                                0.23799541, -0.31750827,  0.67545176 };

    floatVector microDeformation = { 1.03155137,  0.03182759,  0.13440096,
                                     0.34943179,  1.22445532,  0.11102351,
                                     0.22244338, -0.17704109,  0.86178866 };

    floatVector previousMicroDeformation = { 0.72826323, -0.20628595,  0.13097612,
                                            -0.40789506,  0.93370117, -0.06913724,
                                            -0.0063149 , -0.07416971,  0.81226122 };

    floatVector gradientMicroDeformation = {-0.07364869,  0.39338916,  0.44416002,  0.00183668,  0.12395295,
                                            -0.3843816 , -0.18271452, -0.08517379,  0.36630916, -0.24954463,
                                            -0.01696574,  0.48555979,  0.01948512,  0.11289453, -0.37937133,
                                             0.3263408 ,  0.10306013,  0.04506801, -0.15723617, -0.19587921,
                                            -0.08297779,  0.18130077,  0.37545684,  0.01042234,  0.16931378,
                                             0.08593655,  0.1249035 };

    floatVector previousGradientMicroDeformation = { 0.17468905,  0.34234244, -0.41680501,  0.26368284, -0.25633363,
                                                    -0.30577704,  0.07245696, -0.40428748,  0.38532683,  0.12724897,
                                                     0.22341636, -0.48387079,  0.09443188,  0.05678519, -0.34104036,
                                                    -0.34692948,  0.19552953, -0.18123357,  0.1919703 ,  0.05438325,
                                                    -0.11104943,  0.42513249,  0.34167   , -0.14260243, -0.45640854,
                                                    -0.19523193, -0.10181432};

    floatVector previousStateVariables = { 0.20495883,  0.49535848, -0.14408513,  0.26254781,  0.09317692,
                                           0.1917018 , -0.34887255, -0.10112371, -0.2591441 , -0.15654399,
                                           0.01312815,  0.16662455, -0.39409151, -0.36910505, -0.17801939,
                                           0.16156434,  0.34650623,  0.05325734,  0.35445249, -0.11516219,
                                          -0.1832121 , -0.14573532, -0.32891817,  0.32911263, -0.16132915,
                                           0.05237008,  0.07855147,  0.02153306, -0.49731194,  0.48834542,
                                           0.40534158, -0.29236414, -0.20751059,  0.02001015,  0.40191137,
                                           0.48363088, -0.24245794,  0.06435904,  0.30696868, -0.10562995,
                                           0.23107304, -0.33893099,  0.10069857,  0.36586446,  0.48352161,
                                          -0.42063421, -0.07165273, -0.29545714, -0.04936351,  0.04776357,
                                          -0.40667329, -0.20313922,  0.42758424,  0.06900373, -0.042588  ,
                                           0.25352599,  0.24186215, -0.45142097,  0.2086974 ,  0.33924335,
                                          -0.33406212,  0.28099794, -0.21346338, -0.19353025,  0.16526147,
                                          -0.38860783,  0.16487245,  0.38785679,  0.19631127, -0.05967212,
                                          -0.06178562,  0.2650961 ,  0.065642  , -0.41509584,  0.08267109,
                                           0.3148437 , -0.16293362,  0.42757658,  0.250717  ,  0.07406383,
                                           0.25164399, -0.42085104,  0.35938908,  0.32150411,  0.40987166,
                                          -0.3713688 , -0.41821991, -0.36158443, -0.10062129, -0.07569314,
                                           0.06221838, -0.37775645, -0.2986005 ,  0.31164435, -0.03201243,
                                           0.30793821, -0.49257362,  0.05159273,  0.43193215,  0.08217546,
                                          -0.29390427,  0.21775756, -0.12101415,  0.16838395, -0.47068028,
                                           0.13590036, -0.46780207,  0.24478066, -0.027087  , -0.37824564,
                                           0.04263593, -0.43322556,  0.15336487,  0.49608633,  0.26939734,
                                           0.07377411, -0.39736474,  0.19983407,  0.16116787, -0.45090287,
                                           0.2922993 ,  0.01871659, -0.07413231,  0.28818717, -0.08843078,
                                          -0.01897372, -0.31837116, -0.1786811 ,  0.345533  , -0.31309625,
                                          -0.08270894,  0.48903451, -0.26340019,  0.41683233,  0.41839747,
                                          -0.40870366, -0.03634728,  0.00221634, -0.18633105, -0.45266046,
                                           0.24168564,  0.09552964,  0.23824991 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

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

    floatVector answer = { 1.79993357, -0.26247057, -0.16201437,
                          -0.7635743 ,  0.6005345 ,  0.34819881,
                           0.16680516, -0.38411455,  0.34087847 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, hydra.getPreviousSubMicroConfiguration( 0, 2 ) ) );

    floatVector answer2 = { 0.87556098, -0.36781419,  1.20660407,
                            0.21141974,  0.82131808, -0.42325899,
                            0.12210416,  1.0388779 ,  2.07092075 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, hydra.getPreviousSubMicroConfiguration( 1, 3 ) ) );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getPreviousPreceedingMicroConfiguration ){

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 1.19646919, -0.21386067, -0.27314855,
                                        0.05131477,  1.21946897, -0.07689354,
                                        0.4807642 ,  0.18482974,  0.9809319  };

    floatVector previousDeformationGradient = { 0.89211752, -0.15682198,  0.22904971,
                                               -0.06142776,  0.5596779 , -0.10195574,
                                                0.23799541, -0.31750827,  0.67545176 };

    floatVector microDeformation = { 1.03155137,  0.03182759,  0.13440096,
                                     0.34943179,  1.22445532,  0.11102351,
                                     0.22244338, -0.17704109,  0.86178866 };

    floatVector previousMicroDeformation = { 0.72826323, -0.20628595,  0.13097612,
                                            -0.40789506,  0.93370117, -0.06913724,
                                            -0.0063149 , -0.07416971,  0.81226122 };

    floatVector gradientMicroDeformation = {-0.07364869,  0.39338916,  0.44416002,  0.00183668,  0.12395295,
                                            -0.3843816 , -0.18271452, -0.08517379,  0.36630916, -0.24954463,
                                            -0.01696574,  0.48555979,  0.01948512,  0.11289453, -0.37937133,
                                             0.3263408 ,  0.10306013,  0.04506801, -0.15723617, -0.19587921,
                                            -0.08297779,  0.18130077,  0.37545684,  0.01042234,  0.16931378,
                                             0.08593655,  0.1249035 };

    floatVector previousGradientMicroDeformation = { 0.17468905,  0.34234244, -0.41680501,  0.26368284, -0.25633363,
                                                    -0.30577704,  0.07245696, -0.40428748,  0.38532683,  0.12724897,
                                                     0.22341636, -0.48387079,  0.09443188,  0.05678519, -0.34104036,
                                                    -0.34692948,  0.19552953, -0.18123357,  0.1919703 ,  0.05438325,
                                                    -0.11104943,  0.42513249,  0.34167   , -0.14260243, -0.45640854,
                                                    -0.19523193, -0.10181432};

    floatVector previousStateVariables = { 0.20495883,  0.49535848, -0.14408513,  0.26254781,  0.09317692,
                                           0.1917018 , -0.34887255, -0.10112371, -0.2591441 , -0.15654399,
                                           0.01312815,  0.16662455, -0.39409151, -0.36910505, -0.17801939,
                                           0.16156434,  0.34650623,  0.05325734,  0.35445249, -0.11516219,
                                          -0.1832121 , -0.14573532, -0.32891817,  0.32911263, -0.16132915,
                                           0.05237008,  0.07855147,  0.02153306, -0.49731194,  0.48834542,
                                           0.40534158, -0.29236414, -0.20751059,  0.02001015,  0.40191137,
                                           0.48363088, -0.24245794,  0.06435904,  0.30696868, -0.10562995,
                                           0.23107304, -0.33893099,  0.10069857,  0.36586446,  0.48352161,
                                          -0.42063421, -0.07165273, -0.29545714, -0.04936351,  0.04776357,
                                          -0.40667329, -0.20313922,  0.42758424,  0.06900373, -0.042588  ,
                                           0.25352599,  0.24186215, -0.45142097,  0.2086974 ,  0.33924335,
                                          -0.33406212,  0.28099794, -0.21346338, -0.19353025,  0.16526147,
                                          -0.38860783,  0.16487245,  0.38785679,  0.19631127, -0.05967212,
                                          -0.06178562,  0.2650961 ,  0.065642  , -0.41509584,  0.08267109,
                                           0.3148437 , -0.16293362,  0.42757658,  0.250717  ,  0.07406383,
                                           0.25164399, -0.42085104,  0.35938908,  0.32150411,  0.40987166,
                                          -0.3713688 , -0.41821991, -0.36158443, -0.10062129, -0.07569314,
                                           0.06221838, -0.37775645, -0.2986005 ,  0.31164435, -0.03201243,
                                           0.30793821, -0.49257362,  0.05159273,  0.43193215,  0.08217546,
                                          -0.29390427,  0.21775756, -0.12101415,  0.16838395, -0.47068028,
                                           0.13590036, -0.46780207,  0.24478066, -0.027087  , -0.37824564,
                                           0.04263593, -0.43322556,  0.15336487,  0.49608633,  0.26939734,
                                           0.07377411, -0.39736474,  0.19983407,  0.16116787, -0.45090287,
                                           0.2922993 ,  0.01871659, -0.07413231,  0.28818717, -0.08843078,
                                          -0.01897372, -0.31837116, -0.1786811 ,  0.345533  , -0.31309625,
                                          -0.08270894,  0.48903451, -0.26340019,  0.41683233,  0.41839747,
                                          -0.40870366, -0.03634728,  0.00221634, -0.18633105, -0.45266046,
                                           0.24168564,  0.09552964,  0.23824991 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

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

    floatVector answer = { 1.79993357, -0.26247057, -0.16201437,
                          -0.7635743 ,  0.6005345 ,  0.34819881,
                           0.16680516, -0.38411455,  0.34087847 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, hydra.getPreviousPrecedingMicroConfiguration( 2 ) ) );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getPreviousFollowingMicroConfiguration ){

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 1.19646919, -0.21386067, -0.27314855,
                                        0.05131477,  1.21946897, -0.07689354,
                                        0.4807642 ,  0.18482974,  0.9809319  };

    floatVector previousDeformationGradient = { 0.89211752, -0.15682198,  0.22904971,
                                               -0.06142776,  0.5596779 , -0.10195574,
                                                0.23799541, -0.31750827,  0.67545176 };

    floatVector microDeformation = { 1.03155137,  0.03182759,  0.13440096,
                                     0.34943179,  1.22445532,  0.11102351,
                                     0.22244338, -0.17704109,  0.86178866 };

    floatVector previousMicroDeformation = { 0.72826323, -0.20628595,  0.13097612,
                                            -0.40789506,  0.93370117, -0.06913724,
                                            -0.0063149 , -0.07416971,  0.81226122 };

    floatVector gradientMicroDeformation = {-0.07364869,  0.39338916,  0.44416002,  0.00183668,  0.12395295,
                                            -0.3843816 , -0.18271452, -0.08517379,  0.36630916, -0.24954463,
                                            -0.01696574,  0.48555979,  0.01948512,  0.11289453, -0.37937133,
                                             0.3263408 ,  0.10306013,  0.04506801, -0.15723617, -0.19587921,
                                            -0.08297779,  0.18130077,  0.37545684,  0.01042234,  0.16931378,
                                             0.08593655,  0.1249035 };

    floatVector previousGradientMicroDeformation = { 0.17468905,  0.34234244, -0.41680501,  0.26368284, -0.25633363,
                                                    -0.30577704,  0.07245696, -0.40428748,  0.38532683,  0.12724897,
                                                     0.22341636, -0.48387079,  0.09443188,  0.05678519, -0.34104036,
                                                    -0.34692948,  0.19552953, -0.18123357,  0.1919703 ,  0.05438325,
                                                    -0.11104943,  0.42513249,  0.34167   , -0.14260243, -0.45640854,
                                                    -0.19523193, -0.10181432};

    floatVector previousStateVariables = { 0.20495883,  0.49535848, -0.14408513,  0.26254781,  0.09317692,
                                           0.1917018 , -0.34887255, -0.10112371, -0.2591441 , -0.15654399,
                                           0.01312815,  0.16662455, -0.39409151, -0.36910505, -0.17801939,
                                           0.16156434,  0.34650623,  0.05325734,  0.35445249, -0.11516219,
                                          -0.1832121 , -0.14573532, -0.32891817,  0.32911263, -0.16132915,
                                           0.05237008,  0.07855147,  0.02153306, -0.49731194,  0.48834542,
                                           0.40534158, -0.29236414, -0.20751059,  0.02001015,  0.40191137,
                                           0.48363088, -0.24245794,  0.06435904,  0.30696868, -0.10562995,
                                           0.23107304, -0.33893099,  0.10069857,  0.36586446,  0.48352161,
                                          -0.42063421, -0.07165273, -0.29545714, -0.04936351,  0.04776357,
                                          -0.40667329, -0.20313922,  0.42758424,  0.06900373, -0.042588  ,
                                           0.25352599,  0.24186215, -0.45142097,  0.2086974 ,  0.33924335,
                                          -0.33406212,  0.28099794, -0.21346338, -0.19353025,  0.16526147,
                                          -0.38860783,  0.16487245,  0.38785679,  0.19631127, -0.05967212,
                                          -0.06178562,  0.2650961 ,  0.065642  , -0.41509584,  0.08267109,
                                           0.3148437 , -0.16293362,  0.42757658,  0.250717  ,  0.07406383,
                                           0.25164399, -0.42085104,  0.35938908,  0.32150411,  0.40987166,
                                          -0.3713688 , -0.41821991, -0.36158443, -0.10062129, -0.07569314,
                                           0.06221838, -0.37775645, -0.2986005 ,  0.31164435, -0.03201243,
                                           0.30793821, -0.49257362,  0.05159273,  0.43193215,  0.08217546,
                                          -0.29390427,  0.21775756, -0.12101415,  0.16838395, -0.47068028,
                                           0.13590036, -0.46780207,  0.24478066, -0.027087  , -0.37824564,
                                           0.04263593, -0.43322556,  0.15336487,  0.49608633,  0.26939734,
                                           0.07377411, -0.39736474,  0.19983407,  0.16116787, -0.45090287,
                                           0.2922993 ,  0.01871659, -0.07413231,  0.28818717, -0.08843078,
                                          -0.01897372, -0.31837116, -0.1786811 ,  0.345533  , -0.31309625,
                                          -0.08270894,  0.48903451, -0.26340019,  0.41683233,  0.41839747,
                                          -0.40870366, -0.03634728,  0.00221634, -0.18633105, -0.45266046,
                                           0.24168564,  0.09552964,  0.23824991 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

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

    floatVector answer = { 0.37335959,  0.14440808,  0.07815635,
                          -0.05311829,  1.15252061, -0.83175389,
                          -0.26108048,  1.01045459,  1.40735072 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, hydra.getPreviousFollowingMicroConfiguration( 1 ) ) );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getPreviousMicroConfiguration ){

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 1.19646919, -0.21386067, -0.27314855,
                                        0.05131477,  1.21946897, -0.07689354,
                                        0.4807642 ,  0.18482974,  0.9809319  };

    floatVector previousDeformationGradient = { 0.89211752, -0.15682198,  0.22904971,
                                               -0.06142776,  0.5596779 , -0.10195574,
                                                0.23799541, -0.31750827,  0.67545176 };

    floatVector microDeformation = { 1.03155137,  0.03182759,  0.13440096,
                                     0.34943179,  1.22445532,  0.11102351,
                                     0.22244338, -0.17704109,  0.86178866 };

    floatVector previousMicroDeformation = { 0.72826323, -0.20628595,  0.13097612,
                                            -0.40789506,  0.93370117, -0.06913724,
                                            -0.0063149 , -0.07416971,  0.81226122 };

    floatVector gradientMicroDeformation = {-0.07364869,  0.39338916,  0.44416002,  0.00183668,  0.12395295,
                                            -0.3843816 , -0.18271452, -0.08517379,  0.36630916, -0.24954463,
                                            -0.01696574,  0.48555979,  0.01948512,  0.11289453, -0.37937133,
                                             0.3263408 ,  0.10306013,  0.04506801, -0.15723617, -0.19587921,
                                            -0.08297779,  0.18130077,  0.37545684,  0.01042234,  0.16931378,
                                             0.08593655,  0.1249035 };

    floatVector previousGradientMicroDeformation = { 0.17468905,  0.34234244, -0.41680501,  0.26368284, -0.25633363,
                                                    -0.30577704,  0.07245696, -0.40428748,  0.38532683,  0.12724897,
                                                     0.22341636, -0.48387079,  0.09443188,  0.05678519, -0.34104036,
                                                    -0.34692948,  0.19552953, -0.18123357,  0.1919703 ,  0.05438325,
                                                    -0.11104943,  0.42513249,  0.34167   , -0.14260243, -0.45640854,
                                                    -0.19523193, -0.10181432};

    floatVector previousStateVariables = { 0.20495883,  0.49535848, -0.14408513,  0.26254781,  0.09317692,
                                           0.1917018 , -0.34887255, -0.10112371, -0.2591441 , -0.15654399,
                                           0.01312815,  0.16662455, -0.39409151, -0.36910505, -0.17801939,
                                           0.16156434,  0.34650623,  0.05325734,  0.35445249, -0.11516219,
                                          -0.1832121 , -0.14573532, -0.32891817,  0.32911263, -0.16132915,
                                           0.05237008,  0.07855147,  0.02153306, -0.49731194,  0.48834542,
                                           0.40534158, -0.29236414, -0.20751059,  0.02001015,  0.40191137,
                                           0.48363088, -0.24245794,  0.06435904,  0.30696868, -0.10562995,
                                           0.23107304, -0.33893099,  0.10069857,  0.36586446,  0.48352161,
                                          -0.42063421, -0.07165273, -0.29545714, -0.04936351,  0.04776357,
                                          -0.40667329, -0.20313922,  0.42758424,  0.06900373, -0.042588  ,
                                           0.25352599,  0.24186215, -0.45142097,  0.2086974 ,  0.33924335,
                                          -0.33406212,  0.28099794, -0.21346338, -0.19353025,  0.16526147,
                                          -0.38860783,  0.16487245,  0.38785679,  0.19631127, -0.05967212,
                                          -0.06178562,  0.2650961 ,  0.065642  , -0.41509584,  0.08267109,
                                           0.3148437 , -0.16293362,  0.42757658,  0.250717  ,  0.07406383,
                                           0.25164399, -0.42085104,  0.35938908,  0.32150411,  0.40987166,
                                          -0.3713688 , -0.41821991, -0.36158443, -0.10062129, -0.07569314,
                                           0.06221838, -0.37775645, -0.2986005 ,  0.31164435, -0.03201243,
                                           0.30793821, -0.49257362,  0.05159273,  0.43193215,  0.08217546,
                                          -0.29390427,  0.21775756, -0.12101415,  0.16838395, -0.47068028,
                                           0.13590036, -0.46780207,  0.24478066, -0.027087  , -0.37824564,
                                           0.04263593, -0.43322556,  0.15336487,  0.49608633,  0.26939734,
                                           0.07377411, -0.39736474,  0.19983407,  0.16116787, -0.45090287,
                                           0.2922993 ,  0.01871659, -0.07413231,  0.28818717, -0.08843078,
                                          -0.01897372, -0.31837116, -0.1786811 ,  0.345533  , -0.31309625,
                                          -0.08270894,  0.48903451, -0.26340019,  0.41683233,  0.41839747,
                                          -0.40870366, -0.03634728,  0.00221634, -0.18633105, -0.45266046,
                                           0.24168564,  0.09552964,  0.23824991 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

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

    floatVector answer = { 1.41998902,  0.88426858, -0.45291886,
                          -0.76867062,  0.02912989,  0.49178024,
                           0.30842234, -0.36954097,  0.0765542 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, hydra.getPreviousMicroConfiguration( 0 ) ) );

    floatVector answer2 = { 0.75754206,  0.06435904,  0.30696868,
                           -0.10562995,  1.23107304, -0.33893099,
                            0.10069857,  0.36586446,  1.48352161 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, hydra.getPreviousMicroConfiguration( 2 ) ) );

}
