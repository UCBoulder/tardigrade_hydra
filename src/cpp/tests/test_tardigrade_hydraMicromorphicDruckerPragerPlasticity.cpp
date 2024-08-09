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

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

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

bool tolerantCheck( const std::vector< double > &v1, const std::vector< double > &v2, double eps = 1e-6, double tol = 1e-9 ){

    if ( v1.size( ) != v2.size( ) ){

        return false;

    }

    BOOST_CHECK( v1.size( ) == v2.size( ) );

    const unsigned int len = v1.size( );

    for ( unsigned int i = 0; i < len; i++ ){

        if ( ( std::fabs( v1[ i ] ) < tol ) || ( std::fabs( v2[ i ] ) < tol ) ){

            if ( std::fabs( v1[ i ] - v2[ i ] ) > eps ){

                return false;

            }

        }
        else{

            if ( ( ( std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v1[ i ] ) ) > eps ) ||
                 ( ( std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v2[ i ] ) ) > eps ) ){

                return false;

            }

        }

    }

    return true;

}

bool tolerantCheck( const double &v1, const double &v2, double eps = 1e-6, double tol = 1e-9 ){

    std::vector< double > _v1 = { v1 };

    std::vector< double > _v2 = { v2 };

    return tolerantCheck( _v1, _v2, eps, tol );

}

BOOST_AUTO_TEST_CASE( test_maxDeltaPlasticDeformation, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    residualMock R( &hydra, 55, 1, stateVariableIndices, parameters );

    BOOST_TEST( 1 == *R.getMaxMacroPlasticDeltaNorm( ) );

    BOOST_TEST( 1 == *R.getMaxMicroPlasticDeltaNorm( ) );

    BOOST_TEST( 1 == *R.getMaxMicroGradientPlasticDeltaNorm( ) );

    R.setMaxMacroPlasticDeltaNorm( 1.1 );

    R.setMaxMicroPlasticDeltaNorm( 2.2 );

    R.setMaxMicroGradientPlasticDeltaNorm( 3.3 );

    BOOST_TEST( 1.1 == *R.getMaxMacroPlasticDeltaNorm( ) );

    BOOST_TEST( 2.2 == *R.getMaxMicroPlasticDeltaNorm( ) );

    BOOST_TEST( 3.3 == *R.getMaxMicroGradientPlasticDeltaNorm( ) );

}

BOOST_AUTO_TEST_CASE( test_weakMac, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of the weakened Macaulay bracket
     */

    floatType answer = 2.6902334101408023;

    BOOST_TEST( tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( 2.5, 0.75 ) == answer );

    floatType dmacdx;

    BOOST_TEST( tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( 2.5, 0.75, dmacdx ) == answer );

    floatType eps = 1e-6;

    floatType delta = eps * 2.5 + eps;

    floatType jac = ( tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( 2.5 + delta, 0.75 ) - tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( 2.5 - delta, 0.75 ) ) / ( 2 * delta );

    BOOST_TEST( jac == dmacdx );

}

BOOST_AUTO_TEST_CASE( test_extractParameters, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
                                     { }, { },
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

    BOOST_TEST( answer1 == *R.get_macroHardeningParameters( )         , CHECK_PER_ELEMENT );

    BOOST_TEST( answer2 == *R.get_microHardeningParameters( )         , CHECK_PER_ELEMENT );

    BOOST_TEST( answer3 == *R.get_microGradientHardeningParameters( ) , CHECK_PER_ELEMENT );

    BOOST_TEST( answer4 == *R.get_macroFlowParameters( )              , CHECK_PER_ELEMENT );

    BOOST_TEST( answer5 == *R.get_microFlowParameters( )              , CHECK_PER_ELEMENT );

    BOOST_TEST( answer6 == *R.get_microGradientFlowParameters( )      , CHECK_PER_ELEMENT );

    BOOST_TEST( answer7 == *R.get_macroYieldParameters( )             , CHECK_PER_ELEMENT );

    BOOST_TEST( answer8 == *R.get_microYieldParameters( )             , CHECK_PER_ELEMENT );

    BOOST_TEST( answer9 == *R.get_microGradientYieldParameters( )     , CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_set_state_variables, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
                                     { }, { },
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

    BOOST_TEST( stateVariablesAnswer == *R.get_plasticStateVariables( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousStateVariablesAnswer == *R.get_previousPlasticStateVariables( ), CHECK_PER_ELEMENT );

    BOOST_TEST( plasticMultipliersAnswer == *R.get_plasticMultipliers( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousPlasticMultipliersAnswer == *R.get_previousPlasticMultipliers( ), CHECK_PER_ELEMENT );

    BOOST_TEST( plasticStrainLikeISVsAnswer == *R.get_plasticStrainLikeISVs( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousPlasticStrainLikeISVsAnswer == *R.get_previousPlasticStrainLikeISVs( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setDrivingStresses, * boost::unit_test::tolerance( 1e-5 ) ){
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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 55, 1, stateVariableIndices, parameters );

    BOOST_TEST( answerMacroStress               == *R.get_macroDrivingStress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answerMicroStress               == *R.get_symmetricMicroDrivingStress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answerHigherOrderStress         == *R.get_higherOrderDrivingStress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousMacroStress       == *R.get_previousMacroDrivingStress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousMicroStress       == *R.get_previousSymmetricMicroDrivingStress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousHigherOrderStress == *R.get_previousHigherOrderDrivingStress( ), CHECK_PER_ELEMENT );

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
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors(  dMacrodF ) == floatVector(  9 * 9, 0 ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors(  dMicrodF ) == floatVector(  9 * 9, 0 ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dHigherdF ) == floatVector( 27 * 9, 0 ), CHECK_PER_ELEMENT );

    // Check the Jacobians of the micro deformation
    for ( unsigned int i = 0; i < microDeformation.size( ); i++ ){

        floatVector delta( microDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation + delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation - delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors(  dMacrodChi ) == floatVector(  9 * 9, 0 ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors(  dMicrodChi ) == floatVector(  9 * 9, 0 ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dHigherdChi ) == floatVector( 27 * 9, 0 ), CHECK_PER_ELEMENT );

    // Check the Jacobians of the unknown vector
    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

        for ( unsigned int v = 0; v < dimension * dimension; v++ ){

            result_dMacrodX[ i ][ col ] = ( *R.get_dMacroDrivingStressdMacroStress( ) )[ dimension * dimension * i + v ];

            col++;

        }

        col = configuration_unknown_count;

        for ( unsigned int v = 0; v < ( numConfigurations - 1 ) * dimension * dimension; v++ ){

            result_dMacrodX[ i ][ col ] = ( *R.get_dMacroDrivingStressdFn( ) )[ ( numConfigurations - 1 ) * dimension * dimension * i + v ];

            col++;

        }

    }

    for ( unsigned int i = 0; i < 9; i++ ){

        unsigned int col = dimension * dimension;

        for ( unsigned int v = 0; v < dimension * dimension; v++ ){

            result_dMicrodX[ i ][ col ] = ( *R.get_dSymmetricMicroDrivingStressdMicroStress( ) )[ dimension * dimension * i + v ];

            col++;

        }

        col = configuration_unknown_count;

        for ( unsigned int v = 0; v < ( numConfigurations - 1 ) * dimension * dimension; v++ ){

            result_dMicrodX[ i ][ col ] = ( *R.get_dSymmetricMicroDrivingStressdFn( ) )[ ( numConfigurations - 1 ) * dimension * dimension * i + v ];

            col++;

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        unsigned int col = 2 * dimension * dimension;

        for ( unsigned int v = 0; v < dimension * dimension * dimension; v++ ){

            result_dHigherdX[ i ][ col ] = ( *R.get_dHigherOrderDrivingStressdHigherOrderStress( ) )[ dimension * dimension * dimension * i + v ];

            col++;

        }

        col = configuration_unknown_count;

        for ( unsigned int v = 0; v < ( numConfigurations - 1 ) * dimension * dimension; v++ ){

            result_dHigherdX[ i ][ col ] = ( *R.get_dHigherOrderDrivingStressdFn( ) )[ ( numConfigurations - 1 ) * dimension * dimension * i + v ];

            col++;

        }

        for ( unsigned int v = 0; v < ( numConfigurations - 1 ) * dimension * dimension; v++ ){

            result_dHigherdX[ i ][ col ] = ( *R.get_dHigherOrderDrivingStressdChin( ) )[ ( numConfigurations - 1 ) * dimension * dimension * i + v ];

            col++;

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors(  dMacrodX ) == tardigradeVectorTools::appendVectors( result_dMacrodX  ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors(  dMicrodX ) == tardigradeVectorTools::appendVectors( result_dMicrodX  ), CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dHigherdX ), tardigradeVectorTools::appendVectors( result_dHigherdX ), 1e-5, 1e-5 ) );

    // Check the Jacobians of the deformation gradient
    floatVector dvec = tardigradeVectorTools::appendVectors( { hydra.elasticity.previousPK2Stress, hydra.elasticity.previousSigma, hydra.elasticity.previousM, previousStateVariables } );

    for ( unsigned int i = 0; i < previousDeformationGradient.size( ); i++ ){

        floatVector delta( previousDeformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors(  previousdMacrodF ) == floatVector(  9 * 9, 0 ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors(  previousdMicrodF ) == floatVector(  9 * 9, 0 ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdHigherdF ) == floatVector( 27 * 9, 0 ), CHECK_PER_ELEMENT );

    // Check the Jacobians of the micro deformation
    for ( unsigned int i = 0; i < previousMicroDeformation.size( ); i++ ){

        floatVector delta( previousMicroDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation + delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation - delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors(  previousdMacrodChi ) == floatVector(  9 * 9, 0 ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors(  previousdMicrodChi ) == floatVector(  9 * 9, 0 ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdHigherdChi ) == floatVector( 27 * 9, 0 ), CHECK_PER_ELEMENT );

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
                                          { }, { },
                                          previousStateVariables + dXi, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

        for ( unsigned int v = 0; v < dimension * dimension; v++ ){

            result_previousdMacrodX[ i ][ col ] = ( *R.get_previousdMacroDrivingStressdMacroStress( ) )[ dimension * dimension * i + v ];

            col++;

        }

        col = configuration_unknown_count;

        for ( unsigned int v = 0; v < ( numConfigurations - 1 ) * dimension * dimension; v++ ){

            result_previousdMacrodX[ i ][ col ] = ( *R.get_previousdMacroDrivingStressdFn( ) )[ ( numConfigurations - 1 ) * dimension * dimension * i + v ];

            col++;

        }

    }

    for ( unsigned int i = 0; i < 9; i++ ){

        unsigned int col = dimension * dimension;

        for ( unsigned int v = 0; v < dimension * dimension; v++ ){

            result_previousdMicrodX[ i ][ col ] = ( *R.get_previousdSymmetricMicroDrivingStressdMicroStress( ) )[ dimension * dimension * i + v ];

            col++;

        }

        col = configuration_unknown_count;

        for ( unsigned int v = 0; v < ( numConfigurations - 1 ) * dimension * dimension; v++ ){

            result_previousdMicrodX[ i ][ col ] = ( *R.get_previousdSymmetricMicroDrivingStressdFn( ) )[ ( numConfigurations - 1 ) * dimension * dimension * i + v ];

            col++;

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        unsigned int col = 2 * dimension * dimension;

        for ( unsigned int v = 0; v < dimension * dimension * dimension; v++ ){

            result_previousdHigherdX[ i ][ col ] = ( *R.get_previousdHigherOrderDrivingStressdHigherOrderStress( ) )[ dimension * dimension * dimension * i + v ];

            col++;

        }

        col = configuration_unknown_count;

        for ( unsigned int v = 0; v < ( numConfigurations - 1 ) * dimension * dimension; v++ ){

            result_previousdHigherdX[ i ][ col ] = ( *R.get_previousdHigherOrderDrivingStressdFn( ) )[ ( numConfigurations - 1 ) * dimension * dimension * i + v ];

            col++;

        }

        for ( unsigned int v = 0; v < ( numConfigurations - 1 ) * dimension * dimension; v++ ){

            result_previousdHigherdX[ i ][ col ] = ( *R.get_previousdHigherOrderDrivingStressdChin( ) )[ ( numConfigurations - 1 ) * dimension * dimension * i + v ];

            col++;

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors(  previousdMacrodX ) == tardigradeVectorTools::appendVectors( result_previousdMacrodX  ), CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors(  previousdMicrodX ), tardigradeVectorTools::appendVectors( result_previousdMicrodX  ), 1e-5, 1e-5 ) );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( previousdHigherdX ), tardigradeVectorTools::appendVectors( result_previousdHigherdX ), 1e-5, 1e-5 ) );

}

BOOST_AUTO_TEST_CASE( test_setCohesion, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
                                     { }, { },
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

    BOOST_TEST( answer1 == *R.get_macroCohesion( )                                     );

    BOOST_TEST( answer2 == *R.get_microCohesion( )                                     );

    BOOST_TEST( answer3 == *R.get_microGradientCohesion( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( answer4 == *R.get_previousMacroCohesion( )                            );

    BOOST_TEST( answer5 == *R.get_previousMicroCohesion( )                            );

    BOOST_TEST( answer6 == *R.get_previousMicroGradientCohesion( ), CHECK_PER_ELEMENT );

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
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < numNonLinearSolveStateVariables; j++ ){

            dMicroGradientCohesiondXAssembled[ i ][ j + numConfigurations * configuration_unknown_count ] = ( *R.get_dMicroGradientCohesiondStateVariables( ) )[ numNonLinearSolveStateVariables * i + j ];

        }

    }

    BOOST_TEST( dMacroCohesiondXAssembled == dMacroCohesiondX, CHECK_PER_ELEMENT );

    BOOST_TEST( dMicroCohesiondXAssembled == dMicroCohesiondX, CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dMicroGradientCohesiondXAssembled ) == tardigradeVectorTools::appendVectors( dMicroGradientCohesiondX ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < previousStateVariables.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    for ( unsigned int i = 0; i < numNonLinearSolveStateVariables; i++ ){

        previousdMacroCohesiondXAssembled[ i + ( numConfigurations - 1 ) * configuration_unknown_count ] = ( *R.get_dMacroCohesiondStateVariables( ) )[ i ];

        previousdMicroCohesiondXAssembled[ i + ( numConfigurations - 1 ) * configuration_unknown_count ] = ( *R.get_dMicroCohesiondStateVariables( ) )[ i ];

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < numNonLinearSolveStateVariables; j++ ){

            previousdMicroGradientCohesiondXAssembled[ i ][ j + ( numConfigurations - 1 ) * configuration_unknown_count ] = ( *R.get_dMicroGradientCohesiondStateVariables( ) )[ numNonLinearSolveStateVariables * i + j ];

        }

    }

    BOOST_TEST( previousdMacroCohesiondXAssembled == previousdMacroCohesiondX, CHECK_PER_ELEMENT );

    BOOST_TEST( previousdMicroCohesiondXAssembled == previousdMicroCohesiondX, CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdMicroGradientCohesiondXAssembled ) == tardigradeVectorTools::appendVectors( previousdMicroGradientCohesiondX ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( testComputeDruckerPragerInternalParameters, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

    BOOST_TEST( Aresult == Aanswer );

    BOOST_TEST( Bresult == Banswer );

}

BOOST_AUTO_TEST_CASE( testComputeSecondOrderDruckerPragerYieldEquation, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

    BOOST_TEST( result == answer );

    //Test the Jacobian
    variableType resultJ, dFdCohesion;
    variableVector dFdS, dFdPrecedingF;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, precedingF,
                                                                                frictionAngle, beta, resultJ,
                                                                                dFdS, dFdCohesion, dFdPrecedingF );

    BOOST_TEST( resultJ == answer );

    //Test the Jacobian
    variableType resultJ2, dFdcohesionJ2;
    variableVector dFdSJ2, dFdPrecedingFJ2;
    variableVector d2FdStress2J2;
    variableVector d2FdStressdElasticRCGJ2;
    variableVector d2FdS2J2, d2FdSdPrecedingFJ2;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion, precedingF,
                                                                                frictionAngle, beta, resultJ2,
                                                                                dFdSJ2, dFdcohesionJ2, dFdPrecedingFJ2, d2FdS2J2,
                                                                                d2FdSdPrecedingFJ2);

    BOOST_TEST( resultJ2 == answer );

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

        BOOST_TEST( gradCol == dFdS[i] );

        BOOST_TEST( gradCol == dFdSJ2[i] );

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

        BOOST_TEST( gradCol == dFdPrecedingF[i] );

        BOOST_TEST( gradCol == dFdPrecedingFJ2[i] );

    }

    //Test dFdcohesion
    for ( unsigned int i = 0; i < 1; i++ ){
        constantType delta = eps * fabs( cohesion ) + eps;

        floatType rp, rm;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion + delta, precedingF,
                                                                                                            frictionAngle, beta, rp );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeSecondOrderDruckerPragerYieldEquation( S, cohesion - delta, precedingF,
                                                                                                                    frictionAngle, beta, rm );

        BOOST_TEST( ( rp - rm ) / ( 2 * delta ) == dFdCohesion );

        BOOST_TEST( ( rp - rm ) / ( 2 * delta ) == dFdcohesionJ2 );

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
            BOOST_TEST( gradCol[j] == d2FdS2J2[9*j+i] );
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
            BOOST_TEST( gradCol[j] == d2FdSdPrecedingFJ2[9*j+i] );
        }
    }

}

BOOST_AUTO_TEST_CASE( testComputeHigherOrderDruckerPragerYieldEquation, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    //Test the Jacobians

    variableVector resultJ;
    variableVector dFdStress, dFdc, dFdPrecedingF;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, precedingF,
                                                                                                        frictionAngle, beta, resultJ,
                                                                                                        dFdStress, dFdc, dFdPrecedingF );

    BOOST_TEST( resultJ == answer, CHECK_PER_ELEMENT );

    variableVector resultJ2;
    variableVector dFdStressJ2, dFdcJ2, dFdPrecedingFJ2;
    variableVector d2FdStress2J2, d2FdStressdPrecedingFJ2;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, precedingF,
                                                                                                        frictionAngle, beta, resultJ2,
                                                                                                        dFdStressJ2, dFdcJ2, dFdPrecedingFJ2,
                                                                                                        d2FdStress2J2, d2FdStressdPrecedingFJ2 );

    BOOST_TEST( resultJ2 == answer, CHECK_PER_ELEMENT );

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
            BOOST_TEST( gradCol[j] == dFdStress[M.size( ) * j + i] );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[j] == dFdStressJ2[M.size( )*j+i] );
        }

        variableVector dFdStressP, dFdcP, dFdPrecedingFP;
        variableVector dFdStressM, dFdcM, dFdPrecedingFM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M + delta, cohesion, precedingF,
                                                                                                            frictionAngle, beta, resultP,
                                                                                                            dFdStressP, dFdcP, dFdPrecedingFP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M - delta, cohesion, precedingF,
                                                                                                            frictionAngle, beta, resultM,
                                                                                                            dFdStressM, dFdcM, dFdPrecedingFM );

        constantVector gradMat = ( dFdStressP - dFdStressM ) / ( 2 * delta[i] );

        unsigned int n, o, p;

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        n = ( int )( i / 9 );
                        o = ( int )( (i - 9 * n ) / 3 );
                        p = ( i - 9 * n - 3 * o ) % 3;
                        BOOST_TEST( gradMat[ 27 * j + 9 * k + 3 * l + m ] == d2FdStress2J2[ 729 * j + 243 * k + 81 * l + 27 * m + 9 * n + 3 * o + p ] );
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
            BOOST_TEST( gradCol[j] == dFdc[cohesion.size( ) * j + i] );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[j] == dFdcJ2[cohesion.size( ) * j + i] );
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
            BOOST_TEST( gradCol[j] == dFdPrecedingF[precedingF.size( ) * j + i] );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[j] == dFdPrecedingF[precedingF.size( ) * j + i] );
        }

        variableVector dFdStressP, dFdcP, dFdPrecedingFP;
        variableVector dFdStressM, dFdcM, dFdPrecedingFM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, precedingF + delta,
                                                                                                            frictionAngle, beta, resultP,
                                                                                                            dFdStressP, dFdcP, dFdPrecedingFP );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::computeHigherOrderDruckerPragerYieldEquation( M, cohesion, precedingF - delta,
                                                                                                            frictionAngle, beta, resultM,
                                                                                                            dFdStressM, dFdcM, dFdPrecedingFM );

        constantVector gradMat = ( dFdStressP - dFdStressM ) / ( 2 * delta[i] );

        unsigned int n, o;

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        n = ( int )( i / 3 );
                        o = ( i % 3 );
                        BOOST_TEST( gradMat[ 27 * j + 9 * k + 3 * l + m ] == d2FdStressdPrecedingFJ2[ 243 * j + 81 * k + 27 * l + 9 * m + 3 * n + o ] );
                    }
                }
            }
        }
    }

}

BOOST_AUTO_TEST_CASE( test_setFlowDerivatives, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

            floatVector dMacroDrivingStressdMacroStress               = initialize(  9 * 9 );

            floatVector dMacroDrivingStressdF                         = initialize(  9 * 9 );

            floatVector dMacroDrivingStressdFn                        = initialize(  9 * 9 );

            floatVector dMicroDrivingStressdMicroStress               = initialize(  9 * 9 );

            floatVector dMicroDrivingStressdF                         = initialize(  9 * 9 );

            floatVector dMicroDrivingStressdFn                        = initialize(  9 * 9 );

            floatVector dMicroGradientDrivingStressdHigherOrderStress = initialize( 27 * 27 );

            floatVector dMicroGradientDrivingStressdF                 = initialize( 27 *  9 );

            floatVector dMicroGradientDrivingStressdFn                = initialize( 27 *  9 );

            floatVector dMicroGradientDrivingStressdChi               = initialize( 27 *  9 );

            floatVector dMicroGradientDrivingStressdChin              = initialize( 27 *  9 );

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 55, 1, stateVariableIndices, parameters );

    residualMock RJ( &hydra, 55, 1, stateVariableIndices, parameters );

    floatType tempYield, dMacroFlowdCohesion, dMicroFlowdCohesion;

    floatVector tempVectorYield, dMacroFlowdDrivingStress, dMicroFlowdDrivingStress, dMacroFlowdPrecedingF, dMicroFlowdPrecedingF;

    floatVector dMicroGradientFlowdDrivingStress, dMicroGradientFlowdCohesion, dMicroGradientFlowdPrecedingF;

    floatType previousdMacroFlowdCohesion, previousdMicroFlowdCohesion;

    floatVector previousdMacroFlowdDrivingStress, previousdMicroFlowdDrivingStress, previousdMacroFlowdPrecedingF, previousdMicroFlowdPrecedingF;

    floatVector previousdMicroGradientFlowdDrivingStress, previousdMicroGradientFlowdCohesion, previousdMicroGradientFlowdPrecedingF;

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

    BOOST_TEST( dMacroFlowdCohesion                      == *R.get_dMacroFlowdc( ) );

    BOOST_TEST( dMicroFlowdCohesion                      == *R.get_dMicroFlowdc( ) );

    BOOST_TEST( dMicroGradientFlowdCohesion              == *R.get_dMicroGradientFlowdc( )                     , CHECK_PER_ELEMENT );

    BOOST_TEST( dMacroFlowdDrivingStress                 == *R.get_dMacroFlowdDrivingStress( )                 , CHECK_PER_ELEMENT );

    BOOST_TEST( dMicroFlowdDrivingStress                 == *R.get_dMicroFlowdDrivingStress( )                 , CHECK_PER_ELEMENT );

    BOOST_TEST( dMicroGradientFlowdDrivingStress         == *R.get_dMicroGradientFlowdDrivingStress( )         , CHECK_PER_ELEMENT );

    BOOST_TEST( previousdMacroFlowdCohesion              == *R.get_previousdMacroFlowdc( ) );

    BOOST_TEST( previousdMicroFlowdCohesion              == *R.get_previousdMicroFlowdc( ) );

    BOOST_TEST( previousdMicroGradientFlowdCohesion      == *R.get_previousdMicroGradientFlowdc( )             , CHECK_PER_ELEMENT );

    BOOST_TEST( previousdMacroFlowdDrivingStress         == *R.get_previousdMacroFlowdDrivingStress( )         , CHECK_PER_ELEMENT );

    BOOST_TEST( previousdMicroFlowdDrivingStress         == *R.get_previousdMicroFlowdDrivingStress( )         , CHECK_PER_ELEMENT );

    BOOST_TEST( previousdMicroGradientFlowdDrivingStress == *R.get_previousdMicroGradientFlowdDrivingStress( ) , CHECK_PER_ELEMENT );

    BOOST_TEST( dMacroFlowdCohesion                      == *RJ.get_dMacroFlowdc( ) );

    BOOST_TEST( dMicroFlowdCohesion                      == *RJ.get_dMicroFlowdc( ) );

    BOOST_TEST( dMicroGradientFlowdCohesion              == *RJ.get_dMicroGradientFlowdc( )                    , CHECK_PER_ELEMENT );

    BOOST_TEST( dMacroFlowdDrivingStress                 == *RJ.get_dMacroFlowdDrivingStress( )                , CHECK_PER_ELEMENT );

    BOOST_TEST( dMicroFlowdDrivingStress                 == *RJ.get_dMicroFlowdDrivingStress( )                , CHECK_PER_ELEMENT );

    BOOST_TEST( dMicroGradientFlowdDrivingStress         == *RJ.get_dMicroGradientFlowdDrivingStress( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( previousdMacroFlowdCohesion              == *RJ.get_previousdMacroFlowdc( ) );

    BOOST_TEST( previousdMicroFlowdCohesion              == *RJ.get_previousdMicroFlowdc( ) );

    BOOST_TEST( previousdMicroGradientFlowdCohesion      == *RJ.get_previousdMicroGradientFlowdc( )            , CHECK_PER_ELEMENT );

    BOOST_TEST( previousdMacroFlowdDrivingStress         == *RJ.get_previousdMacroFlowdDrivingStress( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( previousdMicroFlowdDrivingStress         == *RJ.get_previousdMicroFlowdDrivingStress( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( previousdMicroGradientFlowdDrivingStress == *RJ.get_previousdMicroGradientFlowdDrivingStress( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setFlowDerivatives2, * boost::unit_test::tolerance( 1e-5 ) ){
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
                                     { }, { },
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
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

        vp = *Rp.get_dMicroGradientFlowdDrivingStress( );

        vm = *Rm.get_dMicroGradientFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            d2MicroGradientFlowdDrivingStressdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_d2MacroFlowdDrivingStressdX(              9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix assembled_d2MicroFlowdDrivingStressdX(              9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix assembled_d2MicroGradientFlowdDrivingStressdX( 3 * 27, floatVector( unknownVector.size( ), 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_d2MacroFlowdDrivingStressdX[ i ][ j ] = ( *R.get_d2MacroFlowdDrivingStressdStress( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_d2MacroFlowdDrivingStressdX[ i ][ j + configuration_unknown_count ] = ( *R.get_d2MacroFlowdDrivingStressdFn( ) )[ 9 * i + j ];

        }

    }

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_d2MicroFlowdDrivingStressdX[ i ][ j + 9 ] = ( *R.get_d2MicroFlowdDrivingStressdStress( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_d2MicroFlowdDrivingStressdX[ i ][ j + configuration_unknown_count ] = ( *R.get_d2MicroFlowdDrivingStressdFn( ) )[ 9 * i + j ];

        }

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 27; j++ ){

            for ( unsigned int k = 0; k < 27; k++ ){

                assembled_d2MicroGradientFlowdDrivingStressdX[ 27 * i + j ][ k + 18 ] = ( *R.get_d2MicroGradientFlowdDrivingStressdStress( ) )[ 27 * 27 * i + 27 * j + k ];

            }

            for ( unsigned int k = 0; k < 9; k++ ){

                assembled_d2MicroGradientFlowdDrivingStressdX[ 27 * i + j ][ k + configuration_unknown_count     ] = ( *R.get_d2MicroGradientFlowdDrivingStressdFn( ) )[ 27 * 9 * i + 9 * j + k ];

                assembled_d2MicroGradientFlowdDrivingStressdX[ 27 * i + j ][ k + configuration_unknown_count + 9 ] = ( *R.get_d2MicroGradientFlowdDrivingStressdChin( ) )[ 27 * 9 * i + 9 * j + k ];

            }

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( assembled_d2MacroFlowdDrivingStressdX ) == tardigradeVectorTools::appendVectors( d2MacroFlowdDrivingStressdX ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( assembled_d2MicroFlowdDrivingStressdX ) == tardigradeVectorTools::appendVectors( d2MicroFlowdDrivingStressdX ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( assembled_d2MicroGradientFlowdDrivingStressdX ) == tardigradeVectorTools::appendVectors( d2MicroGradientFlowdDrivingStressdX ), CHECK_PER_ELEMENT );

    // derivatives w.r.t. current F
    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

        vp = *Rp.get_dMicroGradientFlowdDrivingStress( );

        vm = *Rm.get_dMicroGradientFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            d2MicroGradientFlowdDrivingStressdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_d2MicroGradientFlowdDrivingStressdF( 3 * 27, floatVector( 9, 0 ) );

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 27; j++ ){

            for ( unsigned int k = 0; k < 9; k++ ){

                assembled_d2MicroGradientFlowdDrivingStressdF[ 27 * i + j ][ k ] = ( *R.get_d2MicroGradientFlowdDrivingStressdF( ) )[ 27 * 9 * i + 9 * j + k ];

            }

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( d2MacroFlowdDrivingStressdF ) == *R.get_d2MacroFlowdDrivingStressdF( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( d2MicroFlowdDrivingStressdF ) == *R.get_d2MicroFlowdDrivingStressdF( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( d2MicroGradientFlowdDrivingStressdF ) == tardigradeVectorTools::appendVectors(assembled_d2MicroGradientFlowdDrivingStressdF ), CHECK_PER_ELEMENT );

    // derivatives w.r.t. current chi
    for ( unsigned int i = 0; i < microDeformation.size( ); i++ ){

        floatVector delta( microDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation + delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation - delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

        vp = *Rp.get_dMicroGradientFlowdDrivingStress( );

        vm = *Rm.get_dMicroGradientFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            d2MicroGradientFlowdDrivingStressdChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_d2MicroGradientFlowdDrivingStressdChi( 3 * 27, floatVector( 9, 0 ) );

    BOOST_TEST( tardigradeVectorTools::appendVectors( d2MacroFlowdDrivingStressdChi ) == floatVector( 9 * 9, 0 ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( d2MicroFlowdDrivingStressdChi ) == floatVector( 9 * 9, 0 ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( d2MicroGradientFlowdDrivingStressdChi ) == floatVector( 3 * 27 * 9, 0 ), CHECK_PER_ELEMENT );

    // derivatives w.r.t. previous macro stress
    for ( unsigned int i = 0; i < hydra.elasticity.previousPK2Stress.size( ); i++ ){

        floatVector delta( hydra.elasticity.previousPK2Stress.size( ), 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousPK2Stress[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            BOOST_TEST( 0. == ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) );

        }

        vp = *Rp.get_previousdMicroGradientFlowdDrivingStress( );

        vm = *Rm.get_previousdMicroGradientFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_TEST( 0. == ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousd2MacroFlowdDrivingStressdMacroStress ) == *R.get_previousd2MacroFlowdDrivingStressdStress( ), CHECK_PER_ELEMENT );

    // derivatives w.r.t. previous micro stress
    for ( unsigned int i = 0; i < hydra.elasticity.previousSIGMA.size( ); i++ ){

        floatVector delta( hydra.elasticity.previousSIGMA.size( ), 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousSIGMA[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            BOOST_TEST( 0. == ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) );

        }

        vp = *Rp.get_previousdMicroFlowdDrivingStress( );

        vm = *Rm.get_previousdMicroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousd2MicroFlowdDrivingStressdMicroStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_previousdMicroGradientFlowdDrivingStress( );

        vm = *Rm.get_previousdMicroGradientFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_TEST( 0. == ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousd2MicroFlowdDrivingStressdMicroStress ) == *R.get_previousd2MicroFlowdDrivingStressdStress( ), CHECK_PER_ELEMENT );

    // derivatives w.r.t. previous micro gradient stress
    for ( unsigned int i = 0; i < hydra.elasticity.previousM.size( ); i++ ){

        floatVector delta( hydra.elasticity.previousM.size( ), 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousM[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            BOOST_TEST( 0. == ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) );

        }

        vp = *Rp.get_previousdMicroFlowdDrivingStress( );

        vm = *Rm.get_previousdMicroFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_TEST( 0. == ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) );

        }

        vp = *Rp.get_previousdMicroGradientFlowdDrivingStress( );

        vm = *Rm.get_previousdMicroGradientFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousd2MicroGradientFlowdDrivingStressdHigherOrderStress[ j ][ i ] =  ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_previousd2MicroGradientFlowdDrivingStressdStress( 3 * 27, floatVector( 27, 0 ) );

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 27; j++ ){

            for ( unsigned int k = 0; k < 27; k++ ){

                assembled_previousd2MicroGradientFlowdDrivingStressdStress[ 27 * i + j ][ k ] += ( *R.get_previousd2MicroGradientFlowdDrivingStressdStress( ) )[ 27 * 27 * i + 27 * j + k ];

            }

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( previousd2MicroGradientFlowdDrivingStressdHigherOrderStress ), tardigradeVectorTools::appendVectors( assembled_previousd2MicroGradientFlowdDrivingStressdStress ), 1e-6, 1e-5 ) );

    // derivatives w.r.t. previous, Fn, Chin, gradchin, and ISVs
    for ( unsigned int i = 0; i < previousStateVariables.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

        vp = *Rp.get_previousdMicroGradientFlowdDrivingStress( );

        vm = *Rm.get_previousdMicroGradientFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousd2MicroGradientFlowdDrivingStressdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_previousd2MacroFlowdDrivingStressdX(              9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix assembled_previousd2MicroFlowdDrivingStressdX(              9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix assembled_previousd2MicroGradientFlowdDrivingStressdX( 3 * 27, floatVector( previousStateVariables.size( ), 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_previousd2MacroFlowdDrivingStressdX[ i ][ j ] = ( *R.get_previousd2MacroFlowdDrivingStressdFn( ) )[ 9 * i + j ];

        }

    }

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_previousd2MicroFlowdDrivingStressdX[ i ][ j ] = ( *R.get_previousd2MicroFlowdDrivingStressdFn( ) )[ 9 * i + j ];

        }

    }

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 27; j++ ){

            for ( unsigned int k = 0; k < 9; k++ ){

                assembled_previousd2MicroGradientFlowdDrivingStressdX[ 27 * i + j ][ k     ] = ( *R.get_previousd2MicroGradientFlowdDrivingStressdFn( ) )[ 27 * 9 * i + 9 * j + k ];

                assembled_previousd2MicroGradientFlowdDrivingStressdX[ 27 * i + j ][ k + 9 ] = ( *R.get_previousd2MicroGradientFlowdDrivingStressdChin( ) )[ 27 * 9 * i + 9 * j + k ];

            }

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( assembled_previousd2MacroFlowdDrivingStressdX ) == tardigradeVectorTools::appendVectors( previousd2MacroFlowdDrivingStressdX ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( assembled_previousd2MicroFlowdDrivingStressdX ) == tardigradeVectorTools::appendVectors( previousd2MicroFlowdDrivingStressdX ), CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( assembled_previousd2MicroGradientFlowdDrivingStressdX ), tardigradeVectorTools::appendVectors( previousd2MicroGradientFlowdDrivingStressdX ), 1e-5, 1e-5 ) );

    // derivatives w.r.t. previous F
    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

        vp = *Rp.get_previousdMicroGradientFlowdDrivingStress( );

        vm = *Rm.get_previousdMicroGradientFlowdDrivingStress( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousd2MicroGradientFlowdDrivingStressdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_previousd2MicroGradientFlowdDrivingStressdF( 3 * 27, floatVector( 9, 0 ) );

    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 27; j++ ){

            for ( unsigned int k = 0; k < 9; k++ ){

                assembled_previousd2MicroGradientFlowdDrivingStressdF[ 27 * i + j ][ k ] = ( *R.get_previousd2MicroGradientFlowdDrivingStressdF( ) )[ 27 * 9 * i + 9 * j + k ];

            }

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousd2MacroFlowdDrivingStressdF ) == *R.get_previousd2MacroFlowdDrivingStressdF( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousd2MicroFlowdDrivingStressdF ) == *R.get_previousd2MicroFlowdDrivingStressdF( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousd2MicroGradientFlowdDrivingStressdF ) == tardigradeVectorTools::appendVectors( assembled_previousd2MicroGradientFlowdDrivingStressdF ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setPlasticStrainLikeISVEvolutionRates, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

            floatVector dMicroGradientFlowdc = { 0.34, 0.45, 0.56,
                                                 0.67, 0.78, 0.89,
                                                 0.91, 1.01, 1.12 };

            floatType previousdMacroFlowdc = -1.03;

            floatType previousdMicroFlowdc = -2.04;

            floatVector previousdMicroGradientFlowdc = { -3.4,  -4.5,  -5.6,
                                                         -6.7,  -7.8,  -8.9,
                                                         -9.1, -10.1, -11.2 };

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
                                     { }, { },
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

    BOOST_TEST( evolutionRates == *R.get_plasticStrainLikeISVEvolutionRates( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousEvolutionRates == *R.get_previousPlasticStrainLikeISVEvolutionRates( ), CHECK_PER_ELEMENT );

    BOOST_TEST( evolutionRates == *RJ.get_plasticStrainLikeISVEvolutionRates( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousEvolutionRates == *RJ.get_previousPlasticStrainLikeISVEvolutionRates( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setPlasticStrainLikeISVEvolutionRates2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
                                     { }, { },
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
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            assembled_dEvolutionRatesdX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dPlasticStrainLikeISVEvolutionRatesdStateVariables( ) )[ 10 * i + j ];

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dEvolutionRatesdX ) == tardigradeVectorTools::appendVectors( assembled_dEvolutionRatesdX ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < previousStateVariables.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            assembled_previousdEvolutionRatesdX[ i ][ j + configuration_unknown_count ] = ( *R.get_previousdPlasticStrainLikeISVEvolutionRatesdStateVariables( ) )[ 10 * i + j ];

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdEvolutionRatesdX ) == tardigradeVectorTools::appendVectors( assembled_previousdEvolutionRatesdX ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setUpdatedPlasticStrainLikeISVs, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

            floatVector dPlasticStrainLikeISVEvolutionRatesdStateVariables         = initialize( 5 * 10 );

            floatVector previousdPlasticStrainLikeISVEvolutionRatesdStateVariables = initialize( 5 * 10 );

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatType alpha = 0.67;

    residualMock R( &hydra, 55, 1, stateVariableIndices, parameters, alpha );

    residualMock RJ( &hydra, 55, 1, stateVariableIndices, parameters, alpha );

    RJ.get_dUpdatedPlasticStrainLikeISVsdStateVariables( );

    floatVector updatedPlasticISVs = ( alpha * ( *R.get_plasticStrainLikeISVEvolutionRates( ) ) + ( 1 - alpha ) * ( *R.get_previousPlasticStrainLikeISVEvolutionRates( ) ) ) * deltaTime
                                   + floatVector( previousStateVariables.end( ) - 5, previousStateVariables.end( ) );

    BOOST_TEST( updatedPlasticISVs == *R.get_updatedPlasticStrainLikeISVs( ), CHECK_PER_ELEMENT );

    BOOST_TEST( updatedPlasticISVs == *RJ.get_updatedPlasticStrainLikeISVs( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setUpdatedPlasticStrainLikeISVs2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

            floatVector dMacroDrivingStressdMacroStress               = initialize(  9 * 9 );

            floatVector dMacroDrivingStressdF                         = initialize(  9 * 9 );

            floatVector dMacroDrivingStressdFn                        = initialize(  9 * 9 );

            floatVector dMicroDrivingStressdMicroStress               = initialize(  9 * 9 );

            floatVector dMicroDrivingStressdF                         = initialize(  9 * 9 );

            floatVector dMicroDrivingStressdFn                        = initialize(  9 * 9 );

            floatVector dMicroGradientDrivingStressdHigherOrderStress = initialize( 27 * 27 );

            floatVector dMicroGradientDrivingStressdF                 = initialize( 27 *  9 );

            floatVector dMicroGradientDrivingStressdFn                = initialize( 27 *  9 );

            floatVector dMicroGradientDrivingStressdChi               = initialize( 27 *  9 );

            floatVector dMicroGradientDrivingStressdChin              = initialize( 27 *  9 );

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
                                     { }, { },
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
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            assembled_dUpdatedPlasticStrainLikeISVsdX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dUpdatedPlasticStrainLikeISVsdStateVariables( ) )[ 10 * i + j ];

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticStrainLikeISVsdX ) == tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticStrainLikeISVsdX ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < previousStateVariables.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            assembled_dUpdatedPlasticStrainLikeISVsdPreviousISVs[ i ][ j + configuration_unknown_count ] = ( *R.get_dUpdatedPlasticStrainLikeISVsdPreviousStateVariables( ) )[ 10 * i + j ];

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticStrainLikeISVsdPreviousISVs ) == tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticStrainLikeISVsdPreviousISVs ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setYield, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

            floatVector dMacroCohesiondStateVariables                 = initialize( 10 );

            floatVector dMicroCohesiondStateVariables                 = initialize( 10 );

            floatVector dMicroGradientCohesiondStateVariables         = initialize(  3 * 10 );

            floatVector dMacroDrivingStressdMacroStress               = initialize(  9 * 9 );

            floatVector dMacroDrivingStressdF                         = initialize(  9 * 9 );

            floatVector dMacroDrivingStressdFn                        = initialize(  9 * 9 );

            floatVector dMicroDrivingStressdMicroStress               = initialize(  9 * 9 );

            floatVector dMicroDrivingStressdF                         = initialize(  9 * 9 );

            floatVector dMicroDrivingStressdFn                        = initialize(  9 * 9 );

            floatVector dMicroGradientDrivingStressdHigherOrderStress = initialize( 27 * 27 );

            floatVector dMicroGradientDrivingStressdF                 = initialize( 27 *  9 );

            floatVector dMicroGradientDrivingStressdFn                = initialize( 27 *  9 );

            floatVector dMicroGradientDrivingStressdChi               = initialize( 27 *  9 );

            floatVector dMicroGradientDrivingStressdChin              = initialize( 27 *  9 );

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
                                     { }, { },
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

    BOOST_TEST( macroYield                 == *R.get_macroYield( ) );

    BOOST_TEST( microYield                 == *R.get_microYield( ) );

    BOOST_TEST( microGradientYield         == *R.get_microGradientYield( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousMacroYield         == *R.get_previousMacroYield( ) );

    BOOST_TEST( previousMicroYield         == *R.get_previousMicroYield( ) );

    BOOST_TEST( previousMicroGradientYield == *R.get_previousMicroGradientYield( ), CHECK_PER_ELEMENT );

    BOOST_TEST( macroYield                 == *RJ.get_macroYield( ) );

    BOOST_TEST( microYield                 == *RJ.get_microYield( ) );

    BOOST_TEST( microGradientYield         == *RJ.get_microGradientYield( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousMacroYield         == *RJ.get_previousMacroYield( ) );

    BOOST_TEST( previousMicroYield         == *RJ.get_previousMicroYield( ) );

    BOOST_TEST( previousMicroGradientYield == *RJ.get_previousMicroGradientYield( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setYield2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
                                     { }, { },
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
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            assemble_dMicroGradientYielddX[ i ][ j + 18 ] = ( *R.get_dMicroGradientYielddStress( ) )[ 27 * i + j ];

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

            assemble_dMicroGradientYielddX[ i ][ j + configuration_unknown_count ] = ( *R.get_dMicroGradientYielddFn( ) )[ 9 * i + j ];

        }

    }

    // Chin derivatives
    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assemble_dMicroGradientYielddX[ i ][ j + configuration_unknown_count + 9 ] = ( *R.get_dMicroGradientYielddChin( ) )[ 9 * i + j ];

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

            assemble_dMicroGradientYielddX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dMicroGradientYielddStateVariables( ) )[ 10 * i + j ];

        }

    }

    BOOST_TEST( assemble_dMacroYielddX == dMacroYielddX, CHECK_PER_ELEMENT );

    BOOST_TEST( assemble_dMicroYielddX == dMicroYielddX, CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( assemble_dMicroGradientYielddX ), tardigradeVectorTools::appendVectors( dMicroGradientYielddX ), 1e-5, 1e-6 ) );

    // Jacobians w.r.t. F
    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( dMacroYielddF == *R.get_dMacroYielddF( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dMicroYielddF == *R.get_dMicroYielddF( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dMicroGradientYielddF ) == *R.get_dMicroGradientYielddF( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. Chi
    for ( unsigned int i = 0; i < microDeformation.size( ); i++ ){

        floatVector delta( microDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation + delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation - delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

        BOOST_TEST( 0. == ( sp - sm ) / ( 2 * delta[ i ] ) );

        sp = *Rp.get_microYield( );

        sm = *Rm.get_microYield( );

        BOOST_TEST( 0. == ( sp - sm ) / ( 2 * delta[ i ] ) );

        floatVector vp = *Rp.get_microGradientYield( );

        floatVector vm = *Rm.get_microGradientYield( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicroGradientYielddChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dMicroGradientYielddChi ) == *R.get_dMicroGradientYielddChi( ), CHECK_PER_ELEMENT );

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
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            assemble_previousdMicroGradientYielddStress[ i ][ j + 18 ] = ( *R.get_previousdMicroGradientYielddStress( ) )[ 27 * i + j ];

        }

    }

    BOOST_TEST( assemble_previousdMacroYielddStress == previousdMacroYielddStress, CHECK_PER_ELEMENT );

    BOOST_TEST( assemble_previousdMicroYielddStress == previousdMicroYielddStress, CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( assemble_previousdMicroGradientYielddStress ) == tardigradeVectorTools::appendVectors( previousdMicroGradientYielddStress ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the previous Fn, Chin, and state variables
    for ( unsigned int i = 0; i < previousStateVariables.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            assemble_previousdMicroGradientYielddStateVariables[ i ][ j ] = ( *R.get_previousdMicroGradientYielddFn( ) )[ 9 * i + j ];

        }

    }

    // Chin derivatives
    for ( unsigned int i = 0; i < 3; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assemble_previousdMicroGradientYielddStateVariables[ i ][ j + 9 ] = ( *R.get_previousdMicroGradientYielddChin( ) )[ 9 * i + j ];

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

            assemble_previousdMicroGradientYielddStateVariables[ i ][ j + configuration_unknown_count ] = ( *R.get_previousdMicroGradientYielddStateVariables( ) )[ 10 * i + j ];

        }

    }

    BOOST_TEST( assemble_previousdMacroYielddStateVariables == previousdMacroYielddStateVariables, CHECK_PER_ELEMENT );

    BOOST_TEST( assemble_previousdMicroYielddStateVariables == previousdMicroYielddStateVariables, CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( assemble_previousdMicroGradientYielddStateVariables ) == tardigradeVectorTools::appendVectors( previousdMicroGradientYielddStateVariables ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. previous F
    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( previousdMacroYielddF == *R.get_previousdMacroYielddF( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousdMicroYielddF == *R.get_previousdMicroYielddF( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdMicroGradientYielddF ) == *R.get_previousdMicroGradientYielddF( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. Chi
    for ( unsigned int i = 0; i < microDeformation.size( ); i++ ){

        floatVector delta( microDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation + delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation - delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

        BOOST_TEST( 0. == ( sp - sm ) / ( 2 * delta[ i ] ) );

        sp = *Rp.get_microYield( );

        sm = *Rm.get_microYield( );

        BOOST_TEST( 0. == ( sp - sm ) / ( 2 * delta[ i ] ) );

        floatVector vp = *Rp.get_previousMicroGradientYield( );

        floatVector vm = *Rm.get_previousMicroGradientYield( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroGradientYielddChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdMicroGradientYielddChi ) == *R.get_previousdMicroGradientYielddChi( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setPrecedingDeformationGradient, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
                                     { }, { },
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

    BOOST_TEST( answerFe         == *R.get_precedingDeformationGradient( ),          CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousFe == *R.get_previousPrecedingDeformationGradient( ),  CHECK_PER_ELEMENT );

    BOOST_TEST( answerFe         == *RJ.get_precedingDeformationGradient( ),         CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousFe == *RJ.get_previousPrecedingDeformationGradient( ), CHECK_PER_ELEMENT );

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
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPrecedingFdF ) == *R.get_dPrecedingDeformationGradientdF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + configuration_unknown_count ] = eps * std::fabs( unknownVector[ i + configuration_unknown_count ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPrecedingFdFn ), *R.get_dPrecedingDeformationGradientdFn( ), 1e-5, 1e-6 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdPrecedingFdF ) == *R.get_previousdPrecedingDeformationGradientdF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( previousdPrecedingFdFn ), *R.get_previousdPrecedingDeformationGradientdFn( ), 1e-5, 1e-6 ) );

}

BOOST_AUTO_TEST_CASE( test_setPrecedingMicroDeformation, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
                                     { }, { },
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

    BOOST_TEST( answerChie         == *R.get_precedingMicroDeformation( )         , CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousChie == *R.get_previousPrecedingMicroDeformation( ) , CHECK_PER_ELEMENT );

    BOOST_TEST( answerChie         == *RJ.get_precedingMicroDeformation( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousChie == *RJ.get_previousPrecedingMicroDeformation( ), CHECK_PER_ELEMENT );

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
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation - delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPrecedingChidChi ) == *R.get_dPrecedingMicroDeformationdChi( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + configuration_unknown_count + 9 ] = eps * std::fabs( unknownVector[ i + configuration_unknown_count + 9 ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPrecedingChidChin ), *R.get_dPrecedingMicroDeformationdChin( ), 1e-5, 1e-6 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation + delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation - delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdPrecedingChidChi ) == *R.get_previousdPrecedingMicroDeformationdChi( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( previousStateVariables[ i + 9 ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( previousdPrecedingChidChin ), *R.get_previousdPrecedingMicroDeformationdChin( ), 1e-5, 1e-6 ) );

}

BOOST_AUTO_TEST_CASE( test_setPrecedingGradientMicroDeformation, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    residualMock RJ( &hydra, 55, 1, stateVariableIndices, parameters );

    floatVector answerGradChie         = floatVector( hydra.get_gradientMicroConfigurations( )->begin( ),
                                                      hydra.get_gradientMicroConfigurations( )->begin( ) + 27 );

    floatVector answerPreviousGradChie = floatVector( hydra.get_previousGradientMicroConfigurations( )->begin( ),
                                                      hydra.get_previousGradientMicroConfigurations( )->begin( ) + 27 );

    RJ.get_dPrecedingMicroDeformationdChi( );
    RJ.get_previousdPrecedingMicroDeformationdChi( );

    BOOST_TEST( answerGradChie         == *R.get_precedingGradientMicroDeformation( )         , CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousGradChie == *R.get_previousPrecedingGradientMicroDeformation( ) , CHECK_PER_ELEMENT );

    BOOST_TEST( answerGradChie         == *RJ.get_precedingGradientMicroDeformation( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousGradChie == *RJ.get_previousPrecedingGradientMicroDeformation( ), CHECK_PER_ELEMENT );

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
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPrecedingGradChidFn ) == *R.get_dPrecedingGradientMicroDeformationdFn( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation + delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation - delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPrecedingGradChidChi ) == *R.get_dPrecedingGradientMicroDeformationdChi( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + configuration_unknown_count + 9 ] = eps * std::fabs( unknownVector[ i + configuration_unknown_count + 9 ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPrecedingGradChidChin ) == *R.get_dPrecedingGradientMicroDeformationdChin( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 27; i++ ){

        floatVector delta( 27, 0 );

        delta[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation + delta, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation - delta, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPrecedingGradChidGradChi ) == *R.get_dPrecedingGradientMicroDeformationdGradChi( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 27; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + configuration_unknown_count + 18 ] = eps * std::fabs( unknownVector[ i + configuration_unknown_count + 18 ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPrecedingGradChidGradChin ), *R.get_dPrecedingGradientMicroDeformationdGradChin( ), 1e-5, 1e-6 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdPrecedingGradChidFn ) == *R.get_previousdPrecedingGradientMicroDeformationdFn( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation + delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation - delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdPrecedingGradChidChi ) == *R.get_previousdPrecedingGradientMicroDeformationdChi( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( previousStateVariables[ i + 9 ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdPrecedingGradChidChin ) == *R.get_previousdPrecedingGradientMicroDeformationdChin( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 27; i++ ){

        floatVector delta( 27, 0 );

        delta[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation + delta,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation - delta,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdPrecedingGradChidGradChi ) == *R.get_previousdPrecedingGradientMicroDeformationdGradChi( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 27; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + 18 ] = eps * std::fabs( previousStateVariables[ i + 18 ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( previousdPrecedingGradChidGradChin ), *R.get_previousdPrecedingGradientMicroDeformationdGradChin( ), 1e-5, 1e-6 ) );

}

BOOST_AUTO_TEST_CASE( testComputePlasticMacroVelocityGradient, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

    BOOST_TEST( answerMacroLp == resultMacroLp, CHECK_PER_ELEMENT );

    //Tests of the Jacobians
    variableVector resultMacroLpJ;
    variableVector dMacroLdMacroGammaJ, dMacroLdMicroGammaJ;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseCe,
                                                                                               macroFlowDirection, microFlowDirection,
                                                                                               resultMacroLpJ, dMacroLdMacroGammaJ,
                                                                                               dMacroLdMicroGammaJ );

    BOOST_TEST( answerMacroLp == resultMacroLpJ, CHECK_PER_ELEMENT );

    variableVector resultMacroLpJ2;
    variableVector dMacroLdMacroGammaJ2, dMacroLdMicroGammaJ2;
    variableVector dMacroLdElasticRCG, dMacroLdMacroFlowDirection, dMacroLdMicroFlowDirection;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseCe,
                                                                       macroFlowDirection, microFlowDirection,
                                                                       resultMacroLpJ2, dMacroLdMacroGammaJ2,
                                                                       dMacroLdMicroGammaJ2, dMacroLdElasticRCG,
                                                                       dMacroLdMacroFlowDirection,
                                                                       dMacroLdMicroFlowDirection );

    BOOST_TEST( answerMacroLp == resultMacroLpJ2, CHECK_PER_ELEMENT );

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

    BOOST_TEST( gradCol == dMacroLdMacroGammaJ, CHECK_PER_ELEMENT );

    BOOST_TEST( gradCol == dMacroLdMacroGammaJ2, CHECK_PER_ELEMENT );

    //Test Jacobians w.r.t. microGamma
    scalarDelta = eps * fabs( microGamma) + eps;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma + scalarDelta,
                                                                       inverseCe, macroFlowDirection, microFlowDirection,
                                                                       resultMacroLpP );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMacroVelocityGradient( macroGamma, microGamma - scalarDelta,
                                                                       inverseCe, macroFlowDirection, microFlowDirection,
                                                                       resultMacroLpM );

    gradCol = ( resultMacroLpP - resultMacroLpM ) / ( 2 * scalarDelta );

    BOOST_TEST( gradCol == dMacroLdMicroGammaJ, CHECK_PER_ELEMENT );

    BOOST_TEST( gradCol == dMacroLdMicroGammaJ2, CHECK_PER_ELEMENT );

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

            BOOST_TEST( gradCol[ j ] == dMacroLdElasticRCG[ 9 * j + i ] );

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

            BOOST_TEST( gradCol[ j ] == dMacroLdMacroFlowDirection[ 9 * j + i ] );

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

            BOOST_TEST( gradCol[ j ] == dMacroLdMicroFlowDirection[ 9 * j + i ] );

        }

    }

}

BOOST_AUTO_TEST_CASE( testComputePlasticMicroVelocityGradient, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

    BOOST_TEST( answerMicroLp == resultMicroLp, CHECK_PER_ELEMENT );

    //Test the Jacobians
    variableVector resultMicroLpJ;
    variableVector dMicroLpdMicroGamma;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                                               microFlowDirection, resultMicroLpJ,
                                                                                               dMicroLpdMicroGamma );

    BOOST_TEST( answerMicroLp == resultMicroLpJ, CHECK_PER_ELEMENT );

    variableVector resultMicroLpJ2;
    variableVector dMicroLpdMicroGamma2;
    variableVector dMicroLpdMicroRCG, dMicroLpdPsie, dMicroLpdMicroFlowDirection;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma, Ce, Psie, invPsie,
                                                                                               microFlowDirection, resultMicroLpJ2,
                                                                                               dMicroLpdMicroGamma2,
                                                                                               dMicroLpdMicroRCG, dMicroLpdPsie,
                                                                                               dMicroLpdMicroFlowDirection );

    BOOST_TEST( answerMicroLp == resultMicroLpJ , CHECK_PER_ELEMENT );

    BOOST_TEST( answerMicroLp == resultMicroLpJ2, CHECK_PER_ELEMENT );

    constantType eps = 1e-6;
    constantType scalarDelta = eps * fabs( microGamma ) + eps;

    variableVector resultMicroLpP, resultMicroLpM;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma + scalarDelta, Ce, Psie, invPsie,
                                                                                               microFlowDirection, resultMicroLpP );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroVelocityGradient( microGamma - scalarDelta, Ce, Psie, invPsie,
                                                                                               microFlowDirection, resultMicroLpM );

    variableVector gradCol = ( resultMicroLpP - resultMicroLpM ) / ( 2 * scalarDelta );

    BOOST_TEST( gradCol == dMicroLpdMicroGamma , CHECK_PER_ELEMENT );

    BOOST_TEST( gradCol == dMicroLpdMicroGamma2, CHECK_PER_ELEMENT );

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

            BOOST_TEST( gradCol[j] == dMicroLpdMicroRCG[ 9 * j + i ] );

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

            BOOST_TEST( gradCol[j] == dMicroLpdPsie[ 9 * j + i ] );

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

            BOOST_TEST( gradCol[j] == dMicroLpdMicroFlowDirection[ 9 * j + i ] );

        }

    }

}

BOOST_AUTO_TEST_CASE( testComputePlasticMicroGradientVelocityGradient, * boost::unit_test::tolerance( 1e-5 ) ){
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

    BOOST_TEST( answerMicroGradLp == resultMicroGradLp, CHECK_PER_ELEMENT );

    variableVector resultMicroGradLp1;
    variableVector resultSkewTerm;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                                       microLp, resultMicroGradLp1,
                                                                                                       resultSkewTerm );

    BOOST_TEST( answerMicroGradLp == resultMicroGradLp1, CHECK_PER_ELEMENT );

    BOOST_TEST( answerSkewTerm == resultSkewTerm, CHECK_PER_ELEMENT );

    //Test the Jacobians
    variableVector resultMicroGradLpJ;
    variableVector dPlasticMicroGradientLdMicroGradientGamma;
    variableVector dPlasticMicroGradientLdPlasticMicroL;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                                       microLp, resultMicroGradLpJ,
                                                                                                       dPlasticMicroGradientLdMicroGradientGamma,
                                                                                                       dPlasticMicroGradientLdPlasticMicroL );

    BOOST_TEST( answerMicroGradLp == resultMicroGradLpJ, CHECK_PER_ELEMENT );

    variableVector resultMicroGradLpJ2, resultSkewTerm2;
    variableVector dPlasticMicroGradientLdMicroGradientGamma2;
    variableVector dPlasticMicroGradientLdPlasticMicroL2;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                                       microLp, resultMicroGradLpJ2, resultSkewTerm2,
                                                                                                       dPlasticMicroGradientLdMicroGradientGamma2,
                                                                                                       dPlasticMicroGradientLdPlasticMicroL2 );

    BOOST_TEST( answerMicroGradLp == resultMicroGradLpJ2, CHECK_PER_ELEMENT );

    BOOST_TEST( answerSkewTerm == resultSkewTerm2, CHECK_PER_ELEMENT );

    variableVector resultMicroGradLpJ3;
    variableVector dPlasticMicroGradientLdMicroGradientGamma3;
    variableVector dPlasticMicroGradientLdPlasticMicroL3;
    variableVector dPlasticMicroGradientLdElasticPsi;
    variableVector dPlasticMicroGradientLdElasticGamma;
    variableVector dPlasticMicroGradientLdMicroGradientFlowDirection;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::computePlasticMicroGradientVelocityGradient( microGradientGamma, Psie, invPsie,
                                                                                                       elasticGamma, microGradientFlowDirection,
                                                                                                       microLp, resultMicroGradLpJ3,
                                                                                                       dPlasticMicroGradientLdMicroGradientGamma3,
                                                                                                       dPlasticMicroGradientLdPlasticMicroL3,
                                                                                                       dPlasticMicroGradientLdElasticPsi,
                                                                                                       dPlasticMicroGradientLdElasticGamma,
                                                                                                       dPlasticMicroGradientLdMicroGradientFlowDirection );

    BOOST_TEST( answerMicroGradLp == resultMicroGradLpJ3, CHECK_PER_ELEMENT );

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
            BOOST_TEST( gradCol[ j ] == dPlasticMicroGradientLdMicroGradientGamma[ 3 * j + i ] );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dPlasticMicroGradientLdMicroGradientGamma2[ 3 * j + i ] );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dPlasticMicroGradientLdMicroGradientGamma3[ 3 * j + i ] );
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
            BOOST_TEST( gradCol[ j ] == dPlasticMicroGradientLdPlasticMicroL[ 9 * j + i ] );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dPlasticMicroGradientLdPlasticMicroL2[ 9 * j + i ] );
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dPlasticMicroGradientLdPlasticMicroL3[ 9 * j + i ] );
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
            BOOST_TEST( gradCol[ j ] == dPlasticMicroGradientLdElasticPsi[ 9 * j + i ] );
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
            BOOST_TEST( gradCol[ j ] == dPlasticMicroGradientLdElasticGamma[ 27 * j + i ] );
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
            BOOST_TEST( gradCol[ j ] == dPlasticMicroGradientLdMicroGradientFlowDirection[ 3 * 27 * j + i ] );
        }
    }

}

BOOST_AUTO_TEST_CASE( test_setPlasticVelocityGradients, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

            floatVector dPrecedingDeformationGradientdF          = initialize( 9 * 9 );

            floatVector dPrecedingDeformationGradientdFn         = initialize( 9 * 9 );

            floatVector previousdPrecedingDeformationGradientdF  = initialize( 9 * 9 );

            floatVector previousdPrecedingDeformationGradientdFn = initialize( 9 * 9 );

            floatVector precedingMicroDeformation = { 0.8, -0.1, 0.2, 0.3, 1.1, 0.4, -0.2, 0.5, 1.2 };

            floatVector previousPrecedingMicroDeformation = { 0.7, 0.11, -0.2, 0.2, 0.9, 0.2, 0.2, -0.1, 1.1 };

            floatVector precedingGradientMicroDeformation = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7 };

            floatVector previousPrecedingGradientMicroDeformation = { 0.2, 0.2, 0.7, 0.4, 1.5, 0.6, 0.7, 0.2, 0.9, 1.4, 1.1, 1.2, 1.8, 1.4, 1.5, 1.6, 1.7, 1.8, 5.9, 2.0, 2.1, 2.2, 2.3, 1.4, 2.5, 2.6, 2.7 };

            floatVector dPrecedingMicroDeformationdChi                      = initialize(  9 * 9 );

            floatVector dPrecedingMicroDeformationdChin                     = initialize(  9 * 9 );

            floatVector previousdPrecedingMicroDeformationdChi              = initialize(  9 * 9 );

            floatVector previousdPrecedingMicroDeformationdChin             = initialize(  9 * 9 );

            floatVector dPrecedingGradientMicroDeformationdFn               = initialize( 27 * 9 );

            floatVector dPrecedingGradientMicroDeformationdChi              = initialize( 27 * 9 );

            floatVector dPrecedingGradientMicroDeformationdChin             = initialize( 27 * 9 );

            floatVector dPrecedingGradientMicroDeformationdGradChi          = initialize( 27 * 27 );

            floatVector dPrecedingGradientMicroDeformationdGradChin         = initialize( 27 * 27 );

            floatVector previousdPrecedingGradientMicroDeformationdFn       = initialize( 27 * 9 );

            floatVector previousdPrecedingGradientMicroDeformationdChi      = initialize( 27 * 9 );

            floatVector previousdPrecedingGradientMicroDeformationdChin     = initialize( 27 * 9 );

            floatVector previousdPrecedingGradientMicroDeformationdGradChi  = initialize( 27 * 27 );

            floatVector previousdPrecedingGradientMicroDeformationdGradChin = initialize( 27 * 27 );

            floatVector plasticMultipliers = { 1.1, 1.2, 1.3, 1.4, 1.5 };

            floatVector previousPlasticMultipliers = { 1.01, 1.02, 1.03, 1.04, 1.05 };

            floatVector dMacroFlowdDrivingStress = { 0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.09 };

            floatVector dMicroFlowdDrivingStress = { 0.21, 0.32, 0.43, 0.54, 0.15, 0.26, 0.37, 0.18, 0.29 };

            floatVector dMicroGradientFlowdDrivingStress = {  .1,  .2,  .6,  .4,  .5,  .6,  .7,  .8,  .9, 1.0, 1.1, 1.2, 2.3, 1.4, 1.5, 1.6, 1.2, 1.8, 1.9, 2.0, 2.1, 1.2, 2.3, 2.4, 2.5, 2.6, 1.7,
                                                             2.8, 5.9, 2.0, 1.1, 3.2, 4.3, 3.4, 3.6, 2.6, 3.1, 3.6, 3.9, 4.0, 1.1, 4.2, 4.3, 4.4, 5.5, 4.6, 4.7, 8.8, 4.4, 5.0, 2.1, 5.6, 5.3, 5.4,
                                                             1.5, 5.9, 5.7, 5.8, 2.9, 6.0, 3.1, 6.2, 1.3, 5.4, 6.5, 3.6, 6.7, 6.8, 8.9, 7.0, 7.1, 7.2, 1.3, 7.4, 7.5, 7.6, 3.7, 7.8, 7.9, 3.2, 8.1 };

            floatVector previousdMacroFlowdDrivingStress = { 0.11, -0.22, 0.33, 0.44, -0.55, 0.66, 0.77, -0.88, 0.99 };

            floatVector previousdMicroFlowdDrivingStress = { 0.21, 0.32, 0.43, 0.54, -0.15, -0.26, 0.37, -0.18, 0.29 };

            floatVector previousdMicroGradientFlowdDrivingStress = {  .1,  .2,  .3,  .4,  .5,  .6,  .7,  .8,  .9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7,
                                                                     2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4,
                                                                     5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1 };

            floatVector d2MacroFlowdDrivingStressdStress                   = initialize( 9 * 9 );

            floatVector previousd2MacroFlowdDrivingStressdStress           = initialize( 9 * 9 );

            floatVector d2MacroFlowdDrivingStressdF                        = initialize( 9 * 9 );

            floatVector previousd2MacroFlowdDrivingStressdF                = initialize( 9 * 9 );

            floatVector d2MacroFlowdDrivingStressdFn                       = initialize( 9 * 9 );

            floatVector previousd2MacroFlowdDrivingStressdFn               = initialize( 9 * 9 );

            floatVector d2MicroFlowdDrivingStressdStress                   = initialize( 9 * 9 );

            floatVector previousd2MicroFlowdDrivingStressdStress           = initialize( 9 * 9 );

            floatVector d2MicroFlowdDrivingStressdF                        = initialize( 9 * 9 );

            floatVector previousd2MicroFlowdDrivingStressdF                = initialize( 9 * 9 );

            floatVector d2MicroFlowdDrivingStressdFn                       = initialize( 9 * 9 );

            floatVector previousd2MicroFlowdDrivingStressdFn               = initialize( 9 * 9 );

            floatVector d2MicroGradientFlowdDrivingStressdStress           = initialize( 3 * 27 * 27 );

            floatVector d2MicroGradientFlowdDrivingStressdF                = initialize( 3 * 27 *  9 );

            floatVector d2MicroGradientFlowdDrivingStressdFn               = initialize( 3 * 27 *  9 );

            floatVector d2MicroGradientFlowdDrivingStressdChi              = initialize( 3 * 27 *  9 );

            floatVector d2MicroGradientFlowdDrivingStressdChin             = initialize( 3 * 27 *  9 );

            floatVector previousd2MicroGradientFlowdDrivingStressdStress   = initialize( 3 * 27 * 27 );

            floatVector previousd2MicroGradientFlowdDrivingStressdF        = initialize( 3 * 27 *  9 );

            floatVector previousd2MicroGradientFlowdDrivingStressdFn       = initialize( 3 * 27 *  9 );

            floatVector previousd2MicroGradientFlowdDrivingStressdChi      = initialize( 3 * 27 *  9 );

            floatVector previousd2MicroGradientFlowdDrivingStressdChin     = initialize( 3 * 27 *  9 );

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
                                     { }, { },
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

    BOOST_TEST( answerMacroL                 == *R.get_plasticMacroVelocityGradient( )                 , CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousMacroL         == *R.get_previousPlasticMacroVelocityGradient( )         , CHECK_PER_ELEMENT );

    BOOST_TEST( answerMicroL                 == *R.get_plasticMicroVelocityGradient( )                 , CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousMicroL         == *R.get_previousPlasticMicroVelocityGradient( )         , CHECK_PER_ELEMENT );

    BOOST_TEST( answerGradientMicroL         == *R.get_plasticGradientMicroVelocityGradient( )         , CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousGradientMicroL == *R.get_previousPlasticGradientMicroVelocityGradient( ) , CHECK_PER_ELEMENT );

    BOOST_TEST( answerMacroL                 == *RJ.get_plasticMacroVelocityGradient( )                , CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousMacroL         == *RJ.get_previousPlasticMacroVelocityGradient( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( answerMicroL                 == *RJ.get_plasticMicroVelocityGradient( )                , CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousMicroL         == *RJ.get_previousPlasticMicroVelocityGradient( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( answerGradientMicroL         == *RJ.get_plasticGradientMicroVelocityGradient( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( answerPreviousGradientMicroL == *RJ.get_previousPlasticGradientMicroVelocityGradient( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setPlasticVelocityGradients2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
                                     { }, { },
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
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            assembled_dMacroLdX[ i ][ j ]     = ( *R.get_dPlasticMacroVelocityGradientdMacroStress( ) )[ 9 * i + j ];

            assembled_dMacroLdX[ i ][ j + 9 ] = ( *R.get_dPlasticMacroVelocityGradientdMicroStress( ) )[ 9 * i + j ];

            assembled_dMicroLdX[ i ][ j + 9 ] = ( *R.get_dPlasticMicroVelocityGradientdMicroStress( ) )[ 9 * i + j ];

            assembled_dMacroLdX[ i ][ j + configuration_unknown_count ] = ( *R.get_dPlasticMacroVelocityGradientdFn( ) )[ 9 * i + j ];

            assembled_dMicroLdX[ i ][ j + configuration_unknown_count ] = ( *R.get_dPlasticMicroVelocityGradientdFn( ) )[ 9 * i + j ];

            assembled_dMicroLdX[ i ][ j + configuration_unknown_count + 9 ] = ( *R.get_dPlasticMicroVelocityGradientdChin( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dMacroLdX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dPlasticMacroVelocityGradientdStateVariables( ) )[ 10 * i + j ];

            assembled_dMicroLdX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dPlasticMicroVelocityGradientdStateVariables( ) )[ 10 * i + j ];

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dMicroGradientLdX[ i ][ j + 9 ] = ( *R.get_dPlasticGradientMicroVelocityGradientdMicroStress( ) )[ 9 * i + j ];

            assembled_dMicroGradientLdX[ i ][ j + configuration_unknown_count ] = ( *R.get_dPlasticGradientMicroVelocityGradientdFn( ) )[ 9 * i + j ];

            assembled_dMicroGradientLdX[ i ][ j + configuration_unknown_count + 9 ] = ( *R.get_dPlasticGradientMicroVelocityGradientdChin( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            assembled_dMicroGradientLdX[ i ][ j + 18 ] = ( *R.get_dPlasticGradientMicroVelocityGradientdHigherOrderStress( ) )[ 27 * i + j ];

            assembled_dMicroGradientLdX[ i ][ j + configuration_unknown_count + 18 ] = ( *R.get_dPlasticGradientMicroVelocityGradientdGradChin( ) )[ 27 * i +  j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dMicroGradientLdX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dPlasticGradientMicroVelocityGradientdStateVariables( ) )[ 10 * i + j ];

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dMacroLdX )         == tardigradeVectorTools::appendVectors( assembled_dMacroLdX )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dMicroLdX )         == tardigradeVectorTools::appendVectors( assembled_dMicroLdX )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMicroGradientLdX ), tardigradeVectorTools::appendVectors( assembled_dMicroGradientLdX ), 3e-5, 1e-6 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( dMacroLdF )         == *R.get_dPlasticMacroVelocityGradientdF( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dMicroLdF )         == *R.get_dPlasticMicroVelocityGradientdF( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dMicroGradientLdF ) == *R.get_dPlasticGradientMicroVelocityGradientdF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation + delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation - delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            BOOST_TEST( 0. == ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) );

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

    BOOST_TEST( tardigradeVectorTools::appendVectors( dMicroLdChi )         == *R.get_dPlasticMicroVelocityGradientdChi( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMicroGradientLdChi ), *R.get_dPlasticGradientMicroVelocityGradientdChi( ), 1e-5, 1e-6 ) );

    for ( unsigned int i = 0; i < 27; i++ ){

        floatVector delta( 27, 0 );

        delta[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation + delta, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation - delta, previousGradientMicroDeformation,
                                          { }, { },
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

            BOOST_TEST( 0. == ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) );

        }

        vp = *Rp.get_plasticMicroVelocityGradient( );

        vm = *Rm.get_plasticMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_TEST( 0. == ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) );

        }

        vp = *Rp.get_plasticGradientMicroVelocityGradient( );

        vm = *Rm.get_plasticGradientMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dMicroGradientLdGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMicroGradientLdGradChi ), *R.get_dPlasticGradientMicroVelocityGradientdGradChi( ), 1e-5, 1e-6 ) );

    // Check previous Jacobians
    for ( unsigned int i = 0; i < configuration_unknown_count; i++ ){

        floatVector delta( configuration_unknown_count, 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        floatVector dpk2(   delta.begin( ) +  0, delta.begin( ) +  9 );
        floatVector dsigma( delta.begin( ) +  9, delta.begin( ) + 18 );
        floatVector dm(     delta.begin( ) + 18, delta.end( )        );

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydrap._local_deltaPK2   = dpk2;

        hydrap._local_deltaSIGMA = dsigma;

        hydrap._local_deltaM     = dm;

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            assembled_previousdMacroLdStress[ i ][ j ]     = ( *R.get_previousdPlasticMacroVelocityGradientdMacroStress( ) )[ 9 * i + j ];

            assembled_previousdMacroLdStress[ i ][ j + 9 ] = ( *R.get_previousdPlasticMacroVelocityGradientdMicroStress( ) )[ 9 * i + j ];

            assembled_previousdMicroLdStress[ i ][ j + 9 ] = ( *R.get_previousdPlasticMicroVelocityGradientdMicroStress( ) )[ 9 * i + j ];

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_previousdMicroGradientLdStress[ i ][ j + 9 ] = ( *R.get_previousdPlasticGradientMicroVelocityGradientdMicroStress( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            assembled_previousdMicroGradientLdStress[ i ][ j + 18 ] = ( *R.get_previousdPlasticGradientMicroVelocityGradientdHigherOrderStress( ) )[ 27 * i + j ];

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdMacroLdStress )          == tardigradeVectorTools::appendVectors( assembled_previousdMacroLdStress )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( previousdMicroLdStress )        , tardigradeVectorTools::appendVectors( assembled_previousdMicroLdStress )        , 1e-5, 1e-6 ) );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( previousdMicroGradientLdStress ), tardigradeVectorTools::appendVectors( assembled_previousdMicroGradientLdStress ), 3e-5, 1e-6 ) );

    for ( unsigned int i = 0; i < previousStateVariables.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            assembled_previousdMacroLdISVs[ i ][ j ] = ( *R.get_previousdPlasticMacroVelocityGradientdFn( ) )[ 9 * i +  j ];

            assembled_previousdMicroLdISVs[ i ][ j ] = ( *R.get_previousdPlasticMicroVelocityGradientdFn( ) )[ 9 * i + j ];

            assembled_previousdMicroLdISVs[ i ][ j + 9 ] = ( *R.get_previousdPlasticMicroVelocityGradientdChin( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_previousdMacroLdISVs[ i ][ j + configuration_unknown_count ] = ( *R.get_previousdPlasticMacroVelocityGradientdStateVariables( ) )[ 10 * i + j ];

            assembled_previousdMicroLdISVs[ i ][ j + configuration_unknown_count ] = ( *R.get_previousdPlasticMicroVelocityGradientdStateVariables( ) )[ 10 * i + j ];

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_previousdMicroGradientLdISVs[ i ][ j ] = ( *R.get_previousdPlasticGradientMicroVelocityGradientdFn( ) )[ 9 * i + j ];

            assembled_previousdMicroGradientLdISVs[ i ][ j + 9 ] = ( *R.get_previousdPlasticGradientMicroVelocityGradientdChin( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            assembled_previousdMicroGradientLdISVs[ i ][ j + 18 ] = ( *R.get_previousdPlasticGradientMicroVelocityGradientdGradChin( ) )[ 27 * i + j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_previousdMicroGradientLdISVs[ i ][ j + configuration_unknown_count ] = ( *R.get_previousdPlasticGradientMicroVelocityGradientdStateVariables( ) )[ 10 * i + j ];

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdMacroLdISVs ) == tardigradeVectorTools::appendVectors( assembled_previousdMacroLdISVs ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdMicroLdISVs ) == tardigradeVectorTools::appendVectors( assembled_previousdMicroLdISVs ), CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dMicroGradientLdX    ), tardigradeVectorTools::appendVectors( assembled_dMicroGradientLdX    ), 1e-5, 1e-6 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdMacroLdF )         == *R.get_previousdPlasticMacroVelocityGradientdF( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdMicroLdF )         == *R.get_previousdPlasticMicroVelocityGradientdF( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdMicroGradientLdF ) == *R.get_previousdPlasticGradientMicroVelocityGradientdF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation + delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation - delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
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

            BOOST_TEST( 0. == ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) );

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

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdMicroLdChi )         == *R.get_previousdPlasticMicroVelocityGradientdChi( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( previousdMicroGradientLdChi ), *R.get_previousdPlasticGradientMicroVelocityGradientdChi( ), 1e-5, 1e-6 ) );

    for ( unsigned int i = 0; i < 27; i++ ){

        floatVector delta( 27, 0 );

        delta[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation + delta,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation - delta,
                                          { }, { },
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

            BOOST_TEST( 0. == ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) );

        }

        vp = *Rp.get_previousPlasticMicroVelocityGradient( );

        vm = *Rm.get_previousPlasticMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            BOOST_TEST( 0. == ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] ) );

        }

        vp = *Rp.get_previousPlasticGradientMicroVelocityGradient( );

        vm = *Rm.get_previousPlasticGradientMicroVelocityGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            previousdMicroGradientLdGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( previousdMicroGradientLdGradChi ), *R.get_previousdPlasticGradientMicroVelocityGradientdGradChi( ), 3e-5, 1e-6 ) );

}

BOOST_AUTO_TEST_CASE( testEvolvePlasticMicroGradChi, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test the evolution of the plastic micro gradient deformation.
     *
     */

    variableType Dt = 7.888463751831797;

    variableVector currentPlasticMicroDeformation = { -0.49993961, -0.36238477, -0.24525394,
                                                      -0.00566907,  0.66797545,  0.43135092,
                                                      -0.20196189,  0.04922572, -0.32425703 };

    variableVector currentPlasticMacroVelocityGradient = { -0.05489573, -0.01980382, -0.06060589,
                                                            1.1610081 ,  0.4002548 ,  0.86866858,
                                                            0.33607202,  0.12218348,  0.25723268 };

    variableVector currentPlasticMicroVelocityGradient = { -21.41995847,  -4.22959956,  31.83295868,
                                                             0.16258786,   0.03647972,  -0.17997075,
                                                            -9.01448709,  -1.96510413,  13.0843427 };

    variableVector currentPlasticMicroGradientVelocityGradient = { -70.66715632,   28.05094192,  -12.32029375,    9.17810714,
                                                                   -19.71784064,  -49.77471891,   40.06703071,  -18.75257585,
                                                                   -87.64798388,  216.34003077,   11.38331058,   52.06095799,
                                                                    44.55771056,    2.81701824,   11.53413479, -317.46326014,
                                                                   -16.12786945,  -75.01963462,  -38.71407841,   26.25281724,
                                                                   -31.89136164,    2.12480823,   -4.79956514,  -26.12979314,
                                                                    27.25829277,  -30.0572412 ,    1.8995271 };

    variableVector previousPlasticMicroDeformation = { -0.62976501, -0.65543815,  0.07568244,
                                                        0.25069488,  0.04845195,  0.85682533,
                                                       -0.69408256,  0.68092037,  0.06724845 };

    variableVector previousPlasticMicroGradient = { 0.02280625, 0.22512007, 0.85575746, 0.1825644 , 0.97734329,
                                                    0.15265694, 0.0984977 , 0.73878709, 0.14943404, 0.50190449,
                                                    0.1394924 , 0.49979228, 0.69297073, 0.03931335, 0.61839353,
                                                    0.32686324, 0.67661231, 0.54568142, 0.60788606, 0.46255347,
                                                    0.94481369, 0.6944252 , 0.29844176, 0.55288832, 0.33853861,
                                                    0.03626558, 0.56699219 };

    variableVector previousPlasticMacroVelocityGradient = { 0.63999789, 0.35571483, 0.63708287,
                                                            1.82175239, 1.47320642, 2.26307223,
                                                            0.9862676 , 1.15883656, 1.36908512 };

    variableVector previousPlasticMicroVelocityGradient = { -0.10429558, -0.17299031, -0.13991628,
                                                            -0.1431904 ,  0.18911688,  0.20668886,
                                                             0.01882534,  0.20244689,  0.16117068 };

    variableVector previousPlasticMicroGradientVelocityGradient = { -1.28938531, -0.78106843, -0.16682568, -0.98659316, -0.33046864,
                                                                     1.00840067, -0.86124472,  0.16817579,  0.26768819,  2.42835457,
                                                                     2.54123143,  0.66070946,  2.88282966,  3.30378025,  0.12668134,
                                                                     3.52709363,  3.17568822,  0.4061731 ,  2.3619489 ,  1.92411292,
                                                                     0.39428284,  0.84845341,  0.29862915,  2.21825071,  1.58372838,
                                                                     0.01050663,  2.0526354 };

    parameterType alpha = 1. - 0.19639211333133877;

    variableVector answerCurrentPlasticMicroGradient = {  201.71721607,   -8.0663384 ,   63.93960976,  299.19611487,
                                                          -15.75827887,   98.74480552, -319.49770439,   13.5458325 ,
                                                         -102.33587287, -195.29766492,   -5.08537627,  -39.70978333,
                                                         -291.02098952,   -3.24570068,  -64.71587095,  314.13870087,
                                                            5.09285553,   66.08381791,  109.96443867,   -5.79910293,
                                                           37.22865464,  163.04783477,  -10.35151403,   59.19030845,
                                                         -175.10575052,   10.13754922,  -60.7404024 };

    variableVector answerLHS = {
        1.3643808e+02, 7.3598994e+00, 2.1304384e+00, 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0., 0., 0., 0., 0., 0., 0., 0.,
        -1.2554098e-01, 1.3932339e+02, 7.7454940e-01, 0., 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0., 0., 0., 0., 0., 0., 0.,
        -3.8419476e-01, 5.5066914e+00, 1.3841674e+02, 0., 0., 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 1.3643808e+02, 7.3598994e+00, 2.1304384e+00, 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0., 0., 0., 0., 0.,
        0., 0., 0., -1.2554098e-01, 1.3932339e+02, 7.7454940e-01, 0., 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0., 0., 0., 0.,
        0., 0., 0., -3.8419476e-01, 5.5066914e+00, 1.3841674e+02, 0., 0., 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 1.3643808e+02, 7.3598994e+00, 2.1304384e+00, 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0., 0.,
        0., 0., 0., 0., 0., 0., -1.2554098e-01, 1.3932339e+02, 7.7454940e-01, 0., 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02, 0.,
        0., 0., 0., 0., 0., 0., -3.8419476e-01, 5.5066914e+00, 1.3841674e+02, 0., 0., 0., 0., 0., 0., 0., 0., 2.6812412e+01, 0., 0., 0., 0., 0., 0., 0., 0., -2.0179650e+02,
        -1.0306821e+00, 0., 0., 0., 0., 0., 0., 0., 0., 4.2074983e-01, 7.3598994e+00, 2.1304384e+00, 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0., 0., 0., 0., 0., 0., 0., 0.,
        0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., 0., -1.2554098e-01, 3.3060545e+00, 7.7454940e-01, 0., 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0., 0., 0., 0., 0., 0., 0.,
        0., 0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., -3.8419476e-01, 5.5066914e+00, 2.3994041e+00, 0., 0., 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0., 0., 0., 0., 0., 0.,
        0., 0., 0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., 0., 0., 4.2074983e-01, 7.3598994e+00, 2.1304384e+00, 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0., 0., 0., 0., 0.,
        0., 0., 0., 0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., 0., -1.2554098e-01, 3.3060545e+00, 7.7454940e-01, 0., 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0., 0., 0., 0.,
        0., 0., 0., 0., 0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., -3.8419476e-01, 5.5066914e+00, 2.3994041e+00, 0., 0., 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0., 0., 0.,
        0., 0., 0., 0., 0., 0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., 0., 0., 4.2074983e-01, 7.3598994e+00, 2.1304384e+00, 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0., 0.,
        0., 0., 0., 0., 0., 0., 0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., 0., -1.2554098e-01, 3.3060545e+00, 7.7454940e-01, 0., 0., 0., 0., 0., 0., 0., 1.1408763e+00, 0.,
        0., 0., 0., 0., 0., 0., 0., 0., -1.0306821e+00, 0., 0., 0., 0., 0., 0., -3.8419476e-01, 5.5066914e+00, 2.3994041e+00, 0., 0., 0., 0., 0., 0., 0., 0., 1.1408763e+00,
        5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., 0., 0., -8.2292677e+01, 7.3598994e+00, 2.1304384e+00, 0., 0., 0., 0., 0., 0.,
        0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., 0., -1.2554098e-01, -7.9407372e+01, 7.7454940e-01, 0., 0., 0., 0., 0., 0.,
        0., 0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., -3.8419476e-01, 5.5066914e+00, -8.0314022e+01, 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., 0., 0., -8.2292677e+01, 7.3598994e+00, 2.1304384e+00, 0., 0., 0.,
        0., 0., 0., 0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., 0., -1.2554098e-01, -7.9407372e+01, 7.7454940e-01, 0., 0., 0.,
        0., 0., 0., 0., 0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., -3.8419476e-01, 5.5066914e+00, -8.0314022e+01, 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., 0., 0., -8.2292677e+01, 7.3598994e+00, 2.1304384e+00,
        0., 0., 0., 0., 0., 0., 0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., 0., -1.2554098e-01, -7.9407372e+01, 7.7454940e-01,
        0., 0., 0., 0., 0., 0., 0., 0., 5.7144922e+01, 0., 0., 0., 0., 0., 0., 0., 0., 1.2457250e+01, 0., 0., 0., 0., 0., 0., -3.8419476e-01, 5.5066914e+00, -8.0314022e+01
    };

    variableVector resultCurrentPlasticMicroGradient;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                     currentPlasticMacroVelocityGradient,
                                                                                     currentPlasticMicroVelocityGradient,
                                                                                     currentPlasticMicroGradientVelocityGradient,
                                                                                     previousPlasticMicroDeformation,
                                                                                     previousPlasticMicroGradient,
                                                                                     previousPlasticMacroVelocityGradient,
                                                                                     previousPlasticMicroVelocityGradient,
                                                                                     previousPlasticMicroGradientVelocityGradient,
                                                                                     resultCurrentPlasticMicroGradient, alpha );

    BOOST_TEST( answerCurrentPlasticMicroGradient == resultCurrentPlasticMicroGradient, CHECK_PER_ELEMENT );

    variableVector resultCurrentPlasticMicroGradient2;
    variableVector LHS2;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                     currentPlasticMacroVelocityGradient,
                                                                                     currentPlasticMicroVelocityGradient,
                                                                                     currentPlasticMicroGradientVelocityGradient,
                                                                                     previousPlasticMicroDeformation,
                                                                                     previousPlasticMicroGradient,
                                                                                     previousPlasticMacroVelocityGradient,
                                                                                     previousPlasticMicroVelocityGradient,
                                                                                     previousPlasticMicroGradientVelocityGradient,
                                                                                     resultCurrentPlasticMicroGradient2, LHS2, alpha );

    BOOST_TEST( answerCurrentPlasticMicroGradient == resultCurrentPlasticMicroGradient2, CHECK_PER_ELEMENT );

    BOOST_TEST( answerLHS == LHS2, CHECK_PER_ELEMENT );

    //Test the Jacobians
    variableVector resultCurrentPlasticMicroGradientJ;

    variableVector dCurrentPlasticMicroGradientdPlasticMicroDeformation,
                   dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient,
                   dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient,
                   dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                     currentPlasticMacroVelocityGradient,
                                                                                     currentPlasticMicroVelocityGradient,
                                                                                     currentPlasticMicroGradientVelocityGradient,
                                                                                     previousPlasticMicroDeformation,
                                                                                     previousPlasticMicroGradient,
                                                                                     previousPlasticMacroVelocityGradient,
                                                                                     previousPlasticMicroVelocityGradient,
                                                                                     previousPlasticMicroGradientVelocityGradient,
                                                                                     resultCurrentPlasticMicroGradientJ,
                                                                                     dCurrentPlasticMicroGradientdPlasticMicroDeformation,
                                                                                     dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient,
                                                                                     dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient,
                                                                                     dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient,
                                                                                     alpha );

    BOOST_TEST( answerCurrentPlasticMicroGradient == resultCurrentPlasticMicroGradientJ, CHECK_PER_ELEMENT );

    variableVector resultCurrentPlasticMicroGradientJ2;

    variableVector dCurrentPlasticMicroGradientdPlasticMicroDeformationJ2,
                   dCurrentPlasticMicroGradientdPlasticMacroVelocityGradientJ2,
                   dCurrentPlasticMicroGradientdPlasticMicroVelocityGradientJ2,
                   dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradientJ2;

    variableVector dCurrentPlasticMicroGradientdPreviousPlasticMicroDeformation,
                   dCurrentPlasticMicroGradientdPreviousPlasticMicroGradient,
                   dCurrentPlasticMicroGradientdPreviousPlasticMacroVelocityGradient,
                   dCurrentPlasticMicroGradientdPreviousPlasticMicroVelocityGradient,
                   dCurrentPlasticMicroGradientdPreviousPlasticMicroGradientVelocityGradient;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                     currentPlasticMacroVelocityGradient,
                                                                                     currentPlasticMicroVelocityGradient,
                                                                                     currentPlasticMicroGradientVelocityGradient,
                                                                                     previousPlasticMicroDeformation,
                                                                                     previousPlasticMicroGradient,
                                                                                     previousPlasticMacroVelocityGradient,
                                                                                     previousPlasticMicroVelocityGradient,
                                                                                     previousPlasticMicroGradientVelocityGradient,
                                                                                     resultCurrentPlasticMicroGradientJ2,
                                                                                     dCurrentPlasticMicroGradientdPlasticMicroDeformationJ2,
                                                                                     dCurrentPlasticMicroGradientdPlasticMacroVelocityGradientJ2,
                                                                                     dCurrentPlasticMicroGradientdPlasticMicroVelocityGradientJ2,
                                                                                     dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradientJ2,
                                                                                     dCurrentPlasticMicroGradientdPreviousPlasticMicroDeformation,
                                                                                     dCurrentPlasticMicroGradientdPreviousPlasticMicroGradient,
                                                                                     dCurrentPlasticMicroGradientdPreviousPlasticMacroVelocityGradient,
                                                                                     dCurrentPlasticMicroGradientdPreviousPlasticMicroVelocityGradient,
                                                                                     dCurrentPlasticMicroGradientdPreviousPlasticMicroGradientVelocityGradient,
                                                                                     alpha );

    BOOST_TEST( answerCurrentPlasticMicroGradient == resultCurrentPlasticMicroGradientJ2, CHECK_PER_ELEMENT );

    BOOST_TEST( dCurrentPlasticMicroGradientdPlasticMicroDeformation == dCurrentPlasticMicroGradientdPlasticMicroDeformationJ2, CHECK_PER_ELEMENT );

    BOOST_TEST( dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient == dCurrentPlasticMicroGradientdPlasticMacroVelocityGradientJ2, CHECK_PER_ELEMENT );

    BOOST_TEST( dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient == dCurrentPlasticMicroGradientdPlasticMicroVelocityGradientJ2, CHECK_PER_ELEMENT );

    BOOST_TEST( dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient == dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradientJ2, CHECK_PER_ELEMENT );

    //Test the jacobian w.r.t. the current plastic macro deformation
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < currentPlasticMicroDeformation.size(); i++ ){
        constantVector delta( currentPlasticMicroDeformation.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMicroDeformation[ i ] ) + eps;

        variableVector resultP, resultM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation + delta,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultP, alpha );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation - delta,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultM, alpha );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( tolerantCheck( gradCol[ j ], dCurrentPlasticMicroGradientdPlasticMicroDeformation[ 9 * j + i ], 1e-5, 1e-6 ) );
        }
    }

    //Test the jacobian w.r.t. the current plastic macro velocity gradient
    for ( unsigned int i = 0; i < currentPlasticMacroVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMacroVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMacroVelocityGradient[ i ] ) + eps;

        variableVector resultP, resultM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient + delta,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultP, alpha );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient - delta,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultM, alpha );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( tolerantCheck( gradCol[ j ], dCurrentPlasticMicroGradientdPlasticMacroVelocityGradient[ 9 * j + i ], 1e-5, 1e-6 ) );
        }
    }

    //Test the jacobian w.r.t. the current plastic micro velocity gradient
    for ( unsigned int i = 0; i < currentPlasticMicroVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMicroVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMicroVelocityGradient[ i ] ) + eps;

        variableVector resultP, resultM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient + delta,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultP, alpha );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient - delta,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultM, alpha );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dCurrentPlasticMicroGradientdPlasticMicroVelocityGradient[ 9 * j + i ] );
        }
    }

    //Test the jacobian w.r.t. the current plastic micro gradient velocity gradient
    for ( unsigned int i = 0; i < currentPlasticMicroGradientVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMicroGradientVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMicroGradientVelocityGradient[ i ] ) + eps;

        variableVector resultP, resultM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient + delta,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultP, alpha );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient - delta,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultM, alpha );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( tolerantCheck( gradCol[ j ], dCurrentPlasticMicroGradientdPlasticMicroGradientVelocityGradient[ 27 * j + i ], 6e-5, 5e-5 ) );
        }
    }

    //Test the jacobian w.r.t. the previous plastic macro deformation
    for ( unsigned int i = 0; i < currentPlasticMicroDeformation.size(); i++ ){
        constantVector delta( currentPlasticMicroDeformation.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMicroDeformation[ i ] ) + eps;

        variableVector resultP, resultM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation + delta,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultP, alpha );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation - delta,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultM, alpha );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( tolerantCheck( gradCol[ j ], dCurrentPlasticMicroGradientdPreviousPlasticMicroDeformation[ 9 * j + i ], 1e-5, 1e-6 ) );
        }
    }

    //Test the jacobian w.r.t. the previous plastic macro gradient
    for ( unsigned int i = 0; i < previousPlasticMicroGradient.size(); i++ ){
        constantVector delta( previousPlasticMicroGradient.size(), 0 );
        delta[i] = eps * fabs( previousPlasticMicroGradient[ i ] ) + eps;

        variableVector resultP, resultM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient + delta,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultP, alpha );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient - delta,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultM, alpha );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( tolerantCheck( gradCol[ j ], dCurrentPlasticMicroGradientdPreviousPlasticMicroGradient[ 27 * j + i ], 5e-5, 1e-6 ) );
        }
    }

    //Test the jacobian w.r.t. the previous plastic macro velocity gradient
    for ( unsigned int i = 0; i < currentPlasticMacroVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMacroVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMacroVelocityGradient[ i ] ) + eps;

        variableVector resultP, resultM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient + delta,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultP, alpha );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient - delta,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultM, alpha );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( tolerantCheck( gradCol[ j ], dCurrentPlasticMicroGradientdPreviousPlasticMacroVelocityGradient[ 9 * j + i ], 1e-4, 1e-6 ) );
        }
    }

    //Test the jacobian w.r.t. the previous plastic micro velocity gradient
    for ( unsigned int i = 0; i < currentPlasticMicroVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMicroVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMicroVelocityGradient[ i ] ) + eps;

        variableVector resultP, resultM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient + delta,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultP, alpha );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient - delta,
                                                                                         previousPlasticMicroGradientVelocityGradient,
                                                                                         resultM, alpha );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( tolerantCheck( gradCol[ j ], dCurrentPlasticMicroGradientdPreviousPlasticMicroVelocityGradient[ 9 * j + i ], 4e-5, 1e-6 ) );
        }
    }

    //Test the jacobian w.r.t. the previous plastic micro gradient velocity gradient
    for ( unsigned int i = 0; i < currentPlasticMicroGradientVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMicroGradientVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMicroGradientVelocityGradient[ i ] ) + eps;

        variableVector resultP, resultM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient + delta,
                                                                                         resultP, alpha );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticMicroGradChi( Dt, currentPlasticMicroDeformation,
                                                                                         currentPlasticMacroVelocityGradient,
                                                                                         currentPlasticMicroVelocityGradient,
                                                                                         currentPlasticMicroGradientVelocityGradient,
                                                                                         previousPlasticMicroDeformation,
                                                                                         previousPlasticMicroGradient,
                                                                                         previousPlasticMacroVelocityGradient,
                                                                                         previousPlasticMicroVelocityGradient,
                                                                                         previousPlasticMicroGradientVelocityGradient - delta,
                                                                                         resultM, alpha );

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( tolerantCheck( gradCol[ j ], dCurrentPlasticMicroGradientdPreviousPlasticMicroGradientVelocityGradient[ 27 * j + i ], 3e-5, 1e-6 ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( testEvolvePlasticDeformation, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Evolve the plastic deformation.
     *
     */

    variableType Dt = 8.009359239014827;
    variableType alphaMacro = 1 - 0.48581069403548804;
    variableType alphaMicro = 1 - 0.827106668527372;
    variableType alphaMicroGrad = 1 - 0.608865458065455;

    variableVector currentPlasticMacroVelocityGradient = { -0.27014632, -0.35189361,  0.13880036,
                                                           -0.35584915, -0.48046671,  0.0912627 ,
                                                            0.47263309, -0.31133122, -0.2948745 };

    variableVector currentPlasticMicroVelocityGradient = { -0.01186714,  0.06090556,  0.40756696,
                                                            0.32909297,  0.28343008,  0.32993809,
                                                           -0.1899835 ,  0.41541764, -0.04761333 };

    variableVector currentPlasticMicroGradientVelocityGradient = { -0.44741758,  0.19632574,  0.20883732, -0.18740653,  0.30768804,
                                                                    0.46878663,  0.17691702, -0.07494477, -0.29871104, -0.46454855,
                                                                   -0.34339229, -0.24065171, -0.0737629 ,  0.48339896,  0.25954592,
                                                                    0.06579716, -0.43873995,  0.46482816,  0.20699556, -0.2417395 ,
                                                                   -0.46792573,  0.39616741, -0.22954849, -0.097399  ,  0.20759377,
                                                                   -0.35477122,  0.21945037 };

    variableVector previousPlasticDeformationGradient = { -0.45932572, -0.49691079, -0.13107169,
                                                           0.22431024, -0.61998429, -0.62903476,
                                                          -0.90544313,  0.41193656, -0.14749891 };

    variableVector previousPlasticMicroDeformation = { 0.28821738,  0.34407584,  0.06697721,
                                                       0.31626116, -0.62278299,  0.83699009,
                                                       0.16315786, -0.9434692 , -0.09991987 };

    variableVector previousPlasticMicroGradient = { -0.4908299 , -0.09107442,  0.11438472, -0.468681  , -0.09518528,
                                                     0.39306299,  0.30090607, -0.47472135, -0.39873768,  0.2071941 ,
                                                     0.3967767 , -0.40333359, -0.10263381, -0.39238343, -0.21260101,
                                                    -0.36133084, -0.32771357,  0.49025487,  0.01729508, -0.10338882,
                                                     0.15828675,  0.31044611, -0.30362281, -0.03972172, -0.377879  ,
                                                     0.11046526, -0.46997605 };

    variableVector previousPlasticMacroVelocityGradient = { -0.4362157 , -0.37728779,  0.40148895,
                                                            -0.15976616, -0.15839888,  0.26580961,
                                                             0.23215055, -0.02231055,  0.35815526 };

    variableVector previousPlasticMicroVelocityGradient = { 0.07222414, -0.40261712,  0.48556998,
                                                            0.45642476, -0.09409791, -0.49992287,
                                                            0.48252843,  0.48267731, -0.06878613 };

    variableVector previousPlasticMicroGradientVelocityGradient = { -0.15695492,  0.45802473, -0.38830808, -0.29241496, -0.34225716,
                                                                    -0.26127347, -0.41201504, -0.36558117, -0.24714561,  0.09340483,
                                                                     0.40887886, -0.1481247 ,  0.15128499, -0.04882876, -0.30445054,
                                                                     0.11557493, -0.20789811, -0.33913681,  0.3142761 ,  0.09871806,
                                                                     0.11316847, -0.45596559,  0.19030294, -0.33056333, -0.49391146,
                                                                     0.40129161, -0.06836289 };

    variableVector answerMacro = { -1.56924249,  1.75893631,  0.91035498,
                                    0.24241493, -0.43775652, -0.44317337,
                                   -2.6946146 ,  2.06645595,  0.86856673 };

    variableVector answerMicro = { 3.30247398,  2.51469186,  1.61673915,
                                   6.78559987, 11.60984671,  7.25516648,
                                   4.74326427,  4.32671877,  6.17698967 };

    variableVector answerMicroGrad = {  0.18184066,   0.06902067,  -4.65044761,  -4.6051252 ,
                                       -0.67925892,  -5.62411017,  -5.22449727,   0.2149805 ,
                                       -8.35852814, -14.60027145,  11.46879229,   1.1063327 ,
                                      -30.40158864,  24.16079446,  -3.60050806, -20.35119344,
                                       18.13574972,  -5.22922306,   4.00630219,  -2.39707255,
                                       -7.34817085,   3.12285942,  -6.24320576, -15.1020395 ,
                                       -5.11190824,   3.50558822,  -4.01621741 };

//    variableVector answerMacro = { -1.9572735 ,  0.60966488,  0.59517042, -0.44336072,  0.11557451,
//        -0.04102717, -4.31902521,  1.32203177,  1.08289253 };
//
//    variableVector answerMicro = { 8.32601386, -9.56508067,  2.08126161,  1.7434483 , -2.10943576,
//         0.31240494,  7.77866959, -8.85544348,  1.79336953 };
//
//    variableVector answerMicroGrad = {   3.02800309,   2.82050978,  -5.32790487,  -2.92149612,
//        -3.57991191,   5.52402095,   0.03796393,  -0.20153694,
//        -2.20531614,  -4.53803741,  12.87898183,   3.74466913,
//         1.64578078, -13.45516056,  -4.72016453,   0.63316434,
//         3.24138773,   0.69421099,   1.61882075,  -3.61139697,
//        -0.76706456,   0.19881282,   4.22384657,  -1.71230232,
//        -0.7414286 ,   0.85163677,   1.87885778 };

    variableVector resultMacro, resultMicro, resultMicroGrad;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                    currentPlasticMicroVelocityGradient,
                                                                                    currentPlasticMicroGradientVelocityGradient,
                                                                                    previousPlasticDeformationGradient,
                                                                                    previousPlasticMicroDeformation,
                                                                                    previousPlasticMicroGradient,
                                                                                    previousPlasticMacroVelocityGradient,
                                                                                    previousPlasticMicroVelocityGradient,
                                                                                    previousPlasticMicroGradientVelocityGradient,
                                                                                    resultMacro, resultMicro, resultMicroGrad,
                                                                                    alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_TEST( resultMacro == answerMacro, CHECK_PER_ELEMENT );

    BOOST_TEST( resultMicro == answerMicro, CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( resultMicroGrad, answerMicroGrad, 1e-5, 1e-6 ) );

    //Test the computation of the Jacobians
    variableVector resultMacroJ, resultMicroJ, resultMicroGradJ;
    variableVector dFdMacroL, dChidMicroL, dGradChidMacroL, dGradChidMicroL, dGradChidMicroGradL;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                    currentPlasticMicroVelocityGradient,
                                                                                    currentPlasticMicroGradientVelocityGradient,
                                                                                    previousPlasticDeformationGradient,
                                                                                    previousPlasticMicroDeformation,
                                                                                    previousPlasticMicroGradient,
                                                                                    previousPlasticMacroVelocityGradient,
                                                                                    previousPlasticMicroVelocityGradient,
                                                                                    previousPlasticMicroGradientVelocityGradient,
                                                                                    resultMacroJ, resultMicroJ, resultMicroGradJ,
                                                                                    dFdMacroL, dChidMicroL, dGradChidMacroL,
                                                                                    dGradChidMicroL, dGradChidMicroGradL,
                                                                                    alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_TEST( resultMacroJ == answerMacro, CHECK_PER_ELEMENT );

    BOOST_TEST( resultMicroJ == answerMicro, CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( resultMicroGradJ, answerMicroGrad, 1e-5, 1e-6 ) );

    variableVector resultMacroJ2, resultMicroJ2, resultMicroGradJ2;
    variableVector dFdMacroLJ2, dChidMicroLJ2, dGradChidMacroLJ2, dGradChidMicroLJ2, dGradChidMicroGradLJ2,
                   dFdPreviousF, dFdPreviousMacroL, dChidPreviousChi, dChidPreviousMicroL, dGradChidPreviousMacroL,
                   dGradChidPreviousMicroL, dGradChidPreviousMicroGradL, dGradChidPreviousChi, dGradChidPreviousGradChi;

    tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                    currentPlasticMicroVelocityGradient,
                                                                                    currentPlasticMicroGradientVelocityGradient,
                                                                                    previousPlasticDeformationGradient,
                                                                                    previousPlasticMicroDeformation,
                                                                                    previousPlasticMicroGradient,
                                                                                    previousPlasticMacroVelocityGradient,
                                                                                    previousPlasticMicroVelocityGradient,
                                                                                    previousPlasticMicroGradientVelocityGradient,
                                                                                    resultMacroJ2, resultMicroJ2, resultMicroGradJ2,
                                                                                    dFdMacroLJ2, dChidMicroLJ2, dGradChidMacroLJ2,
                                                                                    dGradChidMicroLJ2, dGradChidMicroGradLJ2,
                                                                                    dFdPreviousF, dFdPreviousMacroL,
                                                                                    dChidPreviousChi, dChidPreviousMicroL,
                                                                                    dGradChidPreviousChi, dGradChidPreviousGradChi,
                                                                                    dGradChidPreviousMacroL, dGradChidPreviousMicroL, dGradChidPreviousMicroGradL,
                                                                                    alphaMacro, alphaMicro, alphaMicroGrad );

    BOOST_TEST( resultMacroJ2 == answerMacro, CHECK_PER_ELEMENT );

    BOOST_TEST( resultMicroJ2 == answerMicro, CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( resultMicroGradJ2, answerMicroGrad, 1e-5, 1e-6 ) );

    BOOST_TEST( dFdMacroLJ2 == dFdMacroL, CHECK_PER_ELEMENT );

    BOOST_TEST( dChidMicroLJ2 == dChidMicroL, CHECK_PER_ELEMENT );

    BOOST_TEST( dGradChidMacroLJ2 == dGradChidMacroL, CHECK_PER_ELEMENT );

    BOOST_TEST( dGradChidMicroLJ2 == dGradChidMicroL, CHECK_PER_ELEMENT );

    BOOST_TEST( dGradChidMicroGradLJ2 == dGradChidMicroGradL, CHECK_PER_ELEMENT );

    //Test jacobians w.r.t. the current plastic macro velocity gradient
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < currentPlasticMacroVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMacroVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMacroVelocityGradient[ i ] ) + eps;

        variableVector resultMacroP, resultMicroP, resultMicroGradP;
        variableVector resultMacroM, resultMicroM, resultMicroGradM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient + delta,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroP, resultMicroP, resultMicroGradP,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient - delta,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroM, resultMicroM, resultMicroGradM,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dFdMacroL[ 9 * j + i ] );
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == 0. );
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dGradChidMacroL[ 9 * j + i ] );
        }
    }

    for ( unsigned int i = 0; i < currentPlasticMicroVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMicroVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMicroVelocityGradient[ i ] ) + eps;

        variableVector resultMacroP, resultMicroP, resultMicroGradP;
        variableVector resultMacroM, resultMicroM, resultMicroGradM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient + delta,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroP, resultMicroP, resultMicroGradP,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient - delta,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroM, resultMicroM, resultMicroGradM,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == 0. );
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dChidMicroL[ 9 * j + i ] );
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dGradChidMicroL[ 9 * j + i ] );
        }
    }

    for ( unsigned int i = 0; i < currentPlasticMicroGradientVelocityGradient.size(); i++ ){
        constantVector delta( currentPlasticMicroGradientVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( currentPlasticMicroGradientVelocityGradient[ i ] ) + eps;

        variableVector resultMacroP, resultMicroP, resultMicroGradP;
        variableVector resultMacroM, resultMicroM, resultMicroGradM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient + delta,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroP, resultMicroP, resultMicroGradP,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient - delta,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroM, resultMicroM, resultMicroGradM,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == 0. );
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == 0. );
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dGradChidMicroGradL[ 27 * j + i ] );
        }
    }

    for ( unsigned int i = 0; i < previousPlasticDeformationGradient.size(); i++ ){
        constantVector delta( previousPlasticDeformationGradient.size(), 0 );
        delta[i] = eps * fabs( previousPlasticDeformationGradient[ i ] ) + eps;

        variableVector resultMacroP, resultMicroP, resultMicroGradP;
        variableVector resultMacroM, resultMicroM, resultMicroGradM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient + delta,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroP, resultMicroP, resultMicroGradP,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient - delta,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroM, resultMicroM, resultMicroGradM,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dFdPreviousF[ 9 * j + i ] );
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == 0. );
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == 0. );
        }
    }

    for ( unsigned int i = 0; i < previousPlasticMicroDeformation.size(); i++ ){
        constantVector delta( previousPlasticMicroDeformation.size(), 0 );
        delta[i] = eps * fabs( previousPlasticMicroDeformation[ i ] ) + eps;

        variableVector resultMacroP, resultMicroP, resultMicroGradP;
        variableVector resultMacroM, resultMicroM, resultMicroGradM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation + delta,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroP, resultMicroP, resultMicroGradP,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation - delta,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroM, resultMicroM, resultMicroGradM,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == 0. );
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dChidPreviousChi[ 9 * j + i ] );
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( tolerantCheck( gradCol[ j ], dGradChidPreviousChi[ 9 * j + i ], 1e-6, 1e-6 ) );
        }
    }

    for ( unsigned int i = 0; i < previousPlasticMicroGradient.size(); i++ ){
        constantVector delta( previousPlasticMicroGradient.size(), 0 );
        delta[i] = eps * fabs( previousPlasticMicroGradient[ i ] ) + eps;

        variableVector resultMacroP, resultMicroP, resultMicroGradP;
        variableVector resultMacroM, resultMicroM, resultMicroGradM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient + delta,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroP, resultMicroP, resultMicroGradP,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient - delta,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroM, resultMicroM, resultMicroGradM,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == 0. );
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == 0. );
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( tolerantCheck( gradCol[ j ], dGradChidPreviousGradChi[ 27 * j + i ], 1e-5, 1e-6 ) );
        }
    }

    for ( unsigned int i = 0; i < previousPlasticMacroVelocityGradient.size(); i++ ){
        constantVector delta( previousPlasticMacroVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( previousPlasticMacroVelocityGradient[ i ] ) + eps;

        variableVector resultMacroP, resultMicroP, resultMicroGradP;
        variableVector resultMacroM, resultMicroM, resultMicroGradM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient + delta,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroP, resultMicroP, resultMicroGradP,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient - delta,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroM, resultMicroM, resultMicroGradM,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dFdPreviousMacroL[ 9 * j + i ] );
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == 0. );
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dGradChidPreviousMacroL[ 9 * j + i ] );
        }
    }

    for ( unsigned int i = 0; i < previousPlasticMicroVelocityGradient.size(); i++ ){
        constantVector delta( previousPlasticMicroVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( previousPlasticMicroVelocityGradient[ i ] ) + eps;

        variableVector resultMacroP, resultMicroP, resultMicroGradP;
        variableVector resultMacroM, resultMicroM, resultMicroGradM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient + delta,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroP, resultMicroP, resultMicroGradP,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient - delta,
                                                                                        previousPlasticMicroGradientVelocityGradient,
                                                                                        resultMacroM, resultMicroM, resultMicroGradM,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == 0. );
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dChidPreviousMicroL[ 9 * j + i ] );
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == dGradChidPreviousMicroL[ 9 * j + i ] );
        }
    }

    for ( unsigned int i = 0; i < previousPlasticMicroGradientVelocityGradient.size(); i++ ){
        constantVector delta( previousPlasticMicroGradientVelocityGradient.size(), 0 );
        delta[i] = eps * fabs( previousPlasticMicroGradientVelocityGradient[ i ] ) + eps;

        variableVector resultMacroP, resultMicroP, resultMicroGradP;
        variableVector resultMacroM, resultMicroM, resultMicroGradM;

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient + delta,
                                                                                        resultMacroP, resultMicroP, resultMicroGradP,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( Dt, currentPlasticMacroVelocityGradient,
                                                                                        currentPlasticMicroVelocityGradient,
                                                                                        currentPlasticMicroGradientVelocityGradient,
                                                                                        previousPlasticDeformationGradient,
                                                                                        previousPlasticMicroDeformation,
                                                                                        previousPlasticMicroGradient,
                                                                                        previousPlasticMacroVelocityGradient,
                                                                                        previousPlasticMicroVelocityGradient,
                                                                                        previousPlasticMicroGradientVelocityGradient - delta,
                                                                                        resultMacroM, resultMicroM, resultMicroGradM,
                                                                                        alphaMacro, alphaMicro, alphaMicroGrad );

        variableVector gradCol = ( resultMacroP - resultMacroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == 0. );
        }

        gradCol = ( resultMicroP - resultMicroM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( gradCol[ j ] == 0. );
        }

        gradCol = ( resultMicroGradP - resultMicroGradM ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            BOOST_TEST( tolerantCheck( gradCol[ j ], dGradChidPreviousMicroGradL[ 27 * j + i ], 1e-5, 1e-6 ) );
        }
    }

}

BOOST_AUTO_TEST_CASE( test_setPlasticDeformation, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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

            floatVector macroL = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            floatVector microL = { 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09 };

            floatVector microGradientL = { 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19,
                                           0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28,
                                           0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37 };

            floatVector previousMacroL = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            floatVector previousMicroL = { 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09 };

            floatVector previousMicroGradientL = { 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19,
                                                   0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28,
                                                   0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37 };

            floatVector dMacroLdMacroStress                   = initialize(  9 *  9 );

            floatVector previousdMacroLdMacroStress           = initialize(  9 *  9 );

            floatVector dMacroLdMicroStress                   = initialize(  9 *  9 );

            floatVector previousdMacroLdMicroStress           = initialize(  9 *  9 );

            floatVector dMicroLdMicroStress                   = initialize(  9 *  9 );

            floatVector previousdMicroLdMicroStress           = initialize(  9 *  9 );

            floatVector dGradMicroLdMicroStress               = initialize( 27 *  9 );

            floatVector previousdGradMicroLdMicroStress       = initialize( 27 *  9 );

            floatVector dGradMicroLdHigherOrderStress         = initialize( 27 * 27 );

            floatVector previousdGradMicroLdHigherOrderStress = initialize( 27 * 27 );

            floatVector dMacroLdF                             = initialize(  9 *  9 );

            floatVector previousdMacroLdF                     = initialize(  9 *  9 );

            floatVector dMicroLdF                             = initialize(  9 *  9 );

            floatVector previousdMicroLdF                     = initialize(  9 *  9 );

            floatVector dMicroGradientLdF                     = initialize( 27 *  9 );

            floatVector previousdMicroGradientLdF             = initialize( 27 *  9 );

            floatVector dMacroLdFn                            = initialize(  9 *  9 );

            floatVector previousdMacroLdFn                    = initialize(  9 *  9 );

            floatVector dMicroLdFn                            = initialize(  9 *  9 );

            floatVector previousdMicroLdFn                    = initialize(  9 *  9 );

            floatVector dMicroGradientLdFn                    = initialize( 27 *  9 );

            floatVector previousdMicroGradientLdFn            = initialize( 27 *  9 );

            floatVector dMicroLdChi                           = initialize(  9 *  9 );

            floatVector previousdMicroLdChi                   = initialize(  9 *  9 );

            floatVector dMicroGradientLdChi                   = initialize( 27 *  9 );

            floatVector previousdMicroGradientLdChi           = initialize( 27 *  9 );

            floatVector dMicroLdChin                          = initialize(  9 *  9 );

            floatVector previousdMicroLdChin                  = initialize(  9 *  9 );

            floatVector dMicroGradientLdChin                  = initialize( 27 *  9 );

            floatVector previousdMicroGradientLdChin          = initialize( 27 *  9 );

            floatVector dMicroGradientLdGradChi               = initialize( 27 * 27 );

            floatVector previousdMicroGradientLdGradChi       = initialize( 27 * 27 );

            floatVector dMicroGradientLdGradChin              = initialize( 27 * 27 );

            floatVector previousdMicroGradientLdGradChin      = initialize( 27 * 27 );

            floatVector dMacroLdXi                            = initialize(  9 * 10 );

            floatVector previousdMacroLdXi                    = initialize(  9 * 10 );

            floatVector dMicroLdXi                            = initialize(  9 * 10 );

            floatVector previousdMicroLdXi                    = initialize(  9 * 10 );

            floatVector dMicroGradientLdXi                    = initialize( 27 * 10 );

            floatVector previousdMicroGradientLdXi            = initialize( 27 * 10 );

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

            virtual void setPlasticVelocityGradients( const bool isPrevious ) override{

                if ( isPrevious ){

                    set_previousPlasticMacroVelocityGradient( previousMacroL );

                    set_previousPlasticMicroVelocityGradient( previousMicroL );

                    set_previousPlasticGradientMicroVelocityGradient( previousMicroGradientL );

                }
                else{

                    set_plasticMacroVelocityGradient( macroL );

                    set_plasticMicroVelocityGradient( microL );

                    set_plasticGradientMicroVelocityGradient( microGradientL );

                }

            }

            virtual void setPlasticVelocityGradientsJacobians( const bool isPrevious ) override{

                setPlasticVelocityGradients( isPrevious );

                if ( isPrevious ){

                    set_previousdPlasticMacroVelocityGradientdMacroStress(               previousdMacroLdMacroStress );

                    set_previousdPlasticMacroVelocityGradientdMicroStress(               previousdMacroLdMicroStress );

                    set_previousdPlasticMicroVelocityGradientdMicroStress(               previousdMicroLdMicroStress );

                    set_previousdPlasticGradientMicroVelocityGradientdMicroStress(       previousdGradMicroLdMicroStress );

                    set_previousdPlasticGradientMicroVelocityGradientdHigherOrderStress( previousdGradMicroLdHigherOrderStress );

                    set_previousdPlasticMacroVelocityGradientdF(                         previousdMacroLdF );

                    set_previousdPlasticMicroVelocityGradientdF(                         previousdMicroLdF );

                    set_previousdPlasticGradientMicroVelocityGradientdF(                 previousdMicroGradientLdF );

                    set_previousdPlasticMacroVelocityGradientdFn(                        previousdMacroLdFn );

                    set_previousdPlasticMicroVelocityGradientdFn(                        previousdMicroLdFn );

                    set_previousdPlasticGradientMicroVelocityGradientdFn(                previousdMicroGradientLdFn );

                    set_previousdPlasticMicroVelocityGradientdChi(                       previousdMicroLdChi );

                    set_previousdPlasticGradientMicroVelocityGradientdChi(               previousdMicroGradientLdChi );

                    set_previousdPlasticMicroVelocityGradientdChin(                      previousdMicroLdChin );

                    set_previousdPlasticGradientMicroVelocityGradientdChin(              previousdMicroGradientLdChin );

                    set_previousdPlasticGradientMicroVelocityGradientdGradChi(           previousdMicroGradientLdGradChi );

                    set_previousdPlasticGradientMicroVelocityGradientdGradChin(          previousdMicroGradientLdGradChin );

                    set_previousdPlasticMacroVelocityGradientdStateVariables(            previousdMacroLdXi );

                    set_previousdPlasticMicroVelocityGradientdStateVariables(            previousdMicroLdXi );

                    set_previousdPlasticGradientMicroVelocityGradientdStateVariables(    previousdMicroGradientLdXi );

                }
                else{

                    set_dPlasticMacroVelocityGradientdMacroStress(               dMacroLdMacroStress );

                    set_dPlasticMacroVelocityGradientdMicroStress(               dMacroLdMicroStress );

                    set_dPlasticMicroVelocityGradientdMicroStress(               dMicroLdMicroStress );

                    set_dPlasticGradientMicroVelocityGradientdMicroStress(       dGradMicroLdMicroStress );

                    set_dPlasticGradientMicroVelocityGradientdHigherOrderStress( dGradMicroLdHigherOrderStress );

                    set_dPlasticMacroVelocityGradientdF(                         dMacroLdF );

                    set_dPlasticMicroVelocityGradientdF(                         dMicroLdF );

                    set_dPlasticGradientMicroVelocityGradientdF(                 dMicroGradientLdF );

                    set_dPlasticMacroVelocityGradientdFn(                        dMacroLdFn );

                    set_dPlasticMicroVelocityGradientdFn(                        dMicroLdFn );

                    set_dPlasticGradientMicroVelocityGradientdFn(                dMicroGradientLdFn );

                    set_dPlasticMicroVelocityGradientdChi(                       dMicroLdChi );

                    set_dPlasticGradientMicroVelocityGradientdChi(               dMicroGradientLdChi );

                    set_dPlasticMicroVelocityGradientdChin(                      dMicroLdChin );

                    set_dPlasticGradientMicroVelocityGradientdChin(              dMicroGradientLdChin );

                    set_dPlasticGradientMicroVelocityGradientdGradChi(           dMicroGradientLdGradChi );

                    set_dPlasticGradientMicroVelocityGradientdGradChin(          dMicroGradientLdGradChin );

                    set_dPlasticMacroVelocityGradientdStateVariables(            dMacroLdXi );

                    set_dPlasticMicroVelocityGradientdStateVariables(            dMicroLdXi );

                    set_dPlasticGradientMicroVelocityGradientdStateVariables(    dMicroGradientLdXi );

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    residualMock RJ( &hydra, 55, 1, stateVariableIndices, parameters );

    floatVector plasticF;

    floatVector plasticChi;

    floatVector plasticGradChi;

    floatVector eye( 9 );
    tardigradeVectorTools::eye( eye );

    floatVector previousPlasticF( previousStateVariables.begin( ) + 0,
                                  previousStateVariables.begin( ) + 9 );

    previousPlasticF += eye;

    floatVector previousPlasticChi( previousStateVariables.begin( ) +  9,
                                    previousStateVariables.begin( ) + 18 );

    previousPlasticChi += eye;

    floatVector previousPlasticGradChi( previousStateVariables.begin( ) + 18,
                                        previousStateVariables.begin( ) + 45 );

    tardigradeHydra::micromorphicDruckerPragerPlasticity::evolvePlasticDeformation( deltaTime, R.macroL, R.microL, R.microGradientL,
                                                                                    previousPlasticF, previousPlasticChi, previousPlasticGradChi,
                                                                                    R.previousMacroL, R.previousMicroL, R.previousMicroGradientL,
                                                                                    plasticF, plasticChi, plasticGradChi );

    RJ.get_dUpdatedPlasticDeformationGradientdMacroStress( );
    RJ.get_dUpdatedPlasticDeformationGradientdPreviousMacroStress( );

    BOOST_TEST(       plasticF == *R.get_updatedPlasticDeformationGradient( )      , CHECK_PER_ELEMENT );

    BOOST_TEST(     plasticChi == *R.get_updatedPlasticMicroDeformation( )         , CHECK_PER_ELEMENT );

    BOOST_TEST( plasticGradChi == *R.get_updatedPlasticGradientMicroDeformation( ) , CHECK_PER_ELEMENT );

    BOOST_TEST(       plasticF == *RJ.get_updatedPlasticDeformationGradient( )     , CHECK_PER_ELEMENT );

    BOOST_TEST(     plasticChi == *RJ.get_updatedPlasticMicroDeformation( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( plasticGradChi == *RJ.get_updatedPlasticGradientMicroDeformation( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setPlasticDeformation2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    // Test the jacobians

    floatType eps = 1e-6;

    floatMatrix dUpdatedPlasticFdX(                      9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dUpdatedPlasticChidX(                    9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dUpdatedPlasticGradChidX(               27, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dUpdatedPlasticFdF(                      9, floatVector( 9, 0 ) );

    floatMatrix dUpdatedPlasticChidF(                    9, floatVector( 9, 0 ) );

    floatMatrix dUpdatedPlasticGradChidF(               27, floatVector( 9, 0 ) );

    floatMatrix dUpdatedPlasticFdChi(                    9, floatVector( 9, 0 ) );

    floatMatrix dUpdatedPlasticChidChi(                  9, floatVector( 9, 0 ) );

    floatMatrix dUpdatedPlasticGradChidChi(             27, floatVector( 9, 0 ) );

    floatMatrix dUpdatedPlasticFdGradChi(                9, floatVector( 27, 0 ) );

    floatMatrix dUpdatedPlasticChidGradChi(              9, floatVector( 27, 0 ) );

    floatMatrix dUpdatedPlasticGradChidGradChi(         27, floatVector( 27, 0 ) );

    floatMatrix dUpdatedPlasticFdPreviousStress(         9, floatVector( configuration_unknown_count, 0 ) );

    floatMatrix dUpdatedPlasticChidPreviousStress(       9, floatVector( configuration_unknown_count, 0 ) );

    floatMatrix dUpdatedPlasticGradChidPreviousStress(  27, floatVector( configuration_unknown_count, 0 ) );

    floatMatrix dUpdatedPlasticFdPreviousISVs(           9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix dUpdatedPlasticChidPreviousISVs(         9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix dUpdatedPlasticGradChidPreviousISVs(    27, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix dUpdatedPlasticFdPreviousF(              9, floatVector(  9, 0 ) );

    floatMatrix dUpdatedPlasticChidPreviousF(            9, floatVector(  9, 0 ) );

    floatMatrix dUpdatedPlasticGradChidPreviousF(       27, floatVector(  9, 0 ) );

    floatMatrix dUpdatedPlasticFdPreviousChi(            9, floatVector(  9, 0 ) );

    floatMatrix dUpdatedPlasticChidPreviousChi(          9, floatVector(  9, 0 ) );

    floatMatrix dUpdatedPlasticGradChidPreviousChi(     27, floatVector(  9, 0 ) );

    floatMatrix dUpdatedPlasticFdPreviousGradChi(        9, floatVector( 27, 0 ) );

    floatMatrix dUpdatedPlasticChidPreviousGradChi(      9, floatVector( 27, 0 ) );

    floatMatrix dUpdatedPlasticGradChidPreviousGradChi( 27, floatVector( 27, 0 ) );

    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_dUpdatedPlasticFdX(        9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix assembled_dUpdatedPlasticChidX(      9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix assembled_dUpdatedPlasticGradChidX( 27, floatVector( unknownVector.size( ), 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dUpdatedPlasticFdX[ i ][ j ] = ( *R.get_dUpdatedPlasticDeformationGradientdMacroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticFdX[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticDeformationGradientdMicroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticChidX[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticMicroDeformationdMicroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticFdX[ i ][ j + configuration_unknown_count ] = ( *R.get_dUpdatedPlasticDeformationGradientdFn( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticChidX[ i ][ j + configuration_unknown_count ] = ( *R.get_dUpdatedPlasticMicroDeformationdFn( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticChidX[ i ][ j + configuration_unknown_count + 9 ] = ( *R.get_dUpdatedPlasticMicroDeformationdChin( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dUpdatedPlasticFdX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dUpdatedPlasticDeformationGradientdStateVariables( ) )[ 10 * i + j ];

            assembled_dUpdatedPlasticChidX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dUpdatedPlasticMicroDeformationdStateVariables( ) )[ 10 * i + j ];

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dUpdatedPlasticGradChidX[ i ][ j + 0 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdMacroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticGradChidX[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdMicroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticGradChidX[ i ][ j + configuration_unknown_count ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdFn( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticGradChidX[ i ][ j + configuration_unknown_count + 9 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdChin( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            assembled_dUpdatedPlasticGradChidX[ i ][ j + 18 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdHigherOrderStress( ) )[ 27 * i + j ];

            assembled_dUpdatedPlasticGradChidX[ i ][ j + configuration_unknown_count + 18 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdGradChin( ) )[ 27 * i + j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dUpdatedPlasticGradChidX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdStateVariables( ) )[ 10 * i + j ];

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdX       ), tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticFdX       ), 1e-5, 1e-6 ) );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidX     ) == tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticChidX     ), CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidX ), tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticGradChidX ), 4e-5, 1e-5 ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdF )       == *R.get_dUpdatedPlasticDeformationGradientdF( )     , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidF )     == *R.get_dUpdatedPlasticMicroDeformationdF( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidF ) == *R.get_dUpdatedPlasticGradientMicroDeformationdF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < microDeformation.size( ); i++ ){

        floatVector delta( microDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation + delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation - delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdChi )       == floatVector( 9 * 9, 0 )                              , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidChi )     == *R.get_dUpdatedPlasticMicroDeformationdChi( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidChi ), *R.get_dUpdatedPlasticGradientMicroDeformationdChi( ), 1e-5, 1e-6 ) );

    for ( unsigned int i = 0; i < gradientMicroDeformation.size( ); i++ ){

        floatVector delta( gradientMicroDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation + delta, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation - delta, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdGradChi )       == floatVector( 9 * 27, 0 )                                 , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidGradChi )     == floatVector( 9 * 27, 0 )                                 , CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidGradChi ), *R.get_dUpdatedPlasticGradientMicroDeformationdGradChi( ), 2e-5, 3e-6 ) );

    floatVector previousStress = tardigradeVectorTools::appendVectors( { hydra.elasticity.previousPK2, hydra.elasticity.previousSIGMA, hydra.elasticity.previousM } );

    for ( unsigned int i = 0; i < previousStress.size( ); i++ ){

        floatVector delta( previousStress.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStress[ i ] ) + eps;

        floatVector deltaPK2(   delta.begin( ) +  0, delta.begin( ) +  9 );

        floatVector deltaSIGMA( delta.begin( ) +  9, delta.begin( ) + 18 );

        floatVector deltaM(     delta.begin( ) + 18, delta.begin( ) + 45 );

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydrap._local_deltaPK2   = deltaPK2;

        hydrap._local_deltaSIGMA = deltaSIGMA;

        hydrap._local_deltaM     = deltaM;

        hydram._local_deltaPK2   = -deltaPK2;

        hydram._local_deltaSIGMA = -deltaSIGMA;

        hydram._local_deltaM     = -deltaM;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdPreviousStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidPreviousStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidPreviousStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_dUpdatedPlasticFdPreviousStress(        9, floatVector( previousStress.size( ), 0 ) );

    floatMatrix assembled_dUpdatedPlasticChidPreviousStress(      9, floatVector( previousStress.size( ), 0 ) );

    floatMatrix assembled_dUpdatedPlasticGradChidPreviousStress( 27, floatVector( previousStress.size( ), 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dUpdatedPlasticFdPreviousStress[ i ][ j ] = ( *R.get_dUpdatedPlasticDeformationGradientdPreviousMacroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticFdPreviousStress[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticDeformationGradientdPreviousMicroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticChidPreviousStress[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticMicroDeformationdPreviousMicroStress( ) )[ 9 * i + j ];

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dUpdatedPlasticGradChidPreviousStress[ i ][ j + 0 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousMacroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticGradChidPreviousStress[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousMicroStress( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            assembled_dUpdatedPlasticGradChidPreviousStress[ i ][ j + 18 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousHigherOrderStress( ) )[ 27 * i + j ];

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdPreviousStress       ), tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticFdPreviousStress )      , 1e-5, 1e-6 ) );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidPreviousStress     ) == tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticChidPreviousStress )    , CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidPreviousStress ), tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticGradChidPreviousStress ), 1e-5, 1e-5 ) );

    for ( unsigned int i = 0; i < previousStateVariables.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables - delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdPreviousISVs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidPreviousISVs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidPreviousISVs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_dUpdatedPlasticFdPreviousISVs(        9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix assembled_dUpdatedPlasticChidPreviousISVs(      9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix assembled_dUpdatedPlasticGradChidPreviousISVs( 27, floatVector( previousStateVariables.size( ), 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dUpdatedPlasticFdPreviousISVs[ i ][ j ] = ( *R.get_dUpdatedPlasticDeformationGradientdPreviousFn( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticChidPreviousISVs[ i ][ j ] = ( *R.get_dUpdatedPlasticMicroDeformationdPreviousFn( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticChidPreviousISVs[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticMicroDeformationdPreviousChin( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dUpdatedPlasticFdPreviousISVs[ i ][ j + configuration_unknown_count ] = ( *R.get_dUpdatedPlasticDeformationGradientdPreviousStateVariables( ) )[ 10 * i + j ];

            assembled_dUpdatedPlasticChidPreviousISVs[ i ][ j + configuration_unknown_count ] = ( *R.get_dUpdatedPlasticMicroDeformationdPreviousStateVariables( ) )[ 10 * i +  j ];

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dUpdatedPlasticGradChidPreviousISVs[ i ][ j ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousFn( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticGradChidPreviousISVs[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousChin( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            assembled_dUpdatedPlasticGradChidPreviousISVs[ i ][ j + 18 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousGradChin( ) )[ 27 * i + j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dUpdatedPlasticGradChidPreviousISVs[ i ][ j + configuration_unknown_count ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousStateVariables( ) )[ 10 * i +  j ];

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdPreviousISVs       ) == tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticFdPreviousISVs       ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidPreviousISVs     ) == tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticChidPreviousISVs     ), CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidPreviousISVs ), tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticGradChidPreviousISVs ), 3e-5, 3e-6 ) );

    for ( unsigned int i = 0; i < previousDeformationGradient.size( ); i++ ){

        floatVector delta( previousDeformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdPreviousF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidPreviousF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidPreviousF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdPreviousF       ) == *R.get_dUpdatedPlasticDeformationGradientdPreviousF( )     , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidPreviousF     ) == *R.get_dUpdatedPlasticMicroDeformationdPreviousF( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidPreviousF ) == *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < previousMicroDeformation.size( ); i++ ){

        floatVector delta( previousMicroDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation + delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation - delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdPreviousChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidPreviousChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidPreviousChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdPreviousChi       ) == floatVector( 9 * 9, 0 )                                      , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidPreviousChi     ) == *R.get_dUpdatedPlasticMicroDeformationdPreviousChi( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidPreviousChi ), *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousChi( ), 1e-5, 1e-6 ) );

    for ( unsigned int i = 0; i < previousGradientMicroDeformation.size( ); i++ ){

        floatVector delta( previousGradientMicroDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousGradientMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation + delta,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation - delta,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdPreviousGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidPreviousGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidPreviousGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdPreviousGradChi       ) == floatVector( 9 * 27, 0 )                                         , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidPreviousGradChi     ) == floatVector( 9 * 27, 0 )                                         , CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidPreviousGradChi ), *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousGradChi( ), 3e-5, 3e-6 ) );

}

BOOST_AUTO_TEST_CASE( test_setPlasticDeformation3, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    R.get_dUpdatedPlasticDeformationGradientdPreviousMacroStress( );

    // Test the jacobians

    floatType eps = 1e-6;

    floatMatrix dUpdatedPlasticFdX(                      9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dUpdatedPlasticChidX(                    9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dUpdatedPlasticGradChidX(               27, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dUpdatedPlasticFdF(                      9, floatVector( 9, 0 ) );

    floatMatrix dUpdatedPlasticChidF(                    9, floatVector( 9, 0 ) );

    floatMatrix dUpdatedPlasticGradChidF(               27, floatVector( 9, 0 ) );

    floatMatrix dUpdatedPlasticFdChi(                    9, floatVector( 9, 0 ) );

    floatMatrix dUpdatedPlasticChidChi(                  9, floatVector( 9, 0 ) );

    floatMatrix dUpdatedPlasticGradChidChi(             27, floatVector( 9, 0 ) );

    floatMatrix dUpdatedPlasticFdGradChi(                9, floatVector( 27, 0 ) );

    floatMatrix dUpdatedPlasticChidGradChi(              9, floatVector( 27, 0 ) );

    floatMatrix dUpdatedPlasticGradChidGradChi(         27, floatVector( 27, 0 ) );

    floatMatrix dUpdatedPlasticFdPreviousStress(         9, floatVector( configuration_unknown_count, 0 ) );

    floatMatrix dUpdatedPlasticChidPreviousStress(       9, floatVector( configuration_unknown_count, 0 ) );

    floatMatrix dUpdatedPlasticGradChidPreviousStress(  27, floatVector( configuration_unknown_count, 0 ) );

    floatMatrix dUpdatedPlasticFdPreviousISVs(           9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix dUpdatedPlasticChidPreviousISVs(         9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix dUpdatedPlasticGradChidPreviousISVs(    27, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix dUpdatedPlasticFdPreviousF(              9, floatVector(  9, 0 ) );

    floatMatrix dUpdatedPlasticChidPreviousF(            9, floatVector(  9, 0 ) );

    floatMatrix dUpdatedPlasticGradChidPreviousF(       27, floatVector(  9, 0 ) );

    floatMatrix dUpdatedPlasticFdPreviousChi(            9, floatVector(  9, 0 ) );

    floatMatrix dUpdatedPlasticChidPreviousChi(          9, floatVector(  9, 0 ) );

    floatMatrix dUpdatedPlasticGradChidPreviousChi(     27, floatVector(  9, 0 ) );

    floatMatrix dUpdatedPlasticFdPreviousGradChi(        9, floatVector( 27, 0 ) );

    floatMatrix dUpdatedPlasticChidPreviousGradChi(      9, floatVector( 27, 0 ) );

    floatMatrix dUpdatedPlasticGradChidPreviousGradChi( 27, floatVector( 27, 0 ) );

    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidX[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_dUpdatedPlasticFdX(        9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix assembled_dUpdatedPlasticChidX(      9, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix assembled_dUpdatedPlasticGradChidX( 27, floatVector( unknownVector.size( ), 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dUpdatedPlasticFdX[ i ][ j ] = ( *R.get_dUpdatedPlasticDeformationGradientdMacroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticFdX[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticDeformationGradientdMicroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticChidX[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticMicroDeformationdMicroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticFdX[ i ][ j + configuration_unknown_count ] = ( *R.get_dUpdatedPlasticDeformationGradientdFn( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticChidX[ i ][ j + configuration_unknown_count ] = ( *R.get_dUpdatedPlasticMicroDeformationdFn( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticChidX[ i ][ j + configuration_unknown_count + 9 ] = ( *R.get_dUpdatedPlasticMicroDeformationdChin( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dUpdatedPlasticFdX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dUpdatedPlasticDeformationGradientdStateVariables( ) )[ 10 * i + j ];

            assembled_dUpdatedPlasticChidX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dUpdatedPlasticMicroDeformationdStateVariables( ) )[ 10 * i + j ];

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dUpdatedPlasticGradChidX[ i ][ j + 0 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdMacroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticGradChidX[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdMicroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticGradChidX[ i ][ j + configuration_unknown_count ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdFn( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticGradChidX[ i ][ j + configuration_unknown_count + 9 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdChin( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            assembled_dUpdatedPlasticGradChidX[ i ][ j + 18 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdHigherOrderStress( ) )[ 27 * i + j ];

            assembled_dUpdatedPlasticGradChidX[ i ][ j + configuration_unknown_count + 18 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdGradChin( ) )[ 27 * i + j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dUpdatedPlasticGradChidX[ i ][ j + 2 * configuration_unknown_count ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdStateVariables( ) )[ 10 * i + j ];

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdX       ), tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticFdX       ), 1e-5, 1e-6 ) );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidX     ) == tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticChidX     ), CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidX ), tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticGradChidX ), 4e-5, 1e-5 ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdF )       == *R.get_dUpdatedPlasticDeformationGradientdF( )     , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidF )     == *R.get_dUpdatedPlasticMicroDeformationdF( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidF ) == *R.get_dUpdatedPlasticGradientMicroDeformationdF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < microDeformation.size( ); i++ ){

        floatVector delta( microDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( microDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation + delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation - delta, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdChi       ) == floatVector( 9 * 9, 0 )                              , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidChi     ) == *R.get_dUpdatedPlasticMicroDeformationdChi( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidChi ), *R.get_dUpdatedPlasticGradientMicroDeformationdChi( ), 1e-5, 1e-6 ) );

    for ( unsigned int i = 0; i < gradientMicroDeformation.size( ); i++ ){

        floatVector delta( gradientMicroDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( gradientMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation + delta, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation - delta, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdGradChi       ) == floatVector( 9 * 27, 0 )                                 , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidGradChi     ) == floatVector( 9 * 27, 0 )                                 , CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidGradChi ), *R.get_dUpdatedPlasticGradientMicroDeformationdGradChi( ), 2e-5, 3e-6 ) );

    floatVector previousStress = tardigradeVectorTools::appendVectors( { hydra.elasticity.previousPK2, hydra.elasticity.previousSIGMA, hydra.elasticity.previousM } );

    for ( unsigned int i = 0; i < previousStress.size( ); i++ ){

        floatVector delta( previousStress.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStress[ i ] ) + eps;

        floatVector deltaPK2(   delta.begin( ) +  0, delta.begin( ) +  9 );

        floatVector deltaSIGMA( delta.begin( ) +  9, delta.begin( ) + 18 );

        floatVector deltaM(     delta.begin( ) + 18, delta.begin( ) + 45 );

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydrap._local_deltaPK2   = deltaPK2;

        hydrap._local_deltaSIGMA = deltaSIGMA;

        hydrap._local_deltaM     = deltaM;

        hydram._local_deltaPK2   = -deltaPK2;

        hydram._local_deltaSIGMA = -deltaSIGMA;

        hydram._local_deltaM     = -deltaM;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdPreviousStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidPreviousStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidPreviousStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_dUpdatedPlasticFdPreviousStress(        9, floatVector( previousStress.size( ), 0 ) );

    floatMatrix assembled_dUpdatedPlasticChidPreviousStress(      9, floatVector( previousStress.size( ), 0 ) );

    floatMatrix assembled_dUpdatedPlasticGradChidPreviousStress( 27, floatVector( previousStress.size( ), 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dUpdatedPlasticFdPreviousStress[ i ][ j ] = ( *R.get_dUpdatedPlasticDeformationGradientdPreviousMacroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticFdPreviousStress[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticDeformationGradientdPreviousMicroStress( ) )[ 9 * i  + j ];

            assembled_dUpdatedPlasticChidPreviousStress[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticMicroDeformationdPreviousMicroStress( ) )[ 9 * i + j ];

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dUpdatedPlasticGradChidPreviousStress[ i ][ j + 0 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousMacroStress( ) )[ 9 * i + j ];

            assembled_dUpdatedPlasticGradChidPreviousStress[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousMicroStress( ) )[ 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            assembled_dUpdatedPlasticGradChidPreviousStress[ i ][ j + 18 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousHigherOrderStress( ) )[ 27 * i + j ];

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdPreviousStress       ), tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticFdPreviousStress )      , 1e-5, 1e-6 ) );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidPreviousStress     ) == tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticChidPreviousStress )    , CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidPreviousStress ), tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticGradChidPreviousStress ), 1e-5, 1e-5 ) );

    for ( unsigned int i = 0; i < previousStateVariables.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables - delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdPreviousISVs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidPreviousISVs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidPreviousISVs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    floatMatrix assembled_dUpdatedPlasticFdPreviousISVs(        9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix assembled_dUpdatedPlasticChidPreviousISVs(      9, floatVector( previousStateVariables.size( ), 0 ) );

    floatMatrix assembled_dUpdatedPlasticGradChidPreviousISVs( 27, floatVector( previousStateVariables.size( ), 0 ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dUpdatedPlasticFdPreviousISVs[ i ][ j ] = ( *R.get_dUpdatedPlasticDeformationGradientdPreviousFn( ) )[ ( numConfigurations - 1 ) * 9 * i + j ];

            assembled_dUpdatedPlasticChidPreviousISVs[ i ][ j ] = ( *R.get_dUpdatedPlasticMicroDeformationdPreviousFn( ) )[ ( numConfigurations - 1 ) * 9 * i + j ];

            assembled_dUpdatedPlasticChidPreviousISVs[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticMicroDeformationdPreviousChin( ) )[ ( numConfigurations - 1 ) * 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dUpdatedPlasticFdPreviousISVs[ i ][ j + configuration_unknown_count ] = ( *R.get_dUpdatedPlasticDeformationGradientdPreviousStateVariables( ) )[ 10 * i + j ];

            assembled_dUpdatedPlasticChidPreviousISVs[ i ][ j + configuration_unknown_count ] = ( *R.get_dUpdatedPlasticMicroDeformationdPreviousStateVariables( ) )[ 10 * i + j ];

        }

    }

    for ( unsigned int i = 0; i < 27; i++ ){

        for ( unsigned int j = 0; j < 9; j++ ){

            assembled_dUpdatedPlasticGradChidPreviousISVs[ i ][ j ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousFn( ) )[ ( numConfigurations - 1 ) * 9 * i + j ];

            assembled_dUpdatedPlasticGradChidPreviousISVs[ i ][ j + 9 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousChin( ) )[ ( numConfigurations - 1 ) * 9 * i + j ];

        }

        for ( unsigned int j = 0; j < 27; j++ ){

            assembled_dUpdatedPlasticGradChidPreviousISVs[ i ][ j + 18 ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousGradChin( ) )[ ( numConfigurations - 1 ) * 27 * i + j ];

        }

        for ( unsigned int j = 0; j < 10; j++ ){

            assembled_dUpdatedPlasticGradChidPreviousISVs[ i ][ j + configuration_unknown_count ] = ( *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousStateVariables( ) )[ 10 * i + j ];

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdPreviousISVs       ), tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticFdPreviousISVs )      , 1e-5, 1e-6 ) );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidPreviousISVs     ), tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticChidPreviousISVs )    , 1e-5, 1e-6 ) );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidPreviousISVs ), tardigradeVectorTools::appendVectors( assembled_dUpdatedPlasticGradChidPreviousISVs ), 3e-5, 2e-6 ) );

    for ( unsigned int i = 0; i < previousDeformationGradient.size( ); i++ ){

        floatVector delta( previousDeformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdPreviousF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidPreviousF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidPreviousF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdPreviousF )       == *R.get_dUpdatedPlasticDeformationGradientdPreviousF( )     , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidPreviousF )     == *R.get_dUpdatedPlasticMicroDeformationdPreviousF( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidPreviousF ) == *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < previousMicroDeformation.size( ); i++ ){

        floatVector delta( previousMicroDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation + delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation - delta, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdPreviousChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidPreviousChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidPreviousChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdPreviousChi       ) == floatVector( 9 * 9, 0 )                                      , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidPreviousChi     ) == *R.get_dUpdatedPlasticMicroDeformationdPreviousChi( )        , CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidPreviousChi ), *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousChi( ), 1e-5, 1e-6 ) );

    for ( unsigned int i = 0; i < previousGradientMicroDeformation.size( ); i++ ){

        floatVector delta( previousGradientMicroDeformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousGradientMicroDeformation[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation + delta,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation - delta,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_updatedPlasticDeformationGradient( );

        floatVector vm = *Rm.get_updatedPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticFdPreviousGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticMicroDeformation( );

        vm = *Rm.get_updatedPlasticMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticChidPreviousGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

        vp = *Rp.get_updatedPlasticGradientMicroDeformation( );

        vm = *Rm.get_updatedPlasticGradientMicroDeformation( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dUpdatedPlasticGradChidPreviousGradChi[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticFdPreviousGradChi       ) == floatVector( 9 * 27, 0 )                                         , CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dUpdatedPlasticChidPreviousGradChi     ) == floatVector( 9 * 27, 0 )                                         , CHECK_PER_ELEMENT );

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dUpdatedPlasticGradChidPreviousGradChi ), *R.get_dUpdatedPlasticGradientMicroDeformationdPreviousGradChi( ), 1e-4, 1e-6 ) );

}

BOOST_AUTO_TEST_CASE( test_setStateVariableResiduals, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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

            floatType macroYield = 1.23;

            floatType microYield = 2.34;

            floatVector microGradientYield = { 3.45, 4.56, 5.67 };

            floatVector updatedPlasticStrainLikeISVs = { 0.123, 0.234, 0.345, 0.456, 0.567 };

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

            }

            virtual void setYield( const bool isPrevious ){


                if ( isPrevious ){

                    // We don't need this part

                }
                else{

                    set_macroYield( macroYield );

                    set_microYield( microYield );

                    set_microGradientYield( microGradientYield );

                }

            }

            virtual void setUpdatedPlasticStrainLikeISVs( ) override{

                set_updatedPlasticStrainLikeISVs( updatedPlasticStrainLikeISVs );

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    residualMock R2( &hydra, 55, 1, stateVariableIndices, parameters );

    R2.macroYield = -R.macroYield;

    R2.microYield = -R.microYield;

    R2.microGradientYield = -R.microGradientYield;

    floatVector answer1( 10, 0 );
    answer1[ 0 ] = R.macroYield              + unknownVector[ 2 * configuration_unknown_count + 0 ] * R.macroYield;
    answer1[ 1 ] = R.microYield              + unknownVector[ 2 * configuration_unknown_count + 1 ] * R.microYield;
    answer1[ 2 ] = R.microGradientYield[ 0 ] + unknownVector[ 2 * configuration_unknown_count + 2 ] * R.microGradientYield[ 0 ];
    answer1[ 3 ] = R.microGradientYield[ 1 ] + unknownVector[ 2 * configuration_unknown_count + 3 ] * R.microGradientYield[ 1 ];
    answer1[ 4 ] = R.microGradientYield[ 2 ] + unknownVector[ 2 * configuration_unknown_count + 4 ] * R.microGradientYield[ 2 ];
    answer1[ 5 ] = R.updatedPlasticStrainLikeISVs[ 0 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 0 ];
    answer1[ 6 ] = R.updatedPlasticStrainLikeISVs[ 1 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 1 ];
    answer1[ 7 ] = R.updatedPlasticStrainLikeISVs[ 2 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 2 ];
    answer1[ 8 ] = R.updatedPlasticStrainLikeISVs[ 3 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 3 ];
    answer1[ 9 ] = R.updatedPlasticStrainLikeISVs[ 4 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 4 ];

    floatVector answer2( 10, 0 );
    answer2[ 0 ] = unknownVector[ 2 * configuration_unknown_count + 0 ] * R2.macroYield;
    answer2[ 1 ] = unknownVector[ 2 * configuration_unknown_count + 1 ] * R2.microYield;
    answer2[ 2 ] = unknownVector[ 2 * configuration_unknown_count + 2 ] * R2.microGradientYield[ 0 ];
    answer2[ 3 ] = unknownVector[ 2 * configuration_unknown_count + 3 ] * R2.microGradientYield[ 1 ];
    answer2[ 4 ] = unknownVector[ 2 * configuration_unknown_count + 4 ] * R2.microGradientYield[ 2 ];
    answer2[ 5 ] = R2.updatedPlasticStrainLikeISVs[ 0 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 0 ];
    answer2[ 6 ] = R2.updatedPlasticStrainLikeISVs[ 1 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 1 ];
    answer2[ 7 ] = R2.updatedPlasticStrainLikeISVs[ 2 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 2 ];
    answer2[ 8 ] = R2.updatedPlasticStrainLikeISVs[ 3 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 3 ];
    answer2[ 9 ] = R2.updatedPlasticStrainLikeISVs[ 4 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 4 ];

    BOOST_TEST( answer1 == *R.get_stateVariableResiduals( ) , CHECK_PER_ELEMENT );

    BOOST_TEST( answer2 == *R2.get_stateVariableResiduals( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setStateVariableJacobians, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    floatMatrix jacobian( 10, floatVector( unknownVector.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_stateVariableResiduals( );

        floatVector vm = *Rm.get_stateVariableResiduals( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            jacobian[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( jacobian ), *R.get_stateVariableJacobians( ), 1e-5, 1e-6 ) );

}

BOOST_AUTO_TEST_CASE( test_setStateVariableJacobians2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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
                                  -0.01, -0.02, -0.03, -0.04, -0.05,
                                   0.06,  0.07,  0.08,  0.09,  0.10 };

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

                plasticity.setUseProjection( false );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                setResidualClasses( residuals );

            }

    };

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    R.setUseProjection( false );

    floatMatrix jacobian( 10, floatVector( unknownVector.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        Rp.setUseProjection( false );

        Rm.setUseProjection( false );

        floatVector vp = *Rp.get_stateVariableResiduals( );

        floatVector vm = *Rm.get_stateVariableResiduals( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            jacobian[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( jacobian ), *R.get_stateVariableJacobians( ), 1e-5, 1e-6 ) );

}

BOOST_AUTO_TEST_CASE( test_setdStateVariableResidualsdD, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    floatMatrix dRdD( 10, floatVector( 45, 0 ) );

    floatType eps = 1e-6;

    floatVector deformation = tardigradeVectorTools::appendVectors( { deformationGradient, microDeformation, gradientMicroDeformation } );

    for ( unsigned int i = 0; i < deformation.size( ); i++ ){

        floatVector delta( deformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformation[ i ] ) + eps;

        floatVector delta_F(       delta.begin( ) +  0, delta.begin( ) +  9 );

        floatVector delta_Chi(     delta.begin( ) +  9, delta.begin( ) + 18 );

        floatVector delta_GradChi( delta.begin( ) + 18, delta.begin( ) +  45 );

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta_F, previousDeformationGradient,
                                          microDeformation + delta_Chi, previousMicroDeformation, gradientMicroDeformation + delta_GradChi, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta_F, previousDeformationGradient,
                                          microDeformation - delta_Chi, previousMicroDeformation, gradientMicroDeformation - delta_GradChi, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_stateVariableResiduals( );

        floatVector vm = *Rm.get_stateVariableResiduals( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dRdD[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dRdD ) == *R.get_dStateVariableResidualsdD( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setdStateVariableResidualsdPreviousISVs, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    floatMatrix dRdPreviousISVs( 10, floatVector( previousStateVariables.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < previousStateVariables.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables - delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_stateVariableResiduals( );

        floatVector vm = *Rm.get_stateVariableResiduals( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dRdPreviousISVs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dRdPreviousISVs ) == *R.get_dStateVariableResidualsdPreviousISVs( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setStateVariableResiduals_weak, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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

            floatType macroYield = 1.23;

            floatType microYield = 2.34;

            floatVector microGradientYield = { 3.45, 4.56, 5.67 };

            floatVector updatedPlasticStrainLikeISVs = { 0.123, 0.234, 0.345, 0.456, 0.567 };

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

            }

            virtual void setYield( const bool isPrevious ){


                if ( isPrevious ){

                    // We don't need this part

                }
                else{

                    set_macroYield( macroYield );

                    set_microYield( microYield );

                    set_microGradientYield( microGradientYield );

                }

            }

            virtual void setUpdatedPlasticStrainLikeISVs( ) override{

                set_updatedPlasticStrainLikeISVs( updatedPlasticStrainLikeISVs );

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

                plasticity = residualMock( this, 55, 1, stateVariableIndices, plasticParameters, 0.5, true, 0.75 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                setResidualClasses( residuals );

            }

    };

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatType a = 0.75;

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters, 0.5, true, a );

    residualMock R2( &hydra, 55, 1, stateVariableIndices, parameters, 0.5, true, a );

    R2.macroYield = -R.macroYield;

    R2.microYield = -R.microYield;

    R2.microGradientYield = -R.microGradientYield;

    floatVector answer1( 10, 0 );
    answer1[ 0 ] = tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( R.macroYield, a )
                 + unknownVector[ 2 * configuration_unknown_count + 0 ] * R.macroYield
                 + 1000 * tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( -unknownVector[ 2 * configuration_unknown_count + 0 ], a );
    answer1[ 1 ] = tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( R.microYield, a )
                 + unknownVector[ 2 * configuration_unknown_count + 1 ] * R.microYield
                 + 1000 * tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( -unknownVector[ 2 * configuration_unknown_count + 1 ], a );
    answer1[ 2 ] = tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( R.microGradientYield[ 0 ], a )
                 + unknownVector[ 2 * configuration_unknown_count + 2 ] * R.microGradientYield[ 0 ]
                 + 1000 * tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( -unknownVector[ 2 * configuration_unknown_count + 2 ], a );
    answer1[ 3 ] = tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( R.microGradientYield[ 1 ], a )
                 + unknownVector[ 2 * configuration_unknown_count + 3 ] * R.microGradientYield[ 1 ]
                 + 1000 * tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( -unknownVector[ 2 * configuration_unknown_count + 3 ], a );
    answer1[ 4 ] = tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( R.microGradientYield[ 2 ], a )
                 + unknownVector[ 2 * configuration_unknown_count + 4 ] * R.microGradientYield[ 2 ]
                 + 1000 * tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( -unknownVector[ 2 * configuration_unknown_count + 4 ], a );
    answer1[ 5 ] = R.updatedPlasticStrainLikeISVs[ 0 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 0 ];
    answer1[ 6 ] = R.updatedPlasticStrainLikeISVs[ 1 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 1 ];
    answer1[ 7 ] = R.updatedPlasticStrainLikeISVs[ 2 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 2 ];
    answer1[ 8 ] = R.updatedPlasticStrainLikeISVs[ 3 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 3 ];
    answer1[ 9 ] = R.updatedPlasticStrainLikeISVs[ 4 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 4 ];

    floatVector answer2( 10, 0 );
    answer2[ 0 ] = tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( R2.macroYield, a )
                 + unknownVector[ 2 * configuration_unknown_count + 0 ] * R2.macroYield
                 + 1000 * tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( -unknownVector[ 2 * configuration_unknown_count + 0 ], a );
    answer2[ 1 ] = tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( R2.microYield, a )
                 + unknownVector[ 2 * configuration_unknown_count + 1 ] * R2.microYield
                 + 1000 * tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( -unknownVector[ 2 * configuration_unknown_count + 1 ], a );
    answer2[ 2 ] = tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( R2.microGradientYield[ 0 ], a )
                 + unknownVector[ 2 * configuration_unknown_count + 2 ] * R2.microGradientYield[ 0 ]
                 + 1000 * tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( -unknownVector[ 2 * configuration_unknown_count + 2 ], a );
    answer2[ 3 ] = tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( R2.microGradientYield[ 1 ], a )
                 + unknownVector[ 2 * configuration_unknown_count + 3 ] * R2.microGradientYield[ 1 ]
                 + 1000 * tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( -unknownVector[ 2 * configuration_unknown_count + 3 ], a );
    answer2[ 4 ] = tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( R2.microGradientYield[ 2 ], a )
                 + unknownVector[ 2 * configuration_unknown_count + 4 ] * R2.microGradientYield[ 2 ]
                 + 1000 * tardigradeHydra::micromorphicDruckerPragerPlasticity::weakMac( -unknownVector[ 2 * configuration_unknown_count + 4 ], a );
    answer2[ 5 ] = R2.updatedPlasticStrainLikeISVs[ 0 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 0 ];
    answer2[ 6 ] = R2.updatedPlasticStrainLikeISVs[ 1 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 1 ];
    answer2[ 7 ] = R2.updatedPlasticStrainLikeISVs[ 2 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 2 ];
    answer2[ 8 ] = R2.updatedPlasticStrainLikeISVs[ 3 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 3 ];
    answer2[ 9 ] = R2.updatedPlasticStrainLikeISVs[ 4 ] - unknownVector[ 2 * configuration_unknown_count + 5 + 4 ];

    BOOST_TEST( answer1 == *R.get_stateVariableResiduals( ) , CHECK_PER_ELEMENT );

    BOOST_TEST( answer2 == *R2.get_stateVariableResiduals( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setStateVariableJacobians_weak, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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

            floatType a = 0.75;

        private:

            using tardigradeHydra::hydraBaseMicromorphic::setResidualClasses;

            virtual void setResidualClasses( ) override{

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = stressMock( this, 45 );

                elasticity.previousPK2   += _local_deltaPK2;

                elasticity.previousSIGMA += _local_deltaSIGMA;

                elasticity.previousM     += _local_deltaM;

                plasticity = residualMock( this, 55, 1, stateVariableIndices, plasticParameters, 0.5, true, a );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                setResidualClasses( residuals );

            }

    };

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatType a = 0.75;

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters, 0.5, true, a );

    floatMatrix jacobian( 10, floatVector( unknownVector.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters, 0.5, true, a );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters, 0.5, true, a );

        floatVector vp = *Rp.get_stateVariableResiduals( );

        floatVector vm = *Rm.get_stateVariableResiduals( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            jacobian[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( jacobian ), *R.get_stateVariableJacobians( ), 2e-4, 1e-6 ) );

}

BOOST_AUTO_TEST_CASE( test_setStateVariableJacobians_weak2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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
                                  -0.01, -0.02, -0.03, -0.04, -0.05,
                                   0.06,  0.07,  0.08,  0.09,  0.10 };

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

            floatType a = 0.75;

        private:

            using tardigradeHydra::hydraBaseMicromorphic::setResidualClasses;

            virtual void setResidualClasses( ) override{

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = stressMock( this, 45 );

                elasticity.previousPK2   += _local_deltaPK2;

                elasticity.previousSIGMA += _local_deltaSIGMA;

                elasticity.previousM     += _local_deltaM;

                plasticity = residualMock( this, 55, 1, stateVariableIndices, plasticParameters, 0.5, true, a );

                plasticity.setUseProjection( false );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                setResidualClasses( residuals );

            }

    };

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatType a = 0.75;

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters, 0.5, true, a );

    R.setUseProjection( false );

    floatMatrix jacobian( 10, floatVector( unknownVector.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters, 0.5, true, a );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters, 0.5, true, a );

        Rp.setUseProjection( false );

        Rm.setUseProjection( false );

        floatVector vp = *Rp.get_stateVariableResiduals( );

        floatVector vm = *Rm.get_stateVariableResiduals( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            jacobian[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( jacobian ), *R.get_stateVariableJacobians( ), 1e-3, 3e-6 ) ); //TODO: This seems like a high tolerance

}

BOOST_AUTO_TEST_CASE( test_setdStateVariableResidualsdD_weak, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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

                plasticity = residualMock( this, 55, 1, stateVariableIndices, plasticParameters, 0.5, true, 0.75 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                setResidualClasses( residuals );

            }

    };

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatType a = 0.75;

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters, 0.5, true, a );

    floatMatrix dRdD( 10, floatVector( 45, 0 ) );

    floatType eps = 1e-6;

    floatVector deformation = tardigradeVectorTools::appendVectors( { deformationGradient, microDeformation, gradientMicroDeformation } );

    for ( unsigned int i = 0; i < deformation.size( ); i++ ){

        floatVector delta( deformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformation[ i ] ) + eps;

        floatVector delta_F(       delta.begin( ) +  0, delta.begin( ) +  9 );

        floatVector delta_Chi(     delta.begin( ) +  9, delta.begin( ) + 18 );

        floatVector delta_GradChi( delta.begin( ) + 18, delta.begin( ) +  45 );

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta_F, previousDeformationGradient,
                                          microDeformation + delta_Chi, previousMicroDeformation, gradientMicroDeformation + delta_GradChi, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta_F, previousDeformationGradient,
                                          microDeformation - delta_Chi, previousMicroDeformation, gradientMicroDeformation - delta_GradChi, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters, 0.5, true, a );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters, 0.5, true, a );

        floatVector vp = *Rp.get_stateVariableResiduals( );

        floatVector vm = *Rm.get_stateVariableResiduals( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dRdD[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dRdD ), *R.get_dStateVariableResidualsdD( ), 1e-5, 1e-6 ) );

}

BOOST_AUTO_TEST_CASE( test_setdStateVariableResidualsdPreviousISVs_weak, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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

                plasticity = residualMock( this, 55, 1, stateVariableIndices, plasticParameters, 0.5, true, 0.75 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                setResidualClasses( residuals );

            }

    };

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatType a = 0.75;

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters, 0.5, true, a );

    floatMatrix dRdPreviousISVs( 10, floatVector( previousStateVariables.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < previousStateVariables.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables + delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables - delta, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters, 0.5, true, a );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters, 0.5, true, a );

        floatVector vp = *Rp.get_stateVariableResiduals( );

        floatVector vm = *Rm.get_stateVariableResiduals( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dRdPreviousISVs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dRdPreviousISVs ) == *R.get_dStateVariableResidualsdPreviousISVs( ), CHECK_PER_ELEMENT );

}
BOOST_AUTO_TEST_CASE( test_setResidual, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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

            floatVector updatedPlasticDeformationGradient      = { 1.1, 0.11, 0.12, 0.13, 0.94, 0.15, 0.16, 0.17, 0.88 };

            floatVector updatedPlasticMicroDeformation         = { 0.9, 0.02, 0.03, 0.14, 1.15, 0.06, 0.70, 0.08, 1.99 };

            floatVector updatedPlasticGradientMicroDeformation = { 0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                                                   1.11, 1.22, 1.44, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99,
                                                                   2.11, 2.22, 2.33, 2.44, 2.55, 2.66, 2.77, 2.88, 2.99 }; 

            floatVector stateVariableResiduals = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };

            floatVector initialize( unsigned int nrows ){

                floatVector value( nrows, 0 );

                return value;

            }

            floatMatrix initialize( unsigned int nrows, unsigned int ncols ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                return value;

            }

            virtual void setPlasticDeformation( ) override{

                set_updatedPlasticDeformationGradient(      updatedPlasticDeformationGradient );

                set_updatedPlasticMicroDeformation(         updatedPlasticMicroDeformation    );

                set_updatedPlasticGradientMicroDeformation( updatedPlasticGradientMicroDeformation );

            }

            virtual void setStateVariableResiduals( ) override{

                set_stateVariableResiduals( stateVariableResiduals );

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    floatVector answer = tardigradeVectorTools::appendVectors( { R.updatedPlasticDeformationGradient         - floatVector( unknownVector.begin( ) + configuration_unknown_count +  0,
                                                                                                                            unknownVector.begin( ) + configuration_unknown_count +  9 ),
                                                                 R.updatedPlasticMicroDeformation            - floatVector( unknownVector.begin( ) + configuration_unknown_count +  9,
                                                                                                                            unknownVector.begin( ) + configuration_unknown_count + 18 ),
                                                                 R.updatedPlasticGradientMicroDeformation    - floatVector( unknownVector.begin( ) + configuration_unknown_count + 18,
                                                                                                                            unknownVector.begin( ) + configuration_unknown_count + 45 ),
                                                                 R.stateVariableResiduals } );

    BOOST_TEST( answer == *R.getResidual( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    floatMatrix jacobian( 55, floatVector( unknownVector.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.getResidual( );

        floatVector vm = *Rm.getResidual( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            jacobian[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( jacobian ), *R.getJacobian( ), 2e-5, 5e-6 ) );

}

BOOST_AUTO_TEST_CASE( test_setdRdD, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    floatMatrix dRdD( 55, floatVector( 45, 0 ) );

    floatType eps = 1e-6;

    floatVector deformation = tardigradeVectorTools::appendVectors( { deformationGradient, microDeformation, gradientMicroDeformation } );

    for ( unsigned int i = 0; i < deformation.size( ); i++ ){

        floatVector delta( deformation.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformation[ i ] ) + eps;

        floatVector delta_F(       delta.begin( ) +  0, delta.begin( ) +  9 );

        floatVector delta_Chi(     delta.begin( ) +  9, delta.begin( ) + 18 );

        floatVector delta_GradChi( delta.begin( ) + 18, delta.begin( ) +  45 );

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta_F, previousDeformationGradient,
                                          microDeformation + delta_Chi, previousMicroDeformation, gradientMicroDeformation + delta_GradChi, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta_F, previousDeformationGradient,
                                          microDeformation - delta_Chi, previousMicroDeformation, gradientMicroDeformation - delta_GradChi, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.getResidual( );

        floatVector vm = *Rm.getResidual( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dRdD[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dRdD ), *R.getdRdD( ), 2e-5, 3e-6 ) );

}

BOOST_AUTO_TEST_CASE( test_setdRdT, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the plastic deformation
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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R(  &hydra, 55, 1, stateVariableIndices, parameters );

    floatVector dRdT( 55, 0 );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector delta( 1, 0 );

        delta[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMicromorphicMock hydrap( time, deltaTime, temperature + delta[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        hydraBaseMicromorphicMock hydram( time, deltaTime, temperature - delta[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                                          microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                          { }, { },
                                          previousStateVariables, parameters,
                                          numConfigurations, numNonLinearSolveStateVariables,
                                          dimension, configuration_unknown_count,
                                          tolr, tola, maxIterations, maxLSIterations, lsAlpha );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp(  &hydrap, 55, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 55, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.getResidual( );

        floatVector vm = *Rm.getResidual( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dRdT[ j ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( dRdT == *R.getdRdT( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_projectSuggestedX, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    residualMock R( &hydra, 55, 1, stateVariableIndices, parameters );

    floatVector trialX = { -0.40474133, -1.31581257, -0.24799601, -0.66822701,  1.70799134,
                            0.13089057, -0.17216613, -0.00864461,  0.80232873,  0.13049679,
                            0.3705609 ,  1.04375356,  2.1860939 ,  0.65654621, -0.91634756,
                           -1.54474396, -0.17666782,  1.65885133,  2.29853222,  0.53338233,
                            0.69052704,  0.78619126,  1.18703091,  1.83891347, -1.20737302,
                            0.48385701, -0.75505922, -1.68384194, -0.01766042,  1.91917888,
                            1.34376976, -0.53305317,  0.46286124,  1.48304898,  1.36610329,
                            0.1893202 , -1.42179529, -2.49648211, -0.08692722, -1.2776921 ,
                            2.28801669, -2.00114005,  1.51340049, -2.17481968, -0.96377428,
                           -0.70299002, -1.47433365,  3.44466764, -0.6873039 , -1.30451035,
                           -2.0131308 , -1.62654174,  0.68615533, -0.41187571,  4.54269611,
                           -0.28936713,  1.9595505 ,  3.64415161, -0.21041396, -3.00618013,
                           -3.57395263,  4.55948367,  4.62129283,  0.73171675, -2.97715061,
                            2.36464984,  0.04360128, -3.57926391,  1.66380595, -2.76950485,
                           -2.16234317, -3.68731285, -3.51451628,  2.71405279,  4.41692098,
                            4.03053223,  0.12819252,  2.11006335, -0.13984174, -2.25057865,
                           -3.81652212, -3.42465008, -4.11456265,  4.52170914, -3.62536239,
                           -4.50230405,  3.28010857, -0.13044935, -1.88583842,  0.88096962,
                           -0.11414338,  0.07622781, -0.11073456, -0.48607211, -0.02456945,
                            0.46456126, -0.27655936, -0.383391  ,  0.26888644,  0.35932029 };

    floatVector Xp = { 2.37179207,  2.04571891,  1.43340891,  1.56210501,  0.29650219,
                      -1.81653588, -1.40387454, -1.2059766 ,  1.98287363, -1.17985447,
                      -0.31947747,  1.67487249, -2.2389102 , -1.02632147, -0.04017631,
                       1.53355464,  1.77147411, -0.85353032, -0.54185951, -1.91388754,
                       1.35593112,  2.27108171,  1.20009559, -2.32963703, -1.36833298,
                       1.57008187, -0.68035577,  1.62855286,  2.08180914, -0.22511169,
                       2.23598071,  1.13591538,  1.8935963 ,  1.02065205, -1.52995918,
                       1.86641792, -2.26013904, -1.3023638 , -1.95985108, -2.43277805,
                       1.30597982,  1.15786049, -1.44268755, -2.4201092 , -0.35044801,
                       0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                       0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                       0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                       0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                       0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                       0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                       0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                       0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                       0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                       0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                       0.        ,  0.        ,  0.        ,  0.        ,  0.         };

    floatVector answer = { -0.40474133, -1.31581257, -0.24799601, -0.66822701,  1.70799134,
                            0.13089057, -0.17216613, -0.00864461,  0.80232873,  0.13049679,
                            0.3705609 ,  1.04375356,  2.1860939 ,  0.65654621, -0.91634756,
                           -1.54474396, -0.17666782,  1.65885133,  2.29853222,  0.53338233,
                            0.69052704,  0.78619126,  1.18703091,  1.83891347, -1.20737302,
                            0.48385701, -0.75505922, -1.68384194, -0.01766042,  1.91917888,
                            1.34376976, -0.53305317,  0.46286124,  1.48304898,  1.36610329,
                            0.1893202 , -1.42179529, -2.49648211, -0.08692722, -1.2776921 ,
                            2.28801669, -2.00114005,  1.51340049, -2.17481968, -0.96377428,
                           -0.1433584 , -0.30065593,  0.70245954, -0.14015958, -0.26602443,
                           -0.41053102, -0.33169521,  0.13992536, -0.08399244,  0.45014875,
                           -0.02867422,  0.19417746,  0.3611094 , -0.02085052, -0.29789099,
                           -0.35415319,  0.45181228,  0.45793712,  0.04786126, -0.19473407,
                            0.15467067,  0.00285194, -0.23411803,  0.10882879, -0.18115205,
                           -0.14143788, -0.24118546, -0.22988291,  0.17752496,  0.2889088 ,
                            0.26363528,  0.00838501,  0.13801829, -0.00914698, -0.14720933,
                           -0.24963698, -0.22400481, -0.26913167,  0.29576294, -0.23713331,
                           -0.29449366,  0.21455041, -0.00853263, -0.12335183,  0.05762382,
                            0.        ,  0.07622781,  0.        ,  0.        ,  0.        ,
                            0.46456126,  0.        ,  0.        ,  0.26888644,  0.35932029 };

    R.projectSuggestedX( trialX, Xp );

    BOOST_TEST( trialX == answer, CHECK_PER_ELEMENT );

}
