/**
  * \file test_tardigrade_hydraMicromorphicDruckerPragerPlasticityOptimization.cpp
  *
  * Tests for tardigrade_hydraMicromorphicDruckerPragerPlasticityOptimization
  */

#include<tardigrade_hydraMicromorphicDruckerPragerPlasticityOptimization.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_micromorphic_tools.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_tardigrade_hydraMicromorphicDruckerPragerPlasticityOptimization
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

typedef tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::variableType variableType; //!< Redefinition of the variable type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::variableVector variableVector; //!< Redefinition of the vector of variable types
typedef tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::variableMatrix variableMatrix; //!< Redefinition of the matrix of variable types

typedef tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::parameterType parameterType; //!< Redefinition of the parameter type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::parameterVector parameterVector; //!< Redefinition of the vector of parameters

typedef tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::constantType constantType; //!< Redefinition of the constant type
typedef tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::constantVector constantVector; //!< Redefinition of the vector of constants
typedef tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::constantMatrix constantMatrix; //!< Redefinition of the matrix of constants

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

BOOST_AUTO_TEST_CASE( test_computeStateVariableResidual, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5 };

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
                                   0.06, 0.07, 0.08, 0.09, 0.10,
                                   0.11, 0.12, 0.13, 0.24, 0.15 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 15;

    std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };

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

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::residual::residual;

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

            std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };

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

                plasticity = residualMock( this, 60, 1, stateVariableIndices, plasticParameters );

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

    residualMock R(  &hydra, 60, 1, stateVariableIndices, parameters );

    residualMock R2( &hydra, 60, 1, stateVariableIndices, parameters );

    R2.macroYield = -R.macroYield;

    R2.microYield = -R.microYield;

    R2.microGradientYield = -R.microGradientYield;

    floatVector answer1( 15, 0 );
    answer1[  5 ] = R.updatedPlasticStrainLikeISVs[ 0 ] - unknownVector[ 2 * configuration_unknown_count +  5 + 0 ];
    answer1[  6 ] = R.updatedPlasticStrainLikeISVs[ 1 ] - unknownVector[ 2 * configuration_unknown_count +  5 + 1 ];
    answer1[  7 ] = R.updatedPlasticStrainLikeISVs[ 2 ] - unknownVector[ 2 * configuration_unknown_count +  5 + 2 ];
    answer1[  8 ] = R.updatedPlasticStrainLikeISVs[ 3 ] - unknownVector[ 2 * configuration_unknown_count +  5 + 3 ];
    answer1[  9 ] = R.updatedPlasticStrainLikeISVs[ 4 ] - unknownVector[ 2 * configuration_unknown_count +  5 + 4 ];
    answer1[ 10 ] = -R.macroYield                       - unknownVector[ 2 * configuration_unknown_count + 10 + 0 ];
    answer1[ 11 ] = -R.microYield                       - unknownVector[ 2 * configuration_unknown_count + 10 + 1 ];
    answer1[ 12 ] = -R.microGradientYield[ 0 ]          - unknownVector[ 2 * configuration_unknown_count + 10 + 2 ];
    answer1[ 13 ] = -R.microGradientYield[ 1 ]          - unknownVector[ 2 * configuration_unknown_count + 10 + 3 ];
    answer1[ 14 ] = -R.microGradientYield[ 2 ]          - unknownVector[ 2 * configuration_unknown_count + 10 + 4 ];

    floatVector answer2( 15, 0 );
    answer2[ 5 ] = R2.updatedPlasticStrainLikeISVs[ 0 ] - unknownVector[ 2 * configuration_unknown_count +  5 + 0 ];
    answer2[ 6 ] = R2.updatedPlasticStrainLikeISVs[ 1 ] - unknownVector[ 2 * configuration_unknown_count +  5 + 1 ];
    answer2[ 7 ] = R2.updatedPlasticStrainLikeISVs[ 2 ] - unknownVector[ 2 * configuration_unknown_count +  5 + 2 ];
    answer2[ 8 ] = R2.updatedPlasticStrainLikeISVs[ 3 ] - unknownVector[ 2 * configuration_unknown_count +  5 + 3 ];
    answer2[ 9 ] = R2.updatedPlasticStrainLikeISVs[ 4 ] - unknownVector[ 2 * configuration_unknown_count +  5 + 4 ];
    answer2[ 10 ] = -R2.macroYield                      - unknownVector[ 2 * configuration_unknown_count + 10 + 0 ];
    answer2[ 11 ] = -R2.microYield                      - unknownVector[ 2 * configuration_unknown_count + 10 + 1 ];
    answer2[ 12 ] = -R2.microGradientYield[ 0 ]         - unknownVector[ 2 * configuration_unknown_count + 10 + 2 ];
    answer2[ 13 ] = -R2.microGradientYield[ 1 ]         - unknownVector[ 2 * configuration_unknown_count + 10 + 3 ];
    answer2[ 14 ] = -R2.microGradientYield[ 2 ]         - unknownVector[ 2 * configuration_unknown_count + 10 + 4 ];

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
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5 };

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
                                   0.06, 0.07, 0.08, 0.09, 0.10,
                                   0.11, 0.12, 0.13, 0.14, 0.15 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 15;

    std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };

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

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::residual::residual;

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

            std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };

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

                plasticity = residualMock( this, 60, 1, stateVariableIndices, plasticParameters );

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

    residualMock R(  &hydra, 60, 1, stateVariableIndices, parameters );

    floatMatrix jacobian( 15, floatVector( unknownVector.size( ), 0 ) );

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

        residualMock Rp(  &hydrap, 60, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 60, 1, stateVariableIndices, parameters );

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
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5 };

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
                                   0.06, 0.07, 0.08, 0.09, 0.10,
                                   0.11, 0.12, 0.13, 0.14, 0.15 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 15;

    std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };

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

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::residual::residual;

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

            std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };

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

                plasticity = residualMock( this, 60, 1, stateVariableIndices, plasticParameters );

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

    residualMock R(  &hydra, 60, 1, stateVariableIndices, parameters );

    floatMatrix dRdD( 15, floatVector( 45, 0 ) );

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

        residualMock Rp(  &hydrap, 60, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 60, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.get_stateVariableResiduals( );

        floatVector vm = *Rm.get_stateVariableResiduals( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            dRdD[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dRdD ) == *R.get_dStateVariableResidualsdD( ), CHECK_PER_ELEMENT );

}


BOOST_AUTO_TEST_CASE( test_computeConstraints, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5 };

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
                                   0.06, 0.07, 0.08, 0.09, 0.10,
                                   0.11, 0.12, 0.13, 0.24, 0.15 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 15;

    std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };

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

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::residual::residual;

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

            std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };

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

                plasticity = residualMock( this, 60, 1, stateVariableIndices, plasticParameters );

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

    residualMock R(  &hydra, 60, 1, stateVariableIndices, parameters );

    residualMock R2( &hydra, 60, 1, stateVariableIndices, parameters );

    floatVector answer1( 10, 0 );
    answer1[ 0 ] = unknownVector[ 2 * configuration_unknown_count +  0 + 0 ];
    answer1[ 1 ] = unknownVector[ 2 * configuration_unknown_count +  0 + 1 ];
    answer1[ 2 ] = unknownVector[ 2 * configuration_unknown_count +  0 + 2 ];
    answer1[ 3 ] = unknownVector[ 2 * configuration_unknown_count +  0 + 3 ];
    answer1[ 4 ] = unknownVector[ 2 * configuration_unknown_count +  0 + 4 ];
    answer1[ 5 ] = unknownVector[ 2 * configuration_unknown_count + 10 + 0 ];
    answer1[ 6 ] = unknownVector[ 2 * configuration_unknown_count + 10 + 1 ];
    answer1[ 7 ] = unknownVector[ 2 * configuration_unknown_count + 10 + 2 ];
    answer1[ 8 ] = unknownVector[ 2 * configuration_unknown_count + 10 + 3 ];
    answer1[ 9 ] = unknownVector[ 2 * configuration_unknown_count + 10 + 4 ];

    BOOST_TEST( answer1 == *R.getConstraints( ) , CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setConstraintJacobians, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5 };

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
                                   0.06, 0.07, 0.08, 0.09, 0.10,
                                   0.11, 0.12, 0.13, 0.14, 0.15 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 15;

    std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };

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

    class residualMock : public tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::residual{

        public:

            using tardigradeHydra::micromorphicDruckerPragerPlasticityOptimization::residual::residual;

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

            std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14 };

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

                plasticity = residualMock( this, 60, 1, stateVariableIndices, plasticParameters );

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

    residualMock R(  &hydra, 60, 1, stateVariableIndices, parameters );

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

        residualMock Rp(  &hydrap, 60, 1, stateVariableIndices, parameters );

        residualMock Rm(  &hydram, 60, 1, stateVariableIndices, parameters );

        floatVector vp = *Rp.getConstraints( );

        floatVector vm = *Rm.getConstraints( );

        for ( unsigned int j = 0; j < vp.size( ); j++ ){

            jacobian[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( jacobian ), *R.getConstraintJacobians( ), 1e-5, 1e-6 ) );

}
