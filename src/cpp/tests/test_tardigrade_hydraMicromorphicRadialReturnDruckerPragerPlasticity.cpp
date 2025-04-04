/**
  * \file test_tardigrade_hydraMicromorphicRadialReturnDruckerPragerPlasticity.cpp
  *
  * Tests for tardigrade_hydraMicromorphicRadialReturnDruckerPragerPlasticity
  */

#include<tardigrade_hydraMicromorphicRadialReturnDruckerPragerPlasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_micromorphic_tools.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_tardigrade_hydraMicromorphicRadialReturnDruckerPragerPlasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

typedef tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::variableType variableType; //!< Redefinition of the variable type
typedef tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::variableVector variableVector; //!< Redefinition of the vector of variable types
typedef tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::variableMatrix variableMatrix; //!< Redefinition of the matrix of variable types

typedef tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::parameterType parameterType; //!< Redefinition of the parameter type
typedef tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::parameterVector parameterVector; //!< Redefinition of the vector of parameters

typedef tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::constantType constantType; //!< Redefinition of the constant type
typedef tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::constantVector constantVector; //!< Redefinition of the vector of constants
typedef tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::constantMatrix constantMatrix; //!< Redefinition of the matrix of constants

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

BOOST_AUTO_TEST_CASE( test_computeDeltaIntegratedPlasticMultipliers, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the change in the plastic multipliers
     */

    // Form a hydra class and residual class
    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    variableVector deformationGradient = { -0.50668478, -0.48303615, -1.43907185,
                                            0.24106815,  0.10656166,  1.06851718,
                                           -0.18353482, -1.11676646, -0.68534721 };

    variableVector microDeformation = { -0.36302711,  0.94624033,  0.07877948,
                                        -0.68872716,  0.10276167, -0.81357538,
                                         0.22307023, -0.23788599,  0.67759487 };

    variableVector gradientMicroDeformation = { 0.02197053,  0.0370446 , -0.02739394,  0.06279444, -0.00127596,
                                               -0.01796094,  0.02814145, -0.05906054,  0.02578498, -0.01458269,
                                                0.00048507, -0.03393819, -0.03257968, -0.00138203, -0.04167585,
                                               -0.03382795,  0.01151479, -0.03641219,  0.01271894,  0.04506872,
                                                0.03179861,  0.04033839,  0.0440033 , -0.05450839, -0.05968426,
                                               -0.02910144,  0.0279304 };

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
        0.034167  , -0.01426024, -0.04564085, -0.01952319, -0.01018143 };

    floatVector parameters = {
        2, 0.53895133, -3.7172145,
        2, 0.37773052, -9.2739145,
        2, 0.53186824, -7.5454313,
        2, 0.95338442, 0.74042148,
        2, 0.38093104, 0.49241325,
        2, 0.82121039, 0.90566759,
        2, 0.01166325, 0.05331896,
        2, 0.32982199, 0.60161431,
        2, 0.58881096, 0.11473813
    };

    floatVector unknownVector( 90, 0 );

    floatType val = -1.0;
    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){
        unknownVector[ i ] = val;
        val--;
    }

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    std::vector< unsigned int > stateVariableIndices = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

    floatType tolr = 1e-2;

    floatType tola = 1e-3;

    unsigned int maxIterations = 24;

    unsigned int maxLSIterations = 45;

    floatType lsAlpha = 2.3;

    class residualMock : public tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual {

        public:

            std::vector< double > _gammaDot         = {  1,  2,  3,  4 };
            std::vector< double > _previousGammaDot = { .1, .2, .3, .4 };

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::residual;

            const floatVector * publicGetDeltaIntegratedPlasticMultiplier( ){

                return get_deltaIntegratedPlasticMultipliers( );

            }

        protected:

            using tardigradeHydra::micromorphicDruckerPragerPlasticity::residual::setPlasticMultipliers;

            virtual void setPlasticMultipliers( const bool isPrevious ) override{

                if ( isPrevious ){

                    auto gammaDot = get_setDataStorage_previousPlasticMultipliers( );

                    *gammaDot.value = _previousGammaDot;

                }
                else{

                    auto gammaDot = get_setDataStorage_plasticMultipliers( );

                    *gammaDot.value = _gammaDot;

                }

            }

    };

    std::vector< double > answer = { 0.80262, 1.60524, 2.40786, 3.21048 };

    tardigradeHydra::hydraBaseMicromorphic hydra(
        time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
        microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
        { }, { },
        previousStateVariables, parameters,
        numConfigurations, numNonLinearSolveStateVariables,
        dimension, configuration_unknown_count,
        tolr, tola, maxIterations, maxLSIterations, lsAlpha
    );

    residualMock residual( &hydra, 45, 1, stateVariableIndices, parameters, 0.27 );

    BOOST_TEST( answer == *residual.publicGetDeltaIntegratedPlasticMultiplier( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_computedDeltaIntegratedPlasticMultipliersdPlasticMultipliers, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the change in the plastic multipliers
     */

    // Form a hydra class and residual class
    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    variableVector deformationGradient = { -0.50668478, -0.48303615, -1.43907185,
                                            0.24106815,  0.10656166,  1.06851718,
                                           -0.18353482, -1.11676646, -0.68534721 };

    variableVector microDeformation = { -0.36302711,  0.94624033,  0.07877948,
                                        -0.68872716,  0.10276167, -0.81357538,
                                         0.22307023, -0.23788599,  0.67759487 };

    variableVector gradientMicroDeformation = { 0.02197053,  0.0370446 , -0.02739394,  0.06279444, -0.00127596,
                                               -0.01796094,  0.02814145, -0.05906054,  0.02578498, -0.01458269,
                                                0.00048507, -0.03393819, -0.03257968, -0.00138203, -0.04167585,
                                               -0.03382795,  0.01151479, -0.03641219,  0.01271894,  0.04506872,
                                                0.03179861,  0.04033839,  0.0440033 , -0.05450839, -0.05968426,
                                               -0.02910144,  0.0279304 };

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

    floatVector parameters = { 2, 0.53895133, -3.7172145,
                               2, 0.37773052, -9.2739145,
                               2, 0.53186824, -7.5454313,
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

    class residualMock : public tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual {

        public:

            std::vector< double > _gammaDot         = {  1,  2,  3,  4 };
            std::vector< double > _previousGammaDot = { .1, .2, .3, .4 };

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::residual;

            const floatVector * publicGetDeltaIntegratedPlasticMultiplier( ){

                return get_deltaIntegratedPlasticMultipliers( );

            }

            const floatVector * publicGetdDeltaIntegratedPlasticMultiplierdPlasticMultipliers( ){

                return get_dDeltaIntegratedPlasticMultipliersdPlasticMultipliers( );

            }

        protected:

    };

    std::vector< double > gDot(           std::end( unknownVector ) - 10,          std::end( unknownVector ) - 5 );
    std::vector< double > pGDot( std::end( previousStateVariables ) - 10, std::end( previousStateVariables ) - 5 );

    std::vector< double > answer = deltaTime * ( ( 1 - 0.27 ) * pGDot + 0.27 * gDot );

    tardigradeHydra::hydraBaseMicromorphic hydra(
        time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
        microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
        { }, { },
        previousStateVariables, parameters,
        numConfigurations, numNonLinearSolveStateVariables,
        dimension, configuration_unknown_count,
        tolr, tola, maxIterations, maxLSIterations, lsAlpha
    );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock residual( &hydra, 55, 1, stateVariableIndices, parameters, 0.27 );

    BOOST_TEST( answer == *residual.publicGetDeltaIntegratedPlasticMultiplier( ), CHECK_PER_ELEMENT );

    floatType eps = 1e-5;

    {

        floatVector x = unknownVector;

        constexpr unsigned int VAR_SIZE = 100;

        constexpr unsigned int OUT_SIZE = 5;

        for ( unsigned int i = 0; i < VAR_SIZE; ++i ){

            floatType delta = eps * std::fabs( x[ i ] ) + eps;

            floatVector xp = x;
            floatVector xm = x;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            tardigradeHydra::hydraBaseMicromorphic hydrap(
                time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );

            tardigradeHydra::hydraBaseMicromorphic hydram(
                time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, xp );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, xm );

            residualMock residualp( &hydrap, 55, 1, stateVariableIndices, parameters, 0.27 );

            residualMock residualm( &hydram, 55, 1, stateVariableIndices, parameters, 0.27 );

            floatVector vp = *residualp.publicGetDeltaIntegratedPlasticMultiplier( );

            floatVector vm = *residualm.publicGetDeltaIntegratedPlasticMultiplier( );

            for ( unsigned int j = 0; j < OUT_SIZE; ++j ){

                if ( ( ( i >= 90 ) && ( i < 95 ) ) && ( ( i - 90 ) == j ) ){

                    BOOST_TEST( ( *residual.publicGetdDeltaIntegratedPlasticMultiplierdPlasticMultipliers( ) )[ j ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

                }
                else{

                    BOOST_TEST( 0. == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

                }

            }

        }

    }

}
