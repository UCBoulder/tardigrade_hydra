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

BOOST_AUTO_TEST_CASE( test_setActiveConstraints, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the active constraints
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

            floatType _macYield = 1.3;
            floatType _micYield = 2.4;
            floatVector _micGradYield = { 1.34, 3.4, 0.2 };

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::residual;

            const std::vector< bool > * publicGetActiveConstraints( ){

                return get_activeConstraints( );

            }

        protected:

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setMacroYield;

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setMicroYield;

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setMicroGradientYield;

            virtual void setMacroYield( ) override{

                auto yield = get_setDataStorage_macroYield( );

                *yield.value = _macYield;

            }

            virtual void setMicroYield( ) override{

                auto yield = get_setDataStorage_microYield( );

                *yield.value = _micYield;

            }

            virtual void setMicroGradientYield( ) override{

                auto yield = get_setDataStorage_microGradientYield( );

                *yield.value = _micGradYield;

            }

    };

    std::vector< bool > answer = { true, true, true, true, true };

    tardigradeHydra::hydraBaseMicromorphic hydra(
        time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
        microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
        { }, { },
        previousStateVariables, parameters,
        numConfigurations, numNonLinearSolveStateVariables,
        dimension, configuration_unknown_count,
        tolr, tola, maxIterations, maxLSIterations, lsAlpha
    );

    residualMock residual1( &hydra, 55, 1, stateVariableIndices, parameters, 0.27 );

    BOOST_TEST( answer == *residual1.publicGetActiveConstraints( ), CHECK_PER_ELEMENT );

    residualMock residual2( &hydra, 55, 1, stateVariableIndices, parameters, 0.27 );
    
    residual2._macYield = -3;
    residual2._micYield = -2;
    residual2._micGradYield = { -3, -4, -5 };

    answer = { false, false, false, false, false };

    BOOST_TEST( answer == *residual2.publicGetActiveConstraints( ), CHECK_PER_ELEMENT );

    residualMock residual3( &hydra, 55, 1, stateVariableIndices, parameters, 0.27 );
    
    residual3._macYield = -3;
    residual3._micYield =  2;
    residual3._micGradYield = { -3, 4, -5 };

    answer = { false, true, false, true, false };

    BOOST_TEST( answer == *residual3.publicGetActiveConstraints( ), CHECK_PER_ELEMENT );

    residualMock residual4( &hydra, 55, 1, stateVariableIndices, parameters, 0.27 );
    
    residual4._macYield =  3;
    residual4._micYield =-2;
    residual4._micGradYield = { 3, -4, 5 };

    answer = { true, false, true, false, true };

    BOOST_TEST( answer == *residual4.publicGetActiveConstraints( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_updateActiveConstraints, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test updating the active constraints
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

            floatType _macYield = 1.3;
            floatType _micYield = 2.4;
            floatVector _micGradYield = { 1.34, 3.4, 0.2 };

            floatVector _pMult = { 0.1, 0.2, 0.3, 0.4, 0.5 };

            unsigned int _setActiveConstraints_calls = 0;

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::residual;

            const std::vector< bool > * publicGetActiveConstraints( ){

                return get_activeConstraints( );

            }

            const void publicUpdateActiveConstraints( ){

                updateActiveConstraints( );

            }

        protected:

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setMacroYield;

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setMicroYield;

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setMicroGradientYield;

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setActiveConstraints;

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setPlasticMultipliers;

            virtual void setActiveConstraints( ) override{

                tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setActiveConstraints( );

		_setActiveConstraints_calls++;

	    }

            virtual void setMacroYield( ) override{

                auto yield = get_setDataStorage_macroYield( );

                *yield.value = _macYield;

            }

            virtual void setMicroYield( ) override{

                auto yield = get_setDataStorage_microYield( );

                *yield.value = _micYield;

            }

            virtual void setMicroGradientYield( ) override{

                auto yield = get_setDataStorage_microGradientYield( );

                *yield.value = _micGradYield;

            }

            virtual void setPlasticMultipliers( ) override{

               auto plasticMultipliers = get_setDataStorage_plasticMultipliers( );

	       *plasticMultipliers.value = _pMult;

	    }

    };

    class hydraBaseMicromorphicMock : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

    };

    std::vector< bool > answer = { true, true, true, true, true };

    hydraBaseMicromorphicMock hydra(
        time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
        microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
        { }, { },
        previousStateVariables, parameters,
        numConfigurations, numNonLinearSolveStateVariables,
        dimension, configuration_unknown_count,
        tolr, tola, maxIterations, maxLSIterations, lsAlpha
    );

    residualMock residual( &hydra, 55, 1, stateVariableIndices, parameters, 0.27 );

    BOOST_TEST( answer == *residual.publicGetActiveConstraints( ), CHECK_PER_ELEMENT );
    BOOST_TEST( residual._setActiveConstraints_calls == 1 );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );
    
    for ( auto v = std::begin( residual._pMult ); v != std::end( residual._pMult ); ++v ){    

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );
        *v *= -1;

        residual.publicUpdateActiveConstraints( );

        answer[ v - std::begin( residual._pMult ) ] = false;

        BOOST_TEST( answer == *residual.publicGetActiveConstraints( ), CHECK_PER_ELEMENT );
        BOOST_TEST( residual._setActiveConstraints_calls == 1 );

    }

}

BOOST_AUTO_TEST_CASE( test_setStateVariableResiduals, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the state variable residuals
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

            floatVector _updatedZ = { 1, 2, 3, 4, 5 };

            floatVector _Z = { 6, 7, 8, 9, 10 };

            floatVector _gammaDot = { 0.1, -0.2, 0.3, -0.4, 0.5 };

            floatVector _previousGammaDot = { -1, -2, -3, -4, -5 };

            floatType _macYield = 1.3;

            floatType _micYield = 2.4;

            floatVector _micGradYield = { 1.34, 3.4, 0.2 };

            std::vector< bool > _ac = { true, true, false, false, true };

            floatVector _pISVS = { .01, .02, .03, .04, .05, .06, .07, .08, .09, .10 };

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::residual;

            const floatVector * publicGetStateVariableResiduals( ){

                return get_stateVariableResiduals( );

            }

        protected:

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setPlasticMultipliers;

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

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setUpdatedPlasticStrainLikeISVs;

            virtual void setUpdatedPlasticStrainLikeISVs( ) override{

                auto z = get_setDataStorage_updatedPlasticStrainLikeISVs( );

                *z.value = _updatedZ;

            }

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setPlasticStrainLikeISVs;

            virtual void setPlasticStrainLikeISVs( ) override{

                auto z = get_setDataStorage_plasticStrainLikeISVs( );

                *z.value = _Z;

            }

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setMacroYield;

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setMicroYield;

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setMicroGradientYield;

            virtual void setMacroYield( ) override{

                auto yield = get_setDataStorage_macroYield( );

                *yield.value = _macYield;

            }

            virtual void setMicroYield( ) override{

                auto yield = get_setDataStorage_microYield( );

                *yield.value = _micYield;

            }

            virtual void setMicroGradientYield( ) override{

                auto yield = get_setDataStorage_microGradientYield( );

                *yield.value = _micGradYield;

            }

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setActiveConstraints;

            virtual void setActiveConstraints( ) override{

                auto activeConstraints = get_setDataStorage_activeConstraints( );

                *activeConstraints.value = _ac;

            }

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setPlasticStateVariables;

            virtual void setPlasticStateVariables( ) override{

                auto plasticStateVariables = get_setDataStorage_plasticStateVariables( );

                *plasticStateVariables.value = _pISVS;

            }

    };

    tardigradeHydra::hydraBaseMicromorphic hydra(
        time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
        microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
        { }, { },
        previousStateVariables, parameters,
        numConfigurations, numNonLinearSolveStateVariables,
        dimension, configuration_unknown_count,
        tolr, tola, maxIterations, maxLSIterations, lsAlpha
    );

    residualMock residual( &hydra, 55, 1, stateVariableIndices, parameters, 0.27 );

    std::vector< double > answer( 10, 0 );
    for ( unsigned int i = 0; i < 5; ++i ){
        answer[ i ] = residual._updatedZ[ i ] - residual._Z[ i ];
    }
    answer[ 5 ] = residual._macYield;
    answer[ 6 ] = residual._micYield;
    answer[ 7 ] = residual._gammaDot[ 2 ];
    answer[ 8 ] = residual._gammaDot[ 3 ];
    answer[ 9 ] = residual._micGradYield[ 2 ];

    BOOST_TEST( answer == *residual.publicGetStateVariableResiduals( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setStateVariableJacobians, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the state variable jacobians
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

    floatVector parameters = { 2, 0.53895133, 3.7172145,
                               2, 0.37773052, 9.2739145,
                               2, 0.53186824, 7.5454313,
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

    class residualMock : public tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual {

        public:

            std::vector< double > _gammaDot         = {  1,  2,  3,  4 };
            std::vector< double > _previousGammaDot = { .1, .2, .3, .4 };

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::residual;

            const floatVector * publicGetUpdatedPlasticStrainLikeISVs( ){

                return get_updatedPlasticStrainLikeISVs( );

            }

            const floatVector * publicGetPlasticStrainLikeISVs( ){

                return get_plasticStrainLikeISVs( );

            }

            const std::vector< bool > * publicGetActiveConstraints( ){

                return get_activeConstraints( );

            }

            const floatVector * publicGetPlasticMultipliers( ){

                return get_plasticMultipliers( );

            }

            const floatVector * publicGetStateVariableResiduals( ){

                return get_stateVariableResiduals( );

            }

            const floatVector * publicGetStateVariableJacobians( ){

                return get_stateVariableJacobians( );

            }

            const floatType * publicGetMacroYield( ){

                return get_macroYield( );

            }

            const floatType * publicGetMicroYield( ){

                return get_microYield( );

            }

            const floatVector * publicGetMicroGradientYield( ){

                return get_microGradientYield( );

            }

        protected:

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

    hydraBaseMicromorphicMock hydra(
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

    std::vector< double > answer( 10, 0 );

    for ( unsigned int i = 0; i < 5; ++i ){
       answer[ i ] = ( *residual.publicGetUpdatedPlasticStrainLikeISVs( ) )[ i ] - ( *residual.publicGetPlasticStrainLikeISVs( ) )[ i ];
    }
    answer[ 5 ] = *residual.publicGetMacroYield( );
    answer[ 6 ] = *residual.publicGetMicroYield( );
    answer[ 7 ] = ( *residual.publicGetMicroGradientYield( ) )[ 0 ];
    answer[ 8 ] = ( *residual.publicGetMicroGradientYield( ) )[ 1 ];
    answer[ 9 ] = ( *residual.publicGetMicroGradientYield( ) )[ 2 ];

    BOOST_TEST( answer == *residual.publicGetStateVariableResiduals( ), CHECK_PER_ELEMENT );

    floatType eps = 1e-5;

    {

        floatVector x = unknownVector;

        constexpr unsigned int VAR_SIZE = 100;

        constexpr unsigned int OUT_SIZE = 10;

        for ( unsigned int i = 0; i < VAR_SIZE; ++i ){

            floatType delta = eps * std::fabs( x[ i ] ) + eps;

            floatVector xp = x;
            floatVector xm = x;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            hydraBaseMicromorphicMock hydrap(
                time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );

            hydraBaseMicromorphicMock hydram(
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

            floatVector vp = *residualp.publicGetStateVariableResiduals( );

            floatVector vm = *residualm.publicGetStateVariableResiduals( );

            for ( unsigned int j = 0; j < OUT_SIZE; ++j ){

                BOOST_TEST( ( *residual.publicGetStateVariableJacobians( ) )[ VAR_SIZE * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_setStateVariableJacobians2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the state variable jacobians
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

    floatVector parameters = { 2, 200.53895133, 3.7172145,
                               2, 200.37773052, 9.2739145,
                               2, 200.53186824, 7.5454313,
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

    class residualMock : public tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual {

        public:

            std::vector< double > _gammaDot         = {  1,  2,  3,  4 };
            std::vector< double > _previousGammaDot = { .1, .2, .3, .4 };

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::residual;

            const floatVector * publicGetUpdatedPlasticStrainLikeISVs( ){

                return get_updatedPlasticStrainLikeISVs( );

            }

            const floatVector * publicGetPlasticStrainLikeISVs( ){

                return get_plasticStrainLikeISVs( );

            }

            const std::vector< bool > * publicGetActiveConstraints( ){

                return get_activeConstraints( );

            }

            const floatVector * publicGetPlasticMultipliers( ){

                return get_plasticMultipliers( );

            }

            const floatVector * publicGetStateVariableResiduals( ){

                return get_stateVariableResiduals( );

            }

            const floatVector * publicGetStateVariableJacobians( ){

                return get_stateVariableJacobians( );

            }

            const floatType * publicGetMacroYield( ){

                return get_macroYield( );

            }

            const floatType * publicGetMicroYield( ){

                return get_microYield( );

            }

            const floatVector * publicGetMicroGradientYield( ){

                return get_microGradientYield( );

            }

        protected:

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

    hydraBaseMicromorphicMock hydra(
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

    std::vector< double > answer( 10, 0 );

    for ( unsigned int i = 0; i < 5; ++i ){
       answer[ i ] = ( *residual.publicGetUpdatedPlasticStrainLikeISVs( ) )[ i ] - ( *residual.publicGetPlasticStrainLikeISVs( ) )[ i ];
    }
    answer[ 5 ] = unknownVector[ 90 ];
    answer[ 6 ] = unknownVector[ 91 ];
    answer[ 7 ] = unknownVector[ 92 ];
    answer[ 8 ] = unknownVector[ 93 ];
    answer[ 9 ] = unknownVector[ 94 ];

    BOOST_TEST( answer == *residual.publicGetStateVariableResiduals( ), CHECK_PER_ELEMENT );

    floatType eps = 1e-5;

    {

        floatVector x = unknownVector;

        constexpr unsigned int VAR_SIZE = 100;

        constexpr unsigned int OUT_SIZE = 10;

        for ( unsigned int i = 0; i < VAR_SIZE; ++i ){

            floatType delta = eps * std::fabs( x[ i ] ) + eps;

            floatVector xp = x;
            floatVector xm = x;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            hydraBaseMicromorphicMock hydrap(
                time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );

            hydraBaseMicromorphicMock hydram(
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

            floatVector vp = *residualp.publicGetStateVariableResiduals( );

            floatVector vm = *residualm.publicGetStateVariableResiduals( );

            for ( unsigned int j = 0; j < OUT_SIZE; ++j ){

                BOOST_TEST( ( *residual.publicGetStateVariableJacobians( ) )[ VAR_SIZE * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_setStateVariableJacobians3, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the state variable jacobians
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

    floatVector parameters = { 2, 10.53895133, 3.7172145,
                               2, 10.37773052, 9.2739145,
                               2, 100.53186824, 7.5454313,
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

    class residualMock : public tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual {

        public:

            std::vector< double > _gammaDot         = {  1,  2,  3,  4 };
            std::vector< double > _previousGammaDot = { .1, .2, .3, .4 };

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::residual;

            const floatVector * publicGetUpdatedPlasticStrainLikeISVs( ){

                return get_updatedPlasticStrainLikeISVs( );

            }

            const floatVector * publicGetPlasticStrainLikeISVs( ){

                return get_plasticStrainLikeISVs( );

            }

            const std::vector< bool > * publicGetActiveConstraints( ){

                return get_activeConstraints( );

            }

            const floatVector * publicGetPlasticMultipliers( ){

                return get_plasticMultipliers( );

            }

            const floatVector * publicGetStateVariableResiduals( ){

                return get_stateVariableResiduals( );

            }

            const floatVector * publicGetStateVariableJacobians( ){

                return get_stateVariableJacobians( );

            }

            const floatType * publicGetMacroYield( ){

                return get_macroYield( );

            }

            const floatType * publicGetMicroYield( ){

                return get_microYield( );

            }

            const floatVector * publicGetMicroGradientYield( ){

                return get_microGradientYield( );

            }

        protected:

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

    hydraBaseMicromorphicMock hydra(
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

    std::vector< double > answer( 10, 0 );

    for ( unsigned int i = 0; i < 5; ++i ){
       answer[ i ] = ( *residual.publicGetUpdatedPlasticStrainLikeISVs( ) )[ i ] - ( *residual.publicGetPlasticStrainLikeISVs( ) )[ i ];
    }
    answer[ 5 ] = unknownVector[ 90 ];
    answer[ 6 ] = *residual.publicGetMicroYield( );
    answer[ 7 ] = unknownVector[ 92 ];
    answer[ 8 ] = ( *residual.publicGetMicroGradientYield( ) )[ 1 ];
    answer[ 9 ] = unknownVector[ 94 ];

    BOOST_TEST( answer == *residual.publicGetStateVariableResiduals( ), CHECK_PER_ELEMENT );

    floatType eps = 1e-5;

    {

        floatVector x = unknownVector;

        constexpr unsigned int VAR_SIZE = 100;

        constexpr unsigned int OUT_SIZE = 10;

        for ( unsigned int i = 0; i < VAR_SIZE; ++i ){

            floatType delta = eps * std::fabs( x[ i ] ) + eps;

            floatVector xp = x;
            floatVector xm = x;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            hydraBaseMicromorphicMock hydrap(
                time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );

            hydraBaseMicromorphicMock hydram(
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

            floatVector vp = *residualp.publicGetStateVariableResiduals( );

            floatVector vm = *residualm.publicGetStateVariableResiduals( );

            for ( unsigned int j = 0; j < OUT_SIZE; ++j ){

                BOOST_TEST( ( *residual.publicGetStateVariableJacobians( ) )[ VAR_SIZE * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_setdStateVariableResidualsdD, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the state variable jacobians
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

    floatVector parameters = { 2, 0.53895133, 3.7172145,
                               2, 0.37773052, 9.2739145,
                               2, 0.53186824, 7.5454313,
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

    class residualMock : public tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual {

        public:

            std::vector< double > _gammaDot         = {  1,  2,  3,  4 };
            std::vector< double > _previousGammaDot = { .1, .2, .3, .4 };

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::residual;

            const floatVector * publicGetUpdatedPlasticStrainLikeISVs( ){

                return get_updatedPlasticStrainLikeISVs( );

            }

            const floatVector * publicGetPlasticStrainLikeISVs( ){

                return get_plasticStrainLikeISVs( );

            }

            const std::vector< bool > * publicGetActiveConstraints( ){

                return get_activeConstraints( );

            }

            const floatVector * publicGetPlasticMultipliers( ){

                return get_plasticMultipliers( );

            }

            const floatVector * publicGetStateVariableResiduals( ){

                return get_stateVariableResiduals( );

            }

            const floatVector * publicGetStateVariableJacobians( ){

                return get_stateVariableJacobians( );

            }

            const floatVector * publicGetdStateVariableResidualsdD( ){

                return get_dStateVariableResidualsdD( );

            }

            const floatType * publicGetMacroYield( ){

                return get_macroYield( );

            }

            const floatType * publicGetMicroYield( ){

                return get_microYield( );

            }

            const floatVector * publicGetMicroGradientYield( ){

                return get_microGradientYield( );

            }

        protected:

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

    hydraBaseMicromorphicMock hydra(
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

    std::vector< double > answer( 10, 0 );

    for ( unsigned int i = 0; i < 5; ++i ){
       answer[ i ] = ( *residual.publicGetUpdatedPlasticStrainLikeISVs( ) )[ i ] - ( *residual.publicGetPlasticStrainLikeISVs( ) )[ i ];
    }
    answer[ 5 ] = *residual.publicGetMacroYield( );
    answer[ 6 ] = *residual.publicGetMicroYield( );
    answer[ 7 ] = ( *residual.publicGetMicroGradientYield( ) )[ 0 ];
    answer[ 8 ] = ( *residual.publicGetMicroGradientYield( ) )[ 1 ];
    answer[ 9 ] = ( *residual.publicGetMicroGradientYield( ) )[ 2 ];

    BOOST_TEST( answer == *residual.publicGetStateVariableResiduals( ), CHECK_PER_ELEMENT );

    floatType eps = 1e-5;

    {

        floatVector x = tardigradeVectorTools::appendVectors( { deformationGradient, microDeformation, gradientMicroDeformation } );

        constexpr unsigned int VAR_SIZE = 45;

        constexpr unsigned int OUT_SIZE = 10;

        for ( unsigned int i = 0; i < VAR_SIZE; ++i ){

            floatType delta = eps * std::fabs( x[ i ] ) + eps;

            floatVector xp = x;
            floatVector xm = x;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            floatVector Fp( std::begin( xp ) + 0, std::begin( xp ) + 9 );
            floatVector Fm( std::begin( xm ) + 0, std::begin( xm ) + 9 );

            floatVector chip( std::begin( xp ) + 9, std::begin( xp ) + 18 );
            floatVector chim( std::begin( xm ) + 9, std::begin( xm ) + 18 );

            floatVector gradChip( std::begin( xp ) + 18, std::end( xp ) );
            floatVector gradChim( std::begin( xm ) + 18, std::end( xm ) );

            hydraBaseMicromorphicMock hydrap(
                time, deltaTime, temperature, previousTemperature, Fp, previousDeformationGradient,
                chip, previousMicroDeformation, gradChip, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );

            hydraBaseMicromorphicMock hydram(
                time, deltaTime, temperature, previousTemperature, Fm, previousDeformationGradient,
                chim, previousMicroDeformation, gradChim, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

            residualMock residualp( &hydrap, 55, 1, stateVariableIndices, parameters, 0.27 );

            residualMock residualm( &hydram, 55, 1, stateVariableIndices, parameters, 0.27 );

            floatVector vp = *residualp.publicGetStateVariableResiduals( );

            floatVector vm = *residualm.publicGetStateVariableResiduals( );

            for ( unsigned int j = 0; j < OUT_SIZE; ++j ){

                BOOST_TEST( ( *residual.publicGetdStateVariableResidualsdD( ) )[ VAR_SIZE * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_setdStateVariableResidualsdD2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the state variable jacobians
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

    floatVector parameters = { 2, 200.53895133, 3.7172145,
                               2, 200.37773052, 9.2739145,
                               2, 200.53186824, 7.5454313,
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

    class residualMock : public tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual {

        public:

            std::vector< double > _gammaDot         = {  1,  2,  3,  4 };
            std::vector< double > _previousGammaDot = { .1, .2, .3, .4 };

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::residual;

            const floatVector * publicGetUpdatedPlasticStrainLikeISVs( ){

                return get_updatedPlasticStrainLikeISVs( );

            }

            const floatVector * publicGetPlasticStrainLikeISVs( ){

                return get_plasticStrainLikeISVs( );

            }

            const std::vector< bool > * publicGetActiveConstraints( ){

                return get_activeConstraints( );

            }

            const floatVector * publicGetPlasticMultipliers( ){

                return get_plasticMultipliers( );

            }

            const floatVector * publicGetStateVariableResiduals( ){

                return get_stateVariableResiduals( );

            }

            const floatVector * publicGetStateVariableJacobians( ){

                return get_stateVariableJacobians( );

            }

            const floatVector * publicGetdStateVariableResidualsdD( ){

                return get_dStateVariableResidualsdD( );

            }

            const floatType * publicGetMacroYield( ){

                return get_macroYield( );

            }

            const floatType * publicGetMicroYield( ){

                return get_microYield( );

            }

            const floatVector * publicGetMicroGradientYield( ){

                return get_microGradientYield( );

            }

        protected:

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

    hydraBaseMicromorphicMock hydra(
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

    std::vector< double > answer( 10, 0 );

    for ( unsigned int i = 0; i < 5; ++i ){
       answer[ i ] = ( *residual.publicGetUpdatedPlasticStrainLikeISVs( ) )[ i ] - ( *residual.publicGetPlasticStrainLikeISVs( ) )[ i ];
    }
    answer[ 5 ] = unknownVector[ 90 ];
    answer[ 6 ] = unknownVector[ 91 ];
    answer[ 7 ] = unknownVector[ 92 ];
    answer[ 8 ] = unknownVector[ 93 ];
    answer[ 9 ] = unknownVector[ 94 ];

    BOOST_TEST( answer == *residual.publicGetStateVariableResiduals( ), CHECK_PER_ELEMENT );

    floatType eps = 1e-5;

    {

        floatVector x = tardigradeVectorTools::appendVectors( { deformationGradient, microDeformation, gradientMicroDeformation } );

        constexpr unsigned int VAR_SIZE = 45;

        constexpr unsigned int OUT_SIZE = 10;

        for ( unsigned int i = 0; i < VAR_SIZE; ++i ){

            floatType delta = eps * std::fabs( x[ i ] ) + eps;

            floatVector xp = x;
            floatVector xm = x;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            floatVector Fp( std::begin( xp ) + 0, std::begin( xp ) + 9 );
            floatVector Fm( std::begin( xm ) + 0, std::begin( xm ) + 9 );

            floatVector chip( std::begin( xp ) + 9, std::begin( xp ) + 18 );
            floatVector chim( std::begin( xm ) + 9, std::begin( xm ) + 18 );

            floatVector gradChip( std::begin( xp ) + 18, std::end( xp ) );
            floatVector gradChim( std::begin( xm ) + 18, std::end( xm ) );

            hydraBaseMicromorphicMock hydrap(
                time, deltaTime, temperature, previousTemperature, Fp, previousDeformationGradient,
                chip, previousMicroDeformation, gradChip, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );

            hydraBaseMicromorphicMock hydram(
                time, deltaTime, temperature, previousTemperature, Fm, previousDeformationGradient,
                chim, previousMicroDeformation, gradChim, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

            residualMock residualp( &hydrap, 55, 1, stateVariableIndices, parameters, 0.27 );

            residualMock residualm( &hydram, 55, 1, stateVariableIndices, parameters, 0.27 );

            floatVector vp = *residualp.publicGetStateVariableResiduals( );

            floatVector vm = *residualm.publicGetStateVariableResiduals( );

            for ( unsigned int j = 0; j < OUT_SIZE; ++j ){

                BOOST_TEST( ( *residual.publicGetdStateVariableResidualsdD( ) )[ VAR_SIZE * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_setdStateVariableResidualsdD3, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the state variable jacobians
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

    floatVector parameters = { 2, 10.53895133, 3.7172145,
                               2, 10.37773052, 9.2739145,
                               2, 100.53186824, 7.5454313,
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

    class residualMock : public tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual {

        public:

            std::vector< double > _gammaDot         = {  1,  2,  3,  4 };
            std::vector< double > _previousGammaDot = { .1, .2, .3, .4 };

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::residual;

            const floatVector * publicGetUpdatedPlasticStrainLikeISVs( ){

                return get_updatedPlasticStrainLikeISVs( );

            }

            const floatVector * publicGetPlasticStrainLikeISVs( ){

                return get_plasticStrainLikeISVs( );

            }

            const std::vector< bool > * publicGetActiveConstraints( ){

                return get_activeConstraints( );

            }

            const floatVector * publicGetPlasticMultipliers( ){

                return get_plasticMultipliers( );

            }

            const floatVector * publicGetStateVariableResiduals( ){

                return get_stateVariableResiduals( );

            }

            const floatVector * publicGetStateVariableJacobians( ){

                return get_stateVariableJacobians( );

            }

            const floatVector * publicGetdStateVariableResidualsdD( ){

                return get_dStateVariableResidualsdD( );

            }

            const floatType * publicGetMacroYield( ){

                return get_macroYield( );

            }

            const floatType * publicGetMicroYield( ){

                return get_microYield( );

            }

            const floatVector * publicGetMicroGradientYield( ){

                return get_microGradientYield( );

            }

        protected:

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

    hydraBaseMicromorphicMock hydra(
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

    std::vector< double > answer( 10, 0 );

    for ( unsigned int i = 0; i < 5; ++i ){
       answer[ i ] = ( *residual.publicGetUpdatedPlasticStrainLikeISVs( ) )[ i ] - ( *residual.publicGetPlasticStrainLikeISVs( ) )[ i ];
    }
    answer[ 5 ] = unknownVector[ 90 ];
    answer[ 6 ] = *residual.publicGetMicroYield( );
    answer[ 7 ] = unknownVector[ 92 ];
    answer[ 8 ] = ( *residual.publicGetMicroGradientYield( ) )[ 1 ];
    answer[ 9 ] = unknownVector[ 94 ];

    BOOST_TEST( answer == *residual.publicGetStateVariableResiduals( ), CHECK_PER_ELEMENT );

    floatType eps = 1e-5;

    {

        floatVector x = tardigradeVectorTools::appendVectors( { deformationGradient, microDeformation, gradientMicroDeformation } );

        constexpr unsigned int VAR_SIZE = 45;

        constexpr unsigned int OUT_SIZE = 10;

        for ( unsigned int i = 0; i < VAR_SIZE; ++i ){

            floatType delta = eps * std::fabs( x[ i ] ) + eps;

            floatVector xp = x;
            floatVector xm = x;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            floatVector Fp( std::begin( xp ) + 0, std::begin( xp ) + 9 );
            floatVector Fm( std::begin( xm ) + 0, std::begin( xm ) + 9 );

            floatVector chip( std::begin( xp ) + 9, std::begin( xp ) + 18 );
            floatVector chim( std::begin( xm ) + 9, std::begin( xm ) + 18 );

            floatVector gradChip( std::begin( xp ) + 18, std::end( xp ) );
            floatVector gradChim( std::begin( xm ) + 18, std::end( xm ) );

            hydraBaseMicromorphicMock hydrap(
                time, deltaTime, temperature, previousTemperature, Fp, previousDeformationGradient,
                chip, previousMicroDeformation, gradChip, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );

            hydraBaseMicromorphicMock hydram(
                time, deltaTime, temperature, previousTemperature, Fm, previousDeformationGradient,
                chim, previousMicroDeformation, gradChim, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

            residualMock residualp( &hydrap, 55, 1, stateVariableIndices, parameters, 0.27 );

            residualMock residualm( &hydram, 55, 1, stateVariableIndices, parameters, 0.27 );

            floatVector vp = *residualp.publicGetStateVariableResiduals( );

            floatVector vm = *residualm.publicGetStateVariableResiduals( );

            for ( unsigned int j = 0; j < OUT_SIZE; ++j ){

                BOOST_TEST( ( *residual.publicGetdStateVariableResidualsdD( ) )[ VAR_SIZE * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_correctResiduals, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test correcting the global residual after the active set update
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

    floatVector parameters = { 2, 200.53895133, 3.7172145,
                               2, 200.37773052, 9.2739145,
                               2, 200.53186824, 7.5454313,
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

            virtual void setResidual( ) override{

                auto residual = get_setDataStorage_residual( );
                residual.zero( *getNumEquations( ) );
                for ( auto v = residual.begin( ); v != residual.end( ); ++v ){
                    *v = -0.125 * ( v - residual.begin( ) );
                }

            }

    };

    class residualMock : public tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual {

        public:

            std::vector< bool >   _AS               = { true, true, true, true, true };
            std::vector< double > _Z                = { 0.11, 0.22, 0.33, 0.44, 0.55 };
            std::vector< double > _Zupdate          = { 0.22, 0.33, 0.44, 0.55, 0.66 };
            std::vector< double > _yields           = { .1, .2, .3, .4, .5 };
            std::vector< double > _gammaDot         = {  1,  2,  3,  4, 5 };

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::residual;
            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setPlasticMultipliers;
            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setPlasticStrainLikeISVs;
            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::setUpdatedPlasticStrainLikeISVs;

            virtual void setActiveConstraints( ) override{
                *get_setDataStorage_activeConstraints( ).value = _AS;
            }

            virtual void setMacroYield( ) override{
                *get_setDataStorage_macroYield( ).value = _yields[ 0 ];
            }

            virtual void setMicroYield( ) override{
                *get_setDataStorage_microYield( ).value = _yields[ 1 ];
            }

            virtual void setMicroGradientYield( ) override{
                *get_setDataStorage_microGradientYield( ).value = floatVector( _yields.begin( ) + 2, _yields.end( ) );
            }

            virtual void setPlasticStrainLikeISVs( ) override{
                *get_setDataStorage_plasticStrainLikeISVs( ).value = _Z;
            }

            virtual void setUpdatedPlasticStrainLikeISVs( ) override{
                *get_setDataStorage_updatedPlasticStrainLikeISVs( ).value = _Zupdate;
            }

            virtual void setPlasticMultipliers( ) override{
                *get_setDataStorage_plasticMultipliers( ).value = _gammaDot;
            }

            virtual void setResidual( ){
                auto value = get_setDataStorage_residual( );
                value.zero( 55 );
                for ( auto v = value.begin( ); v != value.begin( ) + 45; ++v ){
                    *v += ( v - value.begin( ) ) + 0;
                }

                std::copy(
                    std::begin( *get_stateVariableResiduals( ) ),
                    std::end( *get_stateVariableResiduals( ) ),
                    value.begin( ) + 45
                );

            }

            virtual const floatVector *publicGetStateVariableResiduals( ){
                return get_stateVariableResiduals( );
            }

            virtual void publicCorrectResiduals( ){
                correctResiduals( );
            }

        protected:

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

            virtual const std::vector< double > *publicGetResidual( ){
                return getResidual( );
            }

            void publicSetAllowModifyGlobalResidual( const bool value ){
                setAllowModifyGlobalResidual( value );
            }

            void publicSetPlasticResidualCurrent( ){
                setCurrentResidualIndexMeaningful( true );
                setCurrentResidualIndex( ( std::begin( *getResidualClasses( ) ) + 1 ) - std::begin( *getResidualClasses( ) ) );
            }

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

    hydraBaseMicromorphicMock hydra(
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

    floatVector answer( std::begin( *hydra.getResidual( ) ), std::end( *hydra.getResidual( ) ) );
    answer[ 96 ] = residual._gammaDot[ 1 ];
    answer[ 98 ] = residual._gammaDot[ 3 ];

    hydra.publicSetAllowModifyGlobalResidual( true );
    residual._AS = { true, false, true, false, true };
    hydra.publicSetPlasticResidualCurrent( );
    residual.publicCorrectResiduals( );
    hydra.publicSetAllowModifyGlobalResidual( false );

    BOOST_TEST( *hydra.getResidual( ) == answer, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_successfulNLStep, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test correcting the global residual after the active set update
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

    floatVector parameters = { 2, 200.53895133, 3.7172145,
                               2, 200.37773052, 9.2739145,
                               2, 200.53186824, 7.5454313,
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

            virtual void setResidual( ) override{

                auto residual = get_setDataStorage_residual( );
                residual.zero( *getNumEquations( ) );
                for ( auto v = residual.begin( ); v != residual.end( ); ++v ){
                    *v = -0.125 * ( v - residual.begin( ) );
                }

            }

    };

    class residualMock : public tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual {

        public:

            unsigned int num_updateActiveConstraints_calls = 0;
            unsigned int num_correctResiduals_calls = 0;

            using tardigradeHydra::micromorphicRadialReturnDruckerPragerPlasticity::residual::residual;

            virtual void updateActiveConstraints( ) override{
                num_updateActiveConstraints_calls++;
            }

            virtual void correctResiduals( ) override{
                num_correctResiduals_calls++;
            }


        protected:

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

    };

    hydraBaseMicromorphicMock hydra(
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

    residual.successfulNLStep( );

    BOOST_TEST( residual.num_updateActiveConstraints_calls == 1 );
    BOOST_TEST( residual.num_correctResiduals_calls == 1 );

}
