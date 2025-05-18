/**
  * \file test_tardigrade_hydraMicromorphicRadialReturnLinearElasticity.cpp
  *
  * Tests for tardigrade_hydraMicromorphicRadialReturnLinearElasticity
  */

#include<tardigrade_hydraMicromorphicRadialReturnLinearElasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_micromorphic_tools.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_tardigrade_hydraMicromorphicRadialReturnLinearElasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::micromorphicRadialReturnLinearElasticity::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::micromorphicRadialReturnLinearElasticity::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::micromorphicRadialReturnLinearElasticity::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

typedef tardigradeHydra::micromorphicRadialReturnLinearElasticity::variableType variableType; //!< Redefinition of the variable type
typedef tardigradeHydra::micromorphicRadialReturnLinearElasticity::variableVector variableVector; //!< Redefinition of the vector of variable types
typedef tardigradeHydra::micromorphicRadialReturnLinearElasticity::variableMatrix variableMatrix; //!< Redefinition of the matrix of variable types

typedef tardigradeHydra::micromorphicRadialReturnLinearElasticity::parameterType parameterType; //!< Redefinition of the parameter type
typedef tardigradeHydra::micromorphicRadialReturnLinearElasticity::parameterVector parameterVector; //!< Redefinition of the vector of parameters

typedef tardigradeHydra::micromorphicRadialReturnLinearElasticity::constantType constantType; //!< Redefinition of the constant type
typedef tardigradeHydra::micromorphicRadialReturnLinearElasticity::constantVector constantVector; //!< Redefinition of the vector of constants
typedef tardigradeHydra::micromorphicRadialReturnLinearElasticity::constantMatrix constantMatrix; //!< Redefinition of the matrix of constants

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

            if ( ( std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v1[ i ] ) > eps ) ||
                 ( std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v2[ i ] ) > eps ) ){

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

BOOST_AUTO_TEST_CASE( test_setTrialStress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector parameters = { 2, 0.1, 0.2, 5, 0.3, 0.4, 0.5, 0.6, 0.7, 11, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2, 1.9, 2.0 };

    floatVector unknownVector( 90, 0 );

    floatType val = -1.0;
    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){
        unknownVector[ i ] = val;
        val--;
    }

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

    floatType tolr = 1e-2;

    floatType tola = 1e-3;

    unsigned int maxIterations = 24;

    unsigned int maxLSIterations = 45;

    floatType lsAlpha = 2.3;

    class residualMock : public tardigradeHydra::micromorphicRadialReturnLinearElasticity::residual {

        public:

            std::vector< double > answer = { 1, 2, 3, 4 };

            using tardigradeHydra::micromorphicRadialReturnLinearElasticity::residual::residual;

            const floatVector * publicGetTrialStress( ){

                return get_trialStress( );

            }

        protected:

            virtual void setStress( ) override{

                auto stress = get_setDataStorage_stress( );

                stress.zero( answer.size( ) );

                *stress.value = answer;

            }

            virtual void setdTrialStressdD( ) override{

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

    residualMock residual( &hydra, 45, parameters );

    BOOST_TEST( residual.answer == *residual.publicGetTrialStress( ) );

}

BOOST_AUTO_TEST_CASE( test_setdTrialStressdD, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector parameters = { 2, 0.1, 0.2, 5, 0.3, 0.4, 0.5, 0.6, 0.7, 11, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2, 1.9, 2.0 };

    floatVector unknownVector( 90, 0 );

    floatType val = -1.0;
    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){
        unknownVector[ i ] = val;
        val--;
    }

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

    floatType tolr = 1e-2;

    floatType tola = 1e-3;

    unsigned int maxIterations = 24;

    unsigned int maxLSIterations = 45;

    floatType lsAlpha = 2.3;

    class residualMock : public tardigradeHydra::micromorphicRadialReturnLinearElasticity::residual {

        public:

            using tardigradeHydra::micromorphicRadialReturnLinearElasticity::residual::residual;

            const floatVector * publicGetTrialStress( ){

                return get_trialStress( );

            }

            const floatVector * publicGetdTrialStressdD( ){

                return get_dTrialStressdD( );

            }

        protected:

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

    residualMock residual( &hydra, 45, parameters );

    BOOST_TEST( *residual.getStress( ) == *residual.publicGetTrialStress( ), CHECK_PER_ELEMENT );

    floatType eps = 1e-5;

    {

        floatVector x = tardigradeVectorTools::appendVectors( { deformationGradient, microDeformation, gradientMicroDeformation } );

        constexpr unsigned int VAR_SIZE = 45;

        constexpr unsigned int OUT_SIZE = 45;

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

            tardigradeHydra::hydraBaseMicromorphic hydrap(
                time, deltaTime, temperature, previousTemperature, Fp, previousDeformationGradient,
                chip, previousMicroDeformation, gradChip, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );
            
            tardigradeHydra::hydraBaseMicromorphic hydram(
                time, deltaTime, temperature, previousTemperature, Fm, previousDeformationGradient,
                chim, previousMicroDeformation, gradChim, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );
            
            residualMock Rp( &hydrap, 45, parameters );
            residualMock Rm( &hydram, 45, parameters );

            floatVector vp = *Rp.publicGetTrialStress( );
            floatVector vm = *Rm.publicGetTrialStress( );

            for ( unsigned int j = 0; j < OUT_SIZE; ++j ){

                BOOST_TEST( ( *residual.publicGetdTrialStressdD( ) )[ VAR_SIZE * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_setResidaul, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector parameters = { 2, 0.1, 0.2, 5, 0.3, 0.4, 0.5, 0.6, 0.7, 11, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2, 1.9, 2.0 };

    floatVector unknownVector( 90, 0 );

    floatType val = -1.0;
    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){
        unknownVector[ i ] = val;
        val--;
    }

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

    floatType tolr = 1e-2;

    floatType tola = 1e-3;

    unsigned int maxIterations = 24;

    unsigned int maxLSIterations = 45;

    floatType lsAlpha = 2.3;

    class residualMock : public tardigradeHydra::micromorphicRadialReturnLinearElasticity::residual {

        public:

            std::vector< double > trialStress = buildTrialStress( );

            std::vector< double > answer = buildAnswer( );

            using tardigradeHydra::micromorphicRadialReturnLinearElasticity::residual::residual;

            const floatVector * publicGetTrialStress( ){

                return get_trialStress( );

            }

        protected:

            virtual std::vector< double > buildTrialStress( ){
                std::vector< double > value( 45, 0 );
                for ( unsigned int i = 0; i < 45; ++i ){
                    value[ i ] = 0.25 * ( i + 1 );
                }
                return value;
            }

            virtual std::vector< double > buildAnswer( ){
                std::vector< double > value( 45, 0 );
                for ( unsigned int i = 0; i < 45; ++i ){
                    value[i] = ( *hydra->getStress( ) )[ i ] - trialStress[ i ];
                }
                return value;
            }

            virtual void setStress( ) override{

                auto stress = get_setDataStorage_stress( );

                stress.zero( trialStress.size( ) );

                *stress.value = trialStress;

            }

            virtual void setdTrialStressdD( ) override{

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

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock residual( &hydra, 45, parameters );

    BOOST_TEST( residual.answer == *residual.getResidual( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_jacobians, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector parameters = { 2, 0.1, 0.2, 5, 0.3, 0.4, 0.5, 0.6, 0.7, 11, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2, 1.9, 2.0 };

    floatVector unknownVector( 90, 0 );

    floatType val = -1.0;
    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){
        unknownVector[ i ] = val;
        val--;
    }

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

    floatType tolr = 1e-2;

    floatType tola = 1e-3;

    unsigned int maxIterations = 24;

    unsigned int maxLSIterations = 45;

    floatType lsAlpha = 2.3;

    class residualMock : public tardigradeHydra::micromorphicRadialReturnLinearElasticity::residual {

        public:

            using tardigradeHydra::micromorphicRadialReturnLinearElasticity::residual::residual;

            const floatVector * publicGetTrialStress( ){

                return get_trialStress( );

            }

        protected:

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

    residualMock residual( &hydra, 45, parameters );

    residual.publicGetTrialStress( );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    std::vector< double > answer( 45, 0 );
    std::transform(
        std::begin( *hydra.getStress( ) ),
        std::end(   *hydra.getStress( ) ),
        std::begin( *residual.publicGetTrialStress( ) ),
        std::begin( answer ),
        std::minus<>( )
    );

    BOOST_TEST( answer == *residual.getResidual( ), CHECK_PER_ELEMENT );

    floatType eps = 1e-5;

    {

        floatVector x = unknownVector;

        constexpr unsigned int VAR_SIZE = 90;
        constexpr unsigned int OUT_SIZE = 45;

        for ( unsigned int i = 0; i < VAR_SIZE; ++i ){

            floatType delta = eps * std::fabs( unknownVector[ i ] ) + eps;

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

            residualMock residualp( &hydrap, 45, parameters );

            residualMock residualm( &hydram, 45, parameters );

            residualp.publicGetTrialStress( );

            residualm.publicGetTrialStress( );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, xp );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, xm );

            floatVector vp = *residualp.getResidual( );

            floatVector vm = *residualm.getResidual( );

            for ( unsigned int j = 0; j < OUT_SIZE; ++j ){

                BOOST_TEST( ( *residual.getJacobian( ) )[ VAR_SIZE * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    {

        floatVector x = tardigradeVectorTools::appendVectors( { deformationGradient, microDeformation, gradientMicroDeformation } );

        constexpr unsigned int VAR_SIZE = 45;

        constexpr unsigned int OUT_SIZE = 45;

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

            tardigradeHydra::hydraBaseMicromorphic hydrap(
                time, deltaTime, temperature, previousTemperature, Fp, previousDeformationGradient,
                chip, previousMicroDeformation, gradChip, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );

            tardigradeHydra::hydraBaseMicromorphic hydram(
                time, deltaTime, temperature, previousTemperature, Fm, previousDeformationGradient,
                chim, previousMicroDeformation, gradChim, previousGradientMicroDeformation,
                { }, { },
                previousStateVariables, parameters,
                numConfigurations, numNonLinearSolveStateVariables,
                dimension, configuration_unknown_count,
                tolr, tola, maxIterations, maxLSIterations, lsAlpha
            );

            residualMock residualp( &hydrap, 45, parameters );

            residualMock residualm( &hydram, 45, parameters );

            residualp.publicGetTrialStress( );

            residualm.publicGetTrialStress( );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

            floatVector vp = *residualp.getResidual( );

            floatVector vm = *residualm.getResidual( );

            for ( unsigned int j = 0; j < OUT_SIZE; ++j ){

                BOOST_TEST( ( *residual.getdRdD( ) )[ VAR_SIZE * j + i ] == ( vp[ j ] - vm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}

//BOOST_AUTO_TEST_CASE( test_setResidual, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
//
//    class residualMock : public tardigradeHydra::micromorphicRadialReturnLinearElasticity::residual {
//
//        public:
//
//            std::vector< double > answer = { 1, 2, 3, 4 };
//
//            using tardigradeHydra::micromorphicRadialReturnLinearElasticity::residual::residual;
//
//            residualMock( ){ }
//
//            const floatVector * publicGetTrialStress( ){
//
//                return get_trialStress( );
//
//            }
//
//        protected:
//
//            virtual void setStress( ) override{
//
//                auto stress = get_setDataStorage_stress( );
//
//                stress.zero( answer.size( ) );
//
//                *stress.value = answer;
//
//            }
//
//            virtual void setdTrialStressdD( ) override{
//
//            }
//
//    };
//
//    class hydraMock : public tardigradeHydra::hydraBaseMicromorphic{
//
//        public:
//
//            hydraMock( ){ }
//
//    };
//
//    hydraMock h;
//
//    residualMock residual;
//
//}
