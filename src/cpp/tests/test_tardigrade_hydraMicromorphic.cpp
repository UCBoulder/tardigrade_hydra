/**
  * \file test_tardigrade_hydraMicromorphic.cpp
  *
  * Tests for tardigrade_hydraMicromorphic
  */

#include<tardigrade_hydraMicromorphic.h>

#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_tardigrade_hydraMicromorphic
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

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

                std::cout << "v1: " << v1[ i ] << "\n";
                std::cout << "v2: " << v2[ i ] << "\n";
                std::cout << "r1: " << std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v1[ i ] ) << "\n";
                std::cout << "r2: " << std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v2[ i ] ) << "\n";

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

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_constructor, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    BOOST_TEST( time == *hydra.getTime( ) );

    BOOST_TEST( deltaTime == *hydra.getDeltaTime( ) );

    BOOST_TEST( temperature == *hydra.getTemperature( ) );

    BOOST_TEST( previousTemperature == *hydra.getPreviousTemperature( ) );

    BOOST_TEST( deformationGradient == *hydra.getDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousDeformationGradient == *hydra.getPreviousDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( microDeformation == *hydra.getMicroDeformation( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousMicroDeformation == *hydra.getPreviousMicroDeformation( ), CHECK_PER_ELEMENT );

    BOOST_TEST( gradientMicroDeformation == *hydra.getGradientMicroDeformation( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousGradientMicroDeformation == *hydra.getPreviousGradientMicroDeformation( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( configurationsAnswer ) == *hydra.get_configurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousConfigurationsAnswer ) == *hydra.get_previousConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( inverseConfigurationsAnswer ) == *hydra.get_inverseConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousInverseConfigurationsAnswer ) == *hydra.get_previousInverseConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( microConfigurationsAnswer ) == *hydra.get_microConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousMicroConfigurationsAnswer ) == *hydra.get_previousMicroConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( inverseMicroConfigurationsAnswer ) == *hydra.get_inverseMicroConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousInverseMicroConfigurationsAnswer ) == *hydra.get_previousInverseMicroConfigurations( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_constructor2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    class residualBaseMock : public tardigradeHydra::residualBase{

        using tardigradeHydra::residualBase::residualBase;

    };

    class residualBaseMockStress : public tardigradeHydra::residualBase{

        using tardigradeHydra::residualBase::residualBase;

        using tardigradeHydra::residualBase::setStress;

        floatVector stress = { 1,  2,  3,  4,  5,  6,  7,  8,  9,
                              10, 11, 12, 13, 14, 15, 16, 17, 18,
                              19, 20, 21, 22, 23, 24, 25, 26, 27,
                              28, 29, 30, 31, 32, 33, 34, 35, 36,
                              37, 38, 39, 40, 41, 42, 43, 44, 45 };

        virtual void setStress( ) override{

            setStress( stress );

        }

    };

    class hydraBaseMicromorphicMock : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            unsigned int ndecomp = 0;
    
            unsigned int nsrc = 0;

            residualBaseMockStress r1;
        
            residualBaseMock r2;

            unsigned int s1 = 45;

            unsigned int s2 = 51;

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                r1 = residualBaseMockStress( this, s1 );

                r2 = residualBaseMock( this, s2 );

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                residuals[ 0 ] = &r1;

                residuals[ 1 ] = &r2;

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

    BOOST_TEST( time == *hydra.getTime( ) );

    BOOST_TEST( deltaTime == *hydra.getDeltaTime( ) );

    BOOST_TEST( temperature == *hydra.getTemperature( ) );

    BOOST_TEST( previousTemperature == *hydra.getPreviousTemperature( ) );

    BOOST_TEST( deformationGradient == *hydra.getDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousDeformationGradient == *hydra.getPreviousDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( microDeformation == *hydra.getMicroDeformation( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousMicroDeformation == *hydra.getPreviousMicroDeformation( ), CHECK_PER_ELEMENT );

    BOOST_TEST( gradientMicroDeformation == *hydra.getGradientMicroDeformation( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousGradientMicroDeformation == *hydra.getPreviousGradientMicroDeformation( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( configurationsAnswer ) == *hydra.get_configurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousConfigurationsAnswer ) == *hydra.get_previousConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( inverseConfigurationsAnswer ) == *hydra.get_inverseConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousInverseConfigurationsAnswer ) == *hydra.get_previousInverseConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( microConfigurationsAnswer ) == *hydra.get_microConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousMicroConfigurationsAnswer ) == *hydra.get_previousMicroConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( inverseMicroConfigurationsAnswer ) == *hydra.get_inverseMicroConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousInverseMicroConfigurationsAnswer ) == *hydra.get_previousInverseMicroConfigurations( ), CHECK_PER_ELEMENT );

    floatType scale_factor = 0.68;

    hydra.setScaleFactor( scale_factor );

    floatType scaled_time = ( scale_factor - 1 ) * deltaTime + time;

    floatType scaled_deltaTime = scale_factor * deltaTime;

    floatType scaled_temperature = scale_factor * ( temperature - previousTemperature ) + previousTemperature;

    floatVector scaled_deformationGradient = scale_factor * ( deformationGradient - previousDeformationGradient ) + previousDeformationGradient;

    floatVector scaled_microDeformation = scale_factor * ( microDeformation - previousMicroDeformation ) + previousMicroDeformation;

    floatVector scaled_gradientMicroDeformation = scale_factor * ( gradientMicroDeformation - previousGradientMicroDeformation ) + previousGradientMicroDeformation;

    configurationsAnswer[ 0 ] = { 3.36544504,  4.2448981 ,  4.44435116,  6.01497997,  6.81066016, 7.60634036,  8.71712942,  9.11721303, 10.15729664 };

    inverseConfigurationsAnswer[ 0 ] = { -0.12122531, -1.8436894 ,  1.4336981 ,  3.69897495, -3.23648435,  0.80516088, -3.21617129,  4.48735505, -1.854682 };

    microConfigurationsAnswer[ 0 ] = { 24.54492337, 27.02899316, 28.83306296, 28.66493751, 30.58240981, 32.4998821 , 32.01684508, 34.02808772, 35.07933035 };

    inverseMicroConfigurationsAnswer[ 0 ] = { -1.2869465 ,  1.28214968, -0.13007957,  1.36076467, -2.41554301, 1.11945836, -0.14539197,  1.1729392 , -0.93868068 };

    BOOST_TEST( scaled_time == *hydra.getTime( ) );

    BOOST_TEST( scaled_deltaTime == *hydra.getDeltaTime( ) );

    BOOST_TEST( scaled_temperature == *hydra.getTemperature( ) );

    BOOST_TEST( previousTemperature == *hydra.getPreviousTemperature( ) );

    BOOST_TEST( scaled_deformationGradient == *hydra.getDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousDeformationGradient == *hydra.getPreviousDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( scaled_microDeformation == *hydra.getMicroDeformation( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousMicroDeformation == *hydra.getPreviousMicroDeformation( ), CHECK_PER_ELEMENT );

    BOOST_TEST( scaled_gradientMicroDeformation == *hydra.getGradientMicroDeformation( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousGradientMicroDeformation == *hydra.getPreviousGradientMicroDeformation( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( configurationsAnswer ) == *hydra.get_configurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousConfigurationsAnswer ) == *hydra.get_previousConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( inverseConfigurationsAnswer ) == *hydra.get_inverseConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousInverseConfigurationsAnswer ) == *hydra.get_previousInverseConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( microConfigurationsAnswer ) == *hydra.get_microConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousMicroConfigurationsAnswer ) == *hydra.get_previousMicroConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( inverseMicroConfigurationsAnswer ) == *hydra.get_inverseMicroConfigurations( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousInverseMicroConfigurationsAnswer ) == *hydra.get_previousInverseMicroConfigurations( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getSubMicroConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector answer = { 2.63306967, -0.16982965, -0.15109712,
                           1.26964411,  0.59017992,  0.35717977,
                           0.72297358, -0.4912966 ,  0.28183864 };

    BOOST_TEST( answer == hydra.getSubMicroConfiguration( 0, 2 ), CHECK_PER_ELEMENT );

    floatVector answer2 = { 0.87556098, -0.36781419,  1.20660407,
                            0.21141974,  0.82131808, -0.42325899,
                            0.12210416,  1.0388779 ,  2.07092075 };

    BOOST_TEST( answer2 == hydra.getSubMicroConfiguration( 1, 3 ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getPreceedingMicroConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector answer = { 2.63306967, -0.16982965, -0.15109712,
                           1.26964411,  0.59017992,  0.35717977,
                           0.72297358, -0.4912966 ,  0.28183864 };

    BOOST_TEST( answer == hydra.getPrecedingMicroConfiguration( 2 ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getUnknownVector, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    class stressMock : public tardigradeHydra::residualBaseMicromorphic{

        public:

            using tardigradeHydra::residualBaseMicromorphic::residualBaseMicromorphic;

            floatVector PK2   = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector SIGMA = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

            floatVector M     = { 19, 20, 21, 22, 23, 24, 25, 26, 27,
                                  28, 29, 30, 31, 32, 33, 34, 35, 36,
                                  37, 38, 39, 40, 41, 42, 43, 44, 45 };

        protected:

            using tardigradeHydra::residualBaseMicromorphic::setStress;

            virtual void setStress( ) override{

                setStress( tardigradeVectorTools::appendVectors( { PK2, SIGMA, M } ) );

            }

    };

    class residualMock : public tardigradeHydra::residualBaseMicromorphic{

        public:

            using tardigradeHydra::residualBaseMicromorphic::residualBaseMicromorphic;

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

            stressMock elasticity;

            residualMock additionalResiduals;

            unsigned int ndecomp = 0;
    
            unsigned int nsrc = 0;

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

        private:

            using tardigradeHydra::hydraBaseMicromorphic::setResidualClasses;

            virtual void setResidualClasses( ) override{

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = stressMock( this, 45 );

                additionalResiduals = residualMock( this, 140 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &additionalResiduals;

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

    floatVector nonLinearPreviousStateVariables( previousStateVariables.begin( ),
                                                 previousStateVariables.begin( ) + ( numConfigurations - 1 ) * configuration_unknown_count + numNonLinearSolveStateVariables );

    floatVector answer = tardigradeVectorTools::appendVectors( { hydra.elasticity.PK2, hydra.elasticity.SIGMA, hydra.elasticity.M,
                                                                 nonLinearPreviousStateVariables } );

    for ( unsigned int i = 1; i < numConfigurations; i++ ){

        for ( unsigned int j = 0; j < 3; j++ ){

            answer[ configuration_unknown_count + 9 * ( i - 1 ) + 3 * j + j ] += 1;

            answer[ configuration_unknown_count + ( numConfigurations - 1 ) * 9 + 9 * ( i - 1 ) + 3 * j + j ] += 1;

        }

    }

    BOOST_TEST( answer == *hydra.getUnknownVector( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getFollowingMicroConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector answer = { 0.37335959,  0.14440808,  0.07815635,
                          -0.05311829,  1.15252061, -0.83175389,
                          -0.26108048,  1.01045459,  1.40735072 };

    BOOST_TEST( answer == hydra.getFollowingMicroConfiguration( 1 ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getMicroConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector answer = { 2.00093808,  1.48052569, -0.5533865 ,
                           0.74087198,  1.25662963,  0.17264556,
                           0.75663016, -0.11948556, -0.07579564 };

    BOOST_TEST( answer == hydra.getMicroConfiguration( 0 ), CHECK_PER_ELEMENT );

    floatVector answer2 = { 0.75754206,  0.06435904,  0.30696868,
                           -0.10562995,  1.23107304, -0.33893099,
                            0.10069857,  0.36586446,  1.48352161 };

    BOOST_TEST( answer2 == hydra.getMicroConfiguration( 2 ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getPreviousSubMicroConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector answer = { 1.79993357, -0.26247057, -0.16201437,
                          -0.7635743 ,  0.6005345 ,  0.34819881,
                           0.16680516, -0.38411455,  0.34087847 };

    BOOST_TEST( answer == hydra.getPreviousSubMicroConfiguration( 0, 2 ), CHECK_PER_ELEMENT );

    floatVector answer2 = { 0.87556098, -0.36781419,  1.20660407,
                            0.21141974,  0.82131808, -0.42325899,
                            0.12210416,  1.0388779 ,  2.07092075 };

    BOOST_TEST( answer2 == hydra.getPreviousSubMicroConfiguration( 1, 3 ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getPreviousPreceedingMicroConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector answer = { 1.79993357, -0.26247057, -0.16201437,
                          -0.7635743 ,  0.6005345 ,  0.34819881,
                           0.16680516, -0.38411455,  0.34087847 };

    BOOST_TEST( answer == hydra.getPreviousPrecedingMicroConfiguration( 2 ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getPreviousFollowingMicroConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector answer = { 0.37335959,  0.14440808,  0.07815635,
                          -0.05311829,  1.15252061, -0.83175389,
                          -0.26108048,  1.01045459,  1.40735072 };

    BOOST_TEST( answer == hydra.getPreviousFollowingMicroConfiguration( 1 ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_get_previousMicroConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector answer = { 1.41998902,  0.88426858, -0.45291886,
                          -0.76867062,  0.02912989,  0.49178024,
                           0.30842234, -0.36954097,  0.0765542 };

    BOOST_TEST( answer == hydra.getPreviousMicroConfiguration( 0 ), CHECK_PER_ELEMENT );

    floatVector answer2 = { 0.75754206,  0.06435904,  0.30696868,
                           -0.10562995,  1.23107304, -0.33893099,
                            0.10069857,  0.36586446,  1.48352161 };

    BOOST_TEST( answer2 == hydra.getPreviousMicroConfiguration( 2 ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getSubMicroConfigurationJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector microConfigurations = *hydra.get_microConfigurations( );

    floatVector x = microConfigurations;

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 0;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension );

        delta[ i ] = std::fabs( eps * microConfigurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( microConfigurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( microConfigurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getSubMicroConfigurationJacobian( lower, upper ), CHECK_PER_ELEMENT );

    lower = 1;

    upper = 3;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension );

        delta[ i ] = std::fabs( eps * microConfigurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( microConfigurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( microConfigurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getSubMicroConfigurationJacobian( lower, upper ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getPrecedingMicroConfigurationJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector microConfigurations = *hydra.get_microConfigurations( );

    floatVector x = microConfigurations;

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 0;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension );

        delta[ i ] = std::fabs( eps * microConfigurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( microConfigurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( microConfigurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPrecedingMicroConfigurationJacobian( upper ), CHECK_PER_ELEMENT );

    lower = 0;

    upper = 3;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension );

        delta[ i ] = std::fabs( eps * microConfigurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( microConfigurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( microConfigurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPrecedingMicroConfigurationJacobian( upper ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getFollowingMicroConfigurationJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector microConfigurations = *hydra.get_microConfigurations( );

    floatVector x = microConfigurations;

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 1;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension );

        delta[ i ] = std::fabs( eps * microConfigurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( microConfigurations + delta, lower + 1, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( microConfigurations - delta, lower + 1, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getFollowingMicroConfigurationJacobian( lower ), CHECK_PER_ELEMENT );

    lower = 2;

    upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension );

        delta[ i ] = std::fabs( eps * microConfigurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( microConfigurations + delta, lower + 1, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( microConfigurations - delta, lower + 1, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getFollowingMicroConfigurationJacobian( lower ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getPreviousSubMicroConfigurationJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector microConfigurations = *hydra.get_previousMicroConfigurations( );

    floatVector x = microConfigurations;

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 0;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension );

        delta[ i ] = std::fabs( eps * microConfigurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( microConfigurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( microConfigurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPreviousSubMicroConfigurationJacobian( lower, upper ), CHECK_PER_ELEMENT );

    lower = 1;

    upper = 3;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension );

        delta[ i ] = std::fabs( eps * microConfigurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( microConfigurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( microConfigurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPreviousSubMicroConfigurationJacobian( lower, upper ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getPreviousPrecedingMicroConfigurationJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector microConfigurations = *hydra.get_previousMicroConfigurations( );

    floatVector x = microConfigurations;

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 0;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension );

        delta[ i ] = std::fabs( eps * microConfigurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( microConfigurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( microConfigurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPreviousPrecedingMicroConfigurationJacobian( upper ), CHECK_PER_ELEMENT );

    lower = 0;

    upper = 3;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension );

        delta[ i ] = std::fabs( eps * microConfigurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( microConfigurations + delta, lower, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( microConfigurations - delta, lower, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPreviousPrecedingMicroConfigurationJacobian( upper ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getPreviousFollowingMicroConfigurationJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector microConfigurations = *hydra.get_previousMicroConfigurations( );

    floatVector x = microConfigurations;

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( x.size( ), 0 ) );

    unsigned int lower = 1;

    unsigned int upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension );

        delta[ i ] = std::fabs( eps * microConfigurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( microConfigurations + delta, lower + 1, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( microConfigurations - delta, lower + 1, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPreviousFollowingMicroConfigurationJacobian( lower ), CHECK_PER_ELEMENT );

    lower = 2;

    upper = 4;

    for ( unsigned int i = 0; i < x.size( ); i++ ){

        floatVector delta( numConfigurations * dimension * dimension );

        delta[ i ] = std::fabs( eps * microConfigurations[ i ] ) + eps;

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydra.getSubConfiguration( microConfigurations + delta, lower + 1, upper ) );

        BOOST_CHECK_NO_THROW( Fscm = hydra.getSubConfiguration( microConfigurations - delta, lower + 1, upper ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == hydra.getPreviousFollowingMicroConfigurationJacobian( lower ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_get_dChi1dChi, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( dimension * dimension, 0 ) );

    for ( unsigned int i = 0; i < dimension * dimension; i++ ){

        floatVector delta( dimension * dimension, 0 );

        delta[ i ] = std::fabs( eps * microDeformation[ i ] ) + eps;

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

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydrap.getMicroConfiguration( 0 ) );

        BOOST_CHECK_NO_THROW( Fscm = hydram.getMicroConfiguration( 0 ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *hydra.get_dChi1dChi( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_getPreviousdChi1dChi, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( dimension * dimension, 0 ) );

    for ( unsigned int i = 0; i < dimension * dimension; i++ ){

        floatVector delta( dimension * dimension, 0 );

        delta[ i ] = std::fabs( eps * microDeformation[ i ] ) + eps;

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

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydrap.getPreviousMicroConfiguration( 0 ) );

        BOOST_CHECK_NO_THROW( Fscm = hydram.getPreviousMicroConfiguration( 0 ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *hydra.get_previousdChi1dChi( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_get_dChi1dChin, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( ( numConfigurations - 1 ) * dimension * dimension, 0 ) );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ) * dimension * dimension; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + ( numConfigurations - 1 ) * dimension * dimension ] = std::fabs( eps * previousStateVariables[ i + ( numConfigurations - 1 ) * dimension * dimension ] ) + eps;

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

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydrap.getMicroConfiguration( 0 ) );

        BOOST_CHECK_NO_THROW( Fscm = hydram.getMicroConfiguration( 0 ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i + ( numConfigurations - 1 ) * dimension * dimension ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *hydra.get_dChi1dChin( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_get_previousdChi1dChin, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatType eps = 1e-6;

    floatMatrix gradient( dimension * dimension, floatVector( ( numConfigurations - 1 ) * dimension * dimension, 0 ) );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ) * dimension * dimension; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + ( numConfigurations - 1 ) * dimension * dimension ] = std::fabs( eps * previousStateVariables[ i + ( numConfigurations - 1 ) * dimension * dimension ] ) + eps;

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

        floatVector Fscp;

        floatVector Fscm;

        BOOST_CHECK_NO_THROW( Fscp = hydrap.getPreviousMicroConfiguration( 0 ) );

        BOOST_CHECK_NO_THROW( Fscm = hydram.getPreviousMicroConfiguration( 0 ) );

        for ( unsigned int j = 0; j < ( dimension * dimension ); j++ ){

            gradient[ j ][ i ] = ( Fscp[ j ] - Fscm[ j ] ) / ( 2 * delta[ i + ( numConfigurations - 1 ) * dimension * dimension ] );

        }


    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *hydra.get_previousdChi1dChin( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_get_gradientMicroConfigurations, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 0.99524734, -0.25750656,  0.20796104,
                                        0.15868703,  0.88336002, -0.07056109,
                                        0.74019985, -0.45975665,  0.7056873 };

    floatVector previousDeformationGradient = {1, 0, 0,
                                               0, 1, 0,
                                               0, 0, 1 };

    floatVector microDeformation = { 0.67871774,  0.09730619,  0.90720157,
                                    -0.46243252,  0.93666086, -1.17911794,
                                     0.2516719 , -0.05288825,  1.17052434 };

    floatVector previousMicroDeformation = {1, 0, 0,
                                            0, 1, 0,
                                            0, 0, 1 };

    floatVector gradientMicroDeformation = { -0.07681577, -0.22420916,  0.00238935,  0.27122075,  0.42206819,
                                             -0.08200685, -0.44304236, -0.25073248,  0.38224183,  0.26594001,
                                              0.17706757,  0.15024769,  0.81775042, -0.31040183, -0.04897751,
                                             -0.00541132,  0.10023166,  0.06705184, -0.44877897,  0.61614687,
                                              0.15015164, -0.70225686,  0.41458656, -0.96940055,  0.08560181,
                                              1.35320298,  0.3626042 };

    floatVector previousGradientMicroDeformation( 27, 0 );

    floatVector previousStateVariables = { -0.10788248, -0.15682198,  0.22904971, -0.06142776, -0.4403221 ,
                                           -0.10195574,  0.23799541, -0.31750827, -0.32454824,  0.03155137,
                                            0.03182759,  0.13440096,  0.34943179,  0.22445532,  0.11102351,
                                            0.22244338, -0.17704109, -0.13821134, -0.07364869,  0.39338916,
                                            0.44416002,  0.00183668,  0.12395295, -0.3843816 , -0.18271452,
                                           -0.08517379,  0.36630916, -0.24954463, -0.01696574,  0.48555979,
                                            0.01948512,  0.11289453, -0.37937133,  0.3263408 ,  0.10306013,
                                            0.04506801,  0.1919703 ,  0.05438325, -0.11104943,  0.42513249,
                                            0.34167   , -0.14260243, -0.45640854, -0.19523193, -0.10181432,
                                            0.20495883,  0.49535848, -0.14408513,  0.26254781,  0.09317692,
                                            0.1917018 , -0.34887255, -0.10112371, -0.2591441 , -0.15654399,
                                            0.01312815,  0.16662455, -0.39409151, -0.36910505, -0.17801939,
                                            0.16156434,  0.34650623,  0.05325734,  0.35445249, -0.11516219,
                                           -0.1832121 , -0.14573532, -0.32891817,  0.32911263, -0.16132915,
                                            0.05237008,  0.07855147,  0.02153306, -0.49731194,  0.48834542,
                                            0.40534158, -0.29236414, -0.20751059,  0.02001015,  0.40191137,
                                            0.48363088, -0.24245794,  0.06435904,  0.30696868, -0.10562995,
                                            0.23107304, -0.33893099,  0.10069857,  0.36586446,  0.48352161 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 0;

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector answer1 = { -0.15723617, -0.19587921, -0.08297779,  0.18130077,  0.37545684,
                             0.01042234,  0.16931378,  0.08593655,  0.1249035 ,  0.17468905,
                             0.34234244, -0.41680501,  0.26368284, -0.25633363, -0.30577704,
                             0.07245696, -0.40428748,  0.38532683,  0.12724897,  0.22341636,
                            -0.48387079,  0.09443188,  0.05678519, -0.34104036, -0.34692948,
                             0.19552953, -0.18123357 };

    floatVector answer2 = { 0.1919703 ,  0.05438325, -0.11104943,  0.42513249,  0.34167   ,
                           -0.14260243, -0.45640854, -0.19523193, -0.10181432,  0.20495883,
                            0.49535848, -0.14408513,  0.26254781,  0.09317692,  0.1917018 ,
                           -0.34887255, -0.10112371, -0.2591441 , -0.15654399,  0.01312815,
                            0.16662455, -0.39409151, -0.36910505, -0.17801939,  0.16156434,
                            0.34650623,  0.05325734 };

    floatVector answer3 = { 0.35445249, -0.11516219, -0.1832121 , -0.14573532, -0.32891817,
                            0.32911263, -0.16132915,  0.05237008,  0.07855147,  0.02153306,
                           -0.49731194,  0.48834542,  0.40534158, -0.29236414, -0.20751059,
                            0.02001015,  0.40191137,  0.48363088, -0.24245794,  0.06435904,
                            0.30696868, -0.10562995,  0.23107304, -0.33893099,  0.10069857,
                            0.36586446,  0.48352161 };

    BOOST_TEST( tolerantCheck( answer1, floatVector( hydra.get_gradientMicroConfigurations( )->begin( ) + 0 * 27,
                                                     hydra.get_gradientMicroConfigurations( )->begin( ) + 1 * 27 ), 1e-5, 1e-5 ) );

    BOOST_TEST( answer2 == floatVector( hydra.get_gradientMicroConfigurations( )->begin( ) + 1 * 27,
                                        hydra.get_gradientMicroConfigurations( )->begin( ) + 2 * 27 ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer3 == floatVector( hydra.get_gradientMicroConfigurations( )->begin( ) + 2 * 27,
                                        hydra.get_gradientMicroConfigurations( )->begin( ) + 3 * 27 ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_get_gradientMicroConfigurations_jacobians, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 0.99524734, -0.25750656,  0.20796104,
                                        0.15868703,  0.88336002, -0.07056109,
                                        0.74019985, -0.45975665,  0.7056873 };

    floatVector previousDeformationGradient = {1, 0, 0,
                                               0, 1, 0,
                                               0, 0, 1 };

    floatVector microDeformation = { 0.67871774,  0.09730619,  0.90720157,
                                    -0.46243252,  0.93666086, -1.17911794,
                                     0.2516719 , -0.05288825,  1.17052434 };

    floatVector previousMicroDeformation = {1, 0, 0,
                                            0, 1, 0,
                                            0, 0, 1 };

    floatVector gradientMicroDeformation = { -0.07681577, -0.22420916,  0.00238935,  0.27122075,  0.42206819,
                                             -0.08200685, -0.44304236, -0.25073248,  0.38224183,  0.26594001,
                                              0.17706757,  0.15024769,  0.81775042, -0.31040183, -0.04897751,
                                             -0.00541132,  0.10023166,  0.06705184, -0.44877897,  0.61614687,
                                              0.15015164, -0.70225686,  0.41458656, -0.96940055,  0.08560181,
                                              1.35320298,  0.3626042 };

    floatVector previousGradientMicroDeformation( 27, 0 );

    floatVector previousStateVariables = { -0.10788248, -0.15682198,  0.22904971, -0.06142776, -0.4403221 ,
                                           -0.10195574,  0.23799541, -0.31750827, -0.32454824,  0.03155137,
                                            0.03182759,  0.13440096,  0.34943179,  0.22445532,  0.11102351,
                                            0.22244338, -0.17704109, -0.13821134, -0.07364869,  0.39338916,
                                            0.44416002,  0.00183668,  0.12395295, -0.3843816 , -0.18271452,
                                           -0.08517379,  0.36630916, -0.24954463, -0.01696574,  0.48555979,
                                            0.01948512,  0.11289453, -0.37937133,  0.3263408 ,  0.10306013,
                                            0.04506801,  0.1919703 ,  0.05438325, -0.11104943,  0.42513249,
                                            0.34167   , -0.14260243, -0.45640854, -0.19523193, -0.10181432,
                                            0.20495883,  0.49535848, -0.14408513,  0.26254781,  0.09317692,
                                            0.1917018 , -0.34887255, -0.10112371, -0.2591441 , -0.15654399,
                                            0.01312815,  0.16662455, -0.39409151, -0.36910505, -0.17801939,
                                            0.16156434,  0.34650623,  0.05325734,  0.35445249, -0.11516219,
                                           -0.1832121 , -0.14573532, -0.32891817,  0.32911263, -0.16132915,
                                            0.05237008,  0.07855147,  0.02153306, -0.49731194,  0.48834542,
                                            0.40534158, -0.29236414, -0.20751059,  0.02001015,  0.40191137,
                                            0.48363088, -0.24245794,  0.06435904,  0.30696868, -0.10562995,
                                            0.23107304, -0.33893099,  0.10069857,  0.36586446,  0.48352161 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 0;

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatType eps = 1e-6;

    unsigned int nterms = dimension * dimension * dimension;

    floatMatrix dGradChi1dF( nterms, floatVector( deformationGradient.size( ), 0 ) );

    floatMatrix dGradChi1dFn( nterms, floatVector( ( numConfigurations - 1 ) * deformationGradient.size( ), 0 ) );

    floatMatrix dGradChi1dChi( nterms, floatVector( microDeformation.size( ), 0 ) );

    floatMatrix dGradChi1dChin( nterms, floatVector( ( numConfigurations - 1 ) * microDeformation.size( ), 0 ) );

    floatMatrix dGradChi1dGradChi( nterms, floatVector( gradientMicroDeformation.size( ), 0 ) );

    floatMatrix dGradChi1dGradChin( nterms, floatVector( ( numConfigurations - 1 ) * gradientMicroDeformation.size( ), 0 ) );

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

        floatVector valp = floatVector( hydrap.get_gradientMicroConfigurations( )->begin( ),
                                        hydrap.get_gradientMicroConfigurations( )->begin( ) + 27 );

        floatVector valm = floatVector( hydram.get_gradientMicroConfigurations( )->begin( ),
                                        hydram.get_gradientMicroConfigurations( )->begin( ) + 27 );

        for ( unsigned int j = 0; j < nterms; j++ ){

            dGradChi1dF[ j ][ i ] = ( valp[ j ] - valm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dGradChi1dF ) == floatVector( nterms * deformationGradient.size( ), 0 ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ) * deformationGradient.size( ); i++ ){

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

        floatVector valp = floatVector( hydrap.get_gradientMicroConfigurations( )->begin( ),
                                        hydrap.get_gradientMicroConfigurations( )->begin( ) + 27 );

        floatVector valm = floatVector( hydram.get_gradientMicroConfigurations( )->begin( ),
                                        hydram.get_gradientMicroConfigurations( )->begin( ) + 27 );

        for ( unsigned int j = 0; j < nterms; j++ ){

            dGradChi1dFn[ j ][ i ] = ( valp[ j ] - valm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dGradChi1dFn ) == *hydra.get_dGradChi1dFn( ), CHECK_PER_ELEMENT );

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

        floatVector valp = floatVector( hydrap.get_gradientMicroConfigurations( )->begin( ),
                                        hydrap.get_gradientMicroConfigurations( )->begin( ) + 27 );

        floatVector valm = floatVector( hydram.get_gradientMicroConfigurations( )->begin( ),
                                        hydram.get_gradientMicroConfigurations( )->begin( ) + 27 );

        for ( unsigned int j = 0; j < nterms; j++ ){

            dGradChi1dChi[ j ][ i ] = ( valp[ j ] - valm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dGradChi1dChi ) == *hydra.get_dGradChi1dChi( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ) * microDeformation.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + ( numConfigurations - 1 ) * deformationGradient.size( ) ] = eps * std::fabs( previousStateVariables[ i  + ( numConfigurations - 1 ) * deformationGradient.size( ) ] ) + eps;

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

        floatVector valp = floatVector( hydrap.get_gradientMicroConfigurations( )->begin( ),
                                        hydrap.get_gradientMicroConfigurations( )->begin( ) + 27 );

        floatVector valm = floatVector( hydram.get_gradientMicroConfigurations( )->begin( ),
                                        hydram.get_gradientMicroConfigurations( )->begin( ) + 27 );

        for ( unsigned int j = 0; j < nterms; j++ ){

            dGradChi1dChin[ j ][ i ] = ( valp[ j ] - valm[ j ] ) / ( 2 * delta[ i + ( numConfigurations - 1 ) * deformationGradient.size( ) ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dGradChi1dChin ) == *hydra.get_dGradChi1dChin( ), CHECK_PER_ELEMENT );

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

        floatVector valp = floatVector( hydrap.get_gradientMicroConfigurations( )->begin( ),
                                        hydrap.get_gradientMicroConfigurations( )->begin( ) + 27 );

        floatVector valm = floatVector( hydram.get_gradientMicroConfigurations( )->begin( ),
                                        hydram.get_gradientMicroConfigurations( )->begin( ) + 27 );

        for ( unsigned int j = 0; j < nterms; j++ ){

            dGradChi1dGradChi[ j ][ i ] = ( valp[ j ] - valm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dGradChi1dGradChi ) == *hydra.get_dGradChi1dGradChi( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ) * gradientMicroDeformation.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + ( numConfigurations - 1 ) * ( deformationGradient.size( ) + microDeformation.size( ) ) ]
             = eps * std::fabs( previousStateVariables[ i + ( numConfigurations - 1 ) * ( deformationGradient.size( ) + microDeformation.size( ) ) ] ) + eps;

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

        floatVector valp = floatVector( hydrap.get_gradientMicroConfigurations( )->begin( ),
                                        hydrap.get_gradientMicroConfigurations( )->begin( ) + 27 );

        floatVector valm = floatVector( hydram.get_gradientMicroConfigurations( )->begin( ),
                                        hydram.get_gradientMicroConfigurations( )->begin( ) + 27 );

        for ( unsigned int j = 0; j < nterms; j++ ){

            dGradChi1dGradChin[ j ][ i ] = ( valp[ j ] - valm[ j ] ) / ( 2 * delta[ i + ( numConfigurations - 1 ) * ( deformationGradient.size( ) + microDeformation.size( ) ) ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChi1dGradChin ), *hydra.get_dGradChi1dGradChin( ), 1e-5, 1e-5 ) );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_get_previousGradientMicroConfigurations, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 0.99524734, -0.25750656,  0.20796104,
                                        0.15868703,  0.88336002, -0.07056109,
                                        0.74019985, -0.45975665,  0.7056873 };

    floatVector previousDeformationGradient = {1, 0, 0,
                                               0, 1, 0,
                                               0, 0, 1 };

    floatVector microDeformation = { 0.67871774,  0.09730619,  0.90720157,
                                    -0.46243252,  0.93666086, -1.17911794,
                                     0.2516719 , -0.05288825,  1.17052434 };

    floatVector previousMicroDeformation = {1, 0, 0,
                                            0, 1, 0,
                                            0, 0, 1 };

    floatVector gradientMicroDeformation = { -0.07681577, -0.22420916,  0.00238935,  0.27122075,  0.42206819,
                                             -0.08200685, -0.44304236, -0.25073248,  0.38224183,  0.26594001,
                                              0.17706757,  0.15024769,  0.81775042, -0.31040183, -0.04897751,
                                             -0.00541132,  0.10023166,  0.06705184, -0.44877897,  0.61614687,
                                              0.15015164, -0.70225686,  0.41458656, -0.96940055,  0.08560181,
                                              1.35320298,  0.3626042 };

    floatVector previousGradientMicroDeformation( 27, 0 );

    floatVector previousStateVariables = { -0.10788248, -0.15682198,  0.22904971, -0.06142776, -0.4403221 ,
                                           -0.10195574,  0.23799541, -0.31750827, -0.32454824,  0.03155137,
                                            0.03182759,  0.13440096,  0.34943179,  0.22445532,  0.11102351,
                                            0.22244338, -0.17704109, -0.13821134, -0.07364869,  0.39338916,
                                            0.44416002,  0.00183668,  0.12395295, -0.3843816 , -0.18271452,
                                           -0.08517379,  0.36630916, -0.24954463, -0.01696574,  0.48555979,
                                            0.01948512,  0.11289453, -0.37937133,  0.3263408 ,  0.10306013,
                                            0.04506801,  0.1919703 ,  0.05438325, -0.11104943,  0.42513249,
                                            0.34167   , -0.14260243, -0.45640854, -0.19523193, -0.10181432,
                                            0.20495883,  0.49535848, -0.14408513,  0.26254781,  0.09317692,
                                            0.1917018 , -0.34887255, -0.10112371, -0.2591441 , -0.15654399,
                                            0.01312815,  0.16662455, -0.39409151, -0.36910505, -0.17801939,
                                            0.16156434,  0.34650623,  0.05325734,  0.35445249, -0.11516219,
                                           -0.1832121 , -0.14573532, -0.32891817,  0.32911263, -0.16132915,
                                            0.05237008,  0.07855147,  0.02153306, -0.49731194,  0.48834542,
                                            0.40534158, -0.29236414, -0.20751059,  0.02001015,  0.40191137,
                                            0.48363088, -0.24245794,  0.06435904,  0.30696868, -0.10562995,
                                            0.23107304, -0.33893099,  0.10069857,  0.36586446,  0.48352161 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 0;

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatVector answer1 = { -2.96204445,  2.37087437,  3.42033939,  0.43233628, -2.39932744,
                            -2.27084966,  2.93494409, -0.63956264, -2.94263772,  0.7463599 ,
                            -0.71998485, -1.72887124, -0.81991449,  1.11505955,  1.46320729,
                            -0.5171075 , -0.46245997,  1.00725566,  1.75955058, -1.22811493,
                            -2.0935175 , -0.27832769,  1.41188723,  1.70606015, -1.30208954,
                            -0.31930701,  1.22169671 };

    floatVector answer2 = { 0.1919703 ,  0.05438325, -0.11104943,  0.42513249,  0.34167   ,
                           -0.14260243, -0.45640854, -0.19523193, -0.10181432,  0.20495883,
                            0.49535848, -0.14408513,  0.26254781,  0.09317692,  0.1917018 ,
                           -0.34887255, -0.10112371, -0.2591441 , -0.15654399,  0.01312815,
                            0.16662455, -0.39409151, -0.36910505, -0.17801939,  0.16156434,
                            0.34650623,  0.05325734 };

    floatVector answer3 = { 0.35445249, -0.11516219, -0.1832121 , -0.14573532, -0.32891817,
                            0.32911263, -0.16132915,  0.05237008,  0.07855147,  0.02153306,
                           -0.49731194,  0.48834542,  0.40534158, -0.29236414, -0.20751059,
                            0.02001015,  0.40191137,  0.48363088, -0.24245794,  0.06435904,
                            0.30696868, -0.10562995,  0.23107304, -0.33893099,  0.10069857,
                            0.36586446,  0.48352161 };

    BOOST_TEST( answer1 == floatVector( hydra.get_previousGradientMicroConfigurations( )->begin( ) + 0 * 27,
                                        hydra.get_previousGradientMicroConfigurations( )->begin( ) + 1 * 27 ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer2 == floatVector( hydra.get_previousGradientMicroConfigurations( )->begin( ) + 1 * 27,
                                        hydra.get_previousGradientMicroConfigurations( )->begin( ) + 2 * 27 ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer3 == floatVector( hydra.get_previousGradientMicroConfigurations( )->begin( ) + 2 * 27,
                                        hydra.get_previousGradientMicroConfigurations( )->begin( ) + 3 * 27 ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_tardigrade_hydraBaseMicromorphic_get_previousGradientMicroConfigurations_jacobians, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 0.99524734, -0.25750656,  0.20796104,
                                        0.15868703,  0.88336002, -0.07056109,
                                        0.74019985, -0.45975665,  0.7056873 };

    floatVector previousDeformationGradient = {1, 0, 0,
                                               0, 1, 0,
                                               0, 0, 1 };

    floatVector microDeformation = { 0.67871774,  0.09730619,  0.90720157,
                                    -0.46243252,  0.93666086, -1.17911794,
                                     0.2516719 , -0.05288825,  1.17052434 };

    floatVector previousMicroDeformation = {1, 0, 0,
                                            0, 1, 0,
                                            0, 0, 1 };

    floatVector gradientMicroDeformation = { -0.07681577, -0.22420916,  0.00238935,  0.27122075,  0.42206819,
                                             -0.08200685, -0.44304236, -0.25073248,  0.38224183,  0.26594001,
                                              0.17706757,  0.15024769,  0.81775042, -0.31040183, -0.04897751,
                                             -0.00541132,  0.10023166,  0.06705184, -0.44877897,  0.61614687,
                                              0.15015164, -0.70225686,  0.41458656, -0.96940055,  0.08560181,
                                              1.35320298,  0.3626042 };

    floatVector previousGradientMicroDeformation( 27, 0 );

    floatVector previousStateVariables = { -0.10788248, -0.15682198,  0.22904971, -0.06142776, -0.4403221 ,
                                           -0.10195574,  0.23799541, -0.31750827, -0.32454824,  0.03155137,
                                            0.03182759,  0.13440096,  0.34943179,  0.22445532,  0.11102351,
                                            0.22244338, -0.17704109, -0.13821134, -0.07364869,  0.39338916,
                                            0.44416002,  0.00183668,  0.12395295, -0.3843816 , -0.18271452,
                                           -0.08517379,  0.36630916, -0.24954463, -0.01696574,  0.48555979,
                                            0.01948512,  0.11289453, -0.37937133,  0.3263408 ,  0.10306013,
                                            0.04506801,  0.1919703 ,  0.05438325, -0.11104943,  0.42513249,
                                            0.34167   , -0.14260243, -0.45640854, -0.19523193, -0.10181432,
                                            0.20495883,  0.49535848, -0.14408513,  0.26254781,  0.09317692,
                                            0.1917018 , -0.34887255, -0.10112371, -0.2591441 , -0.15654399,
                                            0.01312815,  0.16662455, -0.39409151, -0.36910505, -0.17801939,
                                            0.16156434,  0.34650623,  0.05325734,  0.35445249, -0.11516219,
                                           -0.1832121 , -0.14573532, -0.32891817,  0.32911263, -0.16132915,
                                            0.05237008,  0.07855147,  0.02153306, -0.49731194,  0.48834542,
                                            0.40534158, -0.29236414, -0.20751059,  0.02001015,  0.40191137,
                                            0.48363088, -0.24245794,  0.06435904,  0.30696868, -0.10562995,
                                            0.23107304, -0.33893099,  0.10069857,  0.36586446,  0.48352161 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 0;

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
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    floatType eps = 1e-6;

    unsigned int nterms = dimension * dimension * dimension;

    floatMatrix dGradChi1dF( nterms, floatVector( deformationGradient.size( ), 0 ) );

    floatMatrix dGradChi1dFn( nterms, floatVector( ( numConfigurations - 1 ) * deformationGradient.size( ), 0 ) );

    floatMatrix dGradChi1dChi( nterms, floatVector( microDeformation.size( ), 0 ) );

    floatMatrix dGradChi1dChin( nterms, floatVector( ( numConfigurations - 1 ) * microDeformation.size( ), 0 ) );

    floatMatrix dGradChi1dGradChi( nterms, floatVector( gradientMicroDeformation.size( ), 0 ) );

    floatMatrix dGradChi1dGradChin( nterms, floatVector( ( numConfigurations - 1 ) * gradientMicroDeformation.size( ), 0 ) );

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

        floatVector valp = floatVector( hydrap.get_previousGradientMicroConfigurations( )->begin( ),
                                        hydrap.get_previousGradientMicroConfigurations( )->begin( ) + 27 );

        floatVector valm = floatVector( hydram.get_previousGradientMicroConfigurations( )->begin( ),
                                        hydram.get_previousGradientMicroConfigurations( )->begin( ) + 27 );

        for ( unsigned int j = 0; j < nterms; j++ ){

            dGradChi1dF[ j ][ i ] = ( valp[ j ] - valm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dGradChi1dF ) == floatVector( nterms * deformationGradient.size( ), 0 ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ) * deformationGradient.size( ); i++ ){

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

        floatVector valp = floatVector( hydrap.get_previousGradientMicroConfigurations( )->begin( ),
                                        hydrap.get_previousGradientMicroConfigurations( )->begin( ) + 27 );

        floatVector valm = floatVector( hydram.get_previousGradientMicroConfigurations( )->begin( ),
                                        hydram.get_previousGradientMicroConfigurations( )->begin( ) + 27 );

        for ( unsigned int j = 0; j < nterms; j++ ){

            dGradChi1dFn[ j ][ i ] = ( valp[ j ] - valm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dGradChi1dFn ) == *hydra.get_previousdGradChi1dFn( ), CHECK_PER_ELEMENT );

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

        floatVector valp = floatVector( hydrap.get_previousGradientMicroConfigurations( )->begin( ),
                                        hydrap.get_previousGradientMicroConfigurations( )->begin( ) + 27 );

        floatVector valm = floatVector( hydram.get_previousGradientMicroConfigurations( )->begin( ),
                                        hydram.get_previousGradientMicroConfigurations( )->begin( ) + 27 );

        for ( unsigned int j = 0; j < nterms; j++ ){

            dGradChi1dChi[ j ][ i ] = ( valp[ j ] - valm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dGradChi1dChi ) == *hydra.get_previousdGradChi1dChi( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ) * microDeformation.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + ( numConfigurations - 1 ) * deformationGradient.size( ) ] = eps * std::fabs( previousStateVariables[ i  + ( numConfigurations - 1 ) * deformationGradient.size( ) ] ) + eps;

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

        floatVector valp = floatVector( hydrap.get_previousGradientMicroConfigurations( )->begin( ),
                                        hydrap.get_previousGradientMicroConfigurations( )->begin( ) + 27 );

        floatVector valm = floatVector( hydram.get_previousGradientMicroConfigurations( )->begin( ),
                                        hydram.get_previousGradientMicroConfigurations( )->begin( ) + 27 );

        for ( unsigned int j = 0; j < nterms; j++ ){

            dGradChi1dChin[ j ][ i ] = ( valp[ j ] - valm[ j ] ) / ( 2 * delta[ i + ( numConfigurations - 1 ) * deformationGradient.size( ) ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dGradChi1dChin ) == *hydra.get_previousdGradChi1dChin( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < gradientMicroDeformation.size( ); i++ ){

        floatVector delta( gradientMicroDeformation.size( ), 0 );

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

        floatVector valp = floatVector( hydrap.get_previousGradientMicroConfigurations( )->begin( ),
                                        hydrap.get_previousGradientMicroConfigurations( )->begin( ) + 27 );

        floatVector valm = floatVector( hydram.get_previousGradientMicroConfigurations( )->begin( ),
                                        hydram.get_previousGradientMicroConfigurations( )->begin( ) + 27 );

        for ( unsigned int j = 0; j < nterms; j++ ){

            dGradChi1dGradChi[ j ][ i ] = ( valp[ j ] - valm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dGradChi1dGradChi ) == *hydra.get_previousdGradChi1dGradChi( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ) * gradientMicroDeformation.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + ( numConfigurations - 1 ) * ( deformationGradient.size( ) + microDeformation.size( ) ) ]
             = eps * std::fabs( previousStateVariables[ i + ( numConfigurations - 1 ) * ( deformationGradient.size( ) + microDeformation.size( ) ) ] ) + eps;

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

        floatVector valp = floatVector( hydrap.get_previousGradientMicroConfigurations( )->begin( ),
                                        hydrap.get_previousGradientMicroConfigurations( )->begin( ) + 27 );

        floatVector valm = floatVector( hydram.get_previousGradientMicroConfigurations( )->begin( ),
                                        hydram.get_previousGradientMicroConfigurations( )->begin( ) + 27 );

        for ( unsigned int j = 0; j < nterms; j++ ){

            dGradChi1dGradChin[ j ][ i ] = ( valp[ j ] - valm[ j ] ) / ( 2 * delta[ i + ( numConfigurations - 1 ) * ( deformationGradient.size( ) + microDeformation.size( ) ) ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dGradChi1dGradChin ), *hydra.get_previousdGradChi1dGradChin( ), 1e-5, 1e-5 ) );

}

BOOST_AUTO_TEST_CASE( test_residualBaseMicromorphic, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::residualBaseMicromorphic{

        public:
            using tardigradeHydra::residualBaseMicromorphic::residualBaseMicromorphic;

    };

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 0.99524734, -0.25750656,  0.20796104,
                                        0.15868703,  0.88336002, -0.07056109,
                                        0.74019985, -0.45975665,  0.7056873 };

    floatVector previousDeformationGradient = {1, 0, 0,
                                               0, 1, 0,
                                               0, 0, 1 };

    floatVector microDeformation = { 0.67871774,  0.09730619,  0.90720157,
                                    -0.46243252,  0.93666086, -1.17911794,
                                     0.2516719 , -0.05288825,  1.17052434 };

    floatVector previousMicroDeformation = {1, 0, 0,
                                            0, 1, 0,
                                            0, 0, 1 };

    floatVector gradientMicroDeformation = { -0.07681577, -0.22420916,  0.00238935,  0.27122075,  0.42206819,
                                             -0.08200685, -0.44304236, -0.25073248,  0.38224183,  0.26594001,
                                              0.17706757,  0.15024769,  0.81775042, -0.31040183, -0.04897751,
                                             -0.00541132,  0.10023166,  0.06705184, -0.44877897,  0.61614687,
                                              0.15015164, -0.70225686,  0.41458656, -0.96940055,  0.08560181,
                                              1.35320298,  0.3626042 };

    floatVector previousGradientMicroDeformation( 27, 0 );

    floatVector previousStateVariables = { -0.10788248, -0.15682198,  0.22904971, -0.06142776, -0.4403221 ,
                                           -0.10195574,  0.23799541, -0.31750827, -0.32454824,  0.03155137,
                                            0.03182759,  0.13440096,  0.34943179,  0.22445532,  0.11102351,
                                            0.22244338, -0.17704109, -0.13821134, -0.07364869,  0.39338916,
                                            0.44416002,  0.00183668,  0.12395295, -0.3843816 , -0.18271452,
                                           -0.08517379,  0.36630916, -0.24954463, -0.01696574,  0.48555979,
                                            0.01948512,  0.11289453, -0.37937133,  0.3263408 ,  0.10306013,
                                            0.04506801,  0.1919703 ,  0.05438325, -0.11104943,  0.42513249,
                                            0.34167   , -0.14260243, -0.45640854, -0.19523193, -0.10181432,
                                            0.20495883,  0.49535848, -0.14408513,  0.26254781,  0.09317692,
                                            0.1917018 , -0.34887255, -0.10112371, -0.2591441 , -0.15654399,
                                            0.01312815,  0.16662455, -0.39409151, -0.36910505, -0.17801939,
                                            0.16156434,  0.34650623,  0.05325734,  0.35445249, -0.11516219,
                                           -0.1832121 , -0.14573532, -0.32891817,  0.32911263, -0.16132915,
                                            0.05237008,  0.07855147,  0.02153306, -0.49731194,  0.48834542,
                                            0.40534158, -0.29236414, -0.20751059,  0.02001015,  0.40191137,
                                            0.48363088, -0.24245794,  0.06435904,  0.30696868, -0.10562995,
                                            0.23107304, -0.33893099,  0.10069857,  0.36586446,  0.48352161 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

    floatType tolr = 1e-2;

    floatType tola = 1e-3;

    unsigned int maxIterations = 24;

    unsigned int maxLSIterations = 45;

    floatType lsAlpha = 2.3;

    class hydraBaseMicromorphicMock : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

    };

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    residualMock residual( &hydra, 10 );

    BOOST_CHECK( &( *residual.hydra ) == &( *residual.residualBase::hydra ) );

}

BOOST_AUTO_TEST_CASE( test_updateUnknownVector, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    floatType time = 1.23;

    floatType deltaTime = 2.34;

    floatType temperature = 3.45;

    floatType previousTemperature = 4.56;

    floatVector deformationGradient = { 0.99524734, -0.25750656,  0.20796104,
                                        0.15868703,  0.88336002, -0.07056109,
                                        0.74019985, -0.45975665,  0.7056873 };

    floatVector previousDeformationGradient = {1, 0, 0,
                                               0, 1, 0,
                                               0, 0, 1 };

    floatVector microDeformation = { 0.67871774,  0.09730619,  0.90720157,
                                    -0.46243252,  0.93666086, -1.17911794,
                                     0.2516719 , -0.05288825,  1.17052434 };

    floatVector previousMicroDeformation = {1, 0, 0,
                                            0, 1, 0,
                                            0, 0, 1 };

    floatVector gradientMicroDeformation = { -0.07681577, -0.22420916,  0.00238935,  0.27122075,  0.42206819,
                                             -0.08200685, -0.44304236, -0.25073248,  0.38224183,  0.26594001,
                                              0.17706757,  0.15024769,  0.81775042, -0.31040183, -0.04897751,
                                             -0.00541132,  0.10023166,  0.06705184, -0.44877897,  0.61614687,
                                              0.15015164, -0.70225686,  0.41458656, -0.96940055,  0.08560181,
                                              1.35320298,  0.3626042 };

    floatVector previousGradientMicroDeformation( 27, 0 );

    floatVector previousStateVariables = { -0.10788248, -0.15682198,  0.22904971, -0.06142776, -0.4403221 ,
                                           -0.10195574,  0.23799541, -0.31750827, -0.32454824,  0.03155137,
                                            0.03182759,  0.13440096,  0.34943179,  0.22445532,  0.11102351,
                                            0.22244338, -0.17704109, -0.13821134, -0.07364869,  0.39338916,
                                            0.44416002,  0.00183668,  0.12395295, -0.3843816 , -0.18271452,
                                           -0.08517379,  0.36630916, -0.24954463, -0.01696574,  0.48555979,
                                            0.01948512,  0.11289453, -0.37937133,  0.3263408 ,  0.10306013,
                                            0.04506801,  0.1919703 ,  0.05438325, -0.11104943,  0.42513249,
                                            0.34167   , -0.14260243, -0.45640854, -0.19523193, -0.10181432,
                                            0.20495883,  0.49535848, -0.14408513,  0.26254781,  0.09317692,
                                            0.1917018 , -0.34887255, -0.10112371, -0.2591441 , -0.15654399,
                                            0.01312815,  0.16662455, -0.39409151, -0.36910505, -0.17801939,
                                            0.16156434,  0.34650623,  0.05325734,  0.35445249, -0.11516219,
                                           -0.1832121 , -0.14573532, -0.32891817,  0.32911263, -0.16132915,
                                            0.05237008,  0.07855147,  0.02153306, -0.49731194,  0.48834542,
                                            0.40534158, -0.29236414, -0.20751059,  0.02001015,  0.40191137,
                                            0.48363088, -0.24245794,  0.06435904,  0.30696868, -0.10562995,
                                            0.23107304, -0.33893099,  0.10069857,  0.36586446,  0.48352161,
                                            1, 2, 3, 4, 5, 6, 7 };

    floatVector parameters = { 0.1, 0.2, 0.3, 0.4 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 7;

    unsigned int dimension = 3;

    unsigned int configuration_unknown_count = 45;

    floatVector unknownVector( 3 * configuration_unknown_count + numNonLinearSolveStateVariables, 0 );

    floatVector previousAnswers = previousStateVariables;

    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        unknownVector[ i ] = 1e-2 * i + 1e-2;

    } 

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ); i++ ){

        // Make the deformation gradients invertable
        unknownVector[ 9 * i + 45     ] += 1;

        unknownVector[ 9 * i + 45 + 4 ] += 1;

        unknownVector[ 9 * i + 45 + 8 ] += 1;

        // Make the micro-deformations invertable
        unknownVector[ ( numConfigurations - 1 ) * 9 + 9 * i + 45     ] += 1;

        unknownVector[ ( numConfigurations - 1 ) * 9 + 9 * i + 45 + 4 ] += 1;

        unknownVector[ ( numConfigurations - 1 ) * 9 + 9 * i + 45 + 8 ] += 1;

        // Set up the previous answers
        previousAnswers[ 9 * i     ] += 1;

        previousAnswers[ 9 * i + 4 ] += 1;

        previousAnswers[ 9 * i + 8 ] += 1;

        previousAnswers[ ( numConfigurations - 1 ) * 9 + 9 * i     ] += 1;

        previousAnswers[ ( numConfigurations - 1 ) * 9 + 9 * i + 4 ] += 1;

        previousAnswers[ ( numConfigurations - 1 ) * 9 + 9 * i + 8 ] += 1;

    }

    floatType tolr = 1e-2;

    floatType tola = 1e-3;

    unsigned int maxIterations = 24;

    unsigned int maxLSIterations = 45;

    floatType lsAlpha = 2.3;

    class hydraBaseMicromorphicMock : public tardigradeHydra::hydraBaseMicromorphic{

        public:

            using tardigradeHydra::hydraBaseMicromorphic::hydraBaseMicromorphic;

            void callUpdateUnknownVector( const floatVector &newUnknownVector ){

                updateUnknownVector( newUnknownVector );

            }

    };

    hydraBaseMicromorphicMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                     microDeformation, previousMicroDeformation, gradientMicroDeformation, previousGradientMicroDeformation,
                                     { }, { },
                                     previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables,
                                     dimension, configuration_unknown_count,
                                     tolr, tola, maxIterations, maxLSIterations, lsAlpha );

    hydra.callUpdateUnknownVector( unknownVector );

    floatVector stressAnswer( unknownVector.begin( ), unknownVector.begin( ) + 45 );

    floatVector stateVariableAnswer( unknownVector.begin( ) + 45 + 2 * ( numConfigurations - 1 ) * 9 + ( numConfigurations - 1 ) * 27,
                                     unknownVector.begin( ) + 45 + 2 * ( numConfigurations - 1 ) * 9 + ( numConfigurations - 1 ) * 27 + numNonLinearSolveStateVariables );

    floatVector previousStateVariableAnswer( previousAnswers.begin( ) + 2 * ( numConfigurations - 1 ) * 9 + ( numConfigurations - 1 ) * 27,
                                             previousAnswers.begin( ) + 2 * ( numConfigurations - 1 ) * 9 + ( numConfigurations - 1 ) * 27 + numNonLinearSolveStateVariables );

    BOOST_TEST( *hydra.getStress( ) == stressAnswer, CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < ( numConfigurations - 1 ); i++ ){

        floatVector FiAnswer( unknownVector.begin( ) + 45 + i * 9, unknownVector.begin( ) + 45 + ( i + 1 ) * 9 );

        floatVector ChiiAnswer( unknownVector.begin( ) + 45 + ( numConfigurations - 1 ) * 9 + i * 9, unknownVector.begin( ) + 45 + ( numConfigurations - 1 ) * 9 + ( i + 1 ) * 9 );

        floatVector GradChiiAnswer( unknownVector.begin( ) + 45 + 2 * ( numConfigurations - 1 ) * 9 + i * 27,
                                    unknownVector.begin( ) + 45 + 2 * ( numConfigurations - 1 ) * 9 + ( i + 1 ) * 27 );

        floatVector previousFiAnswer( previousAnswers.begin( ) + i * 9,
                                      previousAnswers.begin( ) + ( i + 1 ) * 9 );

        floatVector previousChiiAnswer( previousAnswers.begin( ) + ( numConfigurations - 1 ) * 9 + i * 9,
                                        previousAnswers.begin( ) + ( numConfigurations - 1 ) * 9 + ( i + 1 ) * 9 );

        floatVector previousGradChiiAnswer( previousAnswers.begin( ) + 2 * ( numConfigurations - 1 ) * 9 + i * 27,
                                            previousAnswers.begin( ) + 2 * ( numConfigurations - 1 ) * 9 + ( i + 1 ) * 27 );

        BOOST_TEST( FiAnswer == hydra.getConfiguration( i + 1 ), CHECK_PER_ELEMENT );

        BOOST_TEST( FiAnswer == floatVector( hydra.get_configurations( )->begin( ) + 9 * ( i + 1 ),
                                             hydra.get_configurations( )->begin( ) + 9 * ( i + 2 ) ), CHECK_PER_ELEMENT );

        BOOST_TEST( ChiiAnswer == hydra.getMicroConfiguration( i + 1 ), CHECK_PER_ELEMENT );

        BOOST_TEST( ChiiAnswer == floatVector( hydra.get_microConfigurations( )->begin( ) + 9 * ( i + 1 ),
                                               hydra.get_microConfigurations( )->begin( ) + 9 * ( i + 2 ) ), CHECK_PER_ELEMENT );

        BOOST_TEST( GradChiiAnswer == floatVector( hydra.get_gradientMicroConfigurations( )->begin( ) + 27 * ( i + 1 ),
                                                   hydra.get_gradientMicroConfigurations( )->begin( ) + 27 * ( i + 2 ) ), CHECK_PER_ELEMENT );

        BOOST_TEST( previousFiAnswer == hydra.getPreviousConfiguration( i + 1 ), CHECK_PER_ELEMENT );

        BOOST_TEST( previousFiAnswer == floatVector( hydra.get_previousConfigurations( )->begin( ) + 9 * ( i + 1 ),
                                                     hydra.get_previousConfigurations( )->begin( ) + 9 * ( i + 2 ) ), CHECK_PER_ELEMENT );

        BOOST_TEST( previousChiiAnswer == hydra.getPreviousMicroConfiguration( i + 1 ), CHECK_PER_ELEMENT );

        BOOST_TEST( previousChiiAnswer == floatVector( hydra.get_previousMicroConfigurations( )->begin( ) + 9 * ( i + 1 ),
                                                       hydra.get_previousMicroConfigurations( )->begin( ) + 9 * ( i + 2 ) ), CHECK_PER_ELEMENT );

        BOOST_TEST( previousGradChiiAnswer == floatVector( hydra.get_previousGradientMicroConfigurations( )->begin( ) + 27 * ( i + 1 ),
                                                           hydra.get_previousGradientMicroConfigurations( )->begin( ) + 27 * ( i + 2 ) ), CHECK_PER_ELEMENT );

    }

    BOOST_TEST( stateVariableAnswer == *hydra.get_nonLinearSolveStateVariables( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousStateVariableAnswer == *hydra.get_previousNonLinearSolveStateVariables( ), CHECK_PER_ELEMENT );

}
