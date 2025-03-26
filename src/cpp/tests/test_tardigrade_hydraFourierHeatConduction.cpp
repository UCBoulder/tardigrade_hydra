/**
  * \file test_tardigrade_hydraFourierHeatConduction.cpp
  * 
  * Tests for tardigrade_hydraFourierHeatConduction
  */

#include<tardigrade_hydraFourierHeatConduction.h>

#define BOOST_TEST_MODULE test_tardigrade_hydraFourierHeatConduction
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::fourierHeatConduction::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::fourierHeatConduction::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::fourierHeatConduction::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

namespace tardigradeHydra{

    namespace unit_test{

        class hydraBaseTester{

            public:

                static void updateUnknownVector( tardigradeHydra::hydraBase &hydra, const floatVector &value ){

                    BOOST_CHECK_NO_THROW( hydra.updateUnknownVector( value ) );

                }

        };

    }

    namespace fourierHeatConduction{

        namespace unit_test{
    
            class residualTester{
    
                public:
    
                    static void runBasicGetTests( tardigradeHydra::fourierHeatConduction::residual &R ){

                    }

            };
    
        }

    }

}

bool tolerantCheck( const std::vector< double > &v1, const std::vector< double > &v2, double eps = 1e-6, double tol = 1e-9 ){

    if ( v1.size( ) != v2.size( ) ){

        return false;

    }

    BOOST_CHECK( v1.size( ) == v2.size( ) );

    const unsigned int len = v1.size( );

    for ( unsigned int i = 0; i < len; i++ ){

        if ( std::fabs( v1[ i ] ) < tol ){

            if ( std::fabs( v1[ i ] - v2[ i ] ) > eps ){

                return false;

            }

            BOOST_CHECK( std::fabs( v1[ i ] - v2[ i ] ) <= eps );

        }
        else{

            if ( ( std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v1[ i ] ) > eps ) ||
                 ( std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v2[ i ] ) > eps ) ){

                return false;

            }

            BOOST_TEST( v1[ i ] == v2[ i ] );

        }

    }

    return true;

}


BOOST_AUTO_TEST_CASE( test_residual_runBasicGetTests, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 0.82;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector additionalDOF( 155 );

    floatVector previousAdditionalDOF( 155 );

    floatVector previousStateVariables = { };

    floatVector parameters = { 0.245 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                      additionalDOF, previousAdditionalDOF,
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::fourierHeatConduction::residual R( &hydra, 23, parameters, 17 );

    tardigradeHydra::fourierHeatConduction::unit_test::residualTester::runBasicGetTests( R );

    BOOST_TEST(  1 == R.get_expectedParameterVectorSize( ) );

    BOOST_TEST( 17 == R.get_temperatureGradientIndex( ) );

    BOOST_TEST( parameters == *R.get_conductivityParameters( ), CHECK_PER_ELEMENT );

}
