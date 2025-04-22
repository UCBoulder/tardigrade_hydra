/**
  * \file test_tardigrade_hydraLinearInternalEnergy.cpp
  * 
  * Tests for tardigrade_hydraLinearInternalEnergy
  */

#include<tardigrade_hydraLinearInternalEnergy.h>

#define BOOST_TEST_MODULE test_tardigrade_hydraLinearInternalEnergy
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::linearInternalEnergy::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::linearInternalEnergy::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::linearInternalEnergy::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

namespace tardigradeHydra{

    namespace unit_test{

        class hydraBaseTester{

            public:

                static void updateUnknownVector( tardigradeHydra::hydraBase &hydra, const floatVector &value ){

                    BOOST_CHECK_NO_THROW( hydra.updateUnknownVector( value ) );

                }

        };

    }

    namespace linearInternalEnergy{

        namespace unit_test{
    
            class residualTester{
    
                public:
    
                    static void runBasicGetTests( tardigradeHydra::linearInternalEnergy::residual &R ){

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

    tardigradeHydra::linearInternalEnergy::residual R( &hydra, 1, parameters );

    tardigradeHydra::linearInternalEnergy::unit_test::residualTester::runBasicGetTests( R );

    BOOST_TEST( 1 == R.get_expectedParameterVectorSize( ) );

    BOOST_TEST( 9 == R.get_internalEnergyIndex( ) );

    BOOST_TEST( parameters == *R.get_specificHeatParameters( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_getSpecificHeat, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector additionalDOF = {
        +3.929384e-01, -4.277213e-01, -5.462971e-01, +1.026295e-01, +4.389379e-01,
        -1.537871e-01, +9.615284e-01, +3.696595e-01, -3.813620e-02, -2.157650e-01,
        -3.136440e-01, +4.580994e-01, -1.228555e-01, -8.806442e-01, -2.039115e-01,
        +4.759908e-01, -6.350165e-01, -6.490965e-01, +6.310275e-02, +6.365517e-02,
        +2.688019e-01, +6.988636e-01, +4.489106e-01, +2.220470e-01, +4.448868e-01,
        -3.540822e-01, -2.764227e-01, -5.434735e-01, -4.125719e-01, +2.619522e-01,
        -8.157901e-01, -1.325977e-01, -1.382745e-01, -1.262980e-02, -1.483394e-01,
        -3.754776e-01, -1.472974e-01, +7.867783e-01, +8.883200e-01, +3.673352e-03,
        +2.479059e-01, -7.687632e-01, -3.654290e-01, -1.703476e-01, +7.326183e-01,
        -4.990893e-01, -3.393147e-02, +9.711196e-01, +3.897024e-02, +2.257891e-01,
        -7.587427e-01, +6.526816e-01, +2.061203e-01, +9.013601e-02, -3.144723e-01
    };

    floatVector previousAdditionalDOF = {
        -3.917584e-01, -1.659556e-01, +3.626015e-01, +7.509137e-01, +2.084467e-02,
        +3.386276e-01, +1.718731e-01, +2.498070e-01, +3.493781e-01, +6.846849e-01,
        -8.336100e-01, +5.273657e-01, -5.126673e-01, -6.115541e-01, +1.449139e-01,
        -8.085750e-01, +7.706537e-01, +2.544979e-01, +4.468327e-01, -9.677416e-01,
        +1.888638e-01, +1.135704e-01, -6.820807e-01, -6.938590e-01, +3.910591e-01,
        -3.624671e-01, +3.839406e-01, +1.087665e-01, -2.220989e-01, +8.502650e-01,
        +6.833400e-01, -2.852049e-01, -9.128171e-01, -3.904639e-01, -2.036286e-01,
        +4.099177e-01, +9.907170e-01, -2.881703e-01, +5.250956e-01, +1.863538e-01,
        +3.834036e-01, -6.977451e-01, -2.022474e-01, -5.182882e-01, -3.130880e-01,
        +2.625631e-02, +3.332491e-01, -7.881830e-01, -7.382101e-01, -3.560388e-01,
        +3.231287e-01, +6.930125e-01, +1.065147e-01, +7.089050e-01, -2.303244e-01
    };

    floatVector previousStateVariables = { };

    floatVector parameters = { 0.245 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                      additionalDOF, previousAdditionalDOF,
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearInternalEnergy::residual R( &hydra, 1, parameters );

    double answer = 0.245;

    BOOST_TEST( answer == *R.get_specificHeat( )         );

    BOOST_TEST( answer == *R.get_previousSpecificHeat( ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getInternalEnergy, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector additionalDOF = {
        +3.929384e-01, -4.277213e-01, -5.462971e-01, +1.026295e-01, +4.389379e-01,
        -1.537871e-01, +9.615284e-01, +3.696595e-01, -3.813620e-02, -2.157650e-01,
        -3.136440e-01, +4.580994e-01, -1.228555e-01, -8.806442e-01, -2.039115e-01,
        +4.759908e-01, -6.350165e-01, -6.490965e-01, +6.310275e-02, +6.365517e-02,
        +2.688019e-01, +6.988636e-01, +4.489106e-01, +2.220470e-01, +4.448868e-01,
        -3.540822e-01, -2.764227e-01, -5.434735e-01, -4.125719e-01, +2.619522e-01,
        -8.157901e-01, -1.325977e-01, -1.382745e-01, -1.262980e-02, -1.483394e-01,
        -3.754776e-01, -1.472974e-01, +7.867783e-01, +8.883200e-01, +3.673352e-03,
        +2.479059e-01, -7.687632e-01, -3.654290e-01, -1.703476e-01, +7.326183e-01,
        -4.990893e-01, -3.393147e-02, +9.711196e-01, +3.897024e-02, +2.257891e-01,
        -7.587427e-01, +6.526816e-01, +2.061203e-01, +9.013601e-02, -3.144723e-01
    };

    floatVector previousAdditionalDOF = {
        -3.917584e-01, -1.659556e-01, +3.626015e-01, +7.509137e-01, +2.084467e-02,
        +3.386276e-01, +1.718731e-01, +2.498070e-01, +3.493781e-01, +6.846849e-01,
        -8.336100e-01, +5.273657e-01, -5.126673e-01, -6.115541e-01, +1.449139e-01,
        -8.085750e-01, +7.706537e-01, +2.544979e-01, +4.468327e-01, -9.677416e-01,
        +1.888638e-01, +1.135704e-01, -6.820807e-01, -6.938590e-01, +3.910591e-01,
        -3.624671e-01, +3.839406e-01, +1.087665e-01, -2.220989e-01, +8.502650e-01,
        +6.833400e-01, -2.852049e-01, -9.128171e-01, -3.904639e-01, -2.036286e-01,
        +4.099177e-01, +9.907170e-01, -2.881703e-01, +5.250956e-01, +1.863538e-01,
        +3.834036e-01, -6.977451e-01, -2.022474e-01, -5.182882e-01, -3.130880e-01,
        +2.625631e-02, +3.332491e-01, -7.881830e-01, -7.382101e-01, -3.560388e-01,
        +3.231287e-01, +6.930125e-01, +1.065147e-01, +7.089050e-01, -2.303244e-01
    };

    floatVector previousStateVariables = { };

    floatVector parameters = { 0.245 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                      additionalDOF, previousAdditionalDOF,
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearInternalEnergy::residual R( &hydra, 1, parameters );

    double answer = 0.2009;

    double previousAnswer = 5.733;

    BOOST_TEST(         answer == *R.get_internalEnergy( )         );

    BOOST_TEST( previousAnswer == *R.get_previousInternalEnergy( ) );

    // Check the jacobians
    floatType eps = 1e-6;

    {

        double x_base = temperature;

        constexpr unsigned int NUM_VAR = 1;

        constexpr unsigned int NUM_OUT = 1;

        for ( unsigned int i = 0; i < NUM_VAR; ++i ){

            floatType delta = eps * std::fabs( x_base ) + eps;

            double x_p, x_m;

            x_p = x_base;
            x_m = x_base;

            x_p += delta;
            x_m -= delta;

            tardigradeHydra::hydraBase hydrap( time, deltaTime, x_p, previousTemperature, deformationGradient, previousDeformationGradient,
                                               additionalDOF, previousAdditionalDOF,
                                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

            tardigradeHydra::hydraBase hydram( time, deltaTime, x_m, previousTemperature, deformationGradient, previousDeformationGradient,
                                               additionalDOF, previousAdditionalDOF,
                                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

            tardigradeHydra::linearInternalEnergy::residual Rp( &hydrap, 1, parameters );

            tardigradeHydra::linearInternalEnergy::residual Rm( &hydram, 1, parameters );

            floatType rp = *Rp.get_internalEnergy( );
            floatType rm = *Rm.get_internalEnergy( );

            for ( unsigned int j = 0; j < NUM_OUT; ++j ){

                BOOST_TEST( ( *R.get_dInternalEnergydT( ) ) == ( rp - rm ) / ( 2 * delta ) );

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_residual_getResidual, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector additionalDOF = {
        +3.929384e-01, -4.277213e-01, -5.462971e-01, +1.026295e-01, +4.389379e-01,
        -1.537871e-01, +9.615284e-01, +3.696595e-01, -3.813620e-02, -2.157650e-01,
        -3.136440e-01, +4.580994e-01, -1.228555e-01, -8.806442e-01, -2.039115e-01,
        +4.759908e-01, -6.350165e-01, -6.490965e-01, +6.310275e-02, +6.365517e-02,
        +2.688019e-01, +6.988636e-01, +4.489106e-01, +2.220470e-01, +4.448868e-01,
        -3.540822e-01, -2.764227e-01, -5.434735e-01, -4.125719e-01, +2.619522e-01,
        -8.157901e-01, -1.325977e-01, -1.382745e-01, -1.262980e-02, -1.483394e-01,
        -3.754776e-01, -1.472974e-01, +7.867783e-01, +8.883200e-01, +3.673352e-03,
        +2.479059e-01, -7.687632e-01, -3.654290e-01, -1.703476e-01, +7.326183e-01,
        -4.990893e-01, -3.393147e-02, +9.711196e-01, +3.897024e-02, +2.257891e-01,
        -7.587427e-01, +6.526816e-01, +2.061203e-01, +9.013601e-02, -3.144723e-01
    };

    floatVector previousAdditionalDOF = {
        -3.917584e-01, -1.659556e-01, +3.626015e-01, +7.509137e-01, +2.084467e-02,
        +3.386276e-01, +1.718731e-01, +2.498070e-01, +3.493781e-01, +6.846849e-01,
        -8.336100e-01, +5.273657e-01, -5.126673e-01, -6.115541e-01, +1.449139e-01,
        -8.085750e-01, +7.706537e-01, +2.544979e-01, +4.468327e-01, -9.677416e-01,
        +1.888638e-01, +1.135704e-01, -6.820807e-01, -6.938590e-01, +3.910591e-01,
        -3.624671e-01, +3.839406e-01, +1.087665e-01, -2.220989e-01, +8.502650e-01,
        +6.833400e-01, -2.852049e-01, -9.128171e-01, -3.904639e-01, -2.036286e-01,
        +4.099177e-01, +9.907170e-01, -2.881703e-01, +5.250956e-01, +1.863538e-01,
        +3.834036e-01, -6.977451e-01, -2.022474e-01, -5.182882e-01, -3.130880e-01,
        +2.625631e-02, +3.332491e-01, -7.881830e-01, -7.382101e-01, -3.560388e-01,
        +3.231287e-01, +6.930125e-01, +1.065147e-01, +7.089050e-01, -2.303244e-01
    };

    floatVector unknownVector = {
        -3.664242e-01, -2.914706e-01, -6.578363e-01, +6.582253e-01, -3.226583e-01,
        +1.047402e-01, +1.571029e-01, +4.306612e-02, -9.946239e-01, +9.766908e-01,
        +8.106832e-01, -5.847283e-01, -4.150212e-01, +4.002031e-02, +8.038227e-01,
        +9.672618e-01, -4.849159e-01, +1.287181e-01, +6.139374e-01, -2.112599e-01,
        +4.621461e-01, -6.778620e-01, +2.013971e-01, +1.012700e-01, +2.274211e-01,
        -8.212532e-01
    };

    floatVector previousStateVariables( 17, 0 );

    floatVector parameters = { 0.245 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 17;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                      additionalDOF, previousAdditionalDOF,
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearInternalEnergy::residual R( &hydra, 1, parameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    std::vector< double > answer = {
         +7.757908e-01,
    };

    BOOST_TEST(         answer == *R.getResidual( ),         CHECK_PER_ELEMENT );

    // Check the jacobians
    floatType eps = 1e-6;

    // dRdF
    {

        std::vector< double > x_base = deformationGradient;

        constexpr unsigned int NUM_VAR = 9;

        constexpr unsigned int NUM_OUT = 1;

        for ( unsigned int i = 0; i < NUM_VAR; ++i ){

            floatType delta = eps * std::fabs( x_base[ i ] ) + eps;

            std::vector< double > x_p, x_m;

            x_p = x_base;
            x_m = x_base;

            x_p[ i ] += delta;
            x_m[ i ] -= delta;

            tardigradeHydra::hydraBase hydrap( time, deltaTime, temperature, previousTemperature, x_p, previousDeformationGradient,
                                               additionalDOF, previousAdditionalDOF,
                                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

            tardigradeHydra::hydraBase hydram( time, deltaTime, temperature, previousTemperature, x_m, previousDeformationGradient,
                                               additionalDOF, previousAdditionalDOF,
                                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

            tardigradeHydra::linearInternalEnergy::residual Rp( &hydrap, 1, parameters );

            tardigradeHydra::linearInternalEnergy::residual Rm( &hydram, 1, parameters );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

            floatVector rp = *Rp.getResidual( );
            floatVector rm = *Rm.getResidual( );

            for ( unsigned int j = 0; j < NUM_OUT; ++j ){

                BOOST_TEST( ( *R.getdRdF( ) )[ NUM_VAR * j + i ] == ( rp[ j ] - rm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // dRdT
    {

        double x_base = temperature;

        constexpr unsigned int NUM_VAR = 1;

        constexpr unsigned int NUM_OUT = 1;

        for ( unsigned int i = 0; i < NUM_VAR; ++i ){

            floatType delta = eps * std::fabs( temperature ) + eps;

            double x_p, x_m;

            x_p = x_base;
            x_m = x_base;

            x_p += delta;
            x_m -= delta;

            tardigradeHydra::hydraBase hydrap( time, deltaTime, x_p, previousTemperature, deformationGradient, previousDeformationGradient,
                                               additionalDOF, previousAdditionalDOF,
                                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

            tardigradeHydra::hydraBase hydram( time, deltaTime, x_m, previousTemperature, deformationGradient, previousDeformationGradient,
                                               additionalDOF, previousAdditionalDOF,
                                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

            tardigradeHydra::linearInternalEnergy::residual Rp( &hydrap, 1, parameters );

            tardigradeHydra::linearInternalEnergy::residual Rm( &hydram, 1, parameters );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

            floatVector rp = *Rp.getResidual( );
            floatVector rm = *Rm.getResidual( );

            for ( unsigned int j = 0; j < NUM_OUT; ++j ){

                BOOST_TEST( ( *R.getdRdT( ) )[ NUM_VAR * j + i ] == ( rp[ j ] - rm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // dRdX
    {

        std::vector< double > x_base = unknownVector;

        constexpr unsigned int NUM_VAR = 26;

        constexpr unsigned int NUM_OUT = 1;

        for ( unsigned int i = 0; i < NUM_VAR; ++i ){

            floatType delta = eps * std::fabs( x_base[ i ] ) + eps;

            std::vector< double > x_p, x_m;

            x_p = x_base;
            x_m = x_base;

            x_p[ i ] += delta;
            x_m[ i ] -= delta;

            tardigradeHydra::hydraBase hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                               additionalDOF, previousAdditionalDOF,
                                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

            tardigradeHydra::hydraBase hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                               additionalDOF, previousAdditionalDOF,
                                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

            tardigradeHydra::linearInternalEnergy::residual Rp( &hydrap, 1, parameters );

            tardigradeHydra::linearInternalEnergy::residual Rm( &hydram, 1, parameters );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, x_p );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, x_m );

            floatVector rp = *Rp.getResidual( );
            floatVector rm = *Rm.getResidual( );

            for ( unsigned int j = 0; j < NUM_OUT; ++j ){

                BOOST_TEST( ( *R.getJacobian( ) )[ NUM_VAR * j + i ] == ( rp[ j ] - rm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

    // dRdAdditionalDOF
    {

        constexpr unsigned int NUM_VAR = 55;

        constexpr unsigned int NUM_OUT = 1;

        for ( unsigned int i = 0; i < NUM_VAR; ++i ){

            floatType delta = eps * std::fabs( additionalDOF[ i ] ) + eps;

            floatVector x_p = additionalDOF;
            floatVector x_m = additionalDOF;

            x_p[ i ] += delta;
            x_m[ i ] -= delta;

            tardigradeHydra::hydraBase hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                               x_p, previousAdditionalDOF,
                                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

            tardigradeHydra::hydraBase hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                               x_m, previousAdditionalDOF,
                                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

            tardigradeHydra::linearInternalEnergy::residual Rp( &hydrap, 1, parameters );

            tardigradeHydra::linearInternalEnergy::residual Rm( &hydram, 1, parameters );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

            floatVector rp = *Rp.getResidual( );
            floatVector rm = *Rm.getResidual( );

            for ( unsigned int j = 0; j < NUM_OUT; ++j ){

                BOOST_TEST( ( *R.getdRdAdditionalDOF( ) )[ NUM_VAR * j + i ] == ( rp[ j ] - rm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}
