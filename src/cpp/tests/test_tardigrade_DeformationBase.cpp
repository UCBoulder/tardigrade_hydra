/**
  * \file test_tardigrade_DeformationBase.cpp
  *
  * Tests for tardigrade_DeformationBase
  */

#include<tardigrade_DeformationBase.h>

#define BOOST_TEST_MODULE test_tardigrade_DeformationBase
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node

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

BOOST_AUTO_TEST_CASE( test_denseMatrixMultiply, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class DeformationBase_mock : public tardigradeHydra::DeformationBase{

        public:

            DeformationBase_mock( ){ }

            void public_denseMatrixMultiply(
                const std::vector<double> &A,
                const std::vector<double> &B,
                std::vector<double> &C
            ){

                _denseMatrixMultiply<2, 3, 4>(
                    std::begin( A ), std::end( A ),
                    std::begin( B ), std::end( B ),
                    std::begin( C ), std::end( C )
                );

            }

    };

    std::vector< double > A = {
        1, 2, 3,
        4, 5, 6
    };

    std::vector< double > B = {
         4,  5,  6,  7,
         8,  9, 10, 11,
        12, 13, 14, 15
    };

    std::vector< double > answer = {
         56,  62,  68,  74,
         128, 143, 158, 173
    };

    std::vector< double > result( 8, -1 );

    DeformationBase_mock deformation;

    try{
    deformation.public_denseMatrixMultiply(
        A, B, result
    );
    }
    catch(std::exception &e){tardigradeErrorTools::printNestedExceptions(e), throw;}

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_getSubConfiguration, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    std::vector< double > configurations = {
        +3.649834609e-01, -6.490964877e-01, +6.310274768e-02, +6.365517419e-02, +2.688019171e-01,
        +1.698863588e+00, +4.489106497e-01, +2.220470214e-01, +4.448867651e-01, -3.540821723e-01,
        +7.235773112e-01, -5.434735382e-01, -4.125719072e-01, +2.619522477e-01, -8.157901201e-01,
        +8.674023454e-01, +8.617255267e-01, -1.262980470e-02, -1.483394194e-01, -3.754775541e-01,
        -1.472973861e-01, +1.786778326e+00, +8.883200364e-01, +3.673351769e-03, +2.479059036e-01,
        -7.687632098e-01, +6.345709636e-01, -1.703475761e-01, +7.326183158e-01, -4.990892692e-01,
        -3.393147147e-02, +1.971119571e+00, +1.038970239e+00, +2.257890515e-01, -7.587426680e-01,
        +6.526816010e-01, +2.061202568e-01, +1.090136013e+00, -3.144723325e-01, -3.917584219e-01,
        -1.659555780e-01, +3.626015316e-01, +1.750913684e+00, +2.084467496e-02, +3.386275659e-01,
        +1.718731051e-01, +2.498070042e-01, +1.349378102e+00, +1.684684875e+00, -8.336100233e-01,
        +5.273656829e-01, -5.126672509e-01, -6.115540788e-01, +1.144913915e+00, -8.085749668e-01,
        +7.706536526e-01, +2.544979441e-01, +4.468327164e-01, +3.225841339e-02, +1.888637589e-01,
        +1.135703848e-01, -6.820807117e-01, -6.938589698e-01, +1.391059058e+00, +6.375328528e-01,
        +3.839405911e-01, +1.087664994e-01, -2.220988518e-01, +8.502649792e-01, +1.683339994e+00,
        -2.852048666e-01, -9.128170724e-01, -3.904638532e-01, -2.036286362e-01, +1.409917661e+00,
        +9.907169641e-01, -2.881702686e-01, +5.250956276e-01, +1.863538331e-01, +1.383403597e+00
    };

    std::vector< double > answer = {
        -1.854335085e+00, -4.818731479e+00, +2.023409992e+00, +2.560814284e+00,
        +3.876505881e+00, +9.445160110e+00, -4.099758968e+00, -3.346599356e+00,
        -9.017019788e-02, -1.845898419e+00, +1.817802803e+00, -1.782496383e+00,
        -7.111239393e-01, +4.733923547e-01, -1.671751118e+00, +4.093105290e+00
    };

    std::vector< double > result( 16, 0 );

    tardigradeHydra::DeformationBase deformation;

    deformation.getSubConfiguration<4>(
        std::begin( configurations ), std::end( configurations ), std::begin( result ), std::end( result )
    );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_getLeadingSubConfigurationJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    std::vector< double > configurations = {
        +3.649834609e-01, -6.490964877e-01, +6.310274768e-02, +6.365517419e-02, +2.688019171e-01,
        +1.698863588e+00, +4.489106497e-01, +2.220470214e-01, +4.448867651e-01, -3.540821723e-01,
        +7.235773112e-01, -5.434735382e-01, -4.125719072e-01, +2.619522477e-01, -8.157901201e-01,
        +8.674023454e-01, +8.617255267e-01, -1.262980470e-02, -1.483394194e-01, -3.754775541e-01,
        -1.472973861e-01, +1.786778326e+00, +8.883200364e-01, +3.673351769e-03, +2.479059036e-01,
        -7.687632098e-01, +6.345709636e-01, -1.703475761e-01, +7.326183158e-01, -4.990892692e-01,
        -3.393147147e-02, +1.971119571e+00, +1.038970239e+00, +2.257890515e-01, -7.587426680e-01,
        +6.526816010e-01, +2.061202568e-01, +1.090136013e+00, -3.144723325e-01, -3.917584219e-01,
        -1.659555780e-01, +3.626015316e-01, +1.750913684e+00, +2.084467496e-02, +3.386275659e-01,
        +1.718731051e-01, +2.498070042e-01, +1.349378102e+00, +1.684684875e+00, -8.336100233e-01,
        +5.273656829e-01, -5.126672509e-01, -6.115540788e-01, +1.144913915e+00, -8.085749668e-01,
        +7.706536526e-01, +2.544979441e-01, +4.468327164e-01, +3.225841339e-02, +1.888637589e-01,
        +1.135703848e-01, -6.820807117e-01, -6.938589698e-01, +1.391059058e+00, +6.375328528e-01,
        +3.839405911e-01, +1.087664994e-01, -2.220988518e-01, +8.502649792e-01, +1.683339994e+00,
        -2.852048666e-01, -9.128170724e-01, -3.904638532e-01, -2.036286362e-01, +1.409917661e+00,
        +9.907169641e-01, -2.881702686e-01, +5.250956276e-01, +1.863538331e-01, +1.383403597e+00
    };

    double eps = 1e-6;
    unsigned int configuration = 0;

    std::vector< double > jacobian( 256, 0 );

    tardigradeHydra::DeformationBase deformation;

    deformation.getLeadingSubConfigurationJacobian<4>(
        std::begin( configurations ), std::end( configurations ), std::begin( jacobian ), std::end( jacobian )
    );

    {

        constexpr unsigned int NUM_UNKNOWNS = 16;
        constexpr unsigned int NUM_OUTPUTS  = 16;
        std::vector< double > x( std::begin( configurations ), std::end( configurations ) );

        for ( unsigned int i = 0; i < NUM_UNKNOWNS; ++i ){

            double delta = eps * std::fabs( x[ i + 16 * configuration ] ) + eps;

            std::vector< double > xp = x;
            std::vector< double > xm = x;

            xp[ i + 16 * configuration ] += delta;
            xm[ i + 16 * configuration ] -= delta;

            std::vector< double > rp( 16 );
            std::vector< double > rm( 16 );

            deformation.getSubConfiguration<4>(
                std::begin( xp ), std::end( xp ), std::begin( rp ), std::end( rp )
            );
            deformation.getSubConfiguration<4>(
                std::begin( xm ), std::end( xm ), std::begin( rm ), std::end( rm )
            );

            for ( unsigned int j = 0; j < NUM_OUTPUTS; ++j ){

                BOOST_TEST( jacobian[ NUM_OUTPUTS * j + i ] == ( rp[ j ] - rm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_getTrailingSubConfigurationJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    std::vector< double > configurations = {
        +3.649834609e-01, -6.490964877e-01, +6.310274768e-02, +6.365517419e-02, +2.688019171e-01,
        +1.698863588e+00, +4.489106497e-01, +2.220470214e-01, +4.448867651e-01, -3.540821723e-01,
        +7.235773112e-01, -5.434735382e-01, -4.125719072e-01, +2.619522477e-01, -8.157901201e-01,
        +8.674023454e-01, +8.617255267e-01, -1.262980470e-02, -1.483394194e-01, -3.754775541e-01,
        -1.472973861e-01, +1.786778326e+00, +8.883200364e-01, +3.673351769e-03, +2.479059036e-01,
        -7.687632098e-01, +6.345709636e-01, -1.703475761e-01, +7.326183158e-01, -4.990892692e-01,
        -3.393147147e-02, +1.971119571e+00, +1.038970239e+00, +2.257890515e-01, -7.587426680e-01,
        +6.526816010e-01, +2.061202568e-01, +1.090136013e+00, -3.144723325e-01, -3.917584219e-01,
        -1.659555780e-01, +3.626015316e-01, +1.750913684e+00, +2.084467496e-02, +3.386275659e-01,
        +1.718731051e-01, +2.498070042e-01, +1.349378102e+00, +1.684684875e+00, -8.336100233e-01,
        +5.273656829e-01, -5.126672509e-01, -6.115540788e-01, +1.144913915e+00, -8.085749668e-01,
        +7.706536526e-01, +2.544979441e-01, +4.468327164e-01, +3.225841339e-02, +1.888637589e-01,
        +1.135703848e-01, -6.820807117e-01, -6.938589698e-01, +1.391059058e+00, +6.375328528e-01,
        +3.839405911e-01, +1.087664994e-01, -2.220988518e-01, +8.502649792e-01, +1.683339994e+00,
        -2.852048666e-01, -9.128170724e-01, -3.904638532e-01, -2.036286362e-01, +1.409917661e+00,
        +9.907169641e-01, -2.881702686e-01, +5.250956276e-01, +1.863538331e-01, +1.383403597e+00
    };

    double eps = 1e-6;
    unsigned int configuration = 4;

    std::vector< double > jacobian( 256, 0 );

    tardigradeHydra::DeformationBase deformation;

    deformation.getTrailingSubConfigurationJacobian<4>(
        std::begin( configurations ), std::end( configurations ), std::begin( jacobian ), std::end( jacobian )
    );

    {

        constexpr unsigned int NUM_UNKNOWNS = 16;
        constexpr unsigned int NUM_OUTPUTS  = 16;
        std::vector< double > x( std::begin( configurations ), std::end( configurations ) );

        for ( unsigned int i = 0; i < NUM_UNKNOWNS; ++i ){

            double delta = eps * std::fabs( x[ i + 16 * configuration ] ) + eps;

            std::vector< double > xp = x;
            std::vector< double > xm = x;

            xp[ i + 16 * configuration ] += delta;
            xm[ i + 16 * configuration ] -= delta;

            std::vector< double > rp( 16 );
            std::vector< double > rm( 16 );

            deformation.getSubConfiguration<4>(
                std::begin( xp ), std::end( xp ), std::begin( rp ), std::end( rp )
            );
            deformation.getSubConfiguration<4>(
                std::begin( xm ), std::end( xm ), std::begin( rm ), std::end( rm )
            );

            for ( unsigned int j = 0; j < NUM_OUTPUTS; ++j ){

                BOOST_TEST( jacobian[ NUM_OUTPUTS * j + i ] == ( rp[ j ] - rm[ j ] ) / ( 2 * delta ) );

            }

        }

    }

}
