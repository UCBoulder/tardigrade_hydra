
/**
  * \file test_tardigrade_SolverStepBase.cpp
  *
  * Tests for tardigrade_SolverStepBase
  */

#include"tardigrade_SolverStepBase.h"

#define BOOST_TEST_MODULE test_tardigrade_SolverStepBase
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

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


namespace tardigradeHydra{

    namespace unit_test{

        class SolverStepBaseTester{

            public:

                static void checkMuk( SolverStepBase &step ){

                    BOOST_CHECK( step._mu_k == step.getMuk( ) );

                }

                static void checkLMMu( SolverStepBase &step ){

                    BOOST_CHECK( step._lm_mu == step.getLMMu( ) );

                }


        };

    }


}

BOOST_AUTO_TEST_CASE( test_hydraBase_getMuk, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::SolverStepBase step;

    tardigradeHydra::unit_test::SolverStepBaseTester::checkMuk( step );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getLMMu, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::SolverStepBase step;

    tardigradeHydra::unit_test::SolverStepBaseTester::checkLMMu( step );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setMuk, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::SolverStepBase step;

    step.setMuk( 123.4 );

    BOOST_TEST( 123.4 == step.getMuk( ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setLMMu, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::SolverStepBase step;

    step.setLMMu( 123.4 );

    BOOST_TEST( 123.4 == step.getLMMu( ) );

}

