/**
  * \file test_tardigrade_SolverBase.cpp
  *
  * Tests for tardigrade_SolverBase
  */

#include"tardigrade_SolverBase.h"
#include"tardigrade_ResidualBase.h"
#include"tardigrade_hydra.h"
#include"tardigrade_TrialStepBase.h"
#include"tardigrade_ArmijoGradientDamping.h"

#define BOOST_TEST_MODULE test_tardigrade_SolverBase
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

        class SolverBaseTester{

            public:

                static void checkRankDeficientError( SolverBase &solver ){

                    BOOST_CHECK( solver._rank_deficient_error == solver.getRankDeficientError( ) );

                }

        };

        class hydraBaseTester{

            public:

                static void set_residual( hydraBase &hydra, const tardigradeHydra::floatVector &value ){

                    hydra._residual.second = value;
                    hydra._residual.first = true;

                    hydra.addIterationData( &hydra._residual );

                }

                static void set_unknownVector( hydraBase &hydra, const tardigradeHydra::floatVector &value ){

                    hydra._X.second = value;
                    hydra._X.first = true;

                }

                static void set_flatJacobian( hydraBase &hydra, const floatVector &value ){

                    hydra._jacobian.second = value;
                    hydra._jacobian.first = true;

                    hydra.addIterationData( &hydra._jacobian );
                }

                static void resetIterationData( hydraBase &hydra ){

                    hydra.resetIterationData( );

                }

        };

    }

}

BOOST_AUTO_TEST_CASE( test_SolverBase_getRankDeficientError, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::SolverBase solver;

    tardigradeHydra::unit_test::SolverBaseTester::checkRankDeficientError( solver );

}

BOOST_AUTO_TEST_CASE( test_SolverBase_setRankDeficientError, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::SolverBase solver;

    solver.setRankDeficientError( false );

    BOOST_TEST( false == solver.getRankDeficientError( ) );

}
