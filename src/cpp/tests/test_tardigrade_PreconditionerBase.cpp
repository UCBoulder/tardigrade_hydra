/**
  * \file test_tardigrade_PreconditionerBase.cpp
  *
  * Tests for tardigrade_PreconditionerBase
  */

#include<tardigrade_PreconditionerBase.h>

#define BOOST_TEST_MODULE test_tardigrade_PreconditionerBase
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

        class PreconditionerBaseTester{

            public:

                static void checkGetFlatPreconditioner( PreconditionerBase &preconditioner ){

                    BOOST_CHECK( &preconditioner._preconditioner.second == preconditioner.getFlatPreconditioner( ) );

                }

                static void set_preconditioner_nohydra( PreconditionerBase &preconditioner, tardigradeHydra::floatVector &value ){

                    preconditioner._preconditioner.second = value;
                    preconditioner._preconditioner.first  = true;

                }

                static void set_preconditioner( PreconditionerBase &preconditioner, tardigradeHydra::floatVector &value ){

                    set_preconditioner_nohydra( preconditioner, value );

                    preconditioner.addIterationData( &preconditioner._preconditioner );

                }

        };

    }

}

BOOST_AUTO_TEST_CASE( test_PreconditionerBase_getUsePreconditioner, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::PreconditionerBase preconditioner;

    preconditioner._use_preconditioner = true;

    BOOST_TEST( preconditioner.getUsePreconditioner( ) );

    preconditioner._use_preconditioner = false;

    BOOST_TEST( !preconditioner.getUsePreconditioner( ) );

}

BOOST_AUTO_TEST_CASE( test_PreconditionerBase_getPreconditionerType, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::PreconditionerBase preconditioner;

    preconditioner._preconditioner_type = 123;

    BOOST_TEST( preconditioner.getPreconditionerType( ) == 123 );

}


BOOST_AUTO_TEST_CASE( test_PreconditionerBase_getFlatPreconditioner, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class PreconditionerMock : public tardigradeHydra::PreconditionerBase{

        public:

            using tardigradeHydra::PreconditionerBase::PreconditionerBase;

            tardigradeHydra::floatVector expected_value = { 1, 2, 3 };

            virtual void formPreconditioner( ) override{

                tardigradeHydra::unit_test::PreconditionerBaseTester::set_preconditioner_nohydra( *this, expected_value );

            }

    };

    PreconditionerMock preconditioner;

    tardigradeHydra::unit_test::PreconditionerBaseTester::checkGetFlatPreconditioner( preconditioner );

    std::vector< double > answer = { 1, 2, 3 };

    BOOST_TEST( answer == *preconditioner.getFlatPreconditioner( ), CHECK_PER_ELEMENT );

}
