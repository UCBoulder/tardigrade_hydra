/**
  * \file test_tardigrade_hydraMassChange.cpp
  *
  * Tests for tardigrade_hydraMassChange
  */

#include<tardigrade_hydraMassChange.h>
#include<tardigrade_constitutive_tools.h>

#define BOOST_TEST_MODULE test_tardigrade_hydraMassChange
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::massChange::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::massChange::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::massChange::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type


namespace tardigradeHydra{

    namespace unit_test{

        class hydraBaseTester{

            public:

                static void updateUnknownVector( tardigradeHydra::hydraBase &hydra, const floatVector &value ){

                    BOOST_CHECK_NO_THROW( hydra.updateUnknownVector( value ) );

                }

        };

    }

    namespace massChange{

        namespace unit_test{

            class residualTester{

                public:

                    static void runBasicGetTests( tardigradeHydra::massChange::residual &R ){

                        BOOST_CHECK( &R._massChangeConfigurationIndex == R.getMassChangeConfigurationIndex( ) );

                    }

            };

        }

    }

}

BOOST_AUTO_TEST_CASE( test_residual_basicGetTests ){

    BOOST_TEST( true );

}
