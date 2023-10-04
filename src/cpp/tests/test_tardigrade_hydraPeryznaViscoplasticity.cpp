/**
  * \file test_tardigrade-hydraPeryznaViscoplasticity.cpp
  *
  * Tests for tardigrade-hydraPeryznaViscoplasticity
  */

#include<tardigrade_hydraPeryznaViscoplasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_stress_tools.h>

#define BOOST_TEST_MODULE test_tardigrade-hydraPeryznaViscoplasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::peryznaViscoplasticity::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::peryznaViscoplasticity::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::peryznaViscoplasticity::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

namespace tardigradeHydra{

    namespace unit_test{

        class hydraBaseTester{

            public:

                static void updateUnknownVector( tardigradeHydra::hydraBase &hydra, const floatVector &value ){

                    BOOST_CHECK_NO_THROW( hydra.updateUnknownVector( value ) );

                }

        };

    }

    namespace peryznaViscoplasticity{

        namespace unit_test{

            class residualTester{

                static void runBasicGetTests( tardigradeHydra::peryznaViscoplasticity::residual &R ){

                    BOOST_CHECK( &R._drivingStress.second == R.getDrivingStress( ) );

                    BOOST_CHECK( &R._flowDirection.second == R.getFlowDirection( ) );

                    BOOST_CHECK( &R._plasticMultiplier.second == R.getPlasticMultiplier( ) );

                    BOOST_CHECK( &R._velocityGradient.second == R.getVelocityGradient( ) );

                    BOOST_CHECK( &R._plasticDeformationGradient.second == R.getPlasticDeformationGradient( ) );

                    BOOST_CHECK( &R._stateVariables.second == R.getStateVariables( ) );

                    BOOST_CHECK( &R._peryznaParameters.second == R.getPeryznaParameters( ) );

                    BOOST_CHECK( &R._dragStressParameters.second == R.getDragStressParameters( ) );

                    BOOST_CHECK( &R._thermalParameters.second == R.getThermalParameters( ) );

                    BOOST_CHECK( &R._yieldParameters.second == R.getYieldParameters( ) );

                    BOOST_CHECK( &R._flowParameters.second == R.getFlowParameters( ) );

                    BOOST_CHECK( &R._mixingParameters.second == R.getMixingParameters( ) );

                }

            };

        }

    }

}
