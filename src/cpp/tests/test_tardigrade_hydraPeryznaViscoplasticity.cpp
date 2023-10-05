/**
  * \file test_tardigrade-hydraPeryznaViscoplasticity.cpp
  *
  * Tests for tardigrade-hydraPeryznaViscoplasticity
  */

#include<tardigrade_hydraPeryznaViscoplasticity.h>
#include<tardigrade_hydraThermalExpansion.h>
#include<tardigrade_hydraLinearElasticity.h>
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

                public:

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

BOOST_AUTO_TEST_CASE( test_residual_basicGetTests ){

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 9, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = { 1.1, 0.1, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0 };

    floatVector previousDeformationGradient = { 1.0, 0.0, 0.0,
                                                0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0 };

    floatVector previousStateVariables = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           1, 2, 3, 4, 5, 7, 8, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 3, 3, 3, 3, 3, 3, 3, 3,
                                  4, 5, 6 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    try{
        hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );
    
        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );
    
        tardigradeHydra::peryznaViscoplasticity::unit_test::residualTester::runBasicGetTests( R );
    }
    catch( std::exception &e ){
        tardigradeErrorTools::printNestedExceptions( e );
        throw e;
    }

}

BOOST_AUTO_TEST_CASE( test_residual_getDrivingStress ){
    /*!
     * Test of computing the driving stress
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 9, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = { 1.1, 0.1, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0 };

    floatVector previousDeformationGradient = { 1.0, 0.0, 0.0,
                                                0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0 };

    floatVector previousStateVariables = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0.02, 0.03 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 6 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatVector answer = { 0.82636364, 1.0       , 0.75,
                           1.        , 1.21012101, 0.90759076,
                           0.75      , 0.90759076, 0.68069307 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getDrivingStress( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getFlowDirection ){
    /*!
     * Test of computing the driving stress
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

            floatVector drivingStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector flowParameters = { 1.23, 4.56 };

        private:

            virtual void setDrivingStress( ){

                tardigradeHydra::peryznaViscoplasticity::residual::setDrivingStress( drivingStress );

                setFlowParameters( flowParameters );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 9, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = { 1.1, 0.1, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0 };

    floatVector previousDeformationGradient = { 1.0, 0.0, 0.0,
                                                0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0 };

    floatVector previousStateVariables = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0, 0.01, 0.02, 0.03 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 6 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType dpYield;
    floatVector jac;
    floatVector answer;
    tardigradeStressTools::druckerPragerSurface( R.drivingStress, R.flowParameters[1], R.flowParameters[0], dpYield, jac, answer );

    std::cout << "answer: "; tardigradeVectorTools::print( answer );
    std::cout << "R.getFlowDirection: "; tardigradeVectorTools::print( *R.getFlowDirection( ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getFlowDirection( ) ) );

}
