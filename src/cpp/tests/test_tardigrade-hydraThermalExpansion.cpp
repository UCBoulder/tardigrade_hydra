/**
  * \file test_tardigrade-hydraThermalExpansion.cpp
  *
  * Tests for tardigrade-hydraThermalExpansion
  */

#include<tardigrade-hydraThermalExpansion.h>
#include<tardigrade-hydraLinearElasticity.h>
#include<constitutive_tools.h>

#define BOOST_TEST_MODULE test_tardigrade-hydraThermalExpansion
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef errorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::thermalExpansion::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::thermalExpansion::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::thermalExpansion::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type


namespace tardigradeHydra{

    namespace unit_test{

        class hydraBaseTester{

            public:

                static void updateUnknownVector( tardigradeHydra::hydraBase &hydra, const floatVector &value ){

                    BOOST_CHECK_NO_THROW( hydra.updateUnknownVector( value ) );

                }

        };

    }

    namespace thermalExpansion{

        namespace unit_test{

            class residualTester{

                public:

                    static void runBasicGetTests( tardigradeHydra::thermalExpansion::residual &R ){

                        BOOST_CHECK( &R._thermalConfigurationIndex == R.getThermalConfigurationIndex( ) );

                        BOOST_CHECK( &R._referenceTemperature == R.getReferenceTemperature( ) );

                        BOOST_CHECK( &R._linearParameters == R.getLinearParameters( ) );

                        BOOST_CHECK( &R._quadraticParameters == R.getQuadraticParameters( ) );

                        BOOST_CHECK( &R._thermalGreenLagrangeStrain.second == R.getThermalGreenLagrangeStrain( ) );

                        BOOST_CHECK( &R._thermalDeformationGradient.second == R.getThermalDeformationGradient( ) );

                        BOOST_CHECK( &R._dThermalGreenLagrangeStraindT.second == R.getdThermalGreenLagrangeStraindT( ) );

                        BOOST_CHECK( &R._dThermalDeformationGradientdT.second == R.getdThermalDeformationGradientdT( ) );

                    }

            };

        }

    }

}

BOOST_AUTO_TEST_CASE( test_residual_basicGetTests ){

    class residualMock : public tardigradeHydra::thermalExpansion::residual {

        public:

            using tardigradeHydra::thermalExpansion::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock thermalExpansion;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                thermalExpansion = residualMock( this, 9, 1, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &thermalExpansion;

                residuals[ 2 ] = &remainder;

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

    floatVector previousStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.thermalParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::thermalExpansion::unit_test::residualTester::runBasicGetTests( R );

}

BOOST_AUTO_TEST_CASE( test_residual_decomposeParameters ){

    class residualMock : public tardigradeHydra::thermalExpansion::residual {

        public:

            using tardigradeHydra::thermalExpansion::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock thermalExpansion;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                thermalExpansion = residualMock( this, 9, 1, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &thermalExpansion;

                residuals[ 2 ] = &remainder;

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

    floatVector previousStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    unsigned int thermalConfigurationIndex = 1;

    floatType Tref = 293.15;

    floatVector linearParameters = { 1e-5, 2e-5, 3e-5,
                                     2e-5, 4e-5, 5e-5,
                                     3e-5, 5e-5, 6e-5 };

    floatVector quadraticParameters = { 7e-5,  8e-5,  9e-5,
                                        8e-5, 10e-5, 11e-5,
                                        9e-5, 11e-5, 12e-5 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.thermalParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    BOOST_CHECK( vectorTools::fuzzyEquals( thermalConfigurationIndex, *R.getThermalConfigurationIndex( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( Tref, *R.getReferenceTemperature( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( linearParameters, *R.getLinearParameters( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( quadraticParameters, *R.getQuadraticParameters( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setThermalGreenLagrangeStrain ){

    class residualMock : public tardigradeHydra::thermalExpansion::residual {

        public:

            using tardigradeHydra::thermalExpansion::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock thermalExpansion;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                thermalExpansion = residualMock( this, 9, 1, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &thermalExpansion;

                residuals[ 2 ] = &remainder;

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

    floatVector previousStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    floatVector thermalGreenLagrangeStrain = { 0.28448393, 0.3251832 , 0.36588248,
                                               0.3251832 , 0.40658175, 0.44728103,
                                               0.36588248, 0.44728103, 0.4879803 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.thermalParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    BOOST_CHECK( vectorTools::fuzzyEquals( thermalGreenLagrangeStrain, *R.getThermalGreenLagrangeStrain( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setThermalDeformationGradient ){

    class residualMock : public tardigradeHydra::thermalExpansion::residual {

        public:

            using tardigradeHydra::thermalExpansion::residual::residual;

            floatVector thermalGreenLagrangeStrain;

            floatVector dThermalGreenLagrangeStraindT;

        private:

            virtual void setThermalGreenLagrangeStrain( ){

                tardigradeHydra::thermalExpansion::residual::setThermalGreenLagrangeStrain( thermalGreenLagrangeStrain );

            }

            virtual void setdThermalGreenLagrangeStraindT( ){

                tardigradeHydra::thermalExpansion::residual::setdThermalGreenLagrangeStraindT( dThermalGreenLagrangeStraindT );


            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock thermalExpansion;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                thermalExpansion = residualMock( this, 9, 1, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &thermalExpansion;

                residuals[ 2 ] = &remainder;

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

    floatVector previousStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    floatVector thermalDeformationGradient = { 1.56022038, 2.24932153, 2.55073813,
                                               2.24932153, 4.11095852, 4.80005966,
                                               2.55073813, 4.80005966, 6.36028004 };

    floatVector E;
    constitutiveTools::computeGreenLagrangeStrain( thermalDeformationGradient, E );

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.thermalParameters );

    R.thermalGreenLagrangeStrain = E;

    R.dThermalGreenLagrangeStraindT = floatVector( E.size( ), 0 );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    BOOST_CHECK( vectorTools::fuzzyEquals( thermalDeformationGradient, *R.getThermalDeformationGradient( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setResidual ){

    class residualMock : public tardigradeHydra::thermalExpansion::residual {

        public:

            using tardigradeHydra::thermalExpansion::residual::residual;

            floatVector thermalDeformationGradient;

        private:

            virtual void setThermalDeformationGradient( ){

                tardigradeHydra::thermalExpansion::residual::setThermalDeformationGradient( thermalDeformationGradient );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock thermalExpansion;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                thermalExpansion = residualMock( this, 9, 1, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &thermalExpansion;

                residuals[ 2 ] = &remainder;

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

    floatVector previousStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    floatVector thermalDeformationGradient = { 2, 3, 4, 5, 6, 7, 8, 9, 10 };

    floatVector residual = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.thermalParameters );

    R.thermalDeformationGradient = thermalDeformationGradient;

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    BOOST_CHECK( vectorTools::fuzzyEquals( residual, *R.getResidual( ) ) );

}
