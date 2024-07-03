/**
  * \file test_tardigrade_hydraThermalExpansion.cpp
  *
  * Tests for tardigrade_hydraThermalExpansion
  */

#include<tardigrade_hydraThermalExpansion.h>
#include<tardigrade_hydraLinearElasticity.h>
#include<tardigrade_constitutive_tools.h>

#define BOOST_TEST_MODULE test_tardigrade_hydraThermalExpansion
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
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

                        BOOST_CHECK( &R._thermalGreenLagrangeStrain.second == R.get_thermalGreenLagrangeStrain( ) );

                        BOOST_CHECK( &R._thermalDeformationGradient.second == R.get_thermalDeformationGradient( ) );

                        BOOST_CHECK( &R._dThermalGreenLagrangeStraindT.second == R.get_dThermalGreenLagrangeStraindT( ) );

                        BOOST_CHECK( &R._dThermalDeformationGradientdT.second == R.get_dThermalDeformationGradientdT( ) );

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

BOOST_AUTO_TEST_CASE( test_residual_basicGetTests, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

                setResidualClasses( residuals );

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.thermalParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::thermalExpansion::unit_test::residualTester::runBasicGetTests( R );

}

BOOST_AUTO_TEST_CASE( test_residual_decomposeParameters, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

                setResidualClasses( residuals );

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.thermalParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    BOOST_TEST( thermalConfigurationIndex == *R.getThermalConfigurationIndex( ) );

    BOOST_TEST( Tref == *R.getReferenceTemperature( ) );

    BOOST_TEST( linearParameters == *R.getLinearParameters( ), CHECK_PER_ELEMENT );

    BOOST_TEST( quadraticParameters == *R.getQuadraticParameters( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setThermalGreenLagrangeStrain, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

                setResidualClasses( residuals );

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.thermalParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    BOOST_TEST( thermalGreenLagrangeStrain == *R.get_thermalGreenLagrangeStrain( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setThermalGreenLagrangeStrainDerivatives, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

                setResidualClasses( residuals );

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.thermalParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatType eps = 1e-6;

    floatMatrix dThermalGreenLagrangeStraindF( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatVector dThermalGreenLagrangeStraindT( deformationGradient.size( ), 0 );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector deltas( deformationGradient.size( ), 0 );

        deltas[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + deltas, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - deltas, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
    
        residualMock Rp( &hydrap, 9, 1, hydra.thermalParameters );
        
        residualMock Rm( &hydram, 9, 1, hydra.thermalParameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dThermalGreenLagrangeStraindF[ j ][ i ] = ( ( *Rp.get_thermalGreenLagrangeStrain( ) )[ j ] - ( *Rm.get_thermalGreenLagrangeStrain( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector deltas( 1, 0 );

        deltas[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + deltas[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature - deltas[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
    
        residualMock Rp( &hydrap, 9, 1, hydra.thermalParameters );
        
        residualMock Rm( &hydram, 9, 1, hydra.thermalParameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dThermalGreenLagrangeStraindT[ j ] = ( ( *Rp.get_thermalGreenLagrangeStrain( ) )[ j ] - ( *Rm.get_thermalGreenLagrangeStrain( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dThermalGreenLagrangeStraindF ) == floatVector( 9 * 9, 0 ), CHECK_PER_ELEMENT );

    BOOST_TEST( dThermalGreenLagrangeStraindT == *R.get_dThermalGreenLagrangeStraindT( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setThermalDeformationGradient, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::thermalExpansion::residual {

        public:

            using tardigradeHydra::thermalExpansion::residual::residual;

            floatVector thermalGreenLagrangeStrain;

            floatVector dThermalGreenLagrangeStraindT;

        private:

            virtual void setThermalGreenLagrangeStrain( ){

                tardigradeHydra::thermalExpansion::residual::set_thermalGreenLagrangeStrain( thermalGreenLagrangeStrain );

            }

            virtual void setdThermalGreenLagrangeStraindT( ){

                tardigradeHydra::thermalExpansion::residual::set_dThermalGreenLagrangeStraindT( dThermalGreenLagrangeStraindT );


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

                setResidualClasses( residuals );

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
    tardigradeConstitutiveTools::computeGreenLagrangeStrain( thermalDeformationGradient, E );

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.thermalParameters );

    R.thermalGreenLagrangeStrain = E;

    R.dThermalGreenLagrangeStraindT = floatVector( E.size( ), 0 );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    BOOST_TEST( thermalDeformationGradient == *R.get_thermalDeformationGradient( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setThermalDeformationGradientDerivatives, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

                setResidualClasses( residuals );

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.thermalParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatType eps = 1e-6;

    floatMatrix dThermalDeformationGradientdF( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatVector dThermalDeformationGradientdT( deformationGradient.size( ), 0 );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector deltas( deformationGradient.size( ), 0 );

        deltas[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + deltas, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - deltas, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
    
        residualMock Rp( &hydrap, 9, 1, hydra.thermalParameters );
        
        residualMock Rm( &hydram, 9, 1, hydra.thermalParameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dThermalDeformationGradientdF[ j ][ i ] = ( ( *Rp.get_thermalDeformationGradient( ) )[ j ] - ( *Rm.get_thermalDeformationGradient( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector deltas( 1, 0 );

        deltas[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + deltas[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature - deltas[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
    
        residualMock Rp( &hydrap, 9, 1, hydra.thermalParameters );
        
        residualMock Rm( &hydram, 9, 1, hydra.thermalParameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dThermalDeformationGradientdT[ j ] = ( ( *Rp.get_thermalDeformationGradient( ) )[ j ] - ( *Rm.get_thermalDeformationGradient( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dThermalDeformationGradientdF ) == floatVector( 81, 0 ), CHECK_PER_ELEMENT );

    BOOST_TEST( dThermalDeformationGradientdT == *R.get_dThermalDeformationGradientdT( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setResidual, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::thermalExpansion::residual {

        public:

            using tardigradeHydra::thermalExpansion::residual::residual;

            floatVector thermalDeformationGradient;

        private:

            virtual void setThermalDeformationGradient( ){

                tardigradeHydra::thermalExpansion::residual::set_thermalDeformationGradient( thermalDeformationGradient );

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

                setResidualClasses( residuals );

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.thermalParameters );

    R.thermalDeformationGradient = thermalDeformationGradient;

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    BOOST_TEST( residual == *R.getResidual( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setResidualDerivatives, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

                setResidualClasses( residuals );

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.thermalParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatMatrix jacobian( deformationGradient.size( ), floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dRdF( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatVector dRdT( deformationGradient.size( ), 0 );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < hydra.getNumUnknowns( ); i++ ){

        floatVector deltas( hydra.getNumUnknowns( ), 0 );

        deltas[ i ] = eps * std::fabs( deltas[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + deltas );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - deltas );

        residualMock Rp( &hydrap, 9, 1, hydra.thermalParameters );

        residualMock Rm( &hydram, 9, 1, hydra.thermalParameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            jacobian[ j ][ i ] = ( ( *Rp.getResidual( ) )[ j ] - ( *Rm.getResidual( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector deltas( deformationGradient.size( ), 0 );

        deltas[ i ] = eps * std::fabs( deltas[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + deltas, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - deltas, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 1, hydra.thermalParameters );

        residualMock Rm( &hydram, 9, 1, hydra.thermalParameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dRdF[ j ][ i ] = ( ( *Rp.getResidual( ) )[ j ] - ( *Rm.getResidual( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector deltas( 1, 0 );

        deltas[ i ] = eps * std::fabs( deltas[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + deltas[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature - deltas[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 1, hydra.thermalParameters );

        residualMock Rm( &hydram, 9, 1, hydra.thermalParameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dRdT[ j ] = ( ( *Rp.getResidual( ) )[ j ] - ( *Rm.getResidual( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( jacobian ) == *R.getJacobian( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dRdF ) == *R.getdRdF( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dRdT == *R.getdRdT( ), CHECK_PER_ELEMENT );

}
