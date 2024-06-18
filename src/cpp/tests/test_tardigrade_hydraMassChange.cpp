/**
  * \file test_tardigrade_hydraMassChange.cpp
  *
  * Tests for tardigrade_hydraMassChange
  */

#include<tardigrade_hydraLinearElasticity.h>
#include<tardigrade_hydraMassChange.h>
#include<tardigrade_constitutive_tools.h>

#define BOOST_TEST_MODULE test_tardigrade_hydraMassChange
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

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

                        BOOST_CHECK( &R._integrationParameter == R.getIntegrationParameter( ) );

                        BOOST_CHECK( &R._density.second         == R.get_density( ) );

                        BOOST_CHECK( &R._previousDensity.second == R.get_previousDensity( ) );

                        BOOST_CHECK( &R._massChangeRate.second         == R.get_massChangeRate( ) );

                        BOOST_CHECK( &R._previousMassChangeRate.second == R.get_previousMassChangeRate( ) );

                        BOOST_CHECK( &R._directionVector.second         == R.get_directionVector( ) );

                        BOOST_CHECK( &R._previousDirectionVector.second == R.get_previousDirectionVector( ) );

                    }

            };

        }

    }

}

void adaptive_tolerance_test( const floatVector &result, const floatVector &answer, floatType atol = 1e-9 ){
    /*!
     * Test which allows for very small values to be compared absolutely while larger values are compared relatively
     */

    BOOST_TEST( result.size( ) == answer.size( ) );
    for ( unsigned int i = 0; i < result.size( ); i++ ){
        if ( std::abs( answer[ i ] ) < atol ){
            BOOST_REQUIRE_SMALL( std::abs( answer[ i ] - result[ i ] ), atol );
        }
        else{
            BOOST_TEST( answer[ i ] == result[ i ] );
        }
    }
}

BOOST_AUTO_TEST_CASE( test_residual_basicGetTests, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

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

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { 0.11, 0.22, 0.33, 0.44, 0.55 }, { 0.12, 0.23, 0.36, 0.47, 0.58 },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.massChangeParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::massChange::unit_test::residualTester::runBasicGetTests( R );

    floatVector massChangeRateAnswer = { 0.33, 0.44, 0.55 };

    floatVector previousMassChangeRateAnswer = { 0.36, 0.47, 0.58 };

    BOOST_TEST( 0.11 == *R.get_density( ) );

    BOOST_TEST( 0.22 == *R.get_massChangeRate( ) );

    BOOST_TEST( massChangeRateAnswer == *R.get_directionVector( ), CHECK_PER_ELEMENT );

    BOOST_TEST( 0.12 == *R.get_previousDensity( ) );

    BOOST_TEST( 0.23 == *R.get_previousMassChangeRate( ) );

    BOOST_TEST( previousMassChangeRateAnswer == *R.get_previousDirectionVector( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeVelocityGradientTrace, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatType answer = 0.22 / 0.11;

    floatType previousAnswer = 0.23 / 0.12;

    residualMock R( &hydra, 9, 1, hydra.massChangeParameters );

    residualMock Rgrad( &hydra, 9, 1, hydra.massChangeParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    Rgrad.get_dMassChangeVelocityGradientTracedDensity( );
    Rgrad.get_dPreviousMassChangeVelocityGradientTracedPreviousDensity( );

    BOOST_TEST( answer == *R.get_massChangeVelocityGradientTrace( ) );

    BOOST_TEST( previousAnswer == *R.get_previousMassChangeVelocityGradientTrace( ) );

    BOOST_TEST( answer == *Rgrad.get_massChangeVelocityGradientTrace( ) );

    BOOST_TEST( previousAnswer == *Rgrad.get_previousMassChangeVelocityGradientTrace( ) );

    floatType eps = 1e-6;

    floatType dtrLdRho;

    floatType dtrLdC;

    floatType dPrevioustrLdPreviousRho;

    floatType dPrevioustrLdPreviousC;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * additionalDOF[ 0 ] + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;
        
        aDp[ 0 ] += delta;

        aDm[ 0 ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatType vp = *Rp.get_massChangeVelocityGradientTrace( );

        floatType vm = *Rm.get_massChangeVelocityGradientTrace( );

        dtrLdRho = ( vp - vm ) / ( 2 * delta );

    }

    BOOST_TEST( dtrLdRho == *R.get_dMassChangeVelocityGradientTracedDensity( ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * additionalDOF[ 1 ] + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;
        
        aDp[ 1 ] += delta;

        aDm[ 1 ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatType vp = *Rp.get_massChangeVelocityGradientTrace( );

        floatType vm = *Rm.get_massChangeVelocityGradientTrace( );

        dtrLdC = ( vp - vm ) / ( 2 * delta );

    }

    BOOST_TEST( dtrLdC == *R.get_dMassChangeVelocityGradientTracedMassChangeRate( ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * additionalDOF[ 0 ] + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;
        
        aDp[ 0 ] += delta;

        aDm[ 0 ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatType vp = *Rp.get_previousMassChangeVelocityGradientTrace( );

        floatType vm = *Rm.get_previousMassChangeVelocityGradientTrace( );

        dPrevioustrLdPreviousRho = ( vp - vm ) / ( 2 * delta );

    }

    BOOST_TEST( dPrevioustrLdPreviousRho == *R.get_dPreviousMassChangeVelocityGradientTracedPreviousDensity( ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * additionalDOF[ 1 ] + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;
        
        aDp[ 1 ] += delta;

        aDm[ 1 ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatType vp = *Rp.get_previousMassChangeVelocityGradientTrace( );

        floatType vm = *Rm.get_previousMassChangeVelocityGradientTrace( );

        dPrevioustrLdPreviousC = ( vp - vm ) / ( 2 * delta );

    }

    BOOST_TEST( dPrevioustrLdPreviousC == *R.get_dPreviousMassChangeVelocityGradientTracedPreviousMassChangeRate( ) );

}

BOOST_AUTO_TEST_CASE( test_residual_unitDirectionVector, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer = { 0.42426407, 0.56568542, 0.70710678 };

    floatVector previousAnswer = { 0.43436592, 0.56708884, 0.69981176 };

    residualMock R( &hydra, 9, 1, hydra.massChangeParameters );

    residualMock Rgrad( &hydra, 9, 1, hydra.massChangeParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    Rgrad.get_dUnitDirectionVectordDirectionVector( );

    Rgrad.get_dPreviousUnitDirectionVectordPreviousDirectionVector( );

    BOOST_TEST( answer == *R.get_unitDirectionVector( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *R.get_previousUnitDirectionVector( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *Rgrad.get_unitDirectionVector( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *Rgrad.get_previousUnitDirectionVector( ), CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    floatVector dDdGradC( 9, 0 );

    floatVector dPreviousDdPreviousGradC( 9, 0 );

    const unsigned int offset = 2;

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ offset + i ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;
        
        aDp[ offset + i ] += delta;

        aDm[ offset + i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_unitDirectionVector( );

        floatVector vm = *Rm.get_unitDirectionVector( );

        for ( unsigned int j = 0; j < 3; j++ ){

            dDdGradC[ 3 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dDdGradC == *R.get_dUnitDirectionVectordDirectionVector( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( previousAdditionalDOF[ offset + i ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;
        
        aDp[ offset + i ] += delta;

        aDm[ offset + i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousUnitDirectionVector( );

        floatVector vm = *Rm.get_previousUnitDirectionVector( );

        for ( unsigned int j = 0; j < 3; j++ ){

            dPreviousDdPreviousGradC[ 3 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dPreviousDdPreviousGradC == *R.get_dPreviousUnitDirectionVectordPreviousDirectionVector( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_unitDirectionVector2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.0, 0.0, 0.0 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.0, 0.0, 0.0 };

    floatVector previousStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer = { 0.0, 0.0, 0.0 };

    floatVector previousAnswer = { 0.0, 0.0, 0.0 };

    residualMock R( &hydra, 9, 1, hydra.massChangeParameters );

    residualMock Rgrad( &hydra, 9, 1, hydra.massChangeParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    Rgrad.get_dUnitDirectionVectordDirectionVector( );

    Rgrad.get_dPreviousUnitDirectionVectordPreviousDirectionVector( );

    BOOST_TEST( answer == *R.get_unitDirectionVector( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *R.get_previousUnitDirectionVector( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *Rgrad.get_unitDirectionVector( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *Rgrad.get_previousUnitDirectionVector( ), CHECK_PER_ELEMENT );

    floatVector dDdGradC( 9, 0 );

    floatVector dPreviousDdPreviousGradC( 9, 0 );

    BOOST_TEST( dDdGradC == *R.get_dUnitDirectionVectordDirectionVector( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dPreviousDdPreviousGradC == *R.get_dPreviousUnitDirectionVectordPreviousDirectionVector( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeVelocityGradient_1, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

            floatType massChangeVelocityGradientTrace = 0.4;

            floatType previousMassChangeVelocityGradientTrace = 0.9;

            floatType dtrLdRho = 0.123;

            floatType dtrLdC = 0.234;

            floatType previousdtrLdRho = 0.345;

            floatType previousdtrLdC = 0.456;

            floatVector unitDirectionVector = { 0.2, 0.3, 0.4 };

            floatVector previousUnitDirectionVector = { -0.1, 0.4, 0.30 };

            floatVector dDdGradC = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousdDdGradC = { .1, .2, .3, .4, .5, .6, .7, .8, .9 };

        protected:

            virtual void setMassChangeVelocityGradientTrace( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_previousMassChangeVelocityGradientTrace( previousMassChangeVelocityGradientTrace );

                }
                else{

                    set_massChangeVelocityGradientTrace( massChangeVelocityGradientTrace );

                }

            }

            virtual void setUnitDirectionVector( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_previousUnitDirectionVector( previousUnitDirectionVector );

                }
                else{

                    set_unitDirectionVector( unitDirectionVector );

                }

            }

            virtual void setUnitDirectionVectorDerivatives( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_dPreviousUnitDirectionVectordPreviousDirectionVector( previousdDdGradC );

                }
                else{

                    set_dUnitDirectionVectordDirectionVector( dDdGradC );

                }

            }

            virtual void setMassChangeVelocityGradientTraceDerivatives( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_dPreviousMassChangeVelocityGradientTracedPreviousDensity( previousdtrLdRho );

                    set_dPreviousMassChangeVelocityGradientTracedPreviousMassChangeRate( previousdtrLdC );

                }
                else{

                    set_dMassChangeVelocityGradientTracedDensity( dtrLdRho );

                    set_dMassChangeVelocityGradientTracedMassChangeRate( dtrLdC );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer = { 0.06977778, 0.013     , 0.01733333,
                           0.013     , 0.08061111, 0.026     ,
                           0.01733333, 0.026     , 0.09577778 };

    floatVector previousAnswer = { 0.142375, -0.0195  , -0.014625,
                                  -0.0195  ,  0.2155  ,  0.0585  ,
                                  -0.014625,  0.0585  ,  0.181375 };

    residualMock R( &hydra, 9, 1, hydra.massChangeParameters );

    residualMock Rgrad( &hydra, 9, 1, hydra.massChangeParameters );

    Rgrad.get_dMassChangeVelocityGradientdDensity( );

    Rgrad.get_dPreviousMassChangeVelocityGradientdPreviousDensity( );

    BOOST_TEST( answer == *R.get_massChangeVelocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *R.get_previousMassChangeVelocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *Rgrad.get_massChangeVelocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *Rgrad.get_previousMassChangeVelocityGradient( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeVelocityGradient_2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

            floatType massChangeVelocityGradientTrace = 0.4;

            floatType previousMassChangeVelocityGradientTrace = 0.9;

            floatType dtrLdRho = 0.123;

            floatType dtrLdC = 0.234;

            floatType previousdtrLdRho = 0.345;

            floatType previousdtrLdC = 0.456;

            floatVector unitDirectionVector = { 0.2, 0.3, 0.1 };

            floatVector previousUnitDirectionVector = { -0.1, 0.4, 0.30 };

            floatVector dDdGradC = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousdDdGradC = { .1, .2, .3, .4, .5, .6, .7, .8, .9 };

        protected:

            virtual void setMassChangeVelocityGradientTrace( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_previousMassChangeVelocityGradientTrace( previousMassChangeVelocityGradientTrace );

                }
                else{

                    set_massChangeVelocityGradientTrace( massChangeVelocityGradientTrace );

                }

            }

            virtual void setUnitDirectionVector( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_previousUnitDirectionVector( previousUnitDirectionVector );

                }
                else{

                    set_unitDirectionVector( unitDirectionVector );

                }

            }

            virtual void setUnitDirectionVectorDerivatives( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_dPreviousUnitDirectionVectordPreviousDirectionVector( previousdDdGradC );

                }
                else{

                    set_dUnitDirectionVectordDirectionVector( dDdGradC );

                }

            }

            virtual void setMassChangeVelocityGradientTraceDerivatives( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_dPreviousMassChangeVelocityGradientTracedPreviousDensity( previousdtrLdRho );

                    set_dPreviousMassChangeVelocityGradientTracedPreviousMassChangeRate( previousdtrLdC );

                }
                else{

                    set_dMassChangeVelocityGradientTracedDensity( dtrLdRho );

                    set_dMassChangeVelocityGradientTracedMassChangeRate( dtrLdC );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer = { 0.4, 0.0, 0.0,
                           0.0, 0.4, 0.0,
                           0.0, 0.0, 0.4 };

    floatVector previousAnswer = { 0.142375, -0.0195  , -0.014625,
                                  -0.0195  ,  0.2155  ,  0.0585  ,
                                  -0.014625,  0.0585  ,  0.181375 };

    residualMock R( &hydra, 9, 1, hydra.massChangeParameters );

    residualMock Rgrad( &hydra, 9, 1, hydra.massChangeParameters );

    Rgrad.get_dMassChangeVelocityGradientdDensity( );

    Rgrad.get_dPreviousMassChangeVelocityGradientdPreviousDensity( );

    BOOST_TEST( answer == *R.get_massChangeVelocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *R.get_previousMassChangeVelocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *Rgrad.get_massChangeVelocityGradient( ), CHECK_PER_ELEMENT );

    floatVector dtrLdRho = { R.dtrLdRho, 0, 0, 0, R.dtrLdRho, 0, 0, 0, R.dtrLdRho };

    floatVector dtrLdC = { R.dtrLdC, 0, 0, 0, R.dtrLdC, 0, 0, 0, R.dtrLdC };

    BOOST_TEST( dtrLdRho == *Rgrad.get_dMassChangeVelocityGradientdDensity( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dtrLdC == *Rgrad.get_dMassChangeVelocityGradientdMassChangeRate( ), CHECK_PER_ELEMENT );

    BOOST_TEST( floatVector( 27, 0 ) == *Rgrad.get_dMassChangeVelocityGradientdDirectionVector( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *Rgrad.get_previousMassChangeVelocityGradient( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeVelocityGradient_3, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

            floatType massChangeVelocityGradientTrace = 0.4;

            floatType previousMassChangeVelocityGradientTrace = 0.9;

            floatType dtrLdRho = 0.123;

            floatType dtrLdC = 0.234;

            floatType previousdtrLdRho = 0.345;

            floatType previousdtrLdC = 0.456;

            floatVector unitDirectionVector = { 0.2, 0.3, 0.4 };

            floatVector previousUnitDirectionVector = { -0.1, 0.1, 0.30 };

            floatVector dDdGradC = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousdDdGradC = { .1, .2, .3, .4, .5, .6, .7, .8, .9 };

        protected:

            virtual void setMassChangeVelocityGradientTrace( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_previousMassChangeVelocityGradientTrace( previousMassChangeVelocityGradientTrace );

                }
                else{

                    set_massChangeVelocityGradientTrace( massChangeVelocityGradientTrace );

                }

            }

            virtual void setUnitDirectionVector( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_previousUnitDirectionVector( previousUnitDirectionVector );

                }
                else{

                    set_unitDirectionVector( unitDirectionVector );

                }

            }

            virtual void setUnitDirectionVectorDerivatives( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_dPreviousUnitDirectionVectordPreviousDirectionVector( previousdDdGradC );

                }
                else{

                    set_dUnitDirectionVectordDirectionVector( dDdGradC );

                }

            }

            virtual void setMassChangeVelocityGradientTraceDerivatives( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_dPreviousMassChangeVelocityGradientTracedPreviousDensity( previousdtrLdRho );

                    set_dPreviousMassChangeVelocityGradientTracedPreviousMassChangeRate( previousdtrLdC );

                }
                else{

                    set_dMassChangeVelocityGradientTracedDensity( dtrLdRho );

                    set_dMassChangeVelocityGradientTracedMassChangeRate( dtrLdC );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer = { 0.06977778, 0.013     , 0.01733333,
                           0.013     , 0.08061111, 0.026     ,
                           0.01733333, 0.026     , 0.09577778 };

    floatVector previousAnswer = { 0.9, 0.0, 0.0,
                                   0.0, 0.9, 0.0,
                                   0.0, 0.0, 0.9 };

    residualMock R( &hydra, 9, 1, hydra.massChangeParameters );

    residualMock Rgrad( &hydra, 9, 1, hydra.massChangeParameters );

    Rgrad.get_dMassChangeVelocityGradientdDensity( );

    Rgrad.get_dPreviousMassChangeVelocityGradientdPreviousDensity( );

    BOOST_TEST( answer == *R.get_massChangeVelocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *R.get_previousMassChangeVelocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *Rgrad.get_massChangeVelocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *Rgrad.get_previousMassChangeVelocityGradient( ), CHECK_PER_ELEMENT );

    floatVector previousdtrLdRho = { R.previousdtrLdRho, 0, 0, 0, R.previousdtrLdRho, 0, 0, 0, R.previousdtrLdRho };

    floatVector previousdtrLdC = { R.previousdtrLdC, 0, 0, 0, R.previousdtrLdC, 0, 0, 0, R.previousdtrLdC };

    BOOST_TEST( previousdtrLdRho == *Rgrad.get_dPreviousMassChangeVelocityGradientdPreviousDensity( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousdtrLdC == *Rgrad.get_dPreviousMassChangeVelocityGradientdPreviousMassChangeRate( ), CHECK_PER_ELEMENT );

    BOOST_TEST( floatVector( 27, 0 ) == *Rgrad.get_dPreviousMassChangeVelocityGradientdPreviousDirectionVector( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeVelocityGradient_4, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.massChangeParameters );

    floatType eps = 1e-6;

    floatVector dLdRho( 9, 0 );

    floatVector dLdC( 9, 0 );

    floatVector dLdGradC( 27, 0 );

    unsigned int offset = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_massChangeVelocityGradient( );

        floatVector vm = *Rm.get_massChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdRho[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdRho == *R.get_dMassChangeVelocityGradientdDensity( ), CHECK_PER_ELEMENT );

    offset = 1;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_massChangeVelocityGradient( );

        floatVector vm = *Rm.get_massChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdC[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdC == *R.get_dMassChangeVelocityGradientdMassChangeRate( ), CHECK_PER_ELEMENT );

    offset = 2;

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_massChangeVelocityGradient( );

        floatVector vm = *Rm.get_massChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdGradC[ 3 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    adaptive_tolerance_test( *R.get_dMassChangeVelocityGradientdDirectionVector( ), dLdGradC );

    offset = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousMassChangeVelocityGradient( );

        floatVector vm = *Rm.get_previousMassChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdRho[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdRho == *R.get_dPreviousMassChangeVelocityGradientdPreviousDensity( ), CHECK_PER_ELEMENT );

    offset = 1;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousMassChangeVelocityGradient( );

        floatVector vm = *Rm.get_previousMassChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdC[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdC == *R.get_dPreviousMassChangeVelocityGradientdPreviousMassChangeRate( ), CHECK_PER_ELEMENT );

    offset = 2;

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousMassChangeVelocityGradient( );

        floatVector vm = *Rm.get_previousMassChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdGradC[ 3 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdGradC == *R.get_dPreviousMassChangeVelocityGradientdPreviousDirectionVector( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeVelocityGradient_5, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

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

    floatVector additionalDOF = { 0.11, 0.22, 0., 0., 0. };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.massChangeParameters );

    floatType eps = 1e-6;

    floatVector dLdRho( 9, 0 );

    floatVector dLdC( 9, 0 );

    floatVector dLdGradC( 27, 0 );

    unsigned int offset = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_massChangeVelocityGradient( );

        floatVector vm = *Rm.get_massChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdRho[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdRho == *R.get_dMassChangeVelocityGradientdDensity( ), CHECK_PER_ELEMENT );

    offset = 1;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_massChangeVelocityGradient( );

        floatVector vm = *Rm.get_massChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdC[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdC == *R.get_dMassChangeVelocityGradientdMassChangeRate( ), CHECK_PER_ELEMENT );

    offset = 2;

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_massChangeVelocityGradient( );

        floatVector vm = *Rm.get_massChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdGradC[ 3 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdGradC == *R.get_dMassChangeVelocityGradientdDirectionVector( ), CHECK_PER_ELEMENT );

    offset = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousMassChangeVelocityGradient( );

        floatVector vm = *Rm.get_previousMassChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdRho[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdRho == *R.get_dPreviousMassChangeVelocityGradientdPreviousDensity( ), CHECK_PER_ELEMENT );

    offset = 1;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousMassChangeVelocityGradient( );

        floatVector vm = *Rm.get_previousMassChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdC[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdC == *R.get_dPreviousMassChangeVelocityGradientdPreviousMassChangeRate( ), CHECK_PER_ELEMENT );

    offset = 2;

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousMassChangeVelocityGradient( );

        floatVector vm = *Rm.get_previousMassChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdGradC[ 3 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdGradC == *R.get_dPreviousMassChangeVelocityGradientdPreviousDirectionVector( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeVelocityGradient_6, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0., 0., 0. };

    floatVector previousStateVariables = { 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.massChangeParameters );

    floatType eps = 1e-6;

    floatVector dLdRho( 9, 0 );

    floatVector dLdC( 9, 0 );

    floatVector dLdGradC( 27, 0 );

    unsigned int offset = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_massChangeVelocityGradient( );

        floatVector vm = *Rm.get_massChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdRho[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdRho == *R.get_dMassChangeVelocityGradientdDensity( ), CHECK_PER_ELEMENT );

    offset = 1;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_massChangeVelocityGradient( );

        floatVector vm = *Rm.get_massChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdC[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdC == *R.get_dMassChangeVelocityGradientdMassChangeRate( ), CHECK_PER_ELEMENT );

    offset = 2;

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_massChangeVelocityGradient( );

        floatVector vm = *Rm.get_massChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdGradC[ 3 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    adaptive_tolerance_test( *R.get_dMassChangeVelocityGradientdDirectionVector( ), dLdGradC );

    offset = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousMassChangeVelocityGradient( );

        floatVector vm = *Rm.get_previousMassChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdRho[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdRho == *R.get_dPreviousMassChangeVelocityGradientdPreviousDensity( ), CHECK_PER_ELEMENT );

    offset = 1;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousMassChangeVelocityGradient( );

        floatVector vm = *Rm.get_previousMassChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdC[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdC == *R.get_dPreviousMassChangeVelocityGradientdPreviousMassChangeRate( ), CHECK_PER_ELEMENT );

    offset = 2;

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 1, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousMassChangeVelocityGradient( );

        floatVector vm = *Rm.get_previousMassChangeVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dLdGradC[ 3 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dLdGradC == *R.get_dPreviousMassChangeVelocityGradientdPreviousDirectionVector( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangePrecedingDeformationGradient_1, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.2, 0.21, 0.29, 0.23, 0.9, 0.11, 0.30, 0.25, 1.1,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatVector answer = { 0.99004464, -0.04875021, -0.25613675, -0.21008905,  1.17745363,
       -0.06235825, -0.20292692, -0.2263232 ,  0.98522215 };
       
    floatVector previousAnswer = { 1.02155172, -0.06034483, -0.14224138, -0.14655172,  0.81034483,
       -0.23275862, -0.31465517, -0.31896552,  0.67672414 };

    residualMock R( &hydra, 9, 1, hydra.massChangeParameters );

    residualMock Rgrad( &hydra, 9, 1, hydra.massChangeParameters );

    Rgrad.get_dPrecedingDeformationGradientdDeformationGradient( );

    Rgrad.get_dPreviousPrecedingDeformationGradientdPreviousDeformationGradient( );

    BOOST_TEST( answer == *R.get_precedingDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *R.get_previousPrecedingDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *Rgrad.get_precedingDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *Rgrad.get_previousPrecedingDeformationGradient( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangePrecedingDeformationGradient_2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

        protected:

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase plasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                plasticity = tardigradeHydra::residualBase( this, 9 );

                massChange = residualMock( this, 9, 2, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                residuals[ 2 ] = &massChange;

                residuals[ 3 ] = &remainder;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                           0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                           -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.1, 0.12, 0.13, 0.14, 1.3, 0.15, 0.16, 0.17, 1.4,
                                  1.2, 0.21, 0.29, 0.23, 0.9, 0.11, 0.30, 0.25, 1.1,
                                  3, 4, 5 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 9, 2, hydra.massChangeParameters );

    floatType eps = 1e-6;

    floatVector dpFdF( 81, 0 );

    floatVector dpFdFn( 81 * 2, 0 );

    floatVector previousdpFdF( 81, 0 );

    floatVector previousdpFdFn( 81 * 2, 0 );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatVector Fp = deformationGradient;

        floatVector Fm = deformationGradient;

        Fp[ i ] += delta;

        Fm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, Fp, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, Fm, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters );

        floatVector vp = *Rp.get_precedingDeformationGradient( );

        floatVector vm = *Rm.get_precedingDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dpFdF[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dpFdF == *R.get_dPrecedingDeformationGradientdDeformationGradient( ), CHECK_PER_ELEMENT );

    unsigned int offset = 9;

    for ( unsigned int i = 0; i < 18; i++ ){

        floatType delta = eps * std::fabs( unknownVector[ i + offset ] ) + eps;

        floatVector xp = unknownVector;

        floatVector xm = unknownVector;

        xp[ i + offset ] += delta;

        xm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient
                , previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, xp );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, xm );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters );

        floatVector vp = *Rp.get_precedingDeformationGradient( );

        floatVector vm = *Rm.get_precedingDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dpFdFn[ 18 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    adaptive_tolerance_test( *R.get_dPrecedingDeformationGradientdSubDeformationGradients( ), dpFdFn );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatVector Fp = previousDeformationGradient;

        floatVector Fm = previousDeformationGradient;

        Fp[ i ] += delta;

        Fm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, Fp,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, Fm,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousPrecedingDeformationGradient( );

        floatVector vm = *Rm.get_previousPrecedingDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdpFdF[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( previousdpFdF == *R.get_dPreviousPrecedingDeformationGradientdPreviousDeformationGradient( ), CHECK_PER_ELEMENT );

    offset = 0;

    for ( unsigned int i = 0; i < 18; i++ ){

        floatType delta = eps * std::fabs( previousStateVariables[ i + offset ] ) + eps;

        floatVector xp = previousStateVariables;

        floatVector xm = previousStateVariables;

        xp[ i + offset ] += delta;

        xm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              xp, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient
                , previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              xm, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousPrecedingDeformationGradient( );

        floatVector vm = *Rm.get_previousPrecedingDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdpFdFn[ 18 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    adaptive_tolerance_test( *R.get_dPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients( ), previousdpFdFn );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeIntermediateVelocityGradient_1, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

            floatVector massChangeVelocityGradient = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousMassChangeVelocityGradient = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            floatVector dMassChangeVelocityGradientdDensity                                = initializeVector(  9 );

            floatVector dMassChangeVelocityGradientdMassChangeRate                         = initializeVector(  9 );

            floatVector dMassChangeVelocityGradientdDirectionVector                 = initializeVector( 27 );

            floatVector dPreviousMassChangeVelocityGradientdPreviousDensity                = initializeVector(  9 );

            floatVector dPreviousMassChangeVelocityGradientdPreviousMassChangeRate         = initializeVector(  9 );

            floatVector dPreviousMassChangeVelocityGradientdPreviousDirectionVector = initializeVector( 27 );

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

        protected:

            virtual void setMassChangeVelocityGradient( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_previousMassChangeVelocityGradient( previousMassChangeVelocityGradient );

                }
                else{

                    set_massChangeVelocityGradient( massChangeVelocityGradient );

                }

            }

            virtual void setMassChangeVelocityGradientDerivatives( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_dPreviousMassChangeVelocityGradientdPreviousDensity( dPreviousMassChangeVelocityGradientdPreviousDensity );

                    set_dPreviousMassChangeVelocityGradientdPreviousMassChangeRate( dPreviousMassChangeVelocityGradientdPreviousMassChangeRate );

                    set_dPreviousMassChangeVelocityGradientdPreviousDirectionVector( dPreviousMassChangeVelocityGradientdPreviousDirectionVector );

                }
                else{

                    set_dMassChangeVelocityGradientdDensity( dMassChangeVelocityGradientdDensity );

                    set_dMassChangeVelocityGradientdMassChangeRate( dMassChangeVelocityGradientdMassChangeRate );

                    set_dMassChangeVelocityGradientdDirectionVector( dMassChangeVelocityGradientdDirectionVector );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.2, 0.21, 0.29, 0.23, 0.9, 0.11, 0.30, 0.25, 1.1,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatVector answer = { 1.12104841, 4.2545126 , 5.17734672,
                           1.85599459, 4.92506006, 5.28347558,
                           4.13186343, 9.15476705, 8.95389153 };

    floatVector previousAnswer = { 0.1, 0.2, 0.3,
                                   0.4, 0.5, 0.6,
                                   0.7, 0.8, 0.9 };

    residualMock R( &hydra, 9, 1, hydra.massChangeParameters );

    residualMock Rgrad( &hydra, 9, 1, hydra.massChangeParameters );

    Rgrad.get_dMassChangeIntermediateVelocityGradientdDensity( );

    Rgrad.get_dPreviousMassChangeIntermediateVelocityGradientdPreviousDensity( );

    BOOST_TEST( answer == *R.get_massChangeIntermediateVelocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *R.get_previousMassChangeIntermediateVelocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *Rgrad.get_massChangeIntermediateVelocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *Rgrad.get_previousMassChangeIntermediateVelocityGradient( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeIntermediateVelocityGradient_2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

        protected:

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase plasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                plasticity = tardigradeHydra::residualBase ( this, 9 );

                massChange = residualMock( this, 9, 2, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                residuals[ 2 ] = &massChange;

                residuals[ 3 ] = &remainder;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                           0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                           -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.1, 0.12, 0.13, 0.14, 1.3, 0.15, 0.16, 0.17, 1.4,
                                  1.2, 0.21, 0.29, 0.23, 0.9, 0.11, 0.30, 0.25, 1.1,
                                  3, 4, 5 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 9, 2, hydra.massChangeParameters );

    floatType eps = 1e-6;

    floatVector dILdF( 81, 0 );

    floatVector dILdFn( 81 * 2, 0 );

    floatVector dILdRho( 9, 0 );

    floatVector dILdC( 9, 0 );

    floatVector dILdGradC( 27, 0 );

    floatVector previousdILdF( 81, 0 );

    floatVector previousdILdFn( 81 * 2, 0 );

    floatVector previousdILdRho( 9, 0 );

    floatVector previousdILdC( 9, 0 );

    floatVector previousdILdGradC( 27, 0 );

    unsigned int offset = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters );

        floatVector vp = *Rp.get_massChangeIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_massChangeIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dILdRho[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dILdRho == *R.get_dMassChangeIntermediateVelocityGradientdDensity( ), CHECK_PER_ELEMENT );

    offset = 1;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters );

        floatVector vp = *Rp.get_massChangeIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_massChangeIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dILdC[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dILdC == *R.get_dMassChangeIntermediateVelocityGradientdMassChangeRate( ), CHECK_PER_ELEMENT );

    offset = 2;

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters );

        floatVector vp = *Rp.get_massChangeIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_massChangeIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dILdGradC[ 3 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dILdGradC == *R.get_dMassChangeIntermediateVelocityGradientdDirectionVector( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatVector Fp = deformationGradient;

        floatVector Fm = deformationGradient;

        Fp[ i ] += delta;

        Fm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, Fp, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, Fm, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters );

        floatVector vp = *Rp.get_massChangeIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_massChangeIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dILdF[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dILdF == *R.get_dMassChangeIntermediateVelocityGradientdDeformationGradient( ), CHECK_PER_ELEMENT );

    offset = 9;

    for ( unsigned int i = 0; i < 18; i++ ){

        floatType delta = eps * std::fabs( unknownVector[ i + offset ] ) + eps;

        floatVector xp = unknownVector;

        floatVector xm = unknownVector;

        xp[ i + offset ] += delta;

        xm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, xp );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, xm );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters );

        floatVector vp = *Rp.get_massChangeIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_massChangeIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dILdFn[ 18 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    adaptive_tolerance_test( *R.get_dMassChangeIntermediateVelocityGradientdSubDeformationGradients( ), dILdFn );

    offset = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousMassChangeIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_previousMassChangeIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdILdRho[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( previousdILdRho == *R.get_dPreviousMassChangeIntermediateVelocityGradientdPreviousDensity( ), CHECK_PER_ELEMENT );

    offset = 1;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousMassChangeIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_previousMassChangeIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdILdC[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( previousdILdC == *R.get_dPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeRate( ), CHECK_PER_ELEMENT );

    offset = 2;

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousMassChangeIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_previousMassChangeIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdILdGradC[ 3 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( previousdILdGradC == *R.get_dPreviousMassChangeIntermediateVelocityGradientdPreviousDirectionVector( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatVector Fp = previousDeformationGradient;

        floatVector Fm = previousDeformationGradient;

        Fp[ i ] += delta;

        Fm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, Fp,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, Fm,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousMassChangeIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_previousMassChangeIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdILdF[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( previousdILdF == *R.get_dPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient( ), CHECK_PER_ELEMENT );

    offset = 0;

    for ( unsigned int i = 0; i < 18; i++ ){

        floatType delta = eps * std::fabs( previousStateVariables[ i + offset ] ) + eps;

        floatVector xp = previousStateVariables;

        floatVector xm = previousStateVariables;

        xp[ i + offset ] += delta;

        xm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              xp, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              xm, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters );

        floatVector vp = *Rp.get_previousMassChangeIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_previousMassChangeIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdILdFn[ 18 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    adaptive_tolerance_test( *R.get_dPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients( ), previousdILdFn );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeDeformationGradient_1, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

            floatVector massChangeVelocityGradient = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousMassChangeVelocityGradient = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            floatVector dMassChangeIntermediateVelocityGradientdDensity                                 = initializeVector(  9 );

            floatVector dMassChangeIntermediateVelocityGradientdMassChangeRate                          = initializeVector(  9 );

            floatVector dMassChangeIntermediateVelocityGradientdDirectionVector                  = initializeVector( 27 );

            floatVector dMassChangeIntermediateVelocityGradientdDeformationGradient                     = initializeVector( 81 );

            floatVector dMassChangeIntermediateVelocityGradientdSubDeformationGradients                 = initializeVector( 81 );

            floatVector dPreviousMassChangeIntermediateVelocityGradientdPreviousDensity                 = initializeVector(  9 );

            floatVector dPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeRate          = initializeVector(  9 );

            floatVector dPreviousMassChangeIntermediateVelocityGradientdPreviousDirectionVector  = initializeVector( 27 );

            floatVector dPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient     = initializeVector( 81 );

            floatVector dPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients = initializeVector( 81 );

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

        protected:

            virtual void setMassChangeIntermediateVelocityGradient( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_previousMassChangeIntermediateVelocityGradient( previousMassChangeVelocityGradient );

                }
                else{

                    set_massChangeIntermediateVelocityGradient( massChangeVelocityGradient );

                }

            }

            virtual void setMassChangeIntermediateVelocityGradientDerivatives( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_dPreviousMassChangeIntermediateVelocityGradientdPreviousDensity( dPreviousMassChangeIntermediateVelocityGradientdPreviousDensity );

                    set_dPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeRate( dPreviousMassChangeIntermediateVelocityGradientdPreviousMassChangeRate );

                    set_dPreviousMassChangeIntermediateVelocityGradientdPreviousDirectionVector( dPreviousMassChangeIntermediateVelocityGradientdPreviousDirectionVector );

                    set_dPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient( dPreviousMassChangeIntermediateVelocityGradientdPreviousDeformationGradient );

                    set_dPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients( dPreviousMassChangeIntermediateVelocityGradientdPreviousSubDeformationGradients );

                }
                else{

                    set_dMassChangeIntermediateVelocityGradientdDensity( dMassChangeIntermediateVelocityGradientdDensity );

                    set_dMassChangeIntermediateVelocityGradientdMassChangeRate( dMassChangeIntermediateVelocityGradientdMassChangeRate );

                    set_dMassChangeIntermediateVelocityGradientdDirectionVector( dMassChangeIntermediateVelocityGradientdDirectionVector );

                    set_dMassChangeIntermediateVelocityGradientdDeformationGradient( dMassChangeIntermediateVelocityGradientdDeformationGradient );

                    set_dMassChangeIntermediateVelocityGradientdSubDeformationGradients( dMassChangeIntermediateVelocityGradientdSubDeformationGradients );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.2, 0.21, 0.29, 0.23, 0.9, 0.11, 0.30, 0.25, 1.1,
                                  3, 4, 5 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatVector answer = { 1.95711900e+10, 2.40473864e+10, 2.85235828e+10, 4.43210244e+10,
       5.44578434e+10, 6.45946624e+10, 6.90708588e+10, 8.48683004e+10,
       1.00665742e+11 };

    residualMock R( &hydra, 9, 1, hydra.massChangeParameters, 0.67 );

    residualMock Rgrad( &hydra, 9, 1, hydra.massChangeParameters, 0.67 );

    Rgrad.get_dMassChangeDeformationGradientdDensity( );

    Rgrad.get_dMassChangeDeformationGradientdPreviousDensity( );

    BOOST_TEST( answer == *R.get_massChangeDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *Rgrad.get_massChangeDeformationGradient( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeDeformationGradient_2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase plasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                plasticity = tardigradeHydra::residualBase ( this, 9 );

                massChange = residualMock( this, 9, 2, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                residuals[ 2 ] = &massChange;

                residuals[ 3 ] = &remainder;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                           0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                           -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.1, 0.12, 0.13, 0.14, 1.3, 0.15, 0.16, 0.17, 1.4,
                                  1.2, 0.21, 0.29, 0.23, 0.9, 0.11, 0.30, 0.25, 1.1,
                                  3, 4, 5 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatType alpha = 0.67;

    residualMock R( &hydra, 9, 2, hydra.massChangeParameters, alpha );

    floatVector dFmdRho( 9, 0 );

    floatVector dFmdC( 9, 0 );

    floatVector dFmdGradC( 27, 0 );

    floatVector dFmdF( 81, 0 );

    floatVector dFmdFn( 81 * 2, 0 );

    floatVector previousdFmdRho( 9, 0 );

    floatVector previousdFmdC( 9, 0 );

    floatVector previousdFmdGradC( 27, 0 );

    floatVector previousdFmdF( 81, 0 );

    floatVector previousdFmdFn( 81 * 2, 0 );

    floatType eps = 1e-6;

    unsigned int offset = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters, alpha );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters, alpha );

        floatVector vp = *Rp.get_massChangeDeformationGradient( );

        floatVector vm = *Rm.get_massChangeDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dFmdRho[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dFmdRho == *R.get_dMassChangeDeformationGradientdDensity( ), CHECK_PER_ELEMENT );

    offset = 1;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters, alpha );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters, alpha );

        floatVector vp = *Rp.get_massChangeDeformationGradient( );

        floatVector vm = *Rm.get_massChangeDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dFmdC[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dFmdC == *R.get_dMassChangeDeformationGradientdMassChangeRate( ), CHECK_PER_ELEMENT );

    offset = 2;

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = additionalDOF;

        floatVector aDm = additionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              aDm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters, alpha );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters, alpha );

        floatVector vp = *Rp.get_massChangeDeformationGradient( );

        floatVector vm = *Rm.get_massChangeDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dFmdGradC[ 3 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dFmdGradC == *R.get_dMassChangeDeformationGradientdDirectionVector( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatVector Fp = deformationGradient;

        floatVector Fm = deformationGradient;

        Fp[ i ] += delta;

        Fm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, Fp, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, Fm, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters, alpha );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters, alpha );

        floatVector vp = *Rp.get_massChangeDeformationGradient( );

        floatVector vm = *Rm.get_massChangeDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dFmdF[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dFmdF == *R.get_dMassChangeDeformationGradientdDeformationGradient( ), CHECK_PER_ELEMENT );

    offset = 9;

    for ( unsigned int i = 0; i < 18; i++ ){

        floatType delta = eps * std::fabs( unknownVector[ i + offset ] ) + eps;

        floatVector xp = unknownVector;

        floatVector xm = unknownVector;

        xp[ i + offset ] += delta;

        xm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, xp );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, xm );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters, alpha );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters, alpha );

        floatVector vp = *Rp.get_massChangeDeformationGradient( );

        floatVector vm = *Rm.get_massChangeDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dFmdFn[ 18 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    adaptive_tolerance_test( *R.get_dMassChangeDeformationGradientdSubDeformationGradients( ), dFmdFn, 1e-8 );

    offset = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters, alpha );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters, alpha );

        floatVector vp = *Rp.get_massChangeDeformationGradient( );

        floatVector vm = *Rm.get_massChangeDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdFmdRho[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( previousdFmdRho == *R.get_dMassChangeDeformationGradientdPreviousDensity( ), CHECK_PER_ELEMENT );

    offset = 1;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters, alpha );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters, alpha );

        floatVector vp = *Rp.get_massChangeDeformationGradient( );

        floatVector vm = *Rm.get_massChangeDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdFmdC[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( previousdFmdC == *R.get_dMassChangeDeformationGradientdPreviousMassChangeRate( ), CHECK_PER_ELEMENT );

    offset = 2;

    for ( unsigned int i = 0; i < 3; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector aDp = previousAdditionalDOF;

        floatVector aDm = previousAdditionalDOF;

        aDp[ i + offset ] += delta;

        aDm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, aDm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters, alpha );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters, alpha );

        floatVector vp = *Rp.get_massChangeDeformationGradient( );

        floatVector vm = *Rm.get_massChangeDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdFmdGradC[ 3 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( previousdFmdGradC == *R.get_dMassChangeDeformationGradientdPreviousDirectionVector( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatVector Fp = previousDeformationGradient;

        floatVector Fm = previousDeformationGradient;

        Fp[ i ] += delta;

        Fm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, Fp,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, Fm,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters, alpha );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters, alpha );

        floatVector vp = *Rp.get_massChangeDeformationGradient( );

        floatVector vm = *Rm.get_massChangeDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdFmdF[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( previousdFmdF == *R.get_dMassChangeDeformationGradientdPreviousDeformationGradient( ), CHECK_PER_ELEMENT );

    offset = 0;

    for ( unsigned int i = 0; i < 18; i++ ){

        floatType delta = eps * std::fabs( previousStateVariables[ i + offset ] ) + eps;

        floatVector xp = previousStateVariables;

        floatVector xm = previousStateVariables;

        xp[ i + offset ] += delta;

        xm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              xp, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              xm, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters, alpha );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters, alpha );

        floatVector vp = *Rp.get_massChangeDeformationGradient( );

        floatVector vm = *Rm.get_massChangeDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdFmdFn[ 18 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    adaptive_tolerance_test( *R.get_dMassChangeDeformationGradientdPreviousSubDeformationGradients( ), previousdFmdFn, 1e-7 );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeResidual, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

            floatVector massChangeDeformationGradient = { 0.111, 0.222, 0.333, 0.444, 0.555, 0.666, 0.777, 0.888, 0.999 };

        protected:

            virtual void setMassChangeDeformationGradient( ) override{

                set_massChangeDeformationGradient( massChangeDeformationGradient );

            }

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase plasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                plasticity = tardigradeHydra::residualBase ( this, 9 );

                massChange = residualMock( this, 9, 2, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                residuals[ 2 ] = &massChange;

                residuals[ 3 ] = &remainder;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                           0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                           -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.1, 0.12, 0.13, 0.14, 1.3, 0.15, 0.16, 0.17, 1.4,
                                  1.2, 0.21, 0.29, 0.23, 0.9, 0.11, 0.30, 0.25, 1.1,
                                  3, 4, 5 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatType alpha = 0.67;

    residualMock R( &hydra, 9, 2, hydra.massChangeParameters, alpha );

    floatVector answer = { -1.089,  0.012,  0.043,
                            0.214, -0.345,  0.556,
                            0.477,  0.638, -0.101 };

    BOOST_TEST( answer == *R.getResidual( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeResidual_2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

        protected:

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase plasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                plasticity = tardigradeHydra::residualBase ( this, 9 );

                massChange = residualMock( this, 9, 2, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                residuals[ 2 ] = &massChange;

                residuals[ 3 ] = &remainder;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                           0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                           -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.1, 0.12, 0.13, 0.14, 1.3, 0.15, 0.16, 0.17, 1.4,
                                  1.2, 0.21, 0.29, 0.23, 0.9, 0.11, 0.30, 0.25, 1.1,
                                  3, 4, 5 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatType alpha = 0.67;

    residualMock R( &hydra, 9, 2, hydra.massChangeParameters, alpha );

    floatVector jacobian( 9 * unknownVector.size( ), 0 );

    floatVector dRdF( 9 * 9, 0 );

    floatVector dRdT( 9, 0 );

    floatVector dRdAdditionalDOF( 9 * 5, 0 );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatType delta = eps * std::fabs( unknownVector[ i ] ) + eps;

        floatVector xp = unknownVector;

        floatVector xm = unknownVector;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, xp );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, xm );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters, alpha );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters, alpha );

        floatVector vp = *Rp.getResidual( );

        floatVector vm = *Rm.getResidual( );

        for ( unsigned int j = 0; j < 9; j++ ){

            jacobian[ unknownVector.size( ) * j + i ] += ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    adaptive_tolerance_test( *R.getJacobian( ), jacobian, 1e-8 );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatVector xp = deformationGradient;

        floatVector xm = deformationGradient;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, xp, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, xm, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters, alpha );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters, alpha );

        floatVector vp = *Rp.getResidual( );

        floatVector vm = *Rm.getResidual( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dRdF[ 9 * j + i ] += ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dRdF == *R.getdRdF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( temperature ) + eps;

        floatType xp = temperature;

        floatType xm = temperature;

        xp += delta;

        xm -= delta;

        hydraBaseMock hydrap( time, deltaTime, xp, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, xm, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters, alpha );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters, alpha );

        floatVector vp = *Rp.getResidual( );

        floatVector vm = *Rm.getResidual( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dRdT[ 1 * j + i ] += ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dRdT == *R.getdRdT( ), CHECK_PER_ELEMENT );

    eps = 1e-7;

    for ( unsigned int i = 0; i < additionalDOF.size( ); i++ ){

        floatType delta = eps * std::fabs( temperature ) + eps;

        floatVector xp = additionalDOF;

        floatVector xm = additionalDOF;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              xp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              xm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, hydra.massChangeParameters, alpha );

        residualMock Rm( &hydram, 9, 2, hydra.massChangeParameters, alpha );

        floatVector vp = *Rp.getResidual( );

        floatVector vm = *Rm.getResidual( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dRdAdditionalDOF[ 5 * j + i ] += ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dRdAdditionalDOF == *R.getdRdAdditionalDOF( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_exampleModel, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.00, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector dStressdF( 9 * 9, 0 );

    floatVector dStressdT(     9, 0 );

    floatVector dStressdAdditionalDOF( 5 * 9, 0 );

    hydra.evaluate( );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatVector xp = deformationGradient;

        floatVector xm = deformationGradient;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, xp, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, xm, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap.evaluate( );

        hydram.evaluate( );

        floatVector vp = *hydrap.getStress( );

        floatVector vm = *hydram.getStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dStressdF[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    floatVector RdStressdF( hydra.getFlatdXdF( )->begin( ), hydra.getFlatdXdF( )->begin( ) + 81 );

    BOOST_TEST( dStressdF.size( ) == RdStressdF.size( ) );
    for ( unsigned int i = 0; i < dStressdF.size( ); i++ ){
        BOOST_TEST( dStressdF[ i ] == RdStressdF[ i ], boost::test_tools::tolerance( 1e-5 ) );
    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatType xp = temperature;

        floatType xm = temperature;

        xp += delta;

        xm -= delta;

        hydraBaseMock hydrap( time, deltaTime, xp, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, xm, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap.evaluate( );

        hydram.evaluate( );

        floatVector vp = *hydrap.getStress( );

        floatVector vm = *hydram.getStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dStressdT[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    floatVector RdStressdT( hydra.getFlatdXdT( )->begin( ), hydra.getFlatdXdT( )->begin( ) + 9 );

    BOOST_TEST( dStressdT == RdStressdT, CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 5; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i ] ) + eps;

        floatVector xp = additionalDOF;

        floatVector xm = additionalDOF;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              xp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              xm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap.evaluate( );

        hydram.evaluate( );

        floatVector vp = *hydrap.getStress( );

        floatVector vm = *hydram.getStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dStressdAdditionalDOF[ 5 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    floatVector RdStressdAdditionalDOF( hydra.getFlatdXdAdditionalDOF( )->begin( ), hydra.getFlatdXdAdditionalDOF( )->begin( ) + 9 * 5 );

    BOOST_TEST( dStressdAdditionalDOF.size( ) == RdStressdAdditionalDOF.size( ) );
    for ( unsigned int i = 0; i < dStressdAdditionalDOF.size( ); i++ ){
        BOOST_TEST( dStressdAdditionalDOF[ i ] == RdStressdAdditionalDOF[ i ], boost::test_tools::tolerance( 1e-5 ) );
    }

}

BOOST_AUTO_TEST_CASE( test_residual_suggestInitialIterateValues, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

            floatVector massChangeDeformationGradient = { 0.111, 0.222, 0.333, 0.444, 0.555, 0.666, 0.777, 0.888, 0.999 };

        protected:

            virtual void setMassChangeDeformationGradient( ) override{

                set_massChangeDeformationGradient( massChangeDeformationGradient );

            }

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase plasticity;

            residualMock massChange;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                plasticity = tardigradeHydra::residualBase ( this, 9 );

                massChange = residualMock( this, 9, 2, massChangeParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                residuals[ 2 ] = &massChange;

                residuals[ 3 ] = &remainder;

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

    floatVector additionalDOF = { 0.11, 0.22, 0.33, 0.44, 0.55 };

    floatVector previousAdditionalDOF = { 0.12, 0.23, 0.36, 0.47, 0.58 };

    floatVector previousStateVariables = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                           0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,
                                           -1, -2, -3 };

    floatVector parameters = { 123.4, 56.7, 0.78 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.1, 0.12, 0.13, 0.14, 1.3, 0.15, 0.16, 0.17, 1.4,
                                  1.2, 0.21, 0.29, 0.23, 0.9, 0.11, 0.30, 0.25, 1.1,
                                  3, 4, 5 };

    const unsigned int numConfigurations = 3;

    const unsigned int numNonLinearSolveStateVariables = 3;

    const unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 9, 2, hydra.massChangeParameters );

    std::vector< unsigned int > indices = { 18, 19, 20, 21, 22, 23, 24, 25, 26 };

    floatVector values = { 0.111, 0.222, 0.333, 0.444, 0.555, 0.666, 0.777, 0.888, 0.999 };

    std::vector< unsigned int > result_1;
    std::vector< floatType > result_2;

    R.suggestInitialIterateValues( result_1, result_2 );

    BOOST_TEST( indices == result_1, CHECK_PER_ELEMENT );

    BOOST_TEST( values == result_2, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_expectedDeformation, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 1.00 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters, 1.0 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.3862943611198906;

    floatType deltaTime = 1.3862943611198906;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = { 1.1, 0.0, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0 };

    floatVector previousDeformationGradient = { 1.0, 0.0, 0.0,
                                                0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0 };

    floatVector additionalDOF = { 1.00, -0.50, .1, 0, 0 };

    floatVector previousAdditionalDOF = { 1.00, -0.00, 0.00, 0.00, 0.00 };

    floatVector previousStateVariables = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    floatVector parameters = { 123.4, 56.7, 1.00 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer = { 0.5, 0, 0, 0, 1, 0, 0, 0, 1 };

    floatVector dStressdF( 9 * 9, 0 );

    floatVector dStressdT(     9, 0 );

    floatVector dStressdAdditionalDOF( 5 * 9, 0 );

    hydra.evaluate( );

    BOOST_TEST( answer == floatVector( hydra.getUnknownVector( )->begin( ) + 9, hydra.getUnknownVector( )->begin( ) + 18 ), CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatVector xp = deformationGradient;

        floatVector xm = deformationGradient;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, xp, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, xm, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap.evaluate( );

        hydram.evaluate( );

        floatVector vp = *hydrap.getStress( );

        floatVector vm = *hydram.getStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dStressdF[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    floatVector RdStressdF( hydra.getFlatdXdF( )->begin( ), hydra.getFlatdXdF( )->begin( ) + 81 );

    adaptive_tolerance_test( RdStressdF, dStressdF, 2e-7 );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatType xp = temperature;

        floatType xm = temperature;

        xp += delta;

        xm -= delta;

        hydraBaseMock hydrap( time, deltaTime, xp, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, xm, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap.evaluate( );

        hydram.evaluate( );

        floatVector vp = *hydrap.getStress( );

        floatVector vm = *hydram.getStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dStressdT[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    floatVector RdStressdT( hydra.getFlatdXdT( )->begin( ), hydra.getFlatdXdT( )->begin( ) + 9 );

    BOOST_TEST( dStressdT == RdStressdT, CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 5; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i ] ) + eps;

        floatVector xp = additionalDOF;

        floatVector xm = additionalDOF;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              xp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              xm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap.evaluate( );

        hydram.evaluate( );

        floatVector vp = *hydrap.getStress( );

        floatVector vm = *hydram.getStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dStressdAdditionalDOF[ 5 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    floatVector RdStressdAdditionalDOF( hydra.getFlatdXdAdditionalDOF( )->begin( ), hydra.getFlatdXdAdditionalDOF( )->begin( ) + 9 * 5 );

    adaptive_tolerance_test( RdStressdAdditionalDOF, dStressdAdditionalDOF, 1e-7 );

}

BOOST_AUTO_TEST_CASE( test_residual_expectedDeformation2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.00 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters, 1.0 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.3862943611198906;

    floatType deltaTime = 4.158883083359672;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = { 1.1, 0.0, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0 };

    floatVector previousDeformationGradient = { 1.0, 0.0, 0.0,
                                                0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0 };

    floatVector additionalDOF = { 1.00, -0.50, .1, 0, 0 };

    floatVector previousAdditionalDOF = { 1.00, -0.00, 0.00, 0.00, 0.00 };

    floatVector previousStateVariables = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    floatVector parameters = { 123.4, 56.7, 0.00 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer = { 0.5, 0, 0, 0, 0.5, 0, 0, 0, 0.5 };

    floatVector dStressdF( 9 * 9, 0 );

    floatVector dStressdT(     9, 0 );

    floatVector dStressdAdditionalDOF( 5 * 9, 0 );

    hydra.evaluate( );

    BOOST_TEST( answer == floatVector( hydra.getUnknownVector( )->begin( ) + 9, hydra.getUnknownVector( )->begin( ) + 18 ), CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatVector xp = deformationGradient;

        floatVector xm = deformationGradient;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, xp, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, xm, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap.evaluate( );

        hydram.evaluate( );

        floatVector vp = *hydrap.getStress( );

        floatVector vm = *hydram.getStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dStressdF[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    floatVector RdStressdF( hydra.getFlatdXdF( )->begin( ), hydra.getFlatdXdF( )->begin( ) + 81 );

    adaptive_tolerance_test( RdStressdF, dStressdF, 2e-7 );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatType xp = temperature;

        floatType xm = temperature;

        xp += delta;

        xm -= delta;

        hydraBaseMock hydrap( time, deltaTime, xp, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, xm, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap.evaluate( );

        hydram.evaluate( );

        floatVector vp = *hydrap.getStress( );

        floatVector vm = *hydram.getStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dStressdT[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    floatVector RdStressdT( hydra.getFlatdXdT( )->begin( ), hydra.getFlatdXdT( )->begin( ) + 9 );

    BOOST_TEST( dStressdT == RdStressdT, CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 5; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i ] ) + eps;

        floatVector xp = additionalDOF;

        floatVector xm = additionalDOF;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              xp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              xm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap.evaluate( );

        hydram.evaluate( );

        floatVector vp = *hydrap.getStress( );

        floatVector vm = *hydram.getStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dStressdAdditionalDOF[ 5 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    floatVector RdStressdAdditionalDOF( hydra.getFlatdXdAdditionalDOF( )->begin( ), hydra.getFlatdXdAdditionalDOF( )->begin( ) + 9 * 5 );

    adaptive_tolerance_test( RdStressdAdditionalDOF, dStressdAdditionalDOF, 1e-7 );

}

BOOST_AUTO_TEST_CASE( test_residual_expectedDeformation3, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector massChangeParameters = { 0.76 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock massChange;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                massChange = residualMock( this, 9, 1, massChangeParameters, 1.0 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &massChange;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.3862943611198906;

    floatType deltaTime = 4.158883083359672;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = { 1.1, 0.0, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0 };

    floatVector previousDeformationGradient = { 1.0, 0.0, 0.0,
                                                0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0 };

    floatVector additionalDOF = { 1.00, -0.00, 0, 0, 0 };

    floatVector previousAdditionalDOF = { 1.00, -0.00, 0.00, 0.00, 0.00 };

    floatVector previousStateVariables = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    floatVector parameters = { 123.4, 56.7, 0.00 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    floatVector answer = { 1.0, 0, 0, 0, 1.0, 0, 0, 0, 1.0 };

    floatVector dStressdF( 9 * 9, 0 );

    floatVector dStressdT(     9, 0 );

    floatVector dStressdAdditionalDOF( 5 * 9, 0 );

    hydra.evaluate( );

    BOOST_TEST( answer == floatVector( hydra.getUnknownVector( )->begin( ) + 9, hydra.getUnknownVector( )->begin( ) + 18 ), CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatVector xp = deformationGradient;

        floatVector xm = deformationGradient;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, xp, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, xm, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap.evaluate( );

        hydram.evaluate( );

        floatVector vp = *hydrap.getStress( );

        floatVector vm = *hydram.getStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dStressdF[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    floatVector RdStressdF( hydra.getFlatdXdF( )->begin( ), hydra.getFlatdXdF( )->begin( ) + 81 );

    adaptive_tolerance_test( RdStressdF, dStressdF, 2e-7 );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( deformationGradient[ i ] ) + eps;

        floatType xp = temperature;

        floatType xm = temperature;

        xp += delta;

        xm -= delta;

        hydraBaseMock hydrap( time, deltaTime, xp, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, xm, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap.evaluate( );

        hydram.evaluate( );

        floatVector vp = *hydrap.getStress( );

        floatVector vm = *hydram.getStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dStressdT[ 1 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    floatVector RdStressdT( hydra.getFlatdXdT( )->begin( ), hydra.getFlatdXdT( )->begin( ) + 9 );

    BOOST_TEST( dStressdT == RdStressdT, CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 5; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i ] ) + eps;

        floatVector xp = additionalDOF;

        floatVector xm = additionalDOF;

        xp[ i ] += delta;

        xm[ i ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              xp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              xm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap.evaluate( );

        hydram.evaluate( );

        floatVector vp = *hydrap.getStress( );

        floatVector vm = *hydram.getStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dStressdAdditionalDOF[ 5 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    floatVector RdStressdAdditionalDOF( hydra.getFlatdXdAdditionalDOF( )->begin( ), hydra.getFlatdXdAdditionalDOF( )->begin( ) + 9 * 5 );

    adaptive_tolerance_test( RdStressdAdditionalDOF, dStressdAdditionalDOF, 1e-7 );

}
