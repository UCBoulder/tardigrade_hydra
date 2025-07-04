/**
  * \file test_tardigrade_hydraDOFVelocityGradientDeformation.cpp
  *
  * Tests for tardigrade_hydraDOFVelocityGradientDeformation
  */

#include<tardigrade_hydraLinearElasticity.h>
#include<tardigrade_hydraDOFVelocityGradientDeformation.h>
#include<tardigrade_constitutive_tools.h>

#define BOOST_TEST_MODULE test_tardigrade_hydraDOFVelocityGradientDeformation
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::dofVelocityGradientDeformation::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::dofVelocityGradientDeformation::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::dofVelocityGradientDeformation::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type


namespace tardigradeHydra{

    namespace unit_test{

        class hydraBaseTester{

            public:

                static void updateUnknownVector( tardigradeHydra::hydraBase &hydra, const floatVector &value ){

                    BOOST_CHECK_NO_THROW( hydra.updateUnknownVector( value ) );

                }

        };

    }

    namespace dofVelocityGradientDeformation{

        namespace unit_test{

            class residualTester{

                public:

                    static void runBasicGetTests( tardigradeHydra::dofVelocityGradientDeformation::residual &R ){

                        BOOST_CHECK( R._dofConfigurationIndex == R.getDOFConfigurationIndex( ) );

                        BOOST_CHECK( R._dofVelocityGradientIndex == R.getDOFVelocityGradientIndex( ) );

                        BOOST_CHECK( R._integrationParameter == R.getIntegrationParameter( ) );

                        BOOST_CHECK( &R._dofVelocityGradient.second         == R.get_dofVelocityGradient( ) );

                        BOOST_CHECK( &R._previousDOFVelocityGradient.second == R.get_previousDOFVelocityGradient( ) );

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

    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {

        public:

            using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector dofDeformationParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock dofDeformation;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                dofDeformation = residualMock( this, 9, 1, 3, dofDeformationParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &dofDeformation;

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

    floatVector additionalDOF = {
        0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
        1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99
    };

    floatVector previousAdditionalDOF = {
        -0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
        -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99
    };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         additionalDOF, previousAdditionalDOF,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, 3, hydra.dofDeformationParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::dofVelocityGradientDeformation::unit_test::residualTester::runBasicGetTests( R );

    floatVector dofVelocityGradientAnswer = {
        0.44, 0.55, 0.66, 0.77, 0.88, 0.99, 1.11, 1.22, 1.33
    };

    floatVector previousDOFVelocityGradientAnswer = {
        -0.44, -0.55, -0.66, -0.77, -0.88, -0.99, -1.11, -1.22, -1.33
    };

    BOOST_TEST( dofVelocityGradientAnswer         == *R.get_dofVelocityGradient( ),         CHECK_PER_ELEMENT );

    BOOST_TEST( previousDOFVelocityGradientAnswer == *R.get_previousDOFVelocityGradient( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_dofPrecedingDeformationGradient_1, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {

        public:

            using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector dofDeformationParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock dofDeformation;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                dofDeformation = residualMock( this, 9, 1, 3, dofDeformationParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &dofDeformation;

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

    floatVector additionalDOF = {
        0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
        1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99
    };

    floatVector previousAdditionalDOF = {
        -0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
        -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99
    };

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

    residualMock R( &hydra, 9, 1, 3, hydra.dofDeformationParameters );

    residualMock Rgrad( &hydra, 9, 1, 3, hydra.dofDeformationParameters );

    Rgrad.get_dPrecedingDeformationGradientdDeformationGradient( );

    Rgrad.get_dPreviousPrecedingDeformationGradientdPreviousDeformationGradient( );

    BOOST_TEST( answer == *R.get_precedingDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *R.get_previousPrecedingDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *Rgrad.get_precedingDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *Rgrad.get_previousPrecedingDeformationGradient( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_dofPrecedingDeformationGradient_2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {

        public:

            using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

        protected:

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector dofDeformationParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase plasticity;

            residualMock dofDeformation;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                plasticity = tardigradeHydra::residualBase( this, 9 );

                dofDeformation = residualMock( this, 9, 2, 3, dofDeformationParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                residuals[ 2 ] = &dofDeformation;

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

    floatVector additionalDOF = {
        0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
        1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99
    };

    floatVector previousAdditionalDOF = {
        -0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
        -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99
    };

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

    residualMock R( &hydra, 9, 2, 3, hydra.dofDeformationParameters );

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters );

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters );

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters );

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters );

        floatVector vp = *Rp.get_previousPrecedingDeformationGradient( );

        floatVector vm = *Rm.get_previousPrecedingDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdpFdFn[ 18 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    adaptive_tolerance_test( *R.get_dPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients( ), previousdpFdFn );

}

BOOST_AUTO_TEST_CASE( test_residual_dofIntermediateVelocityGradient_1, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {

        public:

            using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        protected:

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector dofDeformationParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock dofDeformation;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                dofDeformation = residualMock( this, 9, 1, 3, dofDeformationParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &dofDeformation;

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

    floatVector additionalDOF = {
        0.11, 0.22, 0.33, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00,
        7.00, 8.00, 9.00, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99
    };

    floatVector previousAdditionalDOF = {
        -0.11, -0.22, -0.33, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60,
        0.70, 0.80, 0.90, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99
    };

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

    residualMock R( &hydra, 9, 1, 3, hydra.dofDeformationParameters );

    residualMock Rgrad( &hydra, 9, 1, 3, hydra.dofDeformationParameters );

    Rgrad.get_dDOFIntermediateVelocityGradientdDOFVelocityGradient( );

    Rgrad.get_dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient( );

    BOOST_TEST( answer == *R.get_dofIntermediateVelocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *R.get_previousDOFIntermediateVelocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *Rgrad.get_dofIntermediateVelocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *Rgrad.get_previousDOFIntermediateVelocityGradient( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_dofIntermediateVelocityGradient_2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {

        public:

            using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

        protected:

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector dofDeformationParameters;

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase plasticity;

            residualMock dofDeformation;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                plasticity = tardigradeHydra::residualBase ( this, 9 );

                dofDeformation = residualMock( this, 9, 2, 3, dofDeformationParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                residuals[ 2 ] = &dofDeformation;

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

    floatVector additionalDOF = {
        0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
        1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99
    };

    floatVector previousAdditionalDOF = {
        -0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
        -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99
    };

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

    residualMock R( &hydra, 9, 2, 3, hydra.dofDeformationParameters );

    floatType eps = 1e-6;

    floatVector dILdF( 81, 0 );

    floatVector dILdFn( 81 * 2, 0 );

    floatVector dILdL( 81, 0 );

    floatVector previousdILdF( 81, 0 );

    floatVector previousdILdFn( 81 * 2, 0 );

    floatVector previousdILdL( 81, 0 );

    unsigned int offset = 0;

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters );

        floatVector vp = *Rp.get_dofIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_dofIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dILdF[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dILdF == *R.get_dDOFIntermediateVelocityGradientdDeformationGradient( ), CHECK_PER_ELEMENT );

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters );

        floatVector vp = *Rp.get_dofIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_dofIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dILdFn[ 18 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    adaptive_tolerance_test( *R.get_dDOFIntermediateVelocityGradientdSubDeformationGradients( ), dILdFn );

    offset = 3;

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector xp = additionalDOF;

        floatVector xm = additionalDOF;

        xp[ i + offset ] += delta;

        xm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              xp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              xm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters );

        floatVector vp = *Rp.get_dofIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_dofIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dILdL[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( *R.get_dDOFIntermediateVelocityGradientdDOFVelocityGradient( ) == dILdL, CHECK_PER_ELEMENT );

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters );

        floatVector vp = *Rp.get_previousDOFIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_previousDOFIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdILdF[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( previousdILdF == *R.get_dPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient( ), CHECK_PER_ELEMENT );

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters );

        floatVector vp = *Rp.get_previousDOFIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_previousDOFIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdILdFn[ 18 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    adaptive_tolerance_test( *R.get_dPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients( ), previousdILdFn );

    offset = 3;

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( previousAdditionalDOF[ i + offset ] ) + eps;

        floatVector xp = previousAdditionalDOF;

        floatVector xm = previousAdditionalDOF;

        xp[ i + offset ] += delta;

        xm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, xp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, xm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters );

        floatVector vp = *Rp.get_previousDOFIntermediateVelocityGradient( );

        floatVector vm = *Rm.get_previousDOFIntermediateVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdILdL[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( *R.get_dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient( ) == previousdILdL, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_dofDeformationGradient_1, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {

        public:

            using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

            floatVector dofIntermediateVelocityGradient = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousDOFIntermediateVelocityGradient = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            floatVector dDOFIntermediateVelocityGradientdDOFVelocityGradient              = initializeVector( 81 );

            floatVector dDOFIntermediateVelocityGradientdDeformationGradient                     = initializeVector( 81 );

            floatVector dDOFIntermediateVelocityGradientdSubDeformationGradients                 = initializeVector( 81 );

            floatVector dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient = initializeVector( 81 );

            floatVector dPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient        = initializeVector( 81 );

            floatVector dPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients    = initializeVector( 81 );

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

        protected:

            using tardigradeHydra::dofVelocityGradientDeformation::residual::setDOFIntermediateVelocityGradient;

            virtual void setDOFIntermediateVelocityGradient( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_previousDOFIntermediateVelocityGradient( previousDOFIntermediateVelocityGradient );

                }
                else{

                    set_dofIntermediateVelocityGradient( dofIntermediateVelocityGradient );

                }

            }

            virtual void setDOFIntermediateVelocityGradientDerivatives( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient( dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient );

                    set_dPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient( dPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient );

                    set_dPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients( dPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients );

                }
                else{

                    set_dDOFIntermediateVelocityGradientdDOFVelocityGradient( dDOFIntermediateVelocityGradientdDOFVelocityGradient );

                    set_dDOFIntermediateVelocityGradientdDeformationGradient( dDOFIntermediateVelocityGradientdDeformationGradient );

                    set_dDOFIntermediateVelocityGradientdSubDeformationGradients( dDOFIntermediateVelocityGradientdSubDeformationGradients );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector dofDeformationParameters;

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock dofDeformation;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                dofDeformation = residualMock( this, 9, 1, 3, dofDeformationParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &dofDeformation;

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

    floatVector additionalDOF = {
        0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
        1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99
    };

    floatVector previousAdditionalDOF = {
        -0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
        -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99
    };

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

    residualMock R( &hydra, 9, 1, 3, hydra.dofDeformationParameters, 0.67 );

    residualMock Rgrad( &hydra, 9, 1, 3, hydra.dofDeformationParameters, 0.67 );

    Rgrad.get_dDOFDeformationGradientdDOFVelocityGradient( );

    Rgrad.get_dDOFDeformationGradientdPreviousDOFVelocityGradient( );

    BOOST_TEST( answer == *R.get_dofDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *Rgrad.get_dofDeformationGradient( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_dofDeformationGradient_2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {

        public:

            using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector dofDeformationParameters;

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase plasticity;

            residualMock dofDeformation;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                plasticity = tardigradeHydra::residualBase ( this, 9 );

                dofDeformation = residualMock( this, 9, 2, 3, dofDeformationParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                residuals[ 2 ] = &dofDeformation;

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

    floatVector additionalDOF = {
        0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
        1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99
    };

    floatVector previousAdditionalDOF = {
        -0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
        -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99
    };

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

    residualMock R( &hydra, 9, 2, 3, hydra.dofDeformationParameters, alpha );

    floatVector dFmdL( 81, 0 );

    floatVector dFmdF( 81, 0 );

    floatVector dFmdFn( 81 * 2, 0 );

    floatVector previousdFmdL( 81, 0 );

    floatVector previousdFmdF( 81, 0 );

    floatVector previousdFmdFn( 81 * 2, 0 );

    floatType eps = 1e-6;

    unsigned int offset = 0;

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        floatVector vp = *Rp.get_dofDeformationGradient( );

        floatVector vm = *Rm.get_dofDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dFmdF[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dFmdF == *R.get_dDOFDeformationGradientdDeformationGradient( ), CHECK_PER_ELEMENT );

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        floatVector vp = *Rp.get_dofDeformationGradient( );

        floatVector vm = *Rm.get_dofDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dFmdFn[ 18 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    adaptive_tolerance_test( *R.get_dDOFDeformationGradientdSubDeformationGradients( ), dFmdFn, 1e-8 );

    offset = 3;

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( additionalDOF[ i + offset ] ) + eps;

        floatVector xp = additionalDOF;

        floatVector xm = additionalDOF;

        xp[ i + offset ] += delta;

        xm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              xp, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              xm, previousAdditionalDOF,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        floatVector vp = *Rp.get_dofDeformationGradient( );

        floatVector vm = *Rm.get_dofDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dFmdL[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( *R.get_dDOFDeformationGradientdDOFVelocityGradient( ) == dFmdL, CHECK_PER_ELEMENT );

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        floatVector vp = *Rp.get_dofDeformationGradient( );

        floatVector vm = *Rm.get_dofDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdFmdF[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( previousdFmdF == *R.get_dDOFDeformationGradientdPreviousDeformationGradient( ), CHECK_PER_ELEMENT );

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        floatVector vp = *Rp.get_dofDeformationGradient( );

        floatVector vm = *Rm.get_dofDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdFmdFn[ 18 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    adaptive_tolerance_test( *R.get_dDOFDeformationGradientdPreviousSubDeformationGradients( ), previousdFmdFn, 1e-7 );

    offset = 3;

    for ( unsigned int i = 0; i < 9; i++ ){

        floatType delta = eps * std::fabs( previousAdditionalDOF[ i + offset ] ) + eps;

        floatVector xp = previousAdditionalDOF;

        floatVector xm = previousAdditionalDOF;

        xp[ i + offset ] += delta;

        xm[ i + offset ] -= delta;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, xp,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              additionalDOF, xm,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        floatVector vp = *Rp.get_dofDeformationGradient( );

        floatVector vm = *Rm.get_dofDeformationGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            previousdFmdL[ 9 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( *R.get_dDOFDeformationGradientdPreviousDOFVelocityGradient( ) == previousdFmdL, CHECK_PER_ELEMENT );
}

BOOST_AUTO_TEST_CASE( test_residual_dofResidual, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {

        public:

            using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

            floatVector dofDeformationGradient = { 0.111, 0.222, 0.333, 0.444, 0.555, 0.666, 0.777, 0.888, 0.999 };

        protected:

            virtual void setDOFDeformationGradient( ) override{

                set_dofDeformationGradient( dofDeformationGradient );

            }

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector dofDeformationParameters;

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase plasticity;

            residualMock dofDeformation;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                plasticity = tardigradeHydra::residualBase ( this, 9 );

                dofDeformation = residualMock( this, 9, 2, 3, dofDeformationParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                residuals[ 2 ] = &dofDeformation;

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

    floatVector additionalDOF = {
        0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
        1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99
    };

    floatVector previousAdditionalDOF = {
        -0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
        -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99
    };

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

    residualMock R( &hydra, 9, 2, 3, hydra.dofDeformationParameters, alpha );

    floatVector answer = { -1.089,  0.012,  0.043,
                            0.214, -0.345,  0.556,
                            0.477,  0.638, -0.101 };

    BOOST_TEST( answer == *R.getResidual( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_dofResidual_2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {

        public:

            using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        protected:

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector dofDeformationParameters;

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase plasticity;

            residualMock dofDeformation;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                plasticity = tardigradeHydra::residualBase ( this, 9 );

                dofDeformation = residualMock( this, 9, 2, 3, dofDeformationParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                residuals[ 2 ] = &dofDeformation;

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

    floatVector additionalDOF = {
        0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
        1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99
    };

    floatVector previousAdditionalDOF = {
        -0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
        -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99
    };

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

    residualMock R( &hydra, 9, 2, 3, hydra.dofDeformationParameters, alpha );

    floatVector jacobian( 9 * unknownVector.size( ), 0 );

    floatVector dRdF( 9 * 9, 0 );

    floatVector dRdT( 9, 0 );

    floatVector dRdAdditionalDOF( 9 * 18, 0 );

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters, alpha );

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters, alpha );

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters, alpha );

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

        residualMock Rp( &hydrap, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        residualMock Rm( &hydram, 9, 2, 3, hydra.dofDeformationParameters, alpha );

        floatVector vp = *Rp.getResidual( );

        floatVector vm = *Rm.getResidual( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dRdAdditionalDOF[ 18 * j + i ] += ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dRdAdditionalDOF == *R.getdRdAdditionalDOF( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_exampleModel, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {

        public:

            using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector dofDeformationParameters;

            tardigradeHydra::linearElasticity::residual elasticity;

            residualMock dofDeformation;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                dofDeformation = residualMock( this, 9, 1, 3, dofDeformationParameters );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &dofDeformation;

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

    floatVector additionalDOF = {
        0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
        1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99
    };

    floatVector previousAdditionalDOF = {
        -0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
        -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99
    };

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

    floatVector dStressdAdditionalDOF( 18 * 9, 0 );

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

    for ( unsigned int i = 0; i < 18; i++ ){

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

            dStressdAdditionalDOF[ 18 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    floatVector RdStressdAdditionalDOF( hydra.getFlatdXdAdditionalDOF( )->begin( ), hydra.getFlatdXdAdditionalDOF( )->begin( ) + 9 * 18 );

    BOOST_TEST( dStressdAdditionalDOF.size( ) == RdStressdAdditionalDOF.size( ) );
    for ( unsigned int i = 0; i < dStressdAdditionalDOF.size( ); i++ ){
        BOOST_TEST( dStressdAdditionalDOF[ i ] == RdStressdAdditionalDOF[ i ], boost::test_tools::tolerance( 1e-5 ) );
    }

}

BOOST_AUTO_TEST_CASE( test_residual_suggestInitialIterateValues, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {

        public:

            using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

            floatVector dofDeformationGradient = { 0.111, 0.222, 0.333, 0.444, 0.555, 0.666, 0.777, 0.888, 0.999 };

        protected:

            virtual void setDOFDeformationGradient( ) override{

                set_dofDeformationGradient( dofDeformationGradient );

            }

            floatVector initializeVector( unsigned int size ){

                return floatVector( size, 0 );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector dofDeformationParameters = { 0.78 };

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase plasticity;

            residualMock dofDeformation;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                plasticity = tardigradeHydra::residualBase ( this, 9 );

                dofDeformation = residualMock( this, 9, 2, 3, dofDeformationParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

                residuals[ 2 ] = &dofDeformation;

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

    floatVector additionalDOF = {
        0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
        1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99
    };

    floatVector previousAdditionalDOF = {
        -0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
        -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99
    };

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

    residualMock R( &hydra, 9, 2, 3, hydra.dofDeformationParameters );

    std::vector< unsigned int > indices = { 18, 19, 20, 21, 22, 23, 24, 25, 26 };

    floatVector values = { 0.111, 0.222, 0.333, 0.444, 0.555, 0.666, 0.777, 0.888, 0.999 };

    std::vector< unsigned int > result_1;
    std::vector< floatType > result_2;

    R.suggestInitialIterateValues( result_1, result_2 );

    BOOST_TEST( indices == result_1, CHECK_PER_ELEMENT );

    BOOST_TEST( values == result_2, CHECK_PER_ELEMENT );

}
