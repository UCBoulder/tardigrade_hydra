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

                        BOOST_CHECK( &R._density.second         == R.get_density( ) );

                        BOOST_CHECK( &R._previousDensity.second == R.get_previousDensity( ) );

                        BOOST_CHECK( &R._massChangeRate.second         == R.get_massChangeRate( ) );

                        BOOST_CHECK( &R._previousMassChangeRate.second == R.get_previousMassChangeRate( ) );

                        BOOST_CHECK( &R._massChangeRateGradient.second         == R.get_massChangeRateGradient( ) );

                        BOOST_CHECK( &R._previousMassChangeRateGradient.second == R.get_previousMassChangeRateGradient( ) );

                    }

            };

        }

    }

}

BOOST_AUTO_TEST_CASE( test_residual_basicGetTests ){

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0.11, *R.get_density( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0.22, *R.get_massChangeRate( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( massChangeRateAnswer, *R.get_massChangeRateGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0.12, *R.get_previousDensity( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( 0.23, *R.get_previousMassChangeRate( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousMassChangeRateAnswer, *R.get_previousMassChangeRateGradient( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeVelocityGradientTrace ){

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.get_massChangeVelocityGradientTrace( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousAnswer, *R.get_previousMassChangeVelocityGradientTrace( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *Rgrad.get_massChangeVelocityGradientTrace( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousAnswer, *Rgrad.get_previousMassChangeVelocityGradientTrace( ) ) );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dtrLdRho, *R.get_dMassChangeVelocityGradientTracedDensity( ) ) );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dtrLdC, *R.get_dMassChangeVelocityGradientTracedMassChangeRate( ) ) );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPrevioustrLdPreviousRho, *R.get_dPreviousMassChangeVelocityGradientTracedPreviousDensity( ) ) );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPrevioustrLdPreviousC, *R.get_dPreviousMassChangeVelocityGradientTracedPreviousMassChangeRate( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_directionVector ){

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

    Rgrad.get_dDirectionVectordMassChangeRateGradient( );

    Rgrad.get_dPreviousDirectionVectordPreviousMassChangeRateGradient( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.get_directionVector( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousAnswer, *R.get_previousDirectionVector( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *Rgrad.get_directionVector( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousAnswer, *Rgrad.get_previousDirectionVector( ) ) );

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

        floatVector vp = *Rp.get_directionVector( );

        floatVector vm = *Rm.get_directionVector( );

        for ( unsigned int j = 0; j < 3; j++ ){

            dDdGradC[ 3 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDdGradC, *R.get_dDirectionVectordMassChangeRateGradient( ) ) );

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

        floatVector vp = *Rp.get_previousDirectionVector( );

        floatVector vm = *Rm.get_previousDirectionVector( );

        for ( unsigned int j = 0; j < 3; j++ ){

            dPreviousDdPreviousGradC[ 3 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousDdPreviousGradC, *R.get_dPreviousDirectionVectordPreviousMassChangeRateGradient( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_massChangeVelocityGradient ){

    class residualMock : public tardigradeHydra::massChange::residual {

        public:

            using tardigradeHydra::massChange::residual::residual;

            floatType massChangeVelocityGradientTrace = 0.4;

            floatType previousMassChangeVelocityGradientTrace = 0.9;

            floatVector directionVector = { 0.2, 0.3, 0.4 };

            floatVector previousDirectionVector = { -0.1, 0.4, 0.30 };

        protected:

            virtual void setMassChangeVelocityGradientTrace( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_previousMassChangeVelocityGradientTrace( previousMassChangeVelocityGradientTrace );

                }
                else{

                    set_massChangeVelocityGradientTrace( massChangeVelocityGradientTrace );

                }

            }

            virtual void setDirectionVector( const bool &isPrevious ) override{

                if ( isPrevious ){

                    set_previousDirectionVector( previousDirectionVector );

                }
                else{

                    set_directionVector( directionVector );

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.get_massChangeVelocityGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( previousAnswer, *R.get_previousMassChangeVelocityGradient( ) ) );

}
