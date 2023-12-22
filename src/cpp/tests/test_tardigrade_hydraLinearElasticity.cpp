/**
  * \file test_tardigrade_hydraLinearElasticity.cpp
  * 
  * Tests for tardigrade_hydraLinearElasticity
  */

#include<tardigrade_hydraLinearElasticity.h>
#include<tardigrade_constitutive_tools.h>

#define BOOST_TEST_MODULE test_tardigrade_hydraLinearElasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::linearElasticity::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::linearElasticity::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::linearElasticity::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

namespace tardigradeHydra{

    namespace unit_test{

        class hydraBaseTester{

            public:

                static void updateUnknownVector( tardigradeHydra::hydraBase &hydra, const floatVector &value ){

                    BOOST_CHECK_NO_THROW( hydra.updateUnknownVector( value ) );

                }

        };

    }

    namespace linearElasticity{

        namespace unit_test{
    
            class residualTester{
    
                public:
    
                    static void runBasicGetTests( tardigradeHydra::linearElasticity::residual &R ){
    
                        BOOST_CHECK( &R._lambda == R.getLambda( ) );
    
                        BOOST_CHECK( &R._mu == R.getMu( ) );
    
                        BOOST_CHECK( &R._Ee.second == R.get_Ee( ) );

                        BOOST_CHECK( &R._dEedFe.second == R.get_dEedFe( ) );

                        BOOST_CHECK( &R._PK2Stress.second == R.get_PK2Stress( ) );
    
                        BOOST_CHECK( &R._dPK2StressdEe.second == R.get_dPK2StressdEe( ) );

                        BOOST_CHECK( &R._dPK2StressdFe.second == R.get_dPK2StressdFe( ) );
    
                        BOOST_CHECK( &R._dCauchyStressdPK2Stress.second == R.get_dCauchyStressdPK2Stress( ) );

                        BOOST_CHECK( &R._dCauchyStressdF.second == R.get_dCauchyStressdF( ) );
    
                        BOOST_CHECK( &R._dCauchyStressdFn.second == R.get_dCauchyStressdFn( ) );
    
                    }
    
            };
    
        }

    }

}

BOOST_AUTO_TEST_CASE( test_residual_runBasicGetTests ){

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { };

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    tardigradeHydra::hydraBase hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    tardigradeHydra::linearElasticity::unit_test::residualTester::runBasicGetTests( R );

}

BOOST_AUTO_TEST_CASE( test_residual_decomposeParameterVector ){

    class hydraMock : public tardigradeHydra::hydraBase{

    };

    hydraMock hydra;

    floatVector parameters = { 123.4, 56.7 };

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( parameters[ 0 ], *R.getLambda( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( parameters[ 1 ], *R.getMu( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setEe ){

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::linearElasticity::residual elasticity;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = tardigradeHydra::linearElasticity::residual( this, elasticitySize, *getParameters( ) );

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                residuals[ 0 ] = &elasticity;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { };

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    floatVector EeAnswer = { 0.04473512,  0.11620898, -0.13355661,
                             0.11620898, -0.24386991,  0.07603126,
                            -0.13355661,  0.07603126, -0.33822733 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( EeAnswer, *R.get_Ee( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setdEedFe ){

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::linearElasticity::residual elasticity;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = tardigradeHydra::linearElasticity::residual( this, elasticitySize, *getParameters( ) );

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                residuals[ 0 ] = &elasticity;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { };

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    floatVector EeAnswer;

    floatMatrix dEedFeAnswer;

    TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeGreenLagrangeStrain( deformationGradient, EeAnswer, dEedFeAnswer ) );

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dEedFeAnswer, *R.get_dEedFe( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setPK2Stress ){

    class residualMock : public tardigradeHydra::linearElasticity::residual {

        public:

            using tardigradeHydra::linearElasticity::residual::residual;

            void setEe( floatVector &Ee ){ tardigradeHydra::linearElasticity::residual::set_Ee( Ee ); }

            floatVector Ee = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setEe( ) override {

                setEe( Ee );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            residualMock elasticity;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = residualMock( this, elasticitySize, *getParameters( ) );

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                residuals[ 0 ] = &elasticity;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { };

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    floatVector PK2StressAnswer = { 1964.4,  226.8,  340.2,
                                     453.6, 2418. ,  680.4,
                                     793.8,  907.2, 2871.6 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( PK2StressAnswer, *R.get_PK2Stress( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setdPK2StressdEe ){

    class residualMock : public tardigradeHydra::linearElasticity::residual {

        public:

            using tardigradeHydra::linearElasticity::residual::residual;

            void setEe( floatVector &Ee ){ tardigradeHydra::linearElasticity::residual::set_Ee( Ee ); }

            floatVector Ee = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setEe( ) override {

                setEe( Ee );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            residualMock elasticity;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = residualMock( this, elasticitySize, *getParameters( ) );

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                residuals[ 0 ] = &elasticity;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { };

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    floatVector PK2StressAnswer = { 1964.4,  226.8,  340.2,
                                     453.6, 2418. ,  680.4,
                                     793.8,  907.2, 2871.6 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters );

    floatType eps = 1e-6;

    floatMatrix gradient( R.Ee.size( ), floatVector( R.Ee.size( ), 0 ) );

    for ( unsigned int i = 0; i < R.Ee.size( ); i++ ){

        floatVector delta( R.Ee.size( ), 0 );

        delta[ i ] = eps * std::fabs( R.Ee[ i ] ) + eps;

        residualMock Rp( &hydra, 9, parameters );

        residualMock Rm( &hydra, 9, parameters );

        Rp.Ee += delta;

        Rm.Ee -= delta;

        for ( unsigned int j = 0; j < R.Ee.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.get_PK2Stress( ) )[ j ] - ( *Rm.get_PK2Stress( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradient, *R.get_dPK2StressdEe( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setdPK2StressdFe ){

    class residualMock : public tardigradeHydra::linearElasticity::residual {

        public:

            using tardigradeHydra::linearElasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            residualMock elasticity;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = residualMock( this, elasticitySize, *getParameters( ) );

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                residuals[ 0 ] = &elasticity;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { };

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    floatVector PK2StressAnswer = { 1964.4,  226.8,  340.2,
                                     453.6, 2418. ,  680.4,
                                     793.8,  907.2, 2871.6 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters );

    floatType eps = 1e-6;

    floatMatrix gradient( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );


        residualMock Rp( &hydrap, 9, parameters );

        residualMock Rm( &hydram, 9, parameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.get_PK2Stress( ) )[ j ] - ( *Rm.get_PK2Stress( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradient, *R.get_dPK2StressdFe( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setStress ){

    class residualMock : public tardigradeHydra::linearElasticity::residual {

        public:

            using tardigradeHydra::linearElasticity::residual::residual;

            void setPK2Stress( floatVector &PK2 ){ tardigradeHydra::linearElasticity::residual::set_PK2Stress( PK2 ); }

            floatVector PK2Stress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPK2Stress( ) override {

                setPK2Stress( PK2Stress );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            residualMock elasticity;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = residualMock( this, elasticitySize, *getParameters( ) );

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                residuals[ 0 ] = &elasticity;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { };

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    floatVector cauchyStressAnswer = {  13.48087078,  -7.20925572, -23.74654442,
                                        -3.63500538,   1.23317213,   3.91668732,
                                       -11.2428648 ,   3.52605015,  11.10630669 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( cauchyStressAnswer, *R.getStress( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setdCauchyStressdPK2Stress ){

    class residualMock : public tardigradeHydra::linearElasticity::residual{

        public:

            using tardigradeHydra::linearElasticity::residual::residual;

            void setPK2Stress( floatVector &value ){

                tardigradeHydra::linearElasticity::residual::set_PK2Stress( value );

            }

            floatVector PK2Stress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPK2Stress( ) override {

                setPK2Stress( PK2Stress );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            residualMock elasticity;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = residualMock( this, elasticitySize, *getParameters( ) );

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                residuals[ 0 ] = &elasticity;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { };

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters );

    floatMatrix gradient( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < R.PK2Stress.size( ); i++ ){

        floatVector delta( R.PK2Stress.size( ), 0 );

        delta[ i ] = eps * std::fabs( R.PK2Stress[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters );

        residualMock Rm( &hydram, 9, parameters );

        Rp.PK2Stress += delta;

        Rm.PK2Stress -= delta;

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.getStress( ) )[ j ] - ( *Rm.getStress( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradient, *R.get_dCauchyStressdPK2Stress( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setdCauchyStressdF ){

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::linearElasticity::residual elasticity;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = tardigradeHydra::linearElasticity::residual( this, elasticitySize, *getParameters( ) );

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                residuals[ 0 ] = &elasticity;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { };

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    floatMatrix gradient( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );

        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.getStress( ) )[ j ] - ( *Rm.getStress( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradient, *R.get_dCauchyStressdF( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setdCauchyStressdFn ){

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::linearElasticity::residual elasticity;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = tardigradeHydra::linearElasticity::residual( this, elasticitySize, *getParameters( ) );

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                residuals[ 0 ] = &elasticity;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { };

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    floatMatrix gradient( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    BOOST_CHECK( R.get_dCauchyStressdFn( )->size( ) == 9 );

    for ( unsigned int i = 0; i < 9; i++ ){

        BOOST_CHECK( ( *R.get_dCauchyStressdFn( ) )[ i ].size( ) == 0 );

    }

}

BOOST_AUTO_TEST_CASE( test_residual_setdCauchyStressdFn2 ){

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::linearElasticity::residual elasticity;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = tardigradeHydra::linearElasticity::residual( this, elasticitySize, *getParameters( ) );

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                residuals[ 0 ] = &elasticity;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.1, 0.2, 0.3,  0.2,  0.4, -0.2, 0.3, 0.2, -0.1,
                                          -0.5, 0.4, 0.2, -0.2, -0.1,  0.1, 0.1, 0.3,  0.4};

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    floatMatrix gradient( deformationGradient.size( ), floatVector( 9 * ( numConfigurations - 1 ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 9 * ( numConfigurations - 1 ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );

        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.getStress( ) )[ j ] - ( *Rm.getStress( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradient, *R.get_dCauchyStressdFn( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setResidual ){

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase remainder;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = tardigradeHydra::linearElasticity::residual( this, elasticitySize, *getParameters( ) );

                remainder = tardigradeHydra::residualBase( this, 18 );

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &remainder;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.1, 0.2, 0.3,  0.2,  0.4, -0.2, 0.3, 0.2, -0.1,
                                          -0.5, 0.4, 0.2, -0.2, -0.1,  0.1, 0.1, 0.3,  0.4};

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    floatVector residualAnswer( 9, 0 );

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( residualAnswer, *R.getResidual( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setResidual2 ){

    class residualMock : public tardigradeHydra::linearElasticity::residual{

        public:

            using tardigradeHydra::linearElasticity::residual::residual;

            void setStress( floatVector &cauchyStress ){ tardigradeHydra::residualBase::setStress( cauchyStress ); }

        private:

            void setStress( ){

                floatVector cauchyStress = { 2, 3, 4, 5, 6, 7, 8, 9, 10 };

                setStress( cauchyStress );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            residualMock elasticity;

            tardigradeHydra::residualBase remainder;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = residualMock( this, elasticitySize, *getParameters( ) );

                remainder = tardigradeHydra::residualBase( this, 18 );

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &remainder;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.1, 0.2, 0.3,  0.2,  0.4, -0.2, 0.3, 0.2, -0.1,
                                          -0.5, 0.4, 0.2, -0.2, -0.1,  0.1, 0.1, 0.3,  0.4};

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    floatVector unknownVector = {   1,  1,  1,  1,  1,  1,  1,  1,  1,
                                   .1, .2, .3, .4, .5, .6, .7, .8, .9,
                                   .10, .11, .12, .13, .14, .15, .16, .17, .18 };

    floatVector residualAnswer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 9, parameters );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( residualAnswer, *R.getResidual( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setJacobian ){

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase remainder;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = tardigradeHydra::linearElasticity::residual( this, elasticitySize, *getParameters( ) );

                remainder = tardigradeHydra::residualBase( this, 18 );

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &remainder;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 1.05, 0, 0,
                                        0.00, 1, 0,
                                        0.00, 1, 1};

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.1, 0.2, 0.3,  0.2,  0.4, -0.2, 0.3, 0.2, -0.1,
                                          -0.5, 0.4, 0.2, -0.2, -0.1,  0.1, 0.1, 0.3,  0.4,
                                           0.11, 0.12, 0.13};

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    floatVector unknownVector = {   1,  1,  1,  1,  1,  1,  1,  1,  1,
                                  1.00, .00, .00, .00, 1.00, .00, .00, .00, 1.00,
                                  1.00, .00, .00, .00, 1.00, .00, .00, .00, 1.00,
                                   0.2, -0.3, 0.7 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    floatMatrix gradient( 9, floatVector( hydra.getUnknownVector( )->size( ) ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < ( 9 * numConfigurations + numNonLinearSolveStateVariables ); i++ ){

        floatVector delta( hydra.getUnknownVector( )->size( ), 0 );

        delta[ i ] = eps * std::fabs( ( *hydra.getUnknownVector( ) )[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, *hydra.getUnknownVector( ) + delta );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, *hydra.getUnknownVector( ) - delta );

        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );

        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.getResidual( ) )[ j ] - ( *Rm.getResidual( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradient, *R.getJacobian( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setdRdT ){

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase remainder;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = tardigradeHydra::linearElasticity::residual( this, elasticitySize, *getParameters( ) );

                remainder = tardigradeHydra::residualBase( this, 18 );

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &remainder;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 1.05, 0, 0,
                                        0.00, 1, 0,
                                        0.00, 1, 1};

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.1, 0.2, 0.3,  0.2,  0.4, -0.2, 0.3, 0.2, -0.1,
                                          -0.5, 0.4, 0.2, -0.2, -0.1,  0.1, 0.1, 0.3,  0.4,
                                           0.11, 0.12, 0.13};

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    floatVector unknownVector = {   1,  1,  1,  1,  1,  1,  1,  1,  1,
                                  1.00, .00, .00, .00, 1.00, .00, .00, .00, 1.00,
                                  1.00, .00, .00, .00, 1.00, .00, .00, .00, 1.00,
                                   0.2, -0.3, 0.7 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    floatVector gradient( 9 );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector delta( 1, 0 );

        delta[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + delta[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, *hydra.getUnknownVector( ) );

        hydraBaseMock hydram( time, deltaTime, temperature - delta[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, *hydra.getUnknownVector( ) );

        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );

        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );

        gradient = ( *Rp.getResidual( ) - *Rm.getResidual( ) ) / ( 2 * delta[ i ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradient, *R.getdRdT( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setdRdF ){

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase remainder;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = tardigradeHydra::linearElasticity::residual( this, elasticitySize, *getParameters( ) );

                remainder = tardigradeHydra::residualBase( this, 18 );

                std::vector< tardigradeHydra::residualBase* > residuals( 2 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &remainder;

                setResidualClasses( residuals );

            }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 5.3;

    floatType previousTemperature = 23.4;

    floatVector deformationGradient = { 1.05, 0, 0,
                                        0.00, 1, 0,
                                        0.00, 1, 1};

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { 0.1, 0.2, 0.3,  0.2,  0.4, -0.2, 0.3, 0.2, -0.1,
                                          -0.5, 0.4, 0.2, -0.2, -0.1,  0.1, 0.1, 0.3,  0.4,
                                           0.11, 0.12, 0.13};

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    floatVector unknownVector = {   1,  1,  1,  1,  1,  1,  1,  1,  1,
                                  1.00, .00, .00, .00, 1.00, .00, .00, .00, 1.00,
                                  1.00, .00, .00, .00, 1.00, .00, .00, .00, 1.00,
                                   0.2, -0.3, 0.7 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    floatMatrix gradient( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, *hydra.getUnknownVector( ) );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, *hydra.getUnknownVector( ) );

        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );

        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.getResidual( ) )[ j ] - ( *Rm.getResidual( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( gradient, *R.getdRdF( ) ) );

}
