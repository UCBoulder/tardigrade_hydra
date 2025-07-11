/**
  * \file test_tardigrade_hydraLinearViscoelasticity.cpp
  * 
  * Tests for tardigrade_hydraLinearViscoelasticity
  */

#include<tardigrade_hydraLinearViscoelasticity.h>

#define BOOST_TEST_MODULE test_tardigrade_hydraLinearViscoelasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include<tardigrade_stress_tools.h>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

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

    namespace linearViscoelasticity{

        namespace unit_test{
    
            class residualTester{
    
                public:
    
                    static void runBasicGetTests( tardigradeHydra::linearViscoelasticity::residual &R ){

                        BOOST_CHECK( &R._viscoelasticISVLowerIndex == R.getViscoelasticISVLowerIndex( ) );

                        BOOST_CHECK( &R._viscoelasticISVUpperIndex == R.getViscoelasticISVUpperIndex( ) );

                        BOOST_CHECK( &R._numVolumetricViscousTerms == R.getNumVolumetricViscousTerms( ) );

                        BOOST_CHECK( &R._numIsochoricViscousTerms == R.getNumIsochoricViscousTerms( ) );

                        BOOST_CHECK( &R._Kinf == R.getKinf( ) );
    
                        BOOST_CHECK( &R._Ginf == R.getGinf( ) );
    
                        BOOST_CHECK( &R._Ks == R.getVolumetricModuli( ) );

                        BOOST_CHECK( &R._Gs == R.getIsochoricModuli( ) );

                        BOOST_CHECK( &R._volumetricTaus == R.getVolumetricTaus( ) );

                        BOOST_CHECK( &R._isochoricTaus == R.getIsochoricTaus( ) );

                        BOOST_CHECK( &R._Je.second == R.get_Je( ) );

                        BOOST_CHECK( &R._Fehat.second == R.get_Fehat( ) );
    
                        BOOST_CHECK( &R._previousJe.second == R.get_previousJe( ) );

                        BOOST_CHECK( &R._previousFehat.second == R.get_previousFehat( ) );

                        BOOST_CHECK( &R._numStateVariables == R.getNumStateVariables( ) );

                        BOOST_CHECK( &R._volumetricTemperatureParameters == R.getVolumetricTemperatureParameters( ) );

                        BOOST_CHECK( &R._isochoricTemperatureParameters == R.getIsochoricTemperatureParameters( ) );

                        BOOST_CHECK( &R._volumetricRateMultiplier.second == R.get_volumetricRateMultiplier( ) );

                        BOOST_CHECK( &R._previousVolumetricRateMultiplier.second == R.get_previousVolumetricRateMultiplier( ) );

                        BOOST_CHECK( &R._isochoricRateMultiplier.second == R.get_isochoricRateMultiplier( ) );

                        BOOST_CHECK( &R._previousIsochoricRateMultiplier.second == R.get_previousIsochoricRateMultiplier( ) );

                        BOOST_CHECK( &R._dVolumetricRateMultiplierdT.second == R.get_dVolumetricRateMultiplierdT( ) );

                        BOOST_CHECK( &R._dPreviousVolumetricRateMultiplierdPreviousT.second == R.get_dPreviousVolumetricRateMultiplierdPreviousT( ) );

                        BOOST_CHECK( &R._dIsochoricRateMultiplierdT.second == R.get_dIsochoricRateMultiplierdT( ) );

                        BOOST_CHECK( &R._dPreviousIsochoricRateMultiplierdPreviousT.second == R.get_dPreviousIsochoricRateMultiplierdPreviousT( ) );

                        BOOST_CHECK( &R._integrationAlpha == R.getIntegrationAlpha( ) );

                        BOOST_CHECK( &R._PK2MeanStress.second == R.get_PK2MeanStress( ) );

                        BOOST_CHECK( &R._PK2IsochoricStress.second == R.get_PK2IsochoricStress( ) );

                        BOOST_CHECK( &R._dPK2MeanStressdFe.second == R.get_dPK2MeanStressdFe( ) );

                        BOOST_CHECK( &R._dPK2MeanStressdT.second == R.get_dPK2MeanStressdT( ) );

                        BOOST_CHECK( &R._dPK2IsochoricStressdFe.second == R.get_dPK2IsochoricStressdFe( ) );

                        BOOST_CHECK( &R._dPK2IsochoricStressdT.second == R.get_dPK2IsochoricStressdT( ) );

                        BOOST_CHECK( &R._dCauchyStressdT.second == R.get_dCauchyStressdT( ) );

                        BOOST_CHECK( &R._volumetricViscoelasticStateVariables.second == R.get_volumetricViscoelasticStateVariables( ) );

                        BOOST_CHECK( &R._isochoricViscoelasticStateVariables.second == R.get_isochoricViscoelasticStateVariables( ) );

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

        if ( ( std::fabs( v1[ i ] ) < tol ) || ( std::fabs( v2[ i ] ) < tol ) ){

            if ( std::fabs( v1[ i ] - v2[ i ] ) > eps ){

                return false;

            }

        }
        else{

            if ( ( std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v1[ i ] ) > eps ) ||
                 ( std::fabs( v1[ i ] - v2[ i ] ) / std::fabs( v2[ i ] ) > eps ) ){

                return false;

            }

        }

    }

    return true;

}

bool tolerantCheck( const double &v1, const double &v2, double eps = 1e-6, double tol = 1e-9 ){

    std::vector< double > _v1 = { v1 };

    std::vector< double > _v2 = { v2 };

    return tolerantCheck( _v1, _v2, eps, tol );

}

BOOST_AUTO_TEST_CASE( test_residual_runBasicGetTests_and_decomposeParameters, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

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

    floatVector previousStateVariables( 2 + 31, 0 );

    floatVector parameters = { 2, 3, 123.4, 56.7, 1, 2, 293.15, 3, 4, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    unsigned int numVolumetricViscousTermsAnswer = 2;

    unsigned int numIsochoricViscousTermsAnswer = 3;

    floatType KinfAnswer = 123.4;

    floatType GinfAnswer = 56.7;

    floatVector volTAnswer = { 1, 2, 293.15 };

    floatVector isoTAnswer = { 3, 4, 293.15 };

    floatVector Ks = { 23.4, 25.6 };

    floatVector KTaus = { 0.1, 0.2 };

    floatVector Gs = { 12.3, 13.4, 14.5 };

    floatVector GTaus = { 0.01, 10.0, 100.0 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual R( &hydra, 9, parameters, ISVlb, ISVub );

    BOOST_CHECK_NO_THROW( tardigradeHydra::linearViscoelasticity::unit_test::residualTester::runBasicGetTests( R ) );

    BOOST_TEST( ISVlb == *R.getViscoelasticISVLowerIndex( ) );

    BOOST_TEST( ISVub == *R.getViscoelasticISVUpperIndex( ) );

    BOOST_TEST( numVolumetricViscousTermsAnswer == *R.getNumVolumetricViscousTerms( ) );

    BOOST_TEST( numIsochoricViscousTermsAnswer == *R.getNumIsochoricViscousTerms( ) );

    BOOST_TEST( KinfAnswer == *R.getKinf( ) );

    BOOST_TEST( GinfAnswer == *R.getGinf( ) );

    BOOST_TEST( Ks == *R.getVolumetricModuli( ), CHECK_PER_ELEMENT );

    BOOST_TEST( Gs == *R.getIsochoricModuli( ) , CHECK_PER_ELEMENT );

    BOOST_TEST( KTaus == *R.getVolumetricTaus( ), CHECK_PER_ELEMENT );

    BOOST_TEST( GTaus == *R.getIsochoricTaus( ) , CHECK_PER_ELEMENT );

    BOOST_TEST( volTAnswer == *R.getVolumetricTemperatureParameters( ), CHECK_PER_ELEMENT );

    BOOST_TEST( isoTAnswer == *R.getIsochoricTemperatureParameters( ) , CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_decomposeElasticDeformation, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

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

    floatVector parameters = { 2, 3, 123.4, 56.7, 1, 2, 293.15, 3, 4, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    floatType JeAnswer = 0.2871738879785901;

    floatVector FehatAnswer = { 0.59558368, -0.64830483, -0.82803222,
                                0.15555742,  0.66530605, -0.23309781,
                                1.45740571,  0.56029945, -0.05780372 };

    floatType previousJeAnswer = 0.18782061270903108;

    floatVector previousFehatAnswer = { -0.37676148, -0.54767451,  0.79991771,
                                        -0.21452614, -1.53775117, -0.35606336,
                                         0.83115907, -1.10884443, -1.13343036 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual R1( &hydra, 9, parameters, ISVlb, ISVub );

    tardigradeHydra::linearViscoelasticity::residual R2( &hydra, 9, parameters, ISVlb, ISVub );

    BOOST_TEST( JeAnswer == *R1.get_Je( ) );

    BOOST_TEST( FehatAnswer == *R2.get_Fehat( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousJeAnswer == *R1.get_previousJe( ) );

    BOOST_TEST( previousFehatAnswer == *R2.get_previousFehat( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_decomposePreviousElasticDeformation, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

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

    floatVector parameters = { 2, 3, 123.4, 56.7, 1, 2, 293.15, 3, 4, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    floatType JeAnswer = 0.1878206127090311;

    floatVector FehatAnswer = { -0.37676148, -0.54767451,  0.79991771,
                                -0.21452614, -1.53775117, -0.35606336,
                                 0.83115907, -1.10884443, -1.13343036 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual R1( &hydra, 9, parameters, ISVlb, ISVub );

    tardigradeHydra::linearViscoelasticity::residual R2( &hydra, 9, parameters, ISVlb, ISVub );

    BOOST_TEST( JeAnswer == *R1.get_previousJe( ) );

    BOOST_TEST( FehatAnswer == *R2.get_previousFehat( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_gradientsOfDecomposedElasticDeformationGradient, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

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

    floatVector parameters = { 2, 3, 123.4, 56.7, 1, 2, 293.15, 3, 4, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual R( &hydra, 9, parameters, ISVlb, ISVub );

    floatVector gradientJe( deformationGradient.size( ), 0 );

    floatMatrix gradientFehat( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
        tardigradeHydra::linearViscoelasticity::residual Rp( &hydrap, 9, parameters, ISVlb, ISVub );

        tardigradeHydra::linearViscoelasticity::residual Rm( &hydram, 9, parameters, ISVlb, ISVub );

        for ( unsigned int j = 0; j < 1; j++ ){

            gradientJe[ i ] = ( ( *Rp.get_Je( ) ) - ( *Rm.get_Je( ) ) ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradientFehat[ j ][ i ] = ( ( *Rp.get_Fehat( ) )[ j ] - ( *Rm.get_Fehat( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( gradientJe == *R.get_dJedFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradientFehat ) == *R.get_dFehatdFe( ), CHECK_PER_ELEMENT );

    floatVector gradientPreviousJe( deformationGradient.size( ), 0 );

    floatMatrix gradientPreviousFehat( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
        tardigradeHydra::linearViscoelasticity::residual Rp( &hydrap, 9, parameters, ISVlb, ISVub );

        tardigradeHydra::linearViscoelasticity::residual Rm( &hydram, 9, parameters, ISVlb, ISVub );

        for ( unsigned int j = 0; j < 1; j++ ){

            gradientPreviousJe[ i ] = ( ( *Rp.get_previousJe( ) ) - ( *Rm.get_previousJe( ) ) ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradientPreviousFehat[ j ][ i ] = ( ( *Rp.get_previousFehat( ) )[ j ] - ( *Rm.get_previousFehat( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( gradientPreviousJe == *R.get_previousdJedFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradientPreviousFehat ) == *R.get_previousdFehatdFe( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_decomposeStateVariableVector, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

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

    floatVector previousStateVariables = { -1, 0, 1, 2,
                                            3,  4,  5,  6,  7,  8,  9, 10, 11,
                                           12, 13, 14, 15, 16, 17, 18, 19, 20,
                                           21, 22, 23, 24, 25, 26, 27, 28, 29 };

    floatVector volumetricStateVariablesAnswer = { 1, 2 };

    floatVector isochoricStateVariablesAnswer = {  3,  4,  5,  6,  7,  8,  9, 10, 11,
                                                  12, 13, 14, 15, 16, 17, 18, 19, 20,
                                                  21, 22, 23, 24, 25, 26, 27, 28, 29 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 1, 2, 293.15, 3, 4, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual R( &hydra, 9, parameters, ISVlb, ISVub );

    floatVector volumetricStateVariables, isochoricStateVariables;

    BOOST_CHECK_NO_THROW( R.decomposeStateVariableVector( volumetricStateVariables,
                                                          isochoricStateVariables ) );

    BOOST_TEST( volumetricStateVariables == volumetricStateVariablesAnswer, CHECK_PER_ELEMENT );

    BOOST_TEST( isochoricStateVariables == isochoricStateVariablesAnswer, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setRateMultipliers, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 310.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { -1, 0, 1, 2,
                                            3,  4,  5,  6,  7,  8,  9, 10, 11,
                                           12, 13, 14, 15, 16, 17, 18, 19, 20,
                                           21, 22, 23, 24, 25, 26, 27, 28, 29 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 10, 5, 293.15, 2, 3, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    floatType volumetricRateMultiplierAnswer = 306.5460721724233;

    floatType previousVolumetricRateMultiplierAnswer = 339.2630810101799;

    floatType isochoricRateMultiplierAnswer = 32.53168221009231;

    floatType previousIsochoricRateMultiplierAnswer = 34.516501010605985;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual R( &hydra, 9, parameters, ISVlb, ISVub );

    BOOST_TEST( volumetricRateMultiplierAnswer == *R.get_volumetricRateMultiplier( ) );

    BOOST_TEST( previousVolumetricRateMultiplierAnswer == *R.get_previousVolumetricRateMultiplier( ) );

    BOOST_TEST( isochoricRateMultiplierAnswer == *R.get_isochoricRateMultiplier( ) );

    BOOST_TEST( previousIsochoricRateMultiplierAnswer == *R.get_previousIsochoricRateMultiplier( ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setdRateMultipliersdT, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 310.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { -1, 0, 1, 2,
                                            3,  4,  5,  6,  7,  8,  9, 10, 11,
                                           12, 13, 14, 15, 16, 17, 18, 19, 20,
                                           21, 22, 23, 24, 25, 26, 27, 28, 29 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 10, 5, 293.15, 2, 3, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    floatType dVolumetricRateMultiplierdT = 0;

    floatType dPreviousVolumetricRateMultiplierdPreviousT = 0;

    floatType dIsochoricRateMultiplierdT = 0;

    floatType dPreviousIsochoricRateMultiplierdPreviousT = 0;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual R( &hydra, 9, parameters, ISVlb, ISVub );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 1; i++ ){

         floatVector deltas( 1, 0 );

        deltas[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + deltas[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature - deltas[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual Rp( &hydrap, 9, parameters, ISVlb, ISVub );

    tardigradeHydra::linearViscoelasticity::residual Rm( &hydram, 9, parameters, ISVlb, ISVub );

        dVolumetricRateMultiplierdT = ( *Rp.get_volumetricRateMultiplier( ) - *Rm.get_volumetricRateMultiplier( ) ) / ( 2 * deltas[ i ] );

        dIsochoricRateMultiplierdT = ( *Rp.get_isochoricRateMultiplier( ) - *Rm.get_isochoricRateMultiplier( ) ) / ( 2 * deltas[ i ] );

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector deltas( 1, 0 );

        deltas[ i ] = eps * std::fabs( previousTemperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature + deltas[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature - deltas[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual Rp( &hydrap, 9, parameters, ISVlb, ISVub );

    tardigradeHydra::linearViscoelasticity::residual Rm( &hydram, 9, parameters, ISVlb, ISVub );

        dPreviousVolumetricRateMultiplierdPreviousT = ( *Rp.get_previousVolumetricRateMultiplier( ) - *Rm.get_previousVolumetricRateMultiplier( ) ) / ( 2 * deltas[ i ] );

        dPreviousIsochoricRateMultiplierdPreviousT = ( *Rp.get_previousIsochoricRateMultiplier( ) - *Rm.get_previousIsochoricRateMultiplier( ) ) / ( 2 * deltas[ i ] );

    }

    BOOST_TEST( dVolumetricRateMultiplierdT == *R.get_dVolumetricRateMultiplierdT( ) );

    BOOST_TEST( dPreviousVolumetricRateMultiplierdPreviousT == *R.get_dPreviousVolumetricRateMultiplierdPreviousT( ) );

    BOOST_TEST( dIsochoricRateMultiplierdT == *R.get_dIsochoricRateMultiplierdT( ) );

    BOOST_TEST( dPreviousIsochoricRateMultiplierdPreviousT == *R.get_dPreviousIsochoricRateMultiplierdPreviousT( ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getViscoelasticParameters, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 310.4;

    floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    floatVector previousStateVariables = { -1, 0, 1, 2,
                                            3,  4,  5,  6,  7,  8,  9, 10, 11,
                                           12, 13, 14, 15, 16, 17, 18, 19, 20,
                                           21, 22, 23, 24, 25, 26, 27, 28, 29 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 1, 100, 293.15, 2, 110, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    floatVector volumetricViscoelasticParametersAnswer = { 123.4, 0.1, 0.2, 23.4, 25.6 };

    floatVector isochoricViscoelasticParametersAnswer = { 113.4, 0.01, 10.0, 100, 24.6, 26.8, 29.0 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual R( &hydra, 9, parameters, ISVlb, ISVub );

    BOOST_TEST( volumetricViscoelasticParametersAnswer == R.getVolumetricViscoelasticParameters( ), CHECK_PER_ELEMENT );

    BOOST_TEST( isochoricViscoelasticParametersAnswer  == R.getIsochoricViscoelasticParameters( ) , CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setPK2MeanStress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearViscoelasticity::residual {

        public:

            floatVector stateVariables = { 0.001, 0.002 };

            floatType Je = 1.1;

            floatType previousJe = 0.8;

            using tardigradeHydra::linearViscoelasticity::residual::residual;

            void setJe( const floatType &Je ){

                tardigradeHydra::linearViscoelasticity::residual::set_Je( Je );

            }

            void setPreviousJe( const floatType &Je ){

                tardigradeHydra::linearViscoelasticity::residual::set_previousJe( Je );

            }

            void decomposeStateVariableVector( floatVector &volumetricISVs,
                                               floatVector &isochoricISVs ){

                volumetricISVs = stateVariables;

            }

        private:

            virtual void decomposeElasticDeformation( ) override {

                setJe( Je );

            }

            virtual void decomposePreviousElasticDeformation( ) override {

                setPreviousJe( previousJe );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = { 1.1, 0.0, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0 };

    floatVector previousDeformationGradient = { 1.0, 0.0, 0.0,
                                                0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0 };

    floatVector previousStateVariables = { -1, 0, 1e-3, 1e-1,
                                            3,  4,  5,  6,  7,  8,  9, 10, 11,
                                           12, 13, 14, 15, 16, 17, 18, 19, 20,
                                           21, 22, 23, 24, 25, 26, 27, 28, 29 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 10.0, 5.0, 293.15, 2, 110, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    floatType currentRateModifier;
    tardigradeConstitutiveTools::WLF( temperature, { parameters[ 4 ], parameters[ 5 ], 293.15 }, currentRateModifier );
    currentRateModifier = 1. / currentRateModifier;

    floatType previousRateModifier;
    tardigradeConstitutiveTools::WLF( previousTemperature, { parameters[ 4 ], parameters[ 5 ], 293.15 }, previousRateModifier );
    previousRateModifier = 1. / previousRateModifier;

    floatVector ISVs = { 1e-3, 2e-3 };
    floatVector params = { 123.4, 0.1, 0.2, 23.4, 25.6 };
    floatVector PK2MeanStressAnswer;
    floatVector updatedISVsAnswer;
    floatVector previousPK2MeanStressAnswer;
    floatVector previousUpdatedISVsAnswer;

    tardigradeStressTools::linearViscoelasticity(  1.1, {  0.1 }, -1.1, { -0.2 }, currentRateModifier, previousRateModifier, ISVs, params, 0.5, PK2MeanStressAnswer, updatedISVsAnswer );

    tardigradeStressTools::linearViscoelasticity( -1.1, { -0.2 }, -1.1, { -0.2 }, previousRateModifier, previousRateModifier, ISVs, params, 0.5, previousPK2MeanStressAnswer, previousUpdatedISVsAnswer );

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters, ISVlb, ISVub, 0.5 );

    BOOST_TEST( PK2MeanStressAnswer[ 0 ] == *R.get_PK2MeanStress( ) );

    BOOST_TEST( updatedISVsAnswer == *R.get_volumetricViscoelasticStateVariables( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousPK2MeanStressAnswer[ 0 ] == *R.get_previousPK2MeanStress( ) );
}

BOOST_AUTO_TEST_CASE( test_residual_setPK2MeanStressDerivatives, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearViscoelasticity::residual {

        public:

            using tardigradeHydra::linearViscoelasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = { 1.1, 0.0, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0 };

    floatVector previousDeformationGradient = { 1.0, 0.0, 0.0,
                                                0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0 };

    floatVector previousStateVariables = { -1, 0, 1e-3, 2e-3,
                                            3,  4,  5,  6,  7,  8,  9, 10, 11,
                                           12, 13, 14, 15, 16, 17, 18, 19, 20,
                                           21, 22, 23, 24, 25, 26, 27, 28, 29 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 10.0, 5.0, 293.15, 2, 110, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters, ISVlb, ISVub, 0.5 );

    floatType eps = 1e-6;

    floatVector dPK2MeanStressdFe( deformationGradient.size( ), 0 );

    floatType   dPK2MeanStressdT;

    floatVector dPK2MeanStressdPreviousFe( deformationGradient.size( ), 0 );

    floatType   dPK2MeanStressdPreviousT;

    floatVector dPK2MeanStressdPreviousISVs( 2 + 3*9, 0 );

    floatMatrix dVolumetricISVsdFe( 2, floatVector( 9, 0 ) );

    floatVector dVolumetricISVsdT( 2, 0 );

    floatMatrix dVolumetricISVsdPreviousFe( 2, floatVector( deformationGradient.size( ), 0 ) );

    floatVector dVolumetricISVsdPreviousT( 2, 0 );

    floatMatrix dVolumetricISVsdPreviousISVs( 2, floatVector( 2 + 3*9, 0 ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector deltas( deformationGradient.size( ), 0 );

        deltas[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + deltas, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - deltas, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < 1; j++ ){

            dPK2MeanStressdFe[ i ] = ( ( *Rp.get_PK2MeanStress( ) ) - ( *Rm.get_PK2MeanStress( ) ) ) / ( 2 * deltas[ i ] );

        }

        for ( unsigned int j = 0; j < 2; j++ ){

            dVolumetricISVsdFe[ j ][ i ] = ( ( *Rp.get_volumetricViscoelasticStateVariables( ) )[ j ] - ( *Rm.get_volumetricViscoelasticStateVariables( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector deltas( 1, 0 );

        deltas[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + deltas[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature - deltas[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < 1; j++ ){

            dPK2MeanStressdT = ( ( *Rp.get_PK2MeanStress( ) ) - ( *Rm.get_PK2MeanStress( ) ) ) / ( 2 * deltas[ i ] );

        }

        for ( unsigned int j = 0; j < 2; j++ ){

            dVolumetricISVsdT[ j ] = ( ( *Rp.get_volumetricViscoelasticStateVariables( ) )[ j ] - ( *Rm.get_volumetricViscoelasticStateVariables( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    BOOST_TEST( dPK2MeanStressdFe == *R.get_dPK2MeanStressdFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dPK2MeanStressdT == *R.get_dPK2MeanStressdT( ) );

    floatVector previousdPK2MeanStressdFe( deformationGradient.size( ), 0 );

    floatType previousdPK2MeanStressdT;

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector deltas( deformationGradient.size( ), 0 );

        deltas[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + deltas,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - deltas,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < 1; j++ ){

            dPK2MeanStressdPreviousFe[ i ] = ( ( *Rp.get_PK2MeanStress( ) ) - ( *Rm.get_PK2MeanStress( ) ) ) / ( 2 * deltas[ i ] );

            previousdPK2MeanStressdFe[ i ] = ( ( *Rp.get_previousPK2MeanStress( ) ) - ( *Rm.get_previousPK2MeanStress( ) ) ) / ( 2 * deltas[ i ] );

        }

        for ( unsigned int j = 0; j < 2; j++ ){

            dVolumetricISVsdPreviousFe[ j ][ i ] = ( ( *Rp.get_volumetricViscoelasticStateVariables( ) )[ j ] - ( *Rm.get_volumetricViscoelasticStateVariables( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector deltas( 1, 0 );

        deltas[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature + deltas[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature - deltas[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < 1; j++ ){

            dPK2MeanStressdPreviousT = ( ( *Rp.get_PK2MeanStress( ) ) - ( *Rm.get_PK2MeanStress( ) ) ) / ( 2 * deltas[ i ] );

            previousdPK2MeanStressdT = ( ( *Rp.get_previousPK2MeanStress( ) ) - ( *Rm.get_previousPK2MeanStress( ) ) ) / ( 2 * deltas[ i ] );

        }

        for ( unsigned int j = 0; j < 2; j++ ){

            dVolumetricISVsdPreviousT[ j ] = ( ( *Rp.get_volumetricViscoelasticStateVariables( ) )[ j ] - ( *Rm.get_volumetricViscoelasticStateVariables( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 2 + 3*9; i++ ){

        floatVector deltas( previousStateVariables.size( ), 0 );

        deltas[ i + ISVlb ] = eps * std::fabs( previousStateVariables[ i + ISVlb ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + deltas, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - deltas, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < 1; j++ ){

            dPK2MeanStressdPreviousISVs[ i ] = ( ( *Rp.get_PK2MeanStress( ) ) - ( *Rm.get_PK2MeanStress( ) ) ) / ( 2 * deltas[ i + ISVlb ] );

        }

        for ( unsigned int j = 0; j < 2; j++ ){

            dVolumetricISVsdPreviousISVs[ j ][ i ] = ( ( *Rp.get_volumetricViscoelasticStateVariables( ) )[ j ] - ( *Rm.get_volumetricViscoelasticStateVariables( ) )[ j ] ) / ( 2 * deltas[ i + ISVlb ] );

        }

    }


    BOOST_TEST( previousdPK2MeanStressdFe == *R.get_previousdPK2MeanStressdFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousdPK2MeanStressdT == *R.get_previousdPK2MeanStressdT( ) );

    BOOST_TEST( dPK2MeanStressdPreviousFe == *R.get_dPK2MeanStressdPreviousFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dPK2MeanStressdPreviousT == *R.get_dPK2MeanStressdPreviousT( ) );

    BOOST_TEST( dPK2MeanStressdPreviousISVs == *R.get_dPK2MeanStressdPreviousISVs( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dVolumetricISVsdFe ) == *R.get_dVolumetricISVsdFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dVolumetricISVsdT == *R.get_dVolumetricISVsdT( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dVolumetricISVsdPreviousFe ) == *R.get_dVolumetricISVsdPreviousFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dVolumetricISVsdPreviousT == *R.get_dVolumetricISVsdPreviousT( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dVolumetricISVsdPreviousISVs ) == *R.get_dVolumetricISVsdPreviousISVs( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setPK2IsochoricStress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearViscoelasticity::residual {

        public:

            using tardigradeHydra::linearViscoelasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = { 1.1, 0.1, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0 };

    floatVector previousDeformationGradient = { 1.0, 0.0, 0.0,
                                                0.1, 1.0, 0.0,
                                                0.0, 0.0, 1.0 };

    floatVector previousStateVariables = { -1, 0, 1e-3, 1e-1,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 10.0, 5.0, 293.15, 2, 3, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.1, 1.0, 10.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    floatVector strain;
    tardigradeConstitutiveTools::computeGreenLagrangeStrain( deformationGradient / std::pow( 1.1, 1/3.), strain );

    floatVector previousStrain;
    tardigradeConstitutiveTools::computeGreenLagrangeStrain( previousDeformationGradient, previousStrain );

    floatType currentRateModifier;
    tardigradeConstitutiveTools::WLF( temperature, { parameters[ 7 ], parameters[ 8 ], 293.15 }, currentRateModifier );
    currentRateModifier = 1. / currentRateModifier;

    floatType previousRateModifier;
    tardigradeConstitutiveTools::WLF( previousTemperature, { parameters[ 7 ], parameters[ 8 ], 293.15 }, previousRateModifier );
    previousRateModifier = 1. / previousRateModifier;

    floatVector ISVs( previousStateVariables.begin( ) + ISVlb + parameters[ 0 ],
                      previousStateVariables.end( ) );
    floatVector params = { 2 * 56.7, 0.1, 1.0, 10.0, 2 * 12.3, 2 * 13.4, 2 * 14.5 };
    floatVector PK2IsochoricStressAnswer;
    floatVector updatedISVsAnswer;
    floatVector previousPK2IsochoricStressAnswer;
    floatVector previousUpdatedISVsAnswer;

    tardigradeStressTools::linearViscoelasticity(  1.1,         strain, -1.1, previousStrain,  currentRateModifier, previousRateModifier, ISVs, params, 0.5,        PK2IsochoricStressAnswer, updatedISVsAnswer );

    tardigradeStressTools::linearViscoelasticity( -1.1, previousStrain, -1.1, previousStrain, previousRateModifier, previousRateModifier, ISVs, params, 0.5, previousPK2IsochoricStressAnswer, previousUpdatedISVsAnswer );

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters, ISVlb, ISVub, 0.5 );

    BOOST_TEST( PK2IsochoricStressAnswer == *R.get_PK2IsochoricStress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( updatedISVsAnswer == *R.get_isochoricViscoelasticStateVariables( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousPK2IsochoricStressAnswer == *R.get_previousPK2IsochoricStress( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setPK2IsochoricStressDerivatives, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearViscoelasticity::residual {

        public:

            using tardigradeHydra::linearViscoelasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

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

    floatVector previousStateVariables = { -1, 0, 1e-3, 1e-1,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 10.0, 5.0, 293.15, 2, 3, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.1, 1.0, 10.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters, ISVlb, ISVub, 0.5 );

    floatType eps = 1e-6;

    floatMatrix dPK2IsochoricStressdFe( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatVector dPK2IsochoricStressdT( deformationGradient.size( ), 0 );

    floatMatrix dPK2IsochoricStressdPreviousFe( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatVector dPK2IsochoricStressdPreviousT( deformationGradient.size( ), 0 );

    floatMatrix dPK2IsochoricStressdPreviousISVs( deformationGradient.size( ), floatVector( 2 + 3*9, 0 ) );

    floatMatrix dIsochoricISVsdFe( 3*9, floatVector( deformationGradient.size( ), 0 ) );

    floatVector dIsochoricISVsdT( 3*9, 0 );

    floatMatrix dIsochoricISVsdPreviousFe( 3*9, floatVector( deformationGradient.size( ), 0 ) );

    floatVector dIsochoricISVsdPreviousT( 3*9, 0 );

    floatMatrix dIsochoricISVsdPreviousISVs( 3*9, floatVector( 2 + 3*9, 0 ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector deltas( deformationGradient.size( ), 0 );

        deltas[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + deltas, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - deltas, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dPK2IsochoricStressdFe[ j ][ i ] = ( ( *Rp.get_PK2IsochoricStress( ) )[ j ] - ( *Rm.get_PK2IsochoricStress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

        for ( unsigned int j = 0; j < 3*9; j++ ){

            dIsochoricISVsdFe[ j ][ i ] = ( ( *Rp.get_isochoricViscoelasticStateVariables( ) )[ j ] - ( *Rm.get_isochoricViscoelasticStateVariables( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector deltas( 1, 0 );

        deltas[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + deltas[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature - deltas[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dPK2IsochoricStressdT[ j ] = ( ( *Rp.get_PK2IsochoricStress( ) )[ j ] - ( *Rm.get_PK2IsochoricStress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

        for ( unsigned int j = 0; j < 3*9; j++ ){

            dIsochoricISVsdT[ j ] = ( ( *Rp.get_isochoricViscoelasticStateVariables( ) )[ j ] - ( *Rm.get_isochoricViscoelasticStateVariables( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPK2IsochoricStressdFe ) == *R.get_dPK2IsochoricStressdFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dPK2IsochoricStressdT == *R.get_dPK2IsochoricStressdT( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dIsochoricISVsdFe ) == *R.get_dIsochoricISVsdFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dIsochoricISVsdT == *R.get_dIsochoricISVsdT( ), CHECK_PER_ELEMENT );

    floatMatrix previousdPK2IsochoricStressdFe( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatVector previousdPK2IsochoricStressdT( deformationGradient.size( ), 0 );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector deltas( deformationGradient.size( ), 0 );

        deltas[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + deltas,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - deltas,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dPK2IsochoricStressdPreviousFe[ j ][ i ] = ( ( *Rp.get_PK2IsochoricStress( ) )[ j ] - ( *Rm.get_PK2IsochoricStress( ) )[ j ] ) / ( 2 * deltas[ i ] );

            previousdPK2IsochoricStressdFe[ j ][ i ] = ( ( *Rp.get_previousPK2IsochoricStress( ) )[ j ] - ( *Rm.get_previousPK2IsochoricStress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

        for ( unsigned int j = 0; j < 3*9; j++ ){

            dIsochoricISVsdPreviousFe[ j ][ i ] = ( ( *Rp.get_isochoricViscoelasticStateVariables( ) )[ j ] - ( *Rm.get_isochoricViscoelasticStateVariables( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector deltas( 1, 0 );

        deltas[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature + deltas[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature - deltas[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dPK2IsochoricStressdPreviousT[ j ] = ( ( *Rp.get_PK2IsochoricStress( ) )[ j ] - ( *Rm.get_PK2IsochoricStress( ) )[ j ] ) / ( 2 * deltas[ i ] );

            previousdPK2IsochoricStressdT[ j ] = ( ( *Rp.get_previousPK2IsochoricStress( ) )[ j ] - ( *Rm.get_previousPK2IsochoricStress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

        for ( unsigned int j = 0; j < 3*9; j++ ){

            dIsochoricISVsdPreviousT[ j ] = ( ( *Rp.get_isochoricViscoelasticStateVariables( ) )[ j ] - ( *Rm.get_isochoricViscoelasticStateVariables( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPK2IsochoricStressdPreviousFe ) == *R.get_dPK2IsochoricStressdPreviousFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dPK2IsochoricStressdPreviousT == *R.get_dPK2IsochoricStressdPreviousT( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dIsochoricISVsdPreviousFe ) == *R.get_dIsochoricISVsdPreviousFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dIsochoricISVsdPreviousT == *R.get_dIsochoricISVsdPreviousT( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdPK2IsochoricStressdFe ) == *R.get_previousdPK2IsochoricStressdFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousdPK2IsochoricStressdT == *R.get_previousdPK2IsochoricStressdT( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 2 + 3*9; i++ ){

        floatVector deltas( previousStateVariables.size( ), 0 );

        deltas[ i + ISVlb ] = eps * std::fabs( previousStateVariables[ i + ISVlb ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + deltas, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - deltas, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dPK2IsochoricStressdPreviousISVs[ j ][ i ] = ( ( *Rp.get_PK2IsochoricStress( ) )[ j ] - ( *Rm.get_PK2IsochoricStress( ) )[ j ] ) / ( 2 * deltas[ i + ISVlb ] );

        }

        for ( unsigned int j = 0; j < 3*9; j++ ){

            dIsochoricISVsdPreviousISVs[ j ][ i ] = ( ( *Rp.get_isochoricViscoelasticStateVariables( ) )[ j ] - ( *Rm.get_isochoricViscoelasticStateVariables( ) )[ j ] ) / ( 2 * deltas[ i + ISVlb ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPK2IsochoricStressdPreviousISVs ) == *R.get_dPK2IsochoricStressdPreviousISVs( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dIsochoricISVsdPreviousISVs ) == *R.get_dIsochoricISVsdPreviousISVs( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setPK2Stress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearViscoelasticity::residual {

        public:

            using tardigradeHydra::linearViscoelasticity::residual::residual;

            floatType PK2MeanStress = 12.3;

            floatVector PK2IsochoricStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatType previousPK2MeanStress = 45.6;

            floatVector previousPK2IsochoricStress = { -1, -2, -3, -4, -5, -6, -7, -8, -9 };

            void setPK2MeanStress( const floatType &value ){

                tardigradeHydra::linearViscoelasticity::residual::set_PK2MeanStress( value );

            }

            void setPK2IsochoricStress( const floatVector &value ){

                tardigradeHydra::linearViscoelasticity::residual::set_PK2IsochoricStress( value );

            }

            void setPreviousPK2MeanStress( const floatType &value ){

                tardigradeHydra::linearViscoelasticity::residual::set_previousPK2MeanStress( value );

            }

            void setPreviousPK2IsochoricStress( const floatVector &value ){

                tardigradeHydra::linearViscoelasticity::residual::set_previousPK2IsochoricStress( value );

            }

        private:

            using tardigradeHydra::linearViscoelasticity::residual::setPK2MeanStress;

            using tardigradeHydra::linearViscoelasticity::residual::setPK2IsochoricStress;

            virtual void setPK2MeanStress( ) override {

                setPK2MeanStress( PK2MeanStress );

            }

            virtual void setPK2IsochoricStress( ) override {

                setPK2IsochoricStress( PK2IsochoricStress );

            }

            virtual void setPreviousPK2MeanStress( ) override {

                setPreviousPK2MeanStress( previousPK2MeanStress );

            }

            virtual void setPreviousPK2IsochoricStress( ) override {

                setPreviousPK2IsochoricStress( previousPK2IsochoricStress );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

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

    floatVector previousStateVariables = { -1, 0, 1e-3, 1e-1,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 10.0, 5.0, 293.15, 2, 3, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.1, 1.0, 10.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    floatVector PK2StressAnswer = { 13.3, 2, 3, 4, 17.3, 6, 7, 8, 21.3 }; 

    floatVector previousPK2StressAnswer = { 44.6, -2. , -3. , -4. , 40.6, -6. , -7. , -8. , 36.6 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters, ISVlb, ISVub, 0.5 );

    BOOST_TEST( PK2StressAnswer == *R.get_PK2Stress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousPK2StressAnswer == *R.get_previousPK2Stress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( { *R.get_volumetricViscoelasticStateVariables( ), *R.get_isochoricViscoelasticStateVariables( ) } ) == *R.getCurrentAdditionalStateVariables( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setPK2StressDerivatives, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearViscoelasticity::residual {

        public:

            using tardigradeHydra::linearViscoelasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

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

    floatVector previousStateVariables = { -1, 0, 1e-3, 1e-1,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 10.0, 5.0, 293.15, 2, 3, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.1, 1.0, 10.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    floatMatrix dPK2StressdFe( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatVector dPK2StressdT( deformationGradient.size( ), 0 );

    floatMatrix dPK2StressdPreviousFe( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatVector dPK2StressdPreviousT( deformationGradient.size( ), 0 );

    floatMatrix dPK2StressdPreviousISVs( deformationGradient.size( ), floatVector( 2+3*9, 0 ) );

    floatMatrix previousdPK2StressdFe( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatVector previousdPK2StressdT( deformationGradient.size( ), 0 );

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters, ISVlb, ISVub, 0.5 );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector deltas( deformationGradient.size( ), 0 );

        deltas[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + deltas, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - deltas, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dPK2StressdFe[ j ][ i ] = ( ( *Rp.get_PK2Stress( ) )[ j ] - ( *Rm.get_PK2Stress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector deltas( 1, 0 );

        deltas[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + deltas[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature - deltas[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dPK2StressdT[ j ] = ( ( *Rp.get_PK2Stress( ) )[ j ] - ( *Rm.get_PK2Stress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPK2StressdFe ) == *R.get_dPK2StressdFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dPK2StressdT == *R.get_dPK2StressdT( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector deltas( deformationGradient.size( ), 0 );

        deltas[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + deltas,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - deltas,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dPK2StressdPreviousFe[ j ][ i ] = ( ( *Rp.get_PK2Stress( ) )[ j ] - ( *Rm.get_PK2Stress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            previousdPK2StressdFe[ j ][ i ] = ( ( *Rp.get_previousPK2Stress( ) )[ j ] - ( *Rm.get_previousPK2Stress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector deltas( 1, 0 );

        deltas[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature + deltas[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature - deltas[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dPK2StressdPreviousT[ j ] = ( ( *Rp.get_PK2Stress( ) )[ j ] - ( *Rm.get_PK2Stress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            previousdPK2StressdT[ j ] = ( ( *Rp.get_previousPK2Stress( ) )[ j ] - ( *Rm.get_previousPK2Stress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 2+3*9; i++ ){

        floatVector deltas( previousStateVariables.size( ), 0 );

        deltas[ i + ISVlb ] = eps * std::fabs( previousStateVariables[ i + ISVlb ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + deltas, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - deltas, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dPK2StressdPreviousISVs[ j ][ i ] = ( ( *Rp.get_PK2Stress( ) )[ j ] - ( *Rm.get_PK2Stress( ) )[ j ] ) / ( 2 * deltas[ i + ISVlb ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdPK2StressdFe ) == *R.get_previousdPK2StressdFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousdPK2StressdT == *R.get_previousdPK2StressdT( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPK2StressdPreviousFe ) == *R.get_dPK2StressdPreviousFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dPK2StressdPreviousT == *R.get_dPK2StressdPreviousT( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPK2StressdPreviousISVs ) == *R.get_dPK2StressdPreviousISVs( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setdCauchyStressdT, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearViscoelasticity::residual {

        public:

            using tardigradeHydra::linearViscoelasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

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

    floatVector previousStateVariables = { -1, 0, 1e-3, 1e-1,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 10.0, 5.0, 293.15, 2, 3, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.1, 1.0, 10.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    floatMatrix dCauchyStressdF( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) ); // NOTE: This isn't strictly required! I'm just doing it for the warm fuzzies.

    floatVector dCauchyStressdT( deformationGradient.size( ), 0 );

    floatMatrix dCauchyStressdPreviousF( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) ); // NOTE: This isn't strictly required! I'm just doing it for the warm fuzzies.

    floatVector dCauchyStressdPreviousT( deformationGradient.size( ), 0 );

    floatMatrix dCauchyStressdPreviousISVs( deformationGradient.size( ), floatVector( 2+3*9, 0 ) );

    floatMatrix previousdCauchyStressdF( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) ); // NOTE: This isn't strictly required! I'm just doing it for the warm fuzzies.

    floatVector previousdCauchyStressdT( deformationGradient.size( ), 0 );

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters, ISVlb, ISVub, 0.5 );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector deltas( deformationGradient.size( ), 0 );

        deltas[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + deltas, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - deltas, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dCauchyStressdF[ j ][ i ] = ( ( *Rp.getStress( ) )[ j ] - ( *Rm.getStress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector deltas( 1, 0 );

        deltas[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + deltas[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature - deltas[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dCauchyStressdT[ j ] = ( ( *Rp.getStress( ) )[ j ] - ( *Rm.getStress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dCauchyStressdF ) == *R.get_dCauchyStressdF( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dCauchyStressdT == *R.get_dCauchyStressdT( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector deltas( deformationGradient.size( ), 0 );

        deltas[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + deltas,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - deltas,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dCauchyStressdPreviousF[ j ][ i ] = ( ( *Rp.getStress( ) )[ j ] - ( *Rm.getStress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }
        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            previousdCauchyStressdF[ j ][ i ] = ( ( *Rp.getPreviousStress( ) )[ j ] - ( *Rm.getPreviousStress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector deltas( 1, 0 );

        deltas[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature + deltas[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature - deltas[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dCauchyStressdPreviousT[ j ] = ( ( *Rp.getStress( ) )[ j ] - ( *Rm.getStress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }
        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dCauchyStressdT[ j ] = ( ( *Rp.getPreviousStress( ) )[ j ] - ( *Rm.getPreviousStress( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 2+3*9; i++ ){

        floatVector deltas( previousStateVariables.size( ), 0 );

        deltas[ i + ISVlb ] = eps * std::fabs( previousStateVariables[ i + ISVlb ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + deltas, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - deltas, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dCauchyStressdPreviousISVs[ j ][ i ] = ( ( *Rp.getStress( ) )[ j ] - ( *Rm.getStress( ) )[ j ] ) / ( 2 * deltas[ i + ISVlb ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdCauchyStressdF ) == *R.get_previousdCauchyStressdF( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousdCauchyStressdT == *R.get_previousdCauchyStressdT( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dCauchyStressdPreviousF ) == *R.get_dCauchyStressdPreviousF( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dCauchyStressdPreviousT == *R.get_dCauchyStressdPreviousT( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dCauchyStressdPreviousISVs ) == *R.get_dCauchyStressdPreviousISVs( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setdRdT, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearViscoelasticity::residual {

        public:

            using tardigradeHydra::linearViscoelasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            virtual void setResidualClasses( ){ }

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

    floatVector previousStateVariables = { -1, 0, 1e-3, 1e-1,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 10.0, 5.0, 293.15, 2, 3, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.1, 1.0, 10.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    floatMatrix dRdF( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) ); // NOTE: This isn't strictly required! I'm just doing it for the warm fuzzies.

    floatVector dRdT( deformationGradient.size( ), 0 );

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 9, parameters, ISVlb, ISVub, 0.5 );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector deltas( deformationGradient.size( ), 0 );

        deltas[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + deltas, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );
        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - deltas, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );
        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dRdF[ j ][ i ] = ( ( *Rp.getResidual( ) )[ j ] - ( *Rm.getResidual( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector deltas( 1, 0 );

        deltas[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + deltas[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );
        residualMock Rp( &hydrap, 9, parameters, ISVlb, ISVub, 0.5 );

        hydraBaseMock hydram( time, deltaTime, temperature - deltas[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );
        residualMock Rm( &hydram, 9, parameters, ISVlb, ISVub, 0.5 );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            dRdT[ j ] = ( ( *Rp.getResidual( ) )[ j ] - ( *Rm.getResidual( ) )[ j ] ) / ( 2 * deltas[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dRdF ) == *R.getdRdF( ), CHECK_PER_ELEMENT );

    BOOST_TEST( dRdT == *R.getdRdT( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_evaluate, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearViscoelasticity::residual {

        public:

            using tardigradeHydra::linearViscoelasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            residualMock elasticity;

            unsigned int ISVlb = 2;

            unsigned int ISVub = 31;

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                elasticity = residualMock( this, 9, *getParameters( ), ISVlb, ISVub, 0.5 );

                residuals[ 0 ] = &elasticity;

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

    floatVector previousStateVariables = { -1, 0, 1e-3, 1e-1,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 10.0, 5.0, 293.15, 2, 3, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.1, 1.0, 10.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    hydra.computeTangents( );

    floatType eps = 1e-6;

    {

        constexpr unsigned int OUT_NUM = 9;
        constexpr unsigned int VAR_NUM = 9;
        floatVector x = deformationGradient;

        for ( unsigned int i = 0; i < VAR_NUM; ++i ){

            floatType delta = eps * std::fabs( x[ i ] ) + eps;

            floatVector xp = x;
            floatVector xm = x;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            hydraBaseMock hydrap(
                time, deltaTime, temperature, previousTemperature, xp, previousDeformationGradient,
                { }, { },
                previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension
            );

            hydraBaseMock hydram(
                time, deltaTime, temperature, previousTemperature, xm, previousDeformationGradient,
                { }, { },
                previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension
            );

            hydrap.evaluate( );
            hydram.evaluate( );

            for ( unsigned int j = 0; j < OUT_NUM; ++j ){

                BOOST_TEST( ( ( *hydrap.getUnknownVector( ) )[ j ] - ( *hydram.getUnknownVector( ) )[ j ] ) / ( 2 * delta ) == ( *hydra.getFlatdXdF( ) )[ VAR_NUM * j + i ] );

            }

        }

    }

    {

        constexpr unsigned int OUT_NUM = 9;
        constexpr unsigned int VAR_NUM = 1;
        floatVector x = { temperature };

        for ( unsigned int i = 0; i < VAR_NUM; ++i ){

            floatType delta = eps * std::fabs( x[ i ] ) + eps;

            floatVector xp = x;
            floatVector xm = x;

            xp[ i ] += delta;
            xm[ i ] -= delta;

            hydraBaseMock hydrap(
                time, deltaTime, xp[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                { }, { },
                previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension
            );

            hydraBaseMock hydram(
                time, deltaTime, xm[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                { }, { },
                previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension
            );

            hydrap.evaluate( );
            hydram.evaluate( );

            for ( unsigned int j = 0; j < OUT_NUM; ++j ){

                BOOST_TEST( ( ( *hydrap.getUnknownVector( ) )[ j ] - ( *hydram.getUnknownVector( ) )[ j ] ) / ( 2 * delta ) == ( *hydra.getFlatdXdT( ) )[ VAR_NUM * j + i ] );

            }

        }

    }

}

BOOST_AUTO_TEST_CASE( test_addParameterizationInfo, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearViscoelasticity::residual {

        public:

            using tardigradeHydra::linearViscoelasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            residualMock elasticity;

            unsigned int ISVlb = 2;

            unsigned int ISVub = 31;

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                elasticity = residualMock( this, 9, *getParameters( ), ISVlb, ISVub, 0.5 );

                residuals[ 0 ] = &elasticity;

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

    floatVector previousStateVariables = { -1, 0, 1e-3, 1e-1,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 10.0, 5.0, 293.15, 2, 3, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.1, 1.0, 10.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters, ISVlb, ISVub, 0.5 );

    std::string output;
    R.addParameterizationInfo( output );

}

BOOST_AUTO_TEST_CASE( test_residual_updateAdditionalStateVariables, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearViscoelasticity::residual {

        public:

            using tardigradeHydra::linearViscoelasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            residualMock elasticity;

            unsigned int ISVlb = 2;

            unsigned int ISVub = 31;

            using tardigradeHydra::hydraBase::hydraBase;

        private:

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 1 );

                elasticity = residualMock( this, 9, *getParameters( ), ISVlb, ISVub, 0.5 );

                residuals[ 0 ] = &elasticity;

                setResidualClasses( residuals );

            }

    };


    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = { 1.1, 0.0, 0.0,
                                        0.0, 1.0, 0.0,
                                        0.0, 0.0, 1.0 };

    floatVector previousDeformationGradient = { 1.0, 0.0, 0.0,
                                                0.0, 1.0, 0.0,
                                                0.0, 0.0, 1.0 };

    floatVector previousStateVariables = { -1, 0, 1e-3, 1e-1,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0 };

    floatVector parameters = { 2, 3, 123.4, 56.7, 10.0, 5.0, 293.15, 2, 3, 293.15, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.1, 1.0, 10.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    floatType currentRateModifier;
    tardigradeConstitutiveTools::WLF( temperature, { parameters[ 4 ], parameters[ 5 ], 293.15 }, currentRateModifier );
    currentRateModifier = 1. / currentRateModifier;

    floatType previousRateModifier;
    tardigradeConstitutiveTools::WLF( previousTemperature, { parameters[ 4 ], parameters[ 5 ], 293.15 }, previousRateModifier );
    previousRateModifier = 1. / previousRateModifier;

    floatVector ISVs = { 1e-3, 1e-1 };
    floatVector params = { 123.4, 0.1, 0.2, 23.4, 25.6 };
    floatVector PK2MeanStressAnswer;
    floatVector updatedVolumetricISVsAnswer;
    floatVector previousPK2MeanStressAnswer;
    floatVector previousUpdatedVolumetricISVsAnswer;

    tardigradeStressTools::linearViscoelasticity(  1.1, {  0.1 }, -1.1, { 0.0 }, currentRateModifier, previousRateModifier, ISVs, params, 0.5, PK2MeanStressAnswer, updatedVolumetricISVsAnswer );

    tardigradeStressTools::linearViscoelasticity( -1.1, {  0.0 }, -1.1, { 0.0 }, previousRateModifier, previousRateModifier, ISVs, params, 0.5, previousPK2MeanStressAnswer, previousUpdatedVolumetricISVsAnswer );

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters, ISVlb, ISVub, 0.5 );

    BOOST_TEST( PK2MeanStressAnswer[ 0 ] == *R.get_PK2MeanStress( ) );

    BOOST_TEST( updatedVolumetricISVsAnswer == *R.get_volumetricViscoelasticStateVariables( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousPK2MeanStressAnswer[ 0 ] == *R.get_previousPK2MeanStress( ) );

    floatVector strain;
    tardigradeConstitutiveTools::computeGreenLagrangeStrain( deformationGradient / std::pow( 1.1, 1/3.), strain );

    floatVector previousStrain;
    tardigradeConstitutiveTools::computeGreenLagrangeStrain( previousDeformationGradient, previousStrain );

    tardigradeConstitutiveTools::WLF( temperature, { parameters[ 7 ], parameters[ 8 ], 293.15 }, currentRateModifier );
    currentRateModifier = 1. / currentRateModifier;

    tardigradeConstitutiveTools::WLF( previousTemperature, { parameters[ 7 ], parameters[ 8 ], 293.15 }, previousRateModifier );
    previousRateModifier = 1. / previousRateModifier;

    ISVs = floatVector( previousStateVariables.begin( ) + ISVlb + parameters[ 0 ],
                        previousStateVariables.end( ) );
    params = { 2 * 56.7, 0.1, 1.0, 10.0, 2 * 12.3, 2 * 13.4, 2 * 14.5 };
    floatVector PK2IsochoricStressAnswer;
    floatVector updatedIsochoricISVsAnswer;
    floatVector previousPK2IsochoricStressAnswer;
    floatVector previousUpdatedIsochoricISVsAnswer;

    tardigradeStressTools::linearViscoelasticity(  1.1,         strain, -1.1, previousStrain,  currentRateModifier, previousRateModifier, ISVs, params, 0.5,        PK2IsochoricStressAnswer, updatedIsochoricISVsAnswer );

    tardigradeStressTools::linearViscoelasticity( -1.1, previousStrain, -1.1, previousStrain, previousRateModifier, previousRateModifier, ISVs, params, 0.5, previousPK2IsochoricStressAnswer, previousUpdatedIsochoricISVsAnswer );

    R.get_PK2IsochoricStress( );
    BOOST_TEST( PK2IsochoricStressAnswer == *R.get_PK2IsochoricStress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( updatedIsochoricISVsAnswer == *R.get_isochoricViscoelasticStateVariables( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousPK2IsochoricStressAnswer == *R.get_previousPK2IsochoricStress( ), CHECK_PER_ELEMENT );

    floatVector otherAnswer( previousStateVariables.begin( ), previousStateVariables.begin( ) + 2 );

    hydra.updateAdditionalStateVariables( );
    BOOST_TEST( tardigradeVectorTools::appendVectors( { otherAnswer, updatedVolumetricISVsAnswer, updatedIsochoricISVsAnswer } ) == *hydra.get_additionalStateVariables( ), CHECK_PER_ELEMENT );

}
