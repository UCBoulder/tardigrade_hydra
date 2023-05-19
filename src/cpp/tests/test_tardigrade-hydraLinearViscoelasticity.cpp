/**
  * \file test_tardigrade-hydraLinearElasticity.cpp
  * 
  * Tests for tardigrade-hydraLinearElasticity
  */

#include<tardigrade-hydraLinearViscoelasticity.h>

#define BOOST_TEST_MODULE test_tardigrade-hydraLinearElasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef errorTools::Node errorNode; //!< Redefinition for the error node
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

                        BOOST_CHECK( &R._Je.second == R.getJe( ) );

                        BOOST_CHECK( &R._Fehat.second == R.getFehat( ) );
    
                        BOOST_CHECK( &R._previousJe.second == R.getPreviousJe( ) );

                        BOOST_CHECK( &R._previousFehat.second == R.getPreviousFehat( ) );

                        BOOST_CHECK( &R._numStateVariables == R.getNumStateVariables( ) );

                    }
    
            };
    
        }

    }

}

BOOST_AUTO_TEST_CASE( test_residual_runBasicGetTests_and_decomposeParameters ){

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

    floatVector parameters = { 2, 3, 123.4, 56.7, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    unsigned int numVolumetricViscousTermsAnswer = 2;

    unsigned int numIsochoricViscousTermsAnswer = 3;

    floatType KinfAnswer = 123.4;

    floatType GinfAnswer = 56.7;

    floatVector Ks = { 23.4, 25.6 };

    floatVector KTaus = { 0.1, 0.2 };

    floatVector Gs = { 12.3, 13.4, 14.5 };

    floatVector GTaus = { 0.01, 10.0, 100.0 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual R( &hydra, 9, parameters, ISVlb, ISVub );

    BOOST_CHECK_NO_THROW( tardigradeHydra::linearViscoelasticity::unit_test::residualTester::runBasicGetTests( R ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( ISVlb, *R.getViscoelasticISVLowerIndex( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( ISVub, *R.getViscoelasticISVUpperIndex( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( numVolumetricViscousTermsAnswer, *R.getNumVolumetricViscousTerms( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( numIsochoricViscousTermsAnswer, *R.getNumIsochoricViscousTerms( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( KinfAnswer, *R.getKinf( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( GinfAnswer, *R.getGinf( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( Ks, *R.getVolumetricModuli( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( Gs, *R.getIsochoricModuli( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( KTaus, *R.getVolumetricTaus( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( GTaus, *R.getIsochoricTaus( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_decomposeElasticDeformation ){

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

    floatVector parameters = { 2, 3, 123.4, 56.7, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    floatType JeAnswer = 0.2871738879785901;

    floatVector FehatAnswer = { 0.59558368, -0.64830483, -0.82803222,
                                0.15555742,  0.66530605, -0.23309781,
                                1.45740571,  0.56029945, -0.05780372 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual R1( &hydra, 9, parameters, ISVlb, ISVub );

    tardigradeHydra::linearViscoelasticity::residual R2( &hydra, 9, parameters, ISVlb, ISVub );

    ERROR_TOOLS_CATCH( vectorTools::fuzzyEquals( JeAnswer, *R1.getJe( ) ) );

    ERROR_TOOLS_CATCH( vectorTools::fuzzyEquals( FehatAnswer, *R2.getFehat( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_decomposePreviousElasticDeformation ){

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

    floatVector parameters = { 2, 3, 123.4, 56.7, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual R1( &hydra, 9, parameters, ISVlb, ISVub );

    tardigradeHydra::linearViscoelasticity::residual R2( &hydra, 9, parameters, ISVlb, ISVub );

    ERROR_TOOLS_CATCH( vectorTools::fuzzyEquals( JeAnswer, *R1.getPreviousJe( ) ) );

    ERROR_TOOLS_CATCH( vectorTools::fuzzyEquals( FehatAnswer, *R2.getPreviousFehat( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_gradientsOfDecomposedElasticDeformationGradient ){

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

    floatVector parameters = { 2, 3, 123.4, 56.7, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual R( &hydra, 9, parameters, ISVlb, ISVub );

    floatVector gradientJe( deformationGradient.size( ), 0 );

    floatMatrix gradientFehat( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
        tardigradeHydra::linearViscoelasticity::residual Rp( &hydrap, 9, parameters, ISVlb, ISVub );

        tardigradeHydra::linearViscoelasticity::residual Rm( &hydram, 9, parameters, ISVlb, ISVub );

        for ( unsigned int j = 0; j < 1; j++ ){

            gradientJe[ i ] = ( ( *Rp.getJe( ) ) - ( *Rm.getJe( ) ) ) / ( 2 * delta[ i ] );

        }

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradientFehat[ j ][ i ] = ( ( *Rp.getFehat( ) )[ j ] - ( *Rm.getFehat( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( vectorTools::fuzzyEquals( gradientJe, *R.getdJedFe( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( gradientFehat, *R.getdFehatdFe( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_decomposeStateVariableVector ){

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

    floatVector parameters = { 2, 3, 123.4, 56.7, 23.4, 25.6, 0.1, 0.2, 12.3, 13.4, 14.5, 0.01, 10.0, 100.0 };

    unsigned int numConfigurations = 1;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    unsigned int ISVlb = 2;

    unsigned int ISVub = 31;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearViscoelasticity::residual R( &hydra, 9, parameters, ISVlb, ISVub );

    floatVector volumetricStateVariables, isochoricStateVariables;

    BOOST_CHECK_NO_THROW( R.decomposeStateVariableVector( volumetricStateVariables,
                                                          isochoricStateVariables ) );

    ERROR_TOOLS_CATCH( vectorTools::fuzzyEquals( volumetricStateVariables,
                                                 volumetricStateVariablesAnswer ) );

    ERROR_TOOLS_CATCH( vectorTools::fuzzyEquals( isochoricStateVariables,
                                                 isochoricStateVariablesAnswer ) );

}
