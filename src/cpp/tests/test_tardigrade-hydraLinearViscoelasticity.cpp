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
    
                    }
    
            };
    
        }

    }

}

BOOST_AUTO_TEST_CASE( test_residual_runBasicGetTests ){

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

    unsigned int ISVub = 10;

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
