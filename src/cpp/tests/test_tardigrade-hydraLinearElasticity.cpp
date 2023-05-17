/**
  * \file test_tardigrade-hydraLinearElasticity.cpp
  * 
  * Tests for tardigrade-hydraLinearElasticity
  */

#include<tardigrade-hydraLinearElasticity.h>
#include<constitutive_tools.h>

#define BOOST_TEST_MODULE test_tardigrade-hydraLinearElasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef errorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::linearElasticity::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::linearElasticity::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::linearElasticity::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

namespace tardigradeHydra{

    namespace linearElasticity{

        namespace unit_test{
    
            class residualTester{
    
                public:
    
                    static void runBasicGetTests( tardigradeHydra::linearElasticity::residual &R ){
    
                        BOOST_CHECK( &R._lambda == R.getLambda( ) );
    
                        BOOST_CHECK( &R._mu == R.getMu( ) );
    
                        BOOST_CHECK( &R._Ee.second == R.getEe( ) );

                        BOOST_CHECK( &R._dEedFe.second == R.getdEedFe( ) );

                        BOOST_CHECK( &R._PK2Stress.second == R.getPK2Stress( ) );
    
                        BOOST_CHECK( &R._dPK2dEe.second == R.getdPK2dEe( ) );
    
                        BOOST_CHECK( &R._dCauchyStressdF.second == R.getdCauchyStressdF( ) );
    
                        BOOST_CHECK( &R._dCauchyStressdFn.second == R.getdCauchyStressdFn( ) );
    
                    }
    
            };
    
        }

    }

}

BOOST_AUTO_TEST_CASE( test_residual_runBasicGetTests ){

    class hydraMock : public tardigradeHydra::hydraBase{

    };

    hydraMock hydra;

    floatVector parameters = { 123.4, 56.7 };

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

}

BOOST_AUTO_TEST_CASE( test_residual_decomposeParameterVector ){

    class hydraMock : public tardigradeHydra::hydraBase{

    };

    hydraMock hydra;

    floatVector parameters = { 123.4, 56.7 };

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    BOOST_CHECK( vectorTools::fuzzyEquals( parameters[ 0 ], *R.getLambda( ) ) );

    BOOST_CHECK( vectorTools::fuzzyEquals( parameters[ 1 ], *R.getMu( ) ) );

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

    BOOST_CHECK( vectorTools::fuzzyEquals( EeAnswer, *R.getEe( ) ) );

}
