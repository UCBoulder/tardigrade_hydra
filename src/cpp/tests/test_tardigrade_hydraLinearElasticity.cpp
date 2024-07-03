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

    namespace linearElasticity{

        namespace unit_test{
    
            class residualTester{
    
                public:
    
                    static void runBasicGetTests( tardigradeHydra::linearElasticity::residual &R ){
    
                        BOOST_CHECK( &R._lambda == R.getLambda( ) );
    
                        BOOST_CHECK( &R._mu == R.getMu( ) );
    
                        BOOST_CHECK( &R._Ee.second == R.get_Ee( ) );

                        BOOST_CHECK( &R._dEedFe.second == R.get_dEedFe( ) );

                        BOOST_CHECK( &R._previousEe.second == R.get_previousEe( ) );

                        BOOST_CHECK( &R._previousdEedFe.second == R.get_previousdEedFe( ) );

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

BOOST_AUTO_TEST_CASE( test_residual_runBasicGetTests, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                                      { }, { },
                                      previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    tardigradeHydra::linearElasticity::unit_test::residualTester::runBasicGetTests( R );

}

BOOST_AUTO_TEST_CASE( test_residual_decomposeParameterVector, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraMock : public tardigradeHydra::hydraBase{

    };

    hydraMock hydra;

    floatVector parameters = { 123.4, 56.7 };

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    BOOST_TEST( parameters[ 0 ] == *R.getLambda( ) );

    BOOST_TEST( parameters[ 1 ] == *R.getMu( ) );

}

BOOST_AUTO_TEST_CASE( test_residual_setFe, * boost::unit_test::tolerance( 1e-5 ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::linearElasticity::residual elasticity;

            tardigradeHydra::residualBase remainder;

            unsigned int elasticitySize = 9;

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            virtual void setResidualClasses( ){

                elasticity = tardigradeHydra::linearElasticity::residual( this, elasticitySize, *getParameters( ) );

                remainder = tardigradeHydra::residualBase( this, 9 );

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

    floatVector deformationGradient = { 9.92294371e-01, -1.32912834e-02,  3.57093896e-02,
                                        2.70156033e-02,  9.77585948e-01,  4.76270457e-04,
                                        4.42875587e-02, -4.95172679e-02,  9.54139963e-01 };

    floatVector previousDeformationGradient = { 1.02150915, -0.01813721, -0.01347062,
                                                0.0406867 ,  1.0450375 ,  0.0038319 ,
                                                0.07107588,  0.0013805 ,  0.98514977 };

    floatVector previousStateVariables = { 0.00315514,  0.00318276,  0.0134401 ,
                                           0.03494318,  0.02244553,  0.01110235,
                                           0.02224434, -0.01770411, -0.01382113 };

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    floatVector FeAnswer = { 0.98921175, -0.0156822 ,  0.02290497,
                            -0.00614278,  0.95596779, -0.01019557,
                             0.02379954, -0.03175083,  0.96754518 };

    floatVector previousFeAnswer = { 1.01964692, -0.02138607, -0.02731485,
                                     0.00513148,  1.0219469 , -0.00768935,
                                     0.04807642,  0.01848297,  0.99809319 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    tardigradeHydra::linearElasticity::residual RJ( &hydra, 9, parameters );

    RJ.get_dFedF( );

    BOOST_TEST( FeAnswer == *R.get_Fe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousFeAnswer == *R.get_previousFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( FeAnswer == *RJ.get_Fe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousFeAnswer == *RJ.get_previousFe( ), CHECK_PER_ELEMENT );

    //Test the Jacobians
    floatType eps = 1e-6;
    floatMatrix jac( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );
    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );
        
        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );
        
        floatVector vp = *Rp.get_Fe( );

        floatVector vm = *Rm.get_Fe( );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            jac[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( jac ) == *R.get_dFedF( ), CHECK_PER_ELEMENT );

    jac = floatMatrix( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );
    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );
        
        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );
        
        floatVector vp = *Rp.get_previousFe( );

        floatVector vm = *Rm.get_previousFe( );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            jac[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( jac ) == *R.get_previousdFedF( ), CHECK_PER_ELEMENT );

    jac = floatMatrix( deformationGradient.size( ), floatVector( previousStateVariables.size( ), 0 ) );
    for ( unsigned int i = 0; i < previousStateVariables.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );
        
        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );
        
        floatVector vp = *Rp.get_Fe( );

        floatVector vm = *Rm.get_Fe( );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            jac[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( jac ) == *R.get_dFedFn( ), CHECK_PER_ELEMENT );

    jac = floatMatrix( deformationGradient.size( ), floatVector( previousStateVariables.size( ), 0 ) );
    for ( unsigned int i = 0; i < previousStateVariables.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );
    
        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );
        
        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );
        
        floatVector vp = *Rp.get_previousFe( );

        floatVector vm = *Rm.get_previousFe( );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            jac[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( jac ) == *R.get_previousdFedFn( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setEe, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector previousEeAnswer = { -0.35589238, -0.06319833, -0.19137706,
                                     -0.06319833,  0.13857639,  0.22404018,
                                     -0.19137706,  0.22404018, -0.16361939 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    BOOST_TEST( EeAnswer == *R.get_Ee( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousEeAnswer == *R.get_previousEe( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setdEedFe, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector previousEeAnswer;

    floatMatrix dEedFeAnswer;

    floatMatrix previousdEedFeAnswer;

    TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeGreenLagrangeStrain( deformationGradient, EeAnswer, dEedFeAnswer ) );

    TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeGreenLagrangeStrain( previousDeformationGradient, previousEeAnswer, previousdEedFeAnswer ) );

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    BOOST_TEST( tardigradeVectorTools::appendVectors( dEedFeAnswer ) == *R.get_dEedFe( ), CHECK_PER_ELEMENT );

    BOOST_TEST( tardigradeVectorTools::appendVectors( previousdEedFeAnswer ) == *R.get_previousdEedFe( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setPK2Stress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearElasticity::residual {

        public:

            using tardigradeHydra::linearElasticity::residual::residual;

            void set_Ee( floatVector &Ee ){ tardigradeHydra::linearElasticity::residual::set_Ee( Ee ); }

            void set_previousEe( floatVector &previousEe ){ tardigradeHydra::linearElasticity::residual::set_previousEe( previousEe ); }

            floatVector Ee = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousEe = { -1, -2, -3, -4, -5, -6, -7, -8, -9 };

        private:

            virtual void setEe( ) override {

                set_Ee( Ee );

            }
            virtual void setPreviousEe( ) override {

                set_previousEe( previousEe );

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters );

    BOOST_TEST( PK2StressAnswer == *R.get_PK2Stress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( -PK2StressAnswer == *R.get_previousPK2Stress( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setdPK2StressdEe, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearElasticity::residual {

        public:

            using tardigradeHydra::linearElasticity::residual::residual;

            void set_Ee( floatVector &Ee ){ tardigradeHydra::linearElasticity::residual::set_Ee( Ee ); }

            void set_previousEe( floatVector &previousEe ){ tardigradeHydra::linearElasticity::residual::set_previousEe( previousEe ); }

            floatVector Ee = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousEe = { -1, -2, -3, -4, -5, -6, -7, -8, -9 };

        private:

            virtual void setEe( ) override {

                set_Ee( Ee );

            }
            virtual void setPreviousEe( ) override {

                set_previousEe( previousEe );

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
                         { }, { },
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

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *R.get_dPK2StressdEe( ), CHECK_PER_ELEMENT );

    gradient = floatMatrix( R.previousEe.size( ), floatVector( R.previousEe.size( ), 0 ) );

    for ( unsigned int i = 0; i < R.previousEe.size( ); i++ ){

        floatVector delta( R.previousEe.size( ), 0 );

        delta[ i ] = eps * std::fabs( R.previousEe[ i ] ) + eps;

        residualMock Rp( &hydra, 9, parameters );

        residualMock Rm( &hydra, 9, parameters );

        Rp.previousEe += delta;

        Rm.previousEe -= delta;

        for ( unsigned int j = 0; j < R.Ee.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.get_previousPK2Stress( ) )[ j ] - ( *Rm.get_previousPK2Stress( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *R.get_previousdPK2StressdEe( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setdPK2StressdFe, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters );

    floatType eps = 1e-6;

    floatMatrix gradient( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );


        residualMock Rp( &hydrap, 9, parameters );

        residualMock Rm( &hydram, 9, parameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.get_PK2Stress( ) )[ j ] - ( *Rm.get_PK2Stress( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *R.get_dPK2StressdFe( ), CHECK_PER_ELEMENT );

    gradient = floatMatrix( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );


        residualMock Rp( &hydrap, 9, parameters );

        residualMock Rm( &hydram, 9, parameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.get_previousPK2Stress( ) )[ j ] - ( *Rm.get_previousPK2Stress( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *R.get_previousdPK2StressdFe( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setCauchyStress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearElasticity::residual {

        public:

            using tardigradeHydra::linearElasticity::residual::residual;

            void set_PK2Stress( floatVector &PK2 ){ tardigradeHydra::linearElasticity::residual::set_PK2Stress( PK2 ); }

            void set_previousPK2Stress( floatVector &previousPK2 ){ tardigradeHydra::linearElasticity::residual::set_previousPK2Stress( previousPK2 ); }

            floatVector PK2Stress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousPK2Stress = { -1, -2, -3, -4, -5, -6, -7, -8, -9 };

        private:

            virtual void setPK2Stress( ) override {

                set_PK2Stress( PK2Stress );

            }

            virtual void setPreviousPK2Stress( ) override {

                set_previousPK2Stress( previousPK2Stress );

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

    floatVector previousCauchyStressAnswer = {  0.88800978,  10.67302743,   6.73691146,
                                                1.94755853, -40.89370509, -34.2541294 ,
                                                0.08381629, -48.02194605, -36.74857461 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters );

    BOOST_TEST( cauchyStressAnswer == *R.get_cauchyStress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousCauchyStressAnswer == *R.get_previousCauchyStress( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setStress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class residualMock : public tardigradeHydra::linearElasticity::residual {

        public:

            using tardigradeHydra::linearElasticity::residual::residual;

            void set_PK2Stress( floatVector &PK2 ){ tardigradeHydra::linearElasticity::residual::set_PK2Stress( PK2 ); }

            void set_previousPK2Stress( floatVector &previousPK2 ){ tardigradeHydra::linearElasticity::residual::set_previousPK2Stress( previousPK2 ); }

            floatVector PK2Stress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousPK2Stress = { -1, -2, -3, -4, -5, -6, -7, -8, -9 };

        private:

            virtual void setPK2Stress( ) override {

                set_PK2Stress( PK2Stress );

            }

            virtual void setPreviousPK2Stress( ) override {

                set_previousPK2Stress( previousPK2Stress );

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

    floatVector previousCauchyStressAnswer = {  0.88800978,  10.67302743,   6.73691146,
                                                1.94755853, -40.89370509, -34.2541294 ,
                                                0.08381629, -48.02194605, -36.74857461 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters );

    BOOST_TEST( cauchyStressAnswer == *R.getStress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousCauchyStressAnswer == *R.getPreviousStress( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setdCauchyStressdPK2Stress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, parameters );

    floatMatrix gradient( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < R.PK2Stress.size( ); i++ ){

        floatVector delta( R.PK2Stress.size( ), 0 );

        delta[ i ] = eps * std::fabs( R.PK2Stress[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, parameters );

        residualMock Rm( &hydram, 9, parameters );

        Rp.PK2Stress += delta;

        Rm.PK2Stress -= delta;

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.get_cauchyStress( ) )[ j ] - ( *Rm.get_cauchyStress( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *R.get_dCauchyStressdPK2Stress( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setdCauchyStressdF, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    floatMatrix gradient( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

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

        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );

        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.get_cauchyStress( ) )[ j ] - ( *Rm.get_cauchyStress( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *R.get_dCauchyStressdF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );

        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.get_cauchyStress( ) )[ j ] - ( *Rm.get_cauchyStress( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *R.get_dCauchyStressdPreviousF( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setdCauchyStressdFn, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    floatMatrix gradient( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    BOOST_CHECK( R.get_dCauchyStressdFn( )->size( ) == 0 );

}

BOOST_AUTO_TEST_CASE( test_residual_setdCauchyStressdFn2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    floatMatrix gradient( deformationGradient.size( ), floatVector( 9 * ( numConfigurations - 1 ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 9 * ( numConfigurations - 1 ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );

        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.get_cauchyStress( ) )[ j ] - ( *Rm.get_cauchyStress( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *R.get_dCauchyStressdFn( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setdCauchyStressdFn3, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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

    floatVector unknownVector = { 0.0, 0.0, 0.0,  0.0,  0.0,  0.0, 0.0, 0.0,  0.0,
                                  1.1, 0.2, 0.3,  0.2,  1.4, -0.2, 0.3, 0.1,  0.9,
                                  0.5, 0.4, 0.2, -0.2,  0.9,  0.1, 0.1, 0.3,  1.4};

    floatVector parameters = { 123.4, 56.7 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    floatMatrix gradient( deformationGradient.size( ), floatVector( 9 * ( numConfigurations - 1 ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 9 * ( numConfigurations - 1 ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );

        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.get_cauchyStress( ) )[ j ] - ( *Rm.get_cauchyStress( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *R.get_dCauchyStressdPreviousFn( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setResidual, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    BOOST_TEST( residualAnswer == *R.getResidual( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setResidual2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 9, parameters );

    BOOST_TEST( residualAnswer == *R.getResidual( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    floatMatrix gradient( 9, floatVector( hydra.getNumUnknowns( ) ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < ( 9 * numConfigurations + numNonLinearSolveStateVariables ); i++ ){

        floatVector delta( hydra.getNumUnknowns( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );

        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.getResidual( ) )[ j ] - ( *Rm.getResidual( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *R.getJacobian( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setdRdT, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    floatVector gradient( 9 );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector delta( 1, 0 );

        delta[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + delta[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        hydraBaseMock hydram( time, deltaTime, temperature - delta[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );

        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );

        gradient = ( *Rp.getResidual( ) - *Rm.getResidual( ) ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( gradient == *R.getdRdT( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_setdRdF, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

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
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::linearElasticity::residual R( &hydra, 9, parameters );

    floatMatrix gradient( deformationGradient.size( ), floatVector( deformationGradient.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::linearElasticity::residual Rp( &hydrap, 9, parameters );

        tardigradeHydra::linearElasticity::residual Rm( &hydram, 9, parameters );

        for ( unsigned int j = 0; j < deformationGradient.size( ); j++ ){

            gradient[ j ][ i ] = ( ( *Rp.getResidual( ) )[ j ] - ( *Rm.getResidual( ) )[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( gradient ) == *R.getdRdF( ), CHECK_PER_ELEMENT );

}
