/**
  * \file test_tardigrade_hydraPerzynaViscodamage.cpp
  *
  * Tests for tardigrade_hydraPerzynaViscodamage
  */

#include<tardigrade_hydraPerzynaViscodamage.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_tardigrade_hydra
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

namespace tardigradeHydra{

    namespace unit_test{

        class hydraBaseTester{

            public:

                static void updateUnknownVector( tardigradeHydra::hydraBase &hydra, const floatVector &value ){

                    BOOST_CHECK_NO_THROW( hydra.updateUnknownVector( value ) );

                }

        };

    }

}

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer)
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

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

BOOST_AUTO_TEST_CASE( test_setStateVariableEvolutionRates, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscodamage::residual{

        public:

            residualMock( ) : tardigradeHydra::perzynaViscodamage::residual( ){ }

            residualMock( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations,
                          const unsigned int &damageConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices,
                          const floatVector &parameters ) : tardigradeHydra::perzynaViscodamage::residual( hydra, numEquations, damageConfigurationIndex, stateVariableIndices, parameters ){ }

            floatVector stateVariableEvolutionRates = { 0.78 };

            floatVector previousStateVariableEvolutionRates = { 0.89 };

            floatType plasticMultiplier = 1.2;

            floatType previousPlasticMultiplier = 3.4;

            floatVector previousdXidotdSigma = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousdXidotdF     = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

            floatVector previousdXidotdSubFs = { 19, 20, 21, 22, 23, 24, 25, 26, 27,
                                                 28, 29, 30, 31, 32, 33, 34, 35, 36 };

            floatVector previousdXidotdT = { 37 };

            floatVector previousdXidotdXi = { 38 };

            floatVector dXidotdSigma = { 39, 40, 41, 42, 43, 44, 45, 46, 47 };

            floatVector dXidotdF     = { 48, 49, 50, 51, 52, 53, 54, 55, 56 };

            floatVector dXidotdSubFs = { 57, 58, 59, 60, 61, 62, 63, 64, 65,
                                         66, 67, 68, 69, 70, 71, 72, 73, 74 };

            floatVector dXidotdT = { 75 };

            floatVector dXidotdXi = { 76 };

            floatVector previousdGammadSigma = { -1, -2, -3, -4, -5, -6, -7, -8, -9 };

            floatVector previousdGammadF = { -10, -11, -12, -13, -14, -15, -16, -17, -18 };

            floatVector previousdGammadSubFs = { -19, -20, -21, -22, -23, -24, -25, -26, -27,
                                                 -28, -29, -30, -31, -32, -33, -34, -35, -36 };

            floatType   previousdGammadT = -37;

            floatVector previousdGammadXi = { -38 };

            floatVector dGammadSigma = { -39, -40, -41, -42, -43, -44, -45, -46, -47 };

            floatVector dGammadF = { -48, -49, -50, -51, -52, -53, -54, -55, -56 };

            floatVector dGammadSubFs = { -57, -58, -59, -60, -61, -62, -63, -64, -65,
                                         -66, -67, -68, -69, -70, -71, -72, -73, -74 };

            floatType   dGammadT = -75;

            floatVector dGammadXi = { -76 };

        protected:

            using tardigradeHydra::perzynaViscodamage::residual::setHardeningFunction;

            using tardigradeHydra::perzynaViscodamage::residual::setStateVariables;
            using tardigradeHydra::perzynaViscodamage::residual::setdStateVariableEvolutionRatesdCauchyStress;
            using tardigradeHydra::perzynaViscodamage::residual::setdStateVariableEvolutionRatesdF;
            using tardigradeHydra::perzynaViscodamage::residual::setdStateVariableEvolutionRatesdSubFs;
            using tardigradeHydra::perzynaViscodamage::residual::setdStateVariableEvolutionRatesdT;
            using tardigradeHydra::perzynaViscodamage::residual::setdStateVariableEvolutionRatesdStateVariables;

            using tardigradeHydra::perzynaViscodamage::residual::setPlasticMultiplier;
            using tardigradeHydra::perzynaViscodamage::residual::setdPlasticMultiplierdCauchyStress;
            using tardigradeHydra::perzynaViscodamage::residual::setdPlasticMultiplierdF;
            using tardigradeHydra::perzynaViscodamage::residual::setdPlasticMultiplierdSubFs;
            using tardigradeHydra::perzynaViscodamage::residual::setdPlasticMultiplierdT;
            using tardigradeHydra::perzynaViscodamage::residual::setdPlasticMultiplierdStateVariables;

            virtual void setHardeningFunction( ) override{

                set_hardeningFunction( stateVariableEvolutionRates / plasticMultiplier );

            }

            virtual void setPreviousHardeningFunction( ) override{

                set_previousHardeningFunction( previousStateVariableEvolutionRates / previousPlasticMultiplier );

            }

            virtual void setPlasticMultiplier( ) override{

                set_plasticMultiplier( plasticMultiplier );

            }

            virtual void setPreviousPlasticMultiplier( ) override{

                set_previousPlasticMultiplier( previousPlasticMultiplier );

            }

            virtual void setdPreviousStateVariableEvolutionRatesdPreviousCauchyStress( ) override{

                set_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( previousdXidotdSigma );

            }

            virtual void setdPreviousStateVariableEvolutionRatesdPreviousF( ) override{

                set_dPreviousStateVariableEvolutionRatesdPreviousF( previousdXidotdF );

            }

            virtual void setdPreviousStateVariableEvolutionRatesdPreviousSubFs( ) override{

                set_dPreviousStateVariableEvolutionRatesdPreviousSubFs( previousdXidotdSubFs );

            }

            virtual void setdPreviousStateVariableEvolutionRatesdPreviousT( ) override{

                set_dPreviousStateVariableEvolutionRatesdPreviousT( previousdXidotdT );

            }

            virtual void setdPreviousStateVariableEvolutionRatesdPreviousStateVariables( ) override{

                set_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( previousdXidotdXi );

            }

            virtual void setdStateVariableEvolutionRatesdCauchyStress( ) override{

                set_dStateVariableEvolutionRatesdCauchyStress( dXidotdSigma );

            }

            virtual void setdStateVariableEvolutionRatesdF( ) override{

                set_dStateVariableEvolutionRatesdF( dXidotdF );

            }

            virtual void setdStateVariableEvolutionRatesdSubFs( ) override{

                set_dStateVariableEvolutionRatesdSubFs( dXidotdSubFs );

            }

            virtual void setdStateVariableEvolutionRatesdT( ) override{

                set_dStateVariableEvolutionRatesdT( dXidotdT );

            }

            virtual void setdStateVariableEvolutionRatesdStateVariables( ) override{

                set_dStateVariableEvolutionRatesdCauchyStress( dXidotdXi );

            }

            virtual void setdPreviousPlasticMultiplierdPreviousCauchyStress( ) override{

                set_dPreviousPlasticMultiplierdPreviousCauchyStress( previousdGammadSigma );

            }

            virtual void setdPreviousPlasticMultiplierdPreviousF( ) override{

                set_dPreviousPlasticMultiplierdPreviousF( previousdGammadF );

            }

            virtual void setdPreviousPlasticMultiplierdPreviousSubFs( ) override{

                set_dPreviousPlasticMultiplierdPreviousSubFs( previousdGammadSubFs );

            }

            virtual void setdPreviousPlasticMultiplierdPreviousT( ) override{

                set_dPreviousPlasticMultiplierdPreviousT( previousdGammadT );

            }

            virtual void setdPreviousPlasticMultiplierdPreviousStateVariables( ) override{

                set_dPreviousPlasticMultiplierdPreviousCauchyStress( previousdGammadXi );

            }

            virtual void setdPlasticMultiplierdCauchyStress( ) override{

                set_dPlasticMultiplierdCauchyStress( dGammadSigma );

            }

            virtual void setdPlasticMultiplierdF( ) override{

                set_dPlasticMultiplierdF( dGammadF );

            }

            virtual void setdPlasticMultiplierdSubFs( ) override{

                set_dPlasticMultiplierdSubFs( dGammadSubFs );

            }

            virtual void setdPlasticMultiplierdT( ) override{

                set_dPlasticMultiplierdT( dGammadT );

            }

            virtual void setdPlasticMultiplierdStateVariables( ) override{

                set_dPlasticMultiplierdCauchyStress( dGammadXi );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector viscoDamageParameters = { 2.0,
                                                  1e1, 1e2,
                                                  10, 200, 293.15,
                                                  1, 0.34,
                                                  0.12,
                                                  13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            residualMock damage;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                damage = residualMock( this, 11, 1, stateVariableIndices, viscoDamageParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &damage;

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

    floatVector previousStateVariables = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0.01, 0.02, 0.03, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 6, 7, 8 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R_ngrad( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    residualMock R_grad1( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    residualMock R_grad2( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    floatVector currentAnswer  = { 0.78, 1.2 };

    floatVector previousAnswer = { 0.89, 3.4 };

    R_grad1.get_dStateVariableEvolutionRatesdCauchyStress( );

    R_grad2.get_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( );

    BOOST_TEST( currentAnswer == *R_ngrad.get_stateVariableEvolutionRates( ), CHECK_PER_ELEMENT );

    BOOST_TEST( currentAnswer == *R_grad1.get_stateVariableEvolutionRates( ), CHECK_PER_ELEMENT );

    BOOST_TEST( currentAnswer == *R_grad2.get_stateVariableEvolutionRates( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *R_ngrad.get_previousStateVariableEvolutionRates( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *R_grad1.get_previousStateVariableEvolutionRates( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previousAnswer == *R_grad2.get_previousStateVariableEvolutionRates( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setStateVariableEvolutionRates2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoDamageParameters = { 2.0,
                                                  1e1, 1e2,
                                                  10, 200, 293.15,
                                                  1, 0.34,
                                                  0.12,
                                                  13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            tardigradeHydra::perzynaViscodamage::residual damage;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                damage = tardigradeHydra::perzynaViscodamage::residual( this, 11, 1, stateVariableIndices, viscoDamageParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &damage;

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

    floatVector previousStateVariables = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0.01, 0.02, 0.03, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 3, 1, 4, 1, 5, 1, 1, 1, 8,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  0.04, 0.05, 0.06, 0.07, 0.08 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::perzynaViscodamage::residual R( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    // Test the computation of the Jacobians
    unsigned int nvals = 2;

    floatMatrix dStateVariableEvolutionRatesdCauchyStress( nvals, floatVector( 9, 0 ) );

    floatMatrix dStateVariableEvolutionRatesdF( nvals, floatVector( 9, 0 ) );

    floatMatrix dStateVariableEvolutionRatesdSubFs( nvals, floatVector( 18, 0 ) );

    floatVector dStateVariableEvolutionRatesdT( nvals, 0 );

    floatMatrix dStateVariableEvolutionRatesdStateVariables( nvals, floatVector( 2, 0 ) );

    floatMatrix dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( nvals, floatVector( 9, 0 ) );

    floatMatrix dPreviousStateVariableEvolutionRatesdPreviousF( nvals, floatVector( 9, 0 ) );

    floatMatrix dPreviousStateVariableEvolutionRatesdPreviousSubFs( nvals, floatVector( 18, 0 ) );

    floatVector dPreviousStateVariableEvolutionRatesdPreviousT( 2, 0 );

    floatMatrix dPreviousStateVariableEvolutionRatesdPreviousStateVariables( nvals, floatVector( 2, 0 ) );

    floatType eps = 1e-6;
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] += eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_stateVariableEvolutionRates( );

        floatVector vm = *Rm.get_stateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dStateVariableEvolutionRatesdCauchyStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dStateVariableEvolutionRatesdCauchyStress ) == *R.get_dStateVariableEvolutionRatesdCauchyStress( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] += eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_stateVariableEvolutionRates( );

        floatVector vm = *Rm.get_stateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dStateVariableEvolutionRatesdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dStateVariableEvolutionRatesdF ) == *R.get_dStateVariableEvolutionRatesdF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 18; i++ ){

        unsigned int offset = 9;

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + offset ] += eps * std::fabs( unknownVector[ i + offset ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_stateVariableEvolutionRates( );

        floatVector vm = *Rm.get_stateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dStateVariableEvolutionRatesdSubFs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + offset ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dStateVariableEvolutionRatesdSubFs ) == *R.get_dStateVariableEvolutionRatesdSubFs( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector delta( 1, 0 );

        delta[ i ] += eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + delta[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature - delta[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_stateVariableEvolutionRates( );

        floatVector vm = *Rm.get_stateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dStateVariableEvolutionRatesdT[ j ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( dStateVariableEvolutionRatesdT == *R.get_dStateVariableEvolutionRatesdT( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 2; i++ ){

        unsigned int offset = 9 + 18 + 1;

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + offset ] += eps * std::fabs( unknownVector[ i + offset ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_stateVariableEvolutionRates( );

        floatVector vm = *Rm.get_stateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dStateVariableEvolutionRatesdStateVariables[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + offset ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dStateVariableEvolutionRatesdStateVariables ) == *R.get_dStateVariableEvolutionRatesdStateVariables( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( hydra.elasticity.previousCauchyStress.size( ), 0 );

        delta[ i ] += eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap._local_deltaPreviousCauchyStress =  delta;

        hydram._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_previousStateVariableEvolutionRates( );

        floatVector vm = *Rm.get_previousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousCauchyStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousStateVariableEvolutionRatesdPreviousCauchyStress ) == *R.get_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( previousDeformationGradient.size( ), 0 );

        delta[ i ] += eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_previousStateVariableEvolutionRates( );

        floatVector vm = *Rm.get_previousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousStateVariableEvolutionRatesdPreviousF ) == *R.get_dPreviousStateVariableEvolutionRatesdPreviousF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 18; i++ ){

        unsigned int offset = 0;

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + offset ] += eps * std::fabs( previousStateVariables[ i + offset ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_previousStateVariableEvolutionRates( );

        floatVector vm = *Rm.get_previousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousSubFs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + offset ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousStateVariableEvolutionRatesdPreviousSubFs ) == *R.get_dPreviousStateVariableEvolutionRatesdPreviousSubFs( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector delta( 1, 0 );

        delta[ i ] += eps * std::fabs( previousTemperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature + delta[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature - delta[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_previousStateVariableEvolutionRates( );

        floatVector vm = *Rm.get_previousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousT[ j ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( dPreviousStateVariableEvolutionRatesdPreviousT == *R.get_dPreviousStateVariableEvolutionRatesdPreviousT( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 2; i++ ){

        unsigned int offset = 18 + 1;

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + offset ] += eps * std::fabs( previousStateVariables[ i + offset ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_previousStateVariableEvolutionRates( );

        floatVector vm = *Rm.get_previousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousStateVariables[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + offset ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousStateVariableEvolutionRatesdPreviousStateVariables ) == *R.get_dPreviousStateVariableEvolutionRatesdPreviousStateVariables( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setDamage, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscodamage::residual{

        public:

            residualMock( ) : tardigradeHydra::perzynaViscodamage::residual( ){ }

            residualMock( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations,
                          const unsigned int &damageConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices,
                          const floatVector &parameters ) : tardigradeHydra::perzynaViscodamage::residual( hydra, numEquations, damageConfigurationIndex, stateVariableIndices, parameters ){ }

            floatType plasticMultiplier = 1.2;

            floatType previousPlasticMultiplier = 3.4;

            floatVector previousStateVariables = { 1, 2 };

            floatVector stateVariables = { 3, 4 };

            floatVector dIsvsdCauchy = initializeVector( 2*9, 1 );

            floatVector dIsvsdF = initializeVector( 2*9, 1.2 );

            floatVector dIsvsdSubF = initializeVector( 2*2*9, 0.5 );

            floatVector dIsvsdT = initializeVector( 2, 5 );

            floatVector dIsvsdIsvs = initializeVector( 2*2, 0.3 );

            floatVector initializeVector( int nvals, floatType val_0 ){

                floatVector value( nvals, 0 );

                for ( int i = 0; i < nvals; i++ ){

                    value[ i ] = i + val_0;

                }

                return value;

            }

        protected:

            using tardigradeHydra::perzynaViscodamage::residual::setStateVariables;
            using tardigradeHydra::perzynaViscodamage::residual::setHardeningFunction;
            using tardigradeHydra::perzynaViscodamage::residual::setPlasticMultiplier;

            virtual void setPlasticMultiplier( ) override{

                set_plasticMultiplier( plasticMultiplier );

            }

            virtual void setPreviousPlasticMultiplier( ) override{

                set_previousPlasticMultiplier( previousPlasticMultiplier );

            }

            virtual void setStateVariables( ) override{

                set_stateVariables( stateVariables );

            }

            virtual void setPreviousStateVariables( ) override{

                set_previousStateVariables( previousStateVariables );

            }

            virtual void setdPlasticStateVariablesdPreviousCauchyStress( ) override{

                set_dPlasticStateVariablesdPreviousCauchyStress( dIsvsdCauchy );

            }

            virtual void setdPlasticStateVariablesdPreviousF( ) override{

                set_dPlasticStateVariablesdPreviousF( dIsvsdF );

            }

            virtual void setdPlasticStateVariablesdPreviousSubFs( ) override{

                set_dPlasticStateVariablesdPreviousSubFs( dIsvsdSubF );

            }

            virtual void setdPlasticStateVariablesdPreviousT( ) override{

                set_dPlasticStateVariablesdPreviousT( dIsvsdT );

            }

            virtual void setdPlasticStateVariablesdPreviousStateVariables( ) override{

                set_dPlasticStateVariablesdPreviousStateVariables( dIsvsdIsvs );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector viscoDamageParameters = { 2.0,
                                                  1e1, 1e2,
                                                  10, 200, 293.15,
                                                  1, 0.34,
                                                  0.12,
                                                  13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            residualMock damage;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                damage = residualMock( this, 11, 1, stateVariableIndices, viscoDamageParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &damage;

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

    floatVector previousStateVariables = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0.01, 0.02, 0.03, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 6, 7, 8 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R_ngrad( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    residualMock R_grad1( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    residualMock R_grad2( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    floatType damageAnswer = 7.06; 

    R_grad1.get_dDamagedCauchyStress( );

    R_grad2.get_dDamagedPreviousCauchyStress( );

    BOOST_TEST( damageAnswer == *R_ngrad.get_damage( ) );

    BOOST_TEST( damageAnswer == *R_grad1.get_damage( ) );

    BOOST_TEST( damageAnswer == *R_grad2.get_damage( ) );

}

BOOST_AUTO_TEST_CASE( test_setDamage2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoDamageParameters = { 2.0,
                                                  1e1, 1e2,
                                                  10, 200, 293.15,
                                                  1, 0.34,
                                                  0.12,
                                                  13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            tardigradeHydra::perzynaViscodamage::residual damage;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                damage = tardigradeHydra::perzynaViscodamage::residual( this, 11, 1, stateVariableIndices, viscoDamageParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &damage;

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

    floatVector previousStateVariables = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0.01, 0.02, 0.03, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 3, 1, 4, 1, 5, 1, 1, 1, 8,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  0.04, 0.05, 0.06, 0.07, 0.08 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::perzynaViscodamage::residual R( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    // Test the computation of the Jacobians
    unsigned int nvals = 1;

    floatVector dDamagedCauchyStress( 9, 0 );

    floatVector dDamagedF( 9, 0 );

    floatVector dDamagedSubFs( 18, 0 );

    floatType   dDamagedT;

    floatVector dDamagedStateVariables( 2, 0 );

    floatVector dDamagedPreviousCauchyStress( 9, 0 );

    floatVector dDamagedPreviousF( 9, 0 );

    floatVector dDamagedPreviousSubFs( 18, 0 );

    floatType   dDamagedPreviousT;

    floatVector dDamagedPreviousStateVariables( 2, 0 );

    floatType eps = 1e-6;
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] += eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedCauchyStress[ i ] = ( vp - vm ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( dDamagedCauchyStress == *R.get_dDamagedCauchyStress( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] += eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedF[ i ] = ( vp - vm ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( dDamagedF == *R.get_dDamagedF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 2*deformationGradient.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] += eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedSubFs[ i ] = ( vp - vm ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_TEST( dDamagedSubFs == *R.get_dDamagedSubFs( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector delta( 1, 0 );

        delta[ i ] += eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + delta[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature - delta[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedT = ( vp - vm ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( dDamagedT == *R.get_dDamagedT( ) );

    for ( unsigned int i = 0; i < 2; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 + 18 + 1 ] += eps * std::fabs( unknownVector[ 9 + 18 + 1 + i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedStateVariables[ i ] = ( vp - vm ) / ( 2 * delta[ i + 9 + 18 + 1 ] );

        }

    }

    BOOST_TEST( dDamagedStateVariables == *R.get_dDamagedStateVariables( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] += eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap._local_deltaPreviousCauchyStress = delta;

        hydram._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedPreviousCauchyStress[ i ] = ( vp - vm ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( dDamagedPreviousCauchyStress == *R.get_dDamagedPreviousCauchyStress( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] += eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedPreviousF[ i ] = ( vp - vm ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( dDamagedPreviousF == *R.get_dDamagedPreviousF( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] += eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedPreviousSubFs[ i ] = ( vp - vm ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( dDamagedPreviousSubFs == *R.get_dDamagedPreviousSubFs( ), CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector delta( 1, 0 );

        delta[ i ] += eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature + delta[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature - delta[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedPreviousT = ( vp - vm ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( dDamagedPreviousT == *R.get_dDamagedPreviousT( ) );

    for ( unsigned int i = 0; i < 2; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + 18 + 1 ] += eps * std::fabs( previousStateVariables[ 18 + 1 + i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedPreviousStateVariables[ i ] = ( vp - vm ) / ( 2 * delta[ i + 18 + 1 ] );

        }

    }

    BOOST_TEST( dDamagedPreviousStateVariables == *R.get_dDamagedPreviousStateVariables( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setDamageDeformationGradient, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector cauchyStress = { 1, 1, 1, 1, 1, 1, 1, 1, 1 };

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setStress( ){

                tardigradeHydra::residualBase::setStress( cauchyStress );

            }

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    }
;
    class residualMock : public tardigradeHydra::perzynaViscodamage::residual{

        public:

            residualMock( ) : tardigradeHydra::perzynaViscodamage::residual( ){ }

            residualMock( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations,
                          const unsigned int &damageConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices,
                          const floatVector &parameters ) : tardigradeHydra::perzynaViscodamage::residual( hydra, numEquations, damageConfigurationIndex, stateVariableIndices, parameters ){ }

            floatType damage = 0.3;

        protected:

            virtual void setDamage( ) override{

                set_damage( damage );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector viscoDamageParameters = { 10.0,
                                                  1e1, 1e2,
                                                  10, 200, 293.15,
                                                  1, 0.34,
                                                  0.12,
                                                  13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            residualMock damage;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                damage = residualMock( this, 11, 1, stateVariableIndices, viscoDamageParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &damage;

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

    floatVector previousStateVariables = { 0.1, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0.0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0.01, 0.02, 0.03, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 6, 7, 8 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R_ngrad( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    residualMock R_grad1( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    residualMock R_grad2( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    floatVector damageDeformationGradientAnswer = { 0.99977077, 0.02141056, 0.        ,
                                                    0.02141056, 1.00191182, 0.        ,
                                                    0.        , 0.        , 1.         };

    R_grad1.get_dDamageDeformationGradientdCauchyStress( );

    R_grad2.get_dDamageDeformationGradientdPreviousCauchyStress( );

    BOOST_TEST( damageDeformationGradientAnswer == *R_ngrad.get_damageDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( damageDeformationGradientAnswer == *R_grad1.get_damageDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( damageDeformationGradientAnswer == *R_grad2.get_damageDeformationGradient( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setDamageDeformationGradient2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoDamageParameters = { 2.0,
                                                  1e1, 1e2,
                                                  10, 200, 293.15,
                                                  1, 0.34,
                                                  0.12,
                                                  13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            tardigradeHydra::perzynaViscodamage::residual damage;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                damage = tardigradeHydra::perzynaViscodamage::residual( this, 11, 1, stateVariableIndices, viscoDamageParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &damage;

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

    floatVector previousStateVariables = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0.01, 0.02, 0.03, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 3, 1, 4, 1, 5, 1, 1, 1, 8,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  0.04, 0.05, 0.06, 0.07, 0.08 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::perzynaViscodamage::residual R( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    // Test the computation of the Jacobians
    unsigned int nvals = 9;

    floatMatrix dDamageDeformationGradientdCauchyStress( 9, floatVector( 9, 0 ) );

    floatMatrix dDamageDeformationGradientdF( 9, floatVector( 9, 0 ) );

    floatMatrix dDamageDeformationGradientdSubFs( 9, floatVector( 18, 0 ) );

    floatVector dDamageDeformationGradientdT( 9, 0 );

    floatMatrix dDamageDeformationGradientdStateVariables( 9, floatVector( 2, 0 ) );

    floatMatrix dDamageDeformationGradientdPreviousCauchyStress( 9, floatVector( 9, 0 ) );

    floatMatrix dDamageDeformationGradientdPreviousF( 9, floatVector( 9, 0 ) );

    floatMatrix dDamageDeformationGradientdPreviousSubFs( 9, floatVector( 18, 0 ) );

    floatVector dDamageDeformationGradientdPreviousT( 9, 0 );

    floatMatrix dDamageDeformationGradientdPreviousStateVariables( 9, floatVector( 2, 0 ) );

    floatType eps = 1e-6;
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] += eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_damageDeformationGradient( );

        floatVector vm = *Rm.get_damageDeformationGradient( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamageDeformationGradientdCauchyStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dDamageDeformationGradientdCauchyStress ), *R.get_dDamageDeformationGradientdCauchyStress( ) ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] += eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_damageDeformationGradient( );

        floatVector vm = *Rm.get_damageDeformationGradient( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamageDeformationGradientdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dDamageDeformationGradientdF ), *R.get_dDamageDeformationGradientdF( ) ) );

    for ( unsigned int i = 0; i < 2*deformationGradient.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] += eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_damageDeformationGradient( );

        floatVector vm = *Rm.get_damageDeformationGradient( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamageDeformationGradientdSubFs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dDamageDeformationGradientdSubFs ), *R.get_dDamageDeformationGradientdSubFs( ) ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector delta( 1, 0 );

        delta[ i ] += eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + delta[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature - delta[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_damageDeformationGradient( );

        floatVector vm = *Rm.get_damageDeformationGradient( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamageDeformationGradientdT[ j ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( dDamageDeformationGradientdT, *R.get_dDamageDeformationGradientdT( ) ) );

    for ( unsigned int i = 0; i < 2; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 + 18 + 1 ] += eps * std::fabs( unknownVector[ 9 + 18 + 1 + i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_damageDeformationGradient( );

        floatVector vm = *Rm.get_damageDeformationGradient( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamageDeformationGradientdStateVariables[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + 9 + 18 + 1 ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dDamageDeformationGradientdStateVariables ), *R.get_dDamageDeformationGradientdStateVariables( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] += eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap._local_deltaPreviousCauchyStress = delta;

        hydram._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_damageDeformationGradient( );

        floatVector vm = *Rm.get_damageDeformationGradient( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamageDeformationGradientdPreviousCauchyStress[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dDamageDeformationGradientdPreviousCauchyStress ), *R.get_dDamageDeformationGradientdPreviousCauchyStress( ) ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] += eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_damageDeformationGradient( );

        floatVector vm = *Rm.get_damageDeformationGradient( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamageDeformationGradientdPreviousF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dDamageDeformationGradientdPreviousF ), *R.get_dDamageDeformationGradientdPreviousF( ) ) );

    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] += eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_damageDeformationGradient( );

        floatVector vm = *Rm.get_damageDeformationGradient( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamageDeformationGradientdPreviousSubFs[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dDamageDeformationGradientdPreviousSubFs ), *R.get_dDamageDeformationGradientdPreviousSubFs( ) ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector delta( 1, 0 );

        delta[ i ] += eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature + delta[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature - delta[ 0 ], deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_damageDeformationGradient( );

        floatVector vm = *Rm.get_damageDeformationGradient( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamageDeformationGradientdPreviousT[ j ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( dDamageDeformationGradientdPreviousT, *R.get_dDamageDeformationGradientdPreviousT( ) ) );

    for ( unsigned int i = 0; i < 2; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + 18 + 1 ] += eps * std::fabs( previousStateVariables[ 18 + 1 + i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.get_damageDeformationGradient( );

        floatVector vm = *Rm.get_damageDeformationGradient( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamageDeformationGradientdPreviousStateVariables[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + 18 + 1 ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dDamageDeformationGradientdPreviousStateVariables ), *R.get_dDamageDeformationGradientdPreviousStateVariables( ) ) );

}

BOOST_AUTO_TEST_CASE( test_setResidual, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector cauchyStress = { 1, 1, 1, 1, 1, 1, 1, 1, 1 };

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setStress( ){

                tardigradeHydra::residualBase::setStress( cauchyStress );

            }

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    }
;
    class residualMock : public tardigradeHydra::perzynaViscodamage::residual{

        public:

            residualMock( ) : tardigradeHydra::perzynaViscodamage::residual( ){ }

            residualMock( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations,
                          const unsigned int &damageConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices,
                          const floatVector &parameters ) : tardigradeHydra::perzynaViscodamage::residual( hydra, numEquations, damageConfigurationIndex, stateVariableIndices, parameters ){ }

            floatVector damageDeformationGradient = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector plasticStateVariables = { 10, 0.3 };

            floatVector stateVariables = { 5, 6 };

        protected:

            using tardigradeHydra::perzynaViscodamage::residual::setStateVariables;

            virtual void setDamageDeformationGradient( ) override{

                set_damageDeformationGradient( damageDeformationGradient );

            }

            virtual void setPlasticStateVariables( ) override{

                set_plasticStateVariables( plasticStateVariables );

            }

            virtual void setStateVariables( ) override{

                set_stateVariables( stateVariables );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector viscoDamageParameters = { 10.0,
                                                  1e1, 1e2,
                                                  10, 200, 293.15,
                                                  1, 0.34,
                                                  0.12,
                                                  13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            residualMock damage;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                damage = residualMock( this, 11, 1, stateVariableIndices, viscoDamageParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &damage;

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

    floatVector previousStateVariables = { 0.1, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0.0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0.01, 0.02, 0.03, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 6, 7, 8 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R_ngrad( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    floatVector residualAnswer = { 0.2, -2. , -3. , -4. , -4. , -6. , -7. , -8. , -8. , -5. ,  5.7 };

    BOOST_TEST( residualAnswer == *R_ngrad.getResidual( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_setResidual2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoDamageParameters = { 2.0,
                                                  1e1, 1e2,
                                                  10, 200, 293.15,
                                                  1, 0.34,
                                                  0.12,
                                                  13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            tardigradeHydra::perzynaViscodamage::residual damage;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                damage = tardigradeHydra::perzynaViscodamage::residual( this, 11, 1, stateVariableIndices, viscoDamageParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &damage;

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

    floatVector previousStateVariables = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0.01, 0.02, 0.03, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 3, 1, 4, 1, 5, 1, 1, 1, 8,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  0.04, 0.05, 0.06, 0.07, 0.08 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::perzynaViscodamage::residual R( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    // Test the computation of the Jacobians
    unsigned int nvals = 11;

    floatMatrix jacobian( nvals, floatVector( unknownVector.size( ), 0 ) );

    floatMatrix dRdF( nvals, floatVector( 9, 0 ) );

    floatVector dRdT( nvals, 0 );

    floatType eps = 1e-6;
    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] += eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.getResidual( );

        floatVector vm = *Rm.getResidual( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            jacobian[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( jacobian ), *R.getJacobian( ) ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] += eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.getResidual( );

        floatVector vm = *Rm.getResidual( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dRdF[ j ][ i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dRdF ), *R.getdRdF( ) ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector delta( 1, 0 );

        delta[ i ] += eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + delta[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature - delta[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::perzynaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::perzynaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatVector vp = *Rp.getResidual( );

        floatVector vm = *Rm.getResidual( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dRdT[ j ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( dRdT, *R.getdRdT( ) ) );

}

BOOST_AUTO_TEST_CASE( test_addParameterizationInfo, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoDamageParameters = { 2.0,
                                                  1e1, 1e2,
                                                  10, 200, 293.15,
                                                  1, 0.34,
                                                  0.12,
                                                  13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            tardigradeHydra::perzynaViscodamage::residual damage;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                damage = tardigradeHydra::perzynaViscodamage::residual( this, 11, 1, stateVariableIndices, viscoDamageParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &damage;

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

    floatVector previousStateVariables = { 0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0, 0, 0, 0, 0, 0, 0, 0, 0,
                                           0.01, 0.02, 0.03, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 3, 1, 4, 1, 5, 1, 1, 1, 8,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  0.04, 0.05, 0.06, 0.07, 0.08 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::perzynaViscodamage::residual R( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    std::string output;
    R.addParameterizationInfo(output);

}
