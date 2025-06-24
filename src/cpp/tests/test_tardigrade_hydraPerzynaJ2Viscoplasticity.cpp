/**
  * \file test_tardigrade_hydraPerzynaJ2Viscoplasticity.cpp
  *
  * Tests for tardigrade_hydraPerzynaJ2Viscoplasticity
  */

#include<tardigrade_hydraPerzynaJ2Viscoplasticity.h>
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

BOOST_AUTO_TEST_CASE( test_sgn, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    BOOST_CHECK( tardigradeHydra::perzynaJ2Viscoplasticity::sgn(  1.4 ) ==  1 );

    BOOST_CHECK( tardigradeHydra::perzynaJ2Viscoplasticity::sgn( -1.4 ) == -1 );

    BOOST_CHECK( tardigradeHydra::perzynaJ2Viscoplasticity::sgn(  0   ) ==  0 );

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

BOOST_AUTO_TEST_CASE( test_get_yieldFunction, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 8 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaJ2Viscoplasticity::residual{

        public:

            residualMock( ) : tardigradeHydra::perzynaJ2Viscoplasticity::residual( ){ }

            residualMock( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations,
                          const unsigned int &plasticConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices,
                          const floatVector &parameters ) : tardigradeHydra::perzynaJ2Viscoplasticity::residual( hydra, numEquations, plasticConfigurationIndex, stateVariableIndices, parameters ){ }

            floatVector driving_stress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previous_driving_stress = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            floatVector dDrivingStressdCauchyStress = fillVector( 81 );

            floatVector dDrivingStressdF = fillVector( 9 * 9 );

            floatVector dDrivingStressdSubFs = fillVector( 9 * 9 * 2 );

        protected:

            using tardigradeHydra::perzynaJ2Viscoplasticity::residual::setDrivingStress;

            floatVector fillVector( const unsigned int npoints ){
                /*!
                 * Fill a vector with values
                 */

                return floatVector( npoints, 0 );

            }

            virtual void setDrivingStress( const bool isPrevious ) override{

                if ( isPrevious ){

                    set_previousDrivingStress( previous_driving_stress );

                }
                else{

                    set_drivingStress( driving_stress );

                }

            }

            virtual void setDrivingStressDerivatives( const bool isPrevious ) override{

                if ( isPrevious ){

                    set_previousDrivingStress( previous_driving_stress );

                    set_dPreviousDrivingStressdPreviousCauchyStress( dDrivingStressdCauchyStress );

                    set_dPreviousDrivingStressdPreviousF( dDrivingStressdF );

                    set_dPreviousDrivingStressdPreviousSubFs( dDrivingStressdSubFs );

                }
                else{

                    set_drivingStress( driving_stress );

                    set_dDrivingStressdCauchyStress( dDrivingStressdCauchyStress );

                    set_dDrivingStressdF( dDrivingStressdF );

                    set_dDrivingStressdSubFs( dDrivingStressdSubFs );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector viscoPlasticParameters = { 2.0, 1e1,
                                                   10, 200, 293.15,
                                                   0.1, 1, 0.34,
                                                   0.12, 13.,
                                                   0.14, 15. };

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            residualMock plasticity;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                plasticity = residualMock( this, 11, 1, stateVariableIndices, viscoPlasticParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

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
                                           0.01, 0.02, 6.00, 0.04, 0.05 };

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

    residualMock R_ngrad( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R_grad1( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R_grad2( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType answer = 10.60823934929885;

    floatType previousAnswer = 0.14517606507011494;

    floatType sign_term_answer = 1;

    floatType previous_sign_term_answer = -1;

    R_grad1.get_dYieldFunctiondStateVariables( );

    R_grad2.get_dPreviousYieldFunctiondPreviousStateVariables( );

    BOOST_TEST( answer == *R_ngrad.get_yieldFunction( ) );

    BOOST_TEST( previousAnswer == *R_ngrad.get_previousYieldFunction( ) );

    BOOST_TEST( sign_term_answer == *R_ngrad.get_signTerm( ) );

    BOOST_TEST( previous_sign_term_answer == *R_ngrad.get_previousSignTerm( ) );

    BOOST_TEST( answer == *R_grad1.get_yieldFunction( ) );

    BOOST_TEST( previousAnswer == *R_grad1.get_previousYieldFunction( ) );

    BOOST_TEST( sign_term_answer == *R_grad1.get_signTerm( ) );

    BOOST_TEST( previous_sign_term_answer == *R_grad1.get_previousSignTerm( ) );

    BOOST_TEST( answer == *R_grad2.get_yieldFunction( ) );

    BOOST_TEST( previousAnswer == *R_grad2.get_previousYieldFunction( ) );

    BOOST_TEST( sign_term_answer == *R_grad2.get_signTerm( ) );

    BOOST_TEST( previous_sign_term_answer == *R_grad2.get_previousSignTerm( ) );

}

BOOST_AUTO_TEST_CASE( test_get_yieldFunction2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { .1, .2, .3, .4, .5, .6, .7, .8, .9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaJ2Viscoplasticity::residual{

        public:

            residualMock( ) : tardigradeHydra::perzynaJ2Viscoplasticity::residual( ){ }

            residualMock( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations,
                          const unsigned int &plasticConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices,
                          const floatVector &parameters ) : tardigradeHydra::perzynaJ2Viscoplasticity::residual( hydra, numEquations, plasticConfigurationIndex, stateVariableIndices, parameters ){ }

        protected:

            floatVector fillVector( const unsigned int npoints ){
                /*!
                 * Fill a vector with values
                 */

                return floatVector( npoints, 0 );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector viscoPlasticParameters = { 2.0, 1e1,
                                                   10, 200, 293.15,
                                                   0.1, 1, 0.34,
                                                   0.12, 13.,
                                                   0.14, 15. };

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            residualMock plasticity;

            tardigradeHydra::residualBase remainder;

            floatVector delta_previous_stress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += delta_previous_stress;

                plasticity = residualMock( this, 11, 1, stateVariableIndices, viscoPlasticParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

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
                                           0.01, 0.02, 6.00, 0.04, 0.05 };

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

    residualMock R_grad( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatVector dYieldFunctiondCauchyStress( 9, 0 );

    floatVector dYieldFunctiondF( 9, 0 );

    floatVector dYieldFunctiondSubFs( 9 * 2, 0 );

    floatVector dYieldFunctiondStateVariables( 2, 0 );

    floatVector dPreviousYieldFunctiondPreviousCauchyStress( 9, 0 );

    floatVector dPreviousYieldFunctiondPreviousF( 9, 0 );

    floatVector dPreviousYieldFunctiondPreviousSubFs( 9 * 2, 0 );

    floatVector dPreviousYieldFunctiondPreviousStateVariables( 2, 0 );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        residualMock Rp( &hydrap, 11, 1, hydrap.stateVariableIndices, hydrap.viscoPlasticParameters );

        residualMock Rm( &hydram, 11, 1, hydram.stateVariableIndices, hydram.viscoPlasticParameters );

        floatType vp = *Rp.get_yieldFunction( );

        floatType vm = *Rm.get_yieldFunction( );

        dYieldFunctiondCauchyStress[ i ] = ( vp - vm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( *R_grad.get_dYieldFunctiondCauchyStress( ) == dYieldFunctiondCauchyStress, CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 11, 1, hydrap.stateVariableIndices, hydrap.viscoPlasticParameters );

        residualMock Rm( &hydram, 11, 1, hydram.stateVariableIndices, hydram.viscoPlasticParameters );

        floatType vp = *Rp.get_yieldFunction( );

        floatType vm = *Rm.get_yieldFunction( );

        dYieldFunctiondF[ i ] = ( vp - vm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( *R_grad.get_dYieldFunctiondF( ) == dYieldFunctiondF, CHECK_PER_ELEMENT );

    unsigned int offset = 9;

    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + offset ] = eps * std::fabs( unknownVector[ i + offset ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        residualMock Rp( &hydrap, 11, 1, hydrap.stateVariableIndices, hydrap.viscoPlasticParameters );

        residualMock Rm( &hydram, 11, 1, hydram.stateVariableIndices, hydram.viscoPlasticParameters );

        floatType vp = *Rp.get_yieldFunction( );

        floatType vm = *Rm.get_yieldFunction( );

        dYieldFunctiondSubFs[ i ] = ( vp - vm ) / ( 2 * delta[ i + offset ] );

    }

    BOOST_TEST( *R_grad.get_dYieldFunctiondSubFs( ) == dYieldFunctiondSubFs, CHECK_PER_ELEMENT );

    offset = 9+18+1;

    for ( unsigned int i = 0; i < 2; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + offset ] = eps * std::fabs( unknownVector[ i + offset ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        residualMock Rp( &hydrap, 11, 1, hydrap.stateVariableIndices, hydrap.viscoPlasticParameters );

        residualMock Rm( &hydram, 11, 1, hydram.stateVariableIndices, hydram.viscoPlasticParameters );

        floatType vp = *Rp.get_yieldFunction( );

        floatType vm = *Rm.get_yieldFunction( );

        dYieldFunctiondStateVariables[ i ] = ( vp - vm ) / ( 2 * delta[ i + offset ] );

    }

    BOOST_TEST( *R_grad.get_dYieldFunctiondStateVariables( ) == dYieldFunctiondStateVariables, CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap.delta_previous_stress = delta;

        hydram.delta_previous_stress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 11, 1, hydrap.stateVariableIndices, hydrap.viscoPlasticParameters );

        residualMock Rm( &hydram, 11, 1, hydram.stateVariableIndices, hydram.viscoPlasticParameters );

        floatType vp = *Rp.get_previousYieldFunction( );

        floatType vm = *Rm.get_previousYieldFunction( );

        dPreviousYieldFunctiondPreviousCauchyStress[ i ] = ( vp - vm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( *R_grad.get_dPreviousYieldFunctiondPreviousCauchyStress( ) == dPreviousYieldFunctiondPreviousCauchyStress, CHECK_PER_ELEMENT );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 11, 1, hydrap.stateVariableIndices, hydrap.viscoPlasticParameters );

        residualMock Rm( &hydram, 11, 1, hydram.stateVariableIndices, hydram.viscoPlasticParameters );

        floatType vp = *Rp.get_previousYieldFunction( );

        floatType vm = *Rm.get_previousYieldFunction( );

        dPreviousYieldFunctiondPreviousF[ i ] = ( vp - vm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( *R_grad.get_dPreviousYieldFunctiondPreviousF( ) == dPreviousYieldFunctiondPreviousF, CHECK_PER_ELEMENT );

    offset = 0;

    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + offset ] = eps * std::fabs( previousStateVariables[ i + offset ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 11, 1, hydrap.stateVariableIndices, hydrap.viscoPlasticParameters );

        residualMock Rm( &hydram, 11, 1, hydram.stateVariableIndices, hydram.viscoPlasticParameters );

        floatType vp = *Rp.get_previousYieldFunction( );

        floatType vm = *Rm.get_previousYieldFunction( );

        dPreviousYieldFunctiondPreviousSubFs[ i ] = ( vp - vm ) / ( 2 * delta[ i + offset ] );

    }

    BOOST_TEST( *R_grad.get_dPreviousYieldFunctiondPreviousSubFs( ) == dPreviousYieldFunctiondPreviousSubFs, CHECK_PER_ELEMENT );

    offset = 18+1;

    for ( unsigned int i = 0; i < 2; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + offset ] = eps * std::fabs( previousStateVariables[ i + offset ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 11, 1, hydrap.stateVariableIndices, hydrap.viscoPlasticParameters );

        residualMock Rm( &hydram, 11, 1, hydram.stateVariableIndices, hydram.viscoPlasticParameters );

        floatType vp = *Rp.get_previousYieldFunction( );

        floatType vm = *Rm.get_previousYieldFunction( );

        dPreviousYieldFunctiondPreviousStateVariables[ i ] = ( vp - vm ) / ( 2 * delta[ i + offset ] );

    }

    BOOST_TEST( *R_grad.get_dPreviousYieldFunctiondPreviousStateVariables( ) == dPreviousYieldFunctiondPreviousStateVariables, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_get_dragStress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { .1, .2, .3, .4, .5, .6, .7, .8, .9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaJ2Viscoplasticity::residual{

        public:

            residualMock( ) : tardigradeHydra::perzynaJ2Viscoplasticity::residual( ){ }

            residualMock( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations,
                          const unsigned int &plasticConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices,
                          const floatVector &parameters ) : tardigradeHydra::perzynaJ2Viscoplasticity::residual( hydra, numEquations, plasticConfigurationIndex, stateVariableIndices, parameters ){ }

        protected:

            floatVector fillVector( const unsigned int npoints ){
                /*!
                 * Fill a vector with values
                 */

                return floatVector( npoints, 0 );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector viscoPlasticParameters = { 2.0, 1e1,
                                                   10, 200, 293.15,
                                                   0.1, 1, 0.34,
                                                   0.12, 13.,
                                                   0.14, 15. };

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            residualMock plasticity;

            tardigradeHydra::residualBase remainder;

            floatVector delta_previous_stress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += delta_previous_stress;

                plasticity = residualMock( this, 11, 1, stateVariableIndices, viscoPlasticParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

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
                                           0.01, 0.02, 6.00, 0.04, 0.05 };

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

    residualMock R( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );
    residualMock R_grad1( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );
    residualMock R_grad2( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType answer = 1e1;

    BOOST_TEST( answer == *R.get_dragStress( ) );
    BOOST_TEST( answer == *R.get_previousDragStress( ) );

    R_grad1.get_dDragStressdStateVariables( );
    R_grad2.get_dPreviousDragStressdPreviousStateVariables( );
    
    BOOST_TEST( answer == *R_grad1.get_dragStress( ) );
    BOOST_TEST( answer == *R_grad1.get_previousDragStress( ) );

    BOOST_TEST( answer == *R_grad2.get_dragStress( ) );
    BOOST_TEST( answer == *R_grad2.get_previousDragStress( ) );

    floatType eps = 1e-6;

    floatVector J1( 2, 0 );
    floatVector J2( 2, 0 );

    unsigned int offset = 28;

    for ( unsigned int i = 0; i < 2; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + offset ] = eps * std::fabs( delta[ i + offset ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        residualMock Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        residualMock Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        floatType rp = *Rp.get_dragStress( );

        floatType rm = *Rm.get_dragStress( );

        J1[ i ] = ( rp - rm ) / ( 2 * delta[ i + offset ] );

    }

    BOOST_TEST( J1 == *R.get_dDragStressdStateVariables( ), CHECK_PER_ELEMENT );

    offset = 19;

    for ( unsigned int i = 0; i < 2; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + offset ] = eps * std::fabs( delta[ i + offset ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        residualMock Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        floatType rp = *Rp.get_dragStress( );

        floatType rm = *Rm.get_dragStress( );

        J2[ i ] = ( rp - rm ) / ( 2 * delta[ i + offset ] );

    }

    BOOST_TEST( J2 == *R.get_dPreviousDragStressdPreviousStateVariables( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_get_hardeningFunction, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 8 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaJ2Viscoplasticity::residual{

        public:

            residualMock( ) : tardigradeHydra::perzynaJ2Viscoplasticity::residual( ){ }

            residualMock( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations,
                          const unsigned int &plasticConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices,
                          const floatVector &parameters ) : tardigradeHydra::perzynaJ2Viscoplasticity::residual( hydra, numEquations, plasticConfigurationIndex, stateVariableIndices, parameters ){ }

            floatType sign_term = 1;

            floatType previous_sign_term = -1;

        protected:

            floatVector fillVector( const unsigned int npoints ){
                /*!
                 * Fill a vector with values
                 */

                return floatVector( npoints, 0 );

            }

            virtual void setSignTerm( ) override{

                set_signTerm( sign_term );

            }

            virtual void setPreviousSignTerm( ) override{

                set_previousSignTerm( previous_sign_term );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector viscoPlasticParameters = { 2.0, 1e1,
                                                   10, 200, 293.15,
                                                   0.1, 1, 0.34,
                                                   0.12, 13.,
                                                   0.14, 15. };

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            residualMock plasticity;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                plasticity = residualMock( this, 11, 1, stateVariableIndices, viscoPlasticParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

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
                                           0.01, 0.02, 6.00, 0.04, 0.05 };

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

    residualMock R_ngrad( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R_grad1( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R_grad2( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatVector answer = { 0.12 + 13. * 5, 0.14 + 15 * 5 };

    floatVector previous_answer = { 0.12 + 13. * 0.02, -(0.14 + 15 * 0.02) };

    R_grad1.get_dHardeningFunctiondStateVariables( );
    R_grad2.get_dPreviousHardeningFunctiondPreviousStateVariables( );

    BOOST_TEST( answer == *R_ngrad.get_hardeningFunction( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previous_answer == *R_ngrad.get_previousHardeningFunction( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *R_grad1.get_hardeningFunction( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previous_answer == *R_grad1.get_previousHardeningFunction( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *R_grad2.get_hardeningFunction( ), CHECK_PER_ELEMENT );

    BOOST_TEST( previous_answer == *R_grad2.get_previousHardeningFunction( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_get_hardeningFunction2, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 8 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaJ2Viscoplasticity::residual{

        public:

            residualMock( ) : tardigradeHydra::perzynaJ2Viscoplasticity::residual( ){ }

            residualMock( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations,
                          const unsigned int &plasticConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices,
                          const floatVector &parameters ) : tardigradeHydra::perzynaJ2Viscoplasticity::residual( hydra, numEquations, plasticConfigurationIndex, stateVariableIndices, parameters ){ }

        protected:

            floatVector fillVector( const unsigned int npoints ){
                /*!
                 * Fill a vector with values
                 */

                return floatVector( npoints, 0 );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector viscoPlasticParameters = { 2.0, 1e1,
                                                   10, 200, 293.15,
                                                   0.1, 1, 0.34,
                                                   0.12, 13.,
                                                   0.14, 15. };

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            residualMock plasticity;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                plasticity = residualMock( this, 11, 1, stateVariableIndices, viscoPlasticParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

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
                                           0.01, 0.02, 6.00, 0.04, 0.05 };

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

    residualMock R_grad( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType eps = 1e-6;

    unsigned int offset = 28;

    floatVector J1( 4, 0 );
    floatVector J2( 4, 0 );

    for ( unsigned int i = 0; i < 2; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + offset ] = eps * std::fabs( unknownVector[ i + offset ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        residualMock Rp( &hydrap, 11, 1, hydrap.stateVariableIndices, hydrap.viscoPlasticParameters );

        residualMock Rm( &hydram, 11, 1, hydram.stateVariableIndices, hydrap.viscoPlasticParameters );

        floatVector vp = *Rp.get_hardeningFunction( );

        floatVector vm = *Rm.get_hardeningFunction( );

        for ( unsigned int j = 0; j < 2; j++ ){
            J1[ 2 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + offset ] );
        }

    }

    BOOST_TEST( J1 == *R_grad.get_dHardeningFunctiondStateVariables( ), CHECK_PER_ELEMENT );

    offset = 19;

    for ( unsigned int i = 0; i < 2; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + offset ] = eps * std::fabs( previousStateVariables[ i + offset ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        residualMock Rp( &hydrap, 11, 1, hydrap.stateVariableIndices, hydrap.viscoPlasticParameters );

        residualMock Rm( &hydram, 11, 1, hydram.stateVariableIndices, hydrap.viscoPlasticParameters );

        floatVector vp = *Rp.get_previousHardeningFunction( );

        floatVector vm = *Rm.get_previousHardeningFunction( );

        for ( unsigned int j = 0; j < 2; j++ ){
            J2[ 2 * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + offset ] );
        }

    }

    BOOST_TEST( J2 == *R_grad.get_dPreviousHardeningFunctiondPreviousStateVariables( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_get_jacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 8 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaJ2Viscoplasticity::residual{

        public:

            residualMock( ) : tardigradeHydra::perzynaJ2Viscoplasticity::residual( ){ }

            residualMock( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations,
                          const unsigned int &plasticConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices,
                          const floatVector &parameters ) : tardigradeHydra::perzynaJ2Viscoplasticity::residual( hydra, numEquations, plasticConfigurationIndex, stateVariableIndices, parameters ){ }

        protected:

            floatVector fillVector( const unsigned int npoints ){
                /*!
                 * Fill a vector with values
                 */

                return floatVector( npoints, 0 );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            floatVector viscoPlasticParameters = { 2.0, 1e1,
                                                   10, 200, 293.15,
                                                   0.1, 1, 0.34,
                                                   0.12, 13.,
                                                   0.14, 15. };

            std::vector< unsigned int > stateVariableIndices = { 1, 2 };

            using tardigradeHydra::hydraBase::hydraBase;

            stressMock elasticity;

            residualMock plasticity;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                plasticity = residualMock( this, 11, 1, stateVariableIndices, viscoPlasticParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &plasticity;

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
                                           0.01, 0.02, 6.00, 0.04, 0.05 };

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

    residualMock R_grad( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    R_grad.getResidual( );

    floatType eps = 1e-6;

    unsigned int offset = 0;

    floatVector J1( 11 * unknownVector.size( ), 0 );

    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + offset ] = eps * std::fabs( unknownVector[ i + offset ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        residualMock Rp( &hydrap, 11, 1, hydrap.stateVariableIndices, hydrap.viscoPlasticParameters );

        residualMock Rm( &hydram, 11, 1, hydram.stateVariableIndices, hydrap.viscoPlasticParameters );

        floatVector vp = *Rp.getResidual( );

        floatVector vm = *Rm.getResidual( );

        for ( unsigned int j = 0; j < 11; j++ ){
            J1[ unknownVector.size( ) * j + i ] = ( vp[ j ] - vm[ j ] ) / ( 2 * delta[ i + offset ] );
        }

    }

}
