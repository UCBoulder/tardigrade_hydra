/**
  * \file test_tardigrade_hydraPeryznaViscodamage.cpp
  *
  * Tests for tardigrade_hydraPeryznaViscodamage
  */

#include<tardigrade_hydraPeryznaViscodamage.h>
#include<sstream>
#include<fstream>

#define BOOST_TEST_MODULE test_tardigrade_hydra
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

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

BOOST_AUTO_TEST_CASE( testSayHello ){
    /*!
     * Test message printed to stdout in sayHello function
     */

    //Setup redirect variables for stdout
    std::stringbuf buffer;
    cout_redirect rd(&buffer);
    boost::test_tools::output_test_stream result;

    //Initialize test variables
    std::string message;
    std::string answer;
    errorOut error = NULL;

    cout_redirect guard( result.rdbuf() );

    //Check normal operation
    message = "World!";
    answer = "Hello World!\n";
    error = tardigradeHydra::sayHello(message);
    BOOST_CHECK( ! error );
    BOOST_CHECK( result.is_equal( answer ) );

    //Reset error code between tests
    error = NULL;

    //Check for "George" error
    message = "George";
    error = tardigradeHydra::sayHello(message);
    BOOST_CHECK( error );

}

BOOST_AUTO_TEST_CASE( test_setDamage ){

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::peryznaViscodamage::residual{

        public:

            residualMock( ) : tardigradeHydra::peryznaViscodamage::residual( ){ }

            residualMock( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations,
                          const unsigned int &damageConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices,
                          const floatVector &parameters ) : tardigradeHydra::peryznaViscodamage::residual( hydra, numEquations, damageConfigurationIndex, stateVariableIndices, parameters ){ }

            floatType plasticMultiplier = 1.2;

            floatType previousPlasticMultiplier = 3.4;

            floatVector previousStateVariables = { 1, 2 };

            floatVector stateVariables = { 3, 4 };

        protected:

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R_ngrad( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    residualMock R_grad1( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    residualMock R_grad2( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    floatType damageAnswer = 7.06; 

    R_grad1.get_dDamagedCauchyStress( );

    R_grad2.get_dDamagedPreviousCauchyStress( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( damageAnswer, *R_ngrad.get_damage( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( damageAnswer, *R_grad1.get_damage( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( damageAnswer, *R_grad2.get_damage( ) ) );

}

BOOST_AUTO_TEST_CASE( test_setDamage2 ){

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

            tardigradeHydra::peryznaViscodamage::residual damage;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                damage = tardigradeHydra::peryznaViscodamage::residual( this, 11, 1, stateVariableIndices, viscoDamageParameters );

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    tardigradeHydra::peryznaViscodamage::residual R( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

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
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        tardigradeHydra::peryznaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::peryznaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedCauchyStress[ i ] = ( vp - vm ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDamagedCauchyStress, *R.get_dDamagedCauchyStress( ) ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] += eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::peryznaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::peryznaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedF[ i ] = ( vp - vm ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDamagedF, *R.get_dDamagedF( ) ) );

    for ( unsigned int i = 0; i < 2*deformationGradient.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] += eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        tardigradeHydra::peryznaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::peryznaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedSubFs[ i ] = ( vp - vm ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDamagedSubFs, *R.get_dDamagedSubFs( ) ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector delta( 1, 0 );

        delta[ i ] += eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + delta[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature - delta[ 0 ], previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::peryznaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::peryznaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedT = ( vp - vm ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDamagedT, *R.get_dDamagedT( ) ) );

    for ( unsigned int i = 0; i < 2; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 + 18 + 1 ] += eps * std::fabs( unknownVector[ 9 + 18 + 1 + i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector - delta );

        tardigradeHydra::peryznaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::peryznaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedStateVariables[ i ] = ( vp - vm ) / ( 2 * delta[ i + 9 + 18 + 1 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDamagedStateVariables, *R.get_dDamagedStateVariables( ) ) );

    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] += eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydrap._local_deltaPreviousCauchyStress = delta;

        hydram._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::peryznaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::peryznaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedPreviousCauchyStress[ i ] = ( vp - vm ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDamagedPreviousCauchyStress, *R.get_dDamagedPreviousCauchyStress( ) ) );

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] += eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::peryznaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::peryznaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedPreviousF[ i ] = ( vp - vm ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDamagedPreviousF, *R.get_dDamagedPreviousF( ) ) );

    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] += eps * std::fabs( previousStateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::peryznaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::peryznaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedPreviousSubFs[ i ] = ( vp - vm ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDamagedPreviousSubFs, *R.get_dDamagedPreviousSubFs( ) ) );

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector delta( 1, 0 );

        delta[ i ] += eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature + delta[ 0 ], deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature - delta[ 0 ], deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::peryznaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::peryznaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedPreviousT = ( vp - vm ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDamagedPreviousT, *R.get_dDamagedPreviousT( ) ) );

    for ( unsigned int i = 0; i < 2; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i + 18 + 1 ] += eps * std::fabs( previousStateVariables[ 18 + 1 + i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydrap, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydram, unknownVector );

        tardigradeHydra::peryznaViscodamage::residual Rp( &hydrap, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        tardigradeHydra::peryznaViscodamage::residual Rm( &hydram, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

        floatType vp = *Rp.get_damage( );

        floatType vm = *Rm.get_damage( );

        for ( unsigned int j = 0; j < nvals; j++ ){

            dDamagedPreviousStateVariables[ i ] = ( vp - vm ) / ( 2 * delta[ i + 18 + 1 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDamagedPreviousStateVariables, *R.get_dDamagedPreviousStateVariables( ) ) );

}

BOOST_AUTO_TEST_CASE( test_setDamageDeformationGradient ){

    class residualMock : public tardigradeHydra::peryznaViscodamage::residual{

        public:

            residualMock( ) : tardigradeHydra::peryznaViscodamage::residual( ){ }

            residualMock( tardigradeHydra::hydraBase* hydra, const unsigned int &numEquations,
                          const unsigned int &damageConfigurationIndex, const std::vector< unsigned int > &stateVariableIndices,
                          const floatVector &parameters ) : tardigradeHydra::peryznaViscodamage::residual( hydra, numEquations, damageConfigurationIndex, stateVariableIndices, parameters ){ }

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

            tardigradeHydra::residualBase elasticity;

            residualMock damage;

            tardigradeHydra::residualBase remainder;

            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){

                tardigradeHydra::hydraBase::setResidualClasses( residuals );

            }

        protected:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 3 );

                elasticity = tardigradeHydra::residualBase( this, 9 );

                damage = residualMock( this, 11, 1, stateVariableIndices, viscoDamageParameters );

                remainder = tardigradeHydra::residualBase( this, 12 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &damage;

                residuals[ 2 ] = &remainder;

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R_ngrad( &hydra, 11, 1, hydra.stateVariableIndices, hydra.viscoDamageParameters );

    floatVector damageDeformationGradientAnswer = { 0.99977077, 0.02141056, 0.        ,
                                                    0.02141056, 1.00191182, 0.        ,
                                                    0.        , 0.        , 1.         };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( damageDeformationGradientAnswer, *R_ngrad.get_damageDeformationGradient( ) ) );

}
