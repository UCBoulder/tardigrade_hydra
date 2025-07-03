/**
  * \file test_tardigrade_hydraPerzynaViscoplasticity.cpp
  *
  * Tests for tardigrade_hydraPerzynaViscoplasticity
  */

#include<tardigrade_hydraPerzynaViscoplasticity.h>
#include<tardigrade_hydraThermalExpansion.h>
#include<tardigrade_hydraLinearElasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_stress_tools.h>

#define BOOST_TEST_MODULE test_tardigrade_hydraPerzynaViscoplasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::perzynaViscoplasticity::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::perzynaViscoplasticity::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::perzynaViscoplasticity::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

namespace tardigradeHydra{

    namespace unit_test{

        class hydraBaseTester{

            public:

                static void updateUnknownVector( tardigradeHydra::hydraBase &hydra, const floatVector &value ){

                    BOOST_CHECK_NO_THROW( hydra.updateUnknownVector( value ) );

                }

        };

    }

    namespace perzynaViscoplasticity{

        namespace unit_test{

            class residualTester{

                public:

                    static void runBasicGetTests( tardigradeHydra::perzynaViscoplasticity::residual &R ){

                        try{
    
                            BOOST_CHECK( &R._drivingStress.second == R.get_drivingStress( ) );

                            BOOST_CHECK( &R._previousDrivingStress.second == R.get_previousDrivingStress( ) );
        
                            BOOST_CHECK( &R._flowDirection.second == R.get_flowDirection( ) );

                            BOOST_CHECK( &R._previousFlowDirection.second == R.get_previousFlowDirection( ) );
    
                            BOOST_CHECK( &R._yieldFunction.second == R.get_yieldFunction( ) );

                            BOOST_CHECK( &R._previousYieldFunction.second == R.get_previousYieldFunction( ) );
    
                            BOOST_CHECK( &R._plasticThermalMultiplier.second == R.get_plasticThermalMultiplier( ) );

                            BOOST_CHECK( &R._previousPlasticThermalMultiplier.second == R.get_previousPlasticThermalMultiplier( ) );

                            BOOST_CHECK( &R._dragStress.second == R.get_dragStress( ) );

                            BOOST_CHECK( &R._previousDragStress.second == R.get_previousDragStress( ) );

                            BOOST_CHECK( &R._hardeningFunction.second == R.get_hardeningFunction( ) );

                            BOOST_CHECK( &R._previousHardeningFunction.second == R.get_previousHardeningFunction( ) );
        
                            BOOST_CHECK( &R._plasticMultiplier.second == R.get_plasticMultiplier( ) );

                            BOOST_CHECK( &R._previousPlasticMultiplier.second == R.get_previousPlasticMultiplier( ) );
        
                            BOOST_CHECK( &R._velocityGradient.second == R.get_velocityGradient( ) );

                            BOOST_CHECK( &R._previousVelocityGradient.second == R.get_previousVelocityGradient( ) );

                            BOOST_CHECK( &R._stateVariableEvolutionRates.second == R.get_stateVariableEvolutionRates( ) );

                            BOOST_CHECK( &R._previousStateVariableEvolutionRates.second == R.get_previousStateVariableEvolutionRates( ) );

                            BOOST_CHECK( &R._plasticDeformationGradient.second == R.get_plasticDeformationGradient( ) );

                            BOOST_CHECK( &R._plasticStateVariables.second == R.get_plasticStateVariables( ) );
        
                            BOOST_CHECK( &R._stateVariables.second == R.get_stateVariables( ) );
        
                            BOOST_CHECK( &R._perzynaParameters.second == R.get_perzynaParameters( ) );
        
                            BOOST_CHECK( &R._dragStressParameters.second == R.get_dragStressParameters( ) );
        
                            BOOST_CHECK( &R._thermalParameters.second == R.get_thermalParameters( ) );
        
                            BOOST_CHECK( &R._yieldParameters.second == R.get_yieldParameters( ) );
        
                            BOOST_CHECK( &R._flowParameters.second == R.get_flowParameters( ) );
        
                            BOOST_CHECK( &R._hardeningParameters.second == R.get_hardeningParameters( ) );

                            BOOST_CHECK( &R._integrationParameter == R.getIntegrationParameter( ) );

                            BOOST_CHECK( &R._dDrivingStressdCauchyStress.second == R.get_dDrivingStressdCauchyStress( ) );

                            BOOST_CHECK( &R._dDrivingStressdF.second == R.get_dDrivingStressdF( ) );

                            BOOST_CHECK( &R._dDrivingStressdSubFs.second == R.get_dDrivingStressdSubFs( ) );

                            BOOST_CHECK( &R._dPreviousDrivingStressdPreviousCauchyStress.second == R.get_dPreviousDrivingStressdPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPreviousDrivingStressdPreviousF.second == R.get_dPreviousDrivingStressdPreviousF( ) );

                            BOOST_CHECK( &R._dPreviousDrivingStressdPreviousSubFs.second == R.get_dPreviousDrivingStressdPreviousSubFs( ) );

                            BOOST_CHECK( &R._dFlowDirectiondCauchyStress.second == R.get_dFlowDirectiondCauchyStress( ) );

                            BOOST_CHECK( &R._dFlowDirectiondF.second == R.get_dFlowDirectiondF( ) );

                            BOOST_CHECK( &R._dFlowDirectiondSubFs.second == R.get_dFlowDirectiondSubFs( ) );

                            BOOST_CHECK( &R._dPreviousFlowDirectiondPreviousCauchyStress.second == R.get_dPreviousFlowDirectiondPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPreviousFlowDirectiondPreviousF.second == R.get_dPreviousFlowDirectiondPreviousF( ) );

                            BOOST_CHECK( &R._dPreviousFlowDirectiondPreviousSubFs.second == R.get_dPreviousFlowDirectiondPreviousSubFs( ) );

                            BOOST_CHECK( &R._dYieldFunctiondCauchyStress.second == R.get_dYieldFunctiondCauchyStress( ) );

                            BOOST_CHECK( &R._dYieldFunctiondF.second == R.get_dYieldFunctiondF( ) );

                            BOOST_CHECK( &R._dYieldFunctiondSubFs.second == R.get_dYieldFunctiondSubFs( ) );

                            BOOST_CHECK( &R._dYieldFunctiondStateVariables.second == R.get_dYieldFunctiondStateVariables( ) );

                            BOOST_CHECK( &R._dPreviousYieldFunctiondPreviousCauchyStress.second == R.get_dPreviousYieldFunctiondPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPreviousYieldFunctiondPreviousF.second == R.get_dPreviousYieldFunctiondPreviousF( ) );

                            BOOST_CHECK( &R._dPreviousYieldFunctiondPreviousSubFs.second == R.get_dPreviousYieldFunctiondPreviousSubFs( ) );

                            BOOST_CHECK( &R._dPreviousYieldFunctiondPreviousStateVariables.second == R.get_dPreviousYieldFunctiondPreviousStateVariables( ) );

                            BOOST_CHECK( &R._dPlasticThermalMultiplierdT.second == R.get_dPlasticThermalMultiplierdT( ) );

                            BOOST_CHECK( &R._dPreviousPlasticThermalMultiplierdPreviousT.second == R.get_dPreviousPlasticThermalMultiplierdPreviousT( ) );

                            BOOST_CHECK( &R._dDragStressdStateVariables.second == R.get_dDragStressdStateVariables( ) );

                            BOOST_CHECK( &R._dPreviousDragStressdPreviousStateVariables.second == R.get_dPreviousDragStressdPreviousStateVariables( ) );

                            BOOST_CHECK( &R._dPlasticMultiplierdCauchyStress.second == R.get_dPlasticMultiplierdCauchyStress( ) );

                            BOOST_CHECK( &R._dPlasticMultiplierdF.second == R.get_dPlasticMultiplierdF( ) );

                            BOOST_CHECK( &R._dPlasticMultiplierdSubFs.second == R.get_dPlasticMultiplierdSubFs( ) );

                            BOOST_CHECK( &R._dPlasticMultiplierdT.second == R.get_dPlasticMultiplierdT( ) );

                            BOOST_CHECK( &R._dPlasticMultiplierdStateVariables.second == R.get_dPlasticMultiplierdStateVariables( ) );

                            BOOST_CHECK( &R._dPreviousPlasticMultiplierdPreviousCauchyStress.second == R.get_dPreviousPlasticMultiplierdPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPreviousPlasticMultiplierdPreviousF.second == R.get_dPreviousPlasticMultiplierdPreviousF( ) );

                            BOOST_CHECK( &R._dPreviousPlasticMultiplierdPreviousSubFs.second == R.get_dPreviousPlasticMultiplierdPreviousSubFs( ) );

                            BOOST_CHECK( &R._dPreviousPlasticMultiplierdPreviousT.second == R.get_dPreviousPlasticMultiplierdPreviousT( ) );

                            BOOST_CHECK( &R._dPreviousPlasticMultiplierdPreviousStateVariables.second == R.get_dPreviousPlasticMultiplierdPreviousStateVariables( ) );

                            BOOST_CHECK( &R._dVelocityGradientdCauchyStress.second   == R.get_dVelocityGradientdCauchyStress( ) );

                            BOOST_CHECK( &R._dVelocityGradientdF.second              == R.get_dVelocityGradientdF( ) );

                            BOOST_CHECK( &R._dVelocityGradientdSubFs.second          == R.get_dVelocityGradientdSubFs( ) );

                            BOOST_CHECK( &R._dVelocityGradientdT.second              == R.get_dVelocityGradientdT( ) );

                            BOOST_CHECK( &R._dVelocityGradientdStateVariables.second == R.get_dVelocityGradientdStateVariables( ) );

                            BOOST_CHECK( &R._dPreviousVelocityGradientdPreviousCauchyStress.second   == R.get_dPreviousVelocityGradientdPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPreviousVelocityGradientdPreviousF.second              == R.get_dPreviousVelocityGradientdPreviousF( ) );

                            BOOST_CHECK( &R._dPreviousVelocityGradientdPreviousSubFs.second          == R.get_dPreviousVelocityGradientdPreviousSubFs( ) );

                            BOOST_CHECK( &R._dPreviousVelocityGradientdPreviousT.second              == R.get_dPreviousVelocityGradientdPreviousT( ) );

                            BOOST_CHECK( &R._dPreviousVelocityGradientdPreviousStateVariables.second == R.get_dPreviousVelocityGradientdPreviousStateVariables( ) );

                            BOOST_CHECK( &R._dStateVariableEvolutionRatesdCauchyStress.second   == R.get_dStateVariableEvolutionRatesdCauchyStress( ) );

                            BOOST_CHECK( &R._dStateVariableEvolutionRatesdF.second              == R.get_dStateVariableEvolutionRatesdF( ) );

                            BOOST_CHECK( &R._dStateVariableEvolutionRatesdSubFs.second          == R.get_dStateVariableEvolutionRatesdSubFs( ) );

                            BOOST_CHECK( &R._dStateVariableEvolutionRatesdT.second              == R.get_dStateVariableEvolutionRatesdT( ) );

                            BOOST_CHECK( &R._dStateVariableEvolutionRatesdStateVariables.second == R.get_dStateVariableEvolutionRatesdStateVariables( ) );

                            BOOST_CHECK( &R._dPreviousStateVariableEvolutionRatesdPreviousCauchyStress.second   == R.get_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPreviousStateVariableEvolutionRatesdPreviousF.second              == R.get_dPreviousStateVariableEvolutionRatesdPreviousF( ) );

                            BOOST_CHECK( &R._dPreviousStateVariableEvolutionRatesdPreviousSubFs.second          == R.get_dPreviousStateVariableEvolutionRatesdPreviousSubFs( ) );

                            BOOST_CHECK( &R._dPreviousStateVariableEvolutionRatesdPreviousT.second              == R.get_dPreviousStateVariableEvolutionRatesdPreviousT( ) );

                            BOOST_CHECK( &R._dPreviousStateVariableEvolutionRatesdPreviousStateVariables.second == R.get_dPreviousStateVariableEvolutionRatesdPreviousStateVariables( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdCauchyStress.second           == R.get_dPlasticDeformationGradientdCauchyStress( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdF.second                      == R.get_dPlasticDeformationGradientdF( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdSubFs.second                  == R.get_dPlasticDeformationGradientdSubFs( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdT.second                      == R.get_dPlasticDeformationGradientdT( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdStateVariables.second         == R.get_dPlasticDeformationGradientdStateVariables( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdPreviousCauchyStress.second   == R.get_dPlasticDeformationGradientdPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdPreviousF.second              == R.get_dPlasticDeformationGradientdPreviousF( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdPreviousSubFs.second          == R.get_dPlasticDeformationGradientdPreviousSubFs( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdPreviousT.second              == R.get_dPlasticDeformationGradientdPreviousT( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdPreviousStateVariables.second == R.get_dPlasticDeformationGradientdPreviousStateVariables( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdCauchyStress.second           == R.get_dPlasticStateVariablesdCauchyStress( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdF.second                      == R.get_dPlasticStateVariablesdF( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdSubFs.second                  == R.get_dPlasticStateVariablesdSubFs( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdT.second                      == R.get_dPlasticStateVariablesdT( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdStateVariables.second         == R.get_dPlasticStateVariablesdStateVariables( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdPreviousCauchyStress.second   == R.get_dPlasticStateVariablesdPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdPreviousF.second              == R.get_dPlasticStateVariablesdPreviousF( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdPreviousSubFs.second          == R.get_dPlasticStateVariablesdPreviousSubFs( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdPreviousT.second              == R.get_dPlasticStateVariablesdPreviousT( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdPreviousStateVariables.second == R.get_dPlasticStateVariablesdPreviousStateVariables( ) );

                        }
                        catch( std::exception &e ){

                            tardigradeErrorTools::printNestedExceptions( e );
                            throw e;

                        }
    
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

BOOST_AUTO_TEST_CASE( test_residual_basicGetTests, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class stressResidualMock : public tardigradeHydra::residualBase {

        public:

            floatVector currentCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousCauchyStress = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

        private:

            using tardigradeHydra::residualBase::residualBase;

            virtual void setStress( ){

                tardigradeHydra::residualBase::setStress( currentCauchyStress );

            }

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14};

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            stressResidualMock elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = stressResidualMock( this, 9 );

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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
                                           1, 2, 3, 4, 5, 7, 8, 8, 9,
                                          -1, -2, -3, -4, -5 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 1, 1, 1, 1, 1, 1, 1, 1, 1,
                                  2, 2, 2, 2, 2, 2, 2, 2, 2,
                                  3, 3, 3, 3, 3, 3, 3, 3, 3,
                                  4, 5, 6, 7, 8 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );
    
    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );
    
    tardigradeHydra::perzynaViscoplasticity::unit_test::residualTester::runBasicGetTests( R );

}

BOOST_AUTO_TEST_CASE( test_residual_get_drivingStress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the driving stress
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
            stressMock elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock Rjac( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatVector answer = { 0.82636364, 1.0       , 0.75,
                           1.        , 1.21012101, 0.90759076,
                           0.75      , 0.90759076, 0.68069307 };

    Rjac.get_dDrivingStressdCauchyStress( );

    Rjac.get_dPreviousDrivingStressdPreviousCauchyStress( );

    BOOST_TEST( answer == *R.get_drivingStress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *Rjac.get_drivingStress( ), CHECK_PER_ELEMENT );

    answer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    BOOST_TEST( answer == *R.get_previousDrivingStress( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *Rjac.get_previousDrivingStress( ), CHECK_PER_ELEMENT );

    // Test the jacobians
    floatType eps = 1e-6;

    floatMatrix dDrivingStressdCauchyStress( 9, floatVector( 9, 0 ) );

    floatMatrix dDrivingStressdF( 9, floatVector( 9, 0 ) );

    floatMatrix dDrivingStressdSubFs( 9, floatVector( 18, 0 ) );

    floatMatrix dPreviousDrivingStressdPreviousCauchyStress( 9, floatVector( 9, 0 ) );

    floatMatrix dPreviousDrivingStressdPreviousF( 9, floatVector( 9, 0 ) );

    floatMatrix dPreviousDrivingStressdPreviousSubFs( 9, floatVector( 18, 0 ) );


    // Jacobians w.r.t. the Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector + delta );

        floatVector dp = *R.get_drivingStress( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatVector dm = *R.get_drivingStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dDrivingStressdCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dDrivingStressdCauchyStress ) == *Rjac.get_dDrivingStressdCauchyStress( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_drivingStress( );

        floatVector dm = *Rm.get_drivingStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dDrivingStressdF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dDrivingStressdF ), *Rjac.get_dDrivingStressdF( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_drivingStress( );

        floatVector dm = *Rm.get_drivingStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dDrivingStressdSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dDrivingStressdSubFs ), *Rjac.get_dDrivingStressdSubFs( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousDrivingStress( );

        floatVector dm = *Rm.get_previousDrivingStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousDrivingStressdPreviousCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousDrivingStressdPreviousCauchyStress ) == *Rjac.get_dPreviousDrivingStressdPreviousCauchyStress( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousDrivingStress( );

        floatVector dm = *Rm.get_previousDrivingStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousDrivingStressdPreviousF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousDrivingStressdPreviousF ) == *Rjac.get_dPreviousDrivingStressdPreviousF( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousDrivingStress( );

        floatVector dm = *Rm.get_previousDrivingStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousDrivingStressdPreviousSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousDrivingStressdPreviousSubFs ) == *Rjac.get_dPreviousDrivingStressdPreviousSubFs( ), CHECK_PER_ELEMENT );
}

BOOST_AUTO_TEST_CASE( test_residual_get_flowDirection, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the driving stress
     */

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

            floatVector drivingStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousDrivingStress = { 1, 0.2, 0, 0.4, 7, -1, 0.5, 3, .2 };

            floatVector flowParameters = { 1.23, 0.46 };

        private:

            using tardigradeHydra::perzynaViscoplasticity::residual::setDrivingStress;

            virtual void setDrivingStress( const bool isPrevious ) override{

                set_flowParameters( flowParameters );

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_previousDrivingStress( previousDrivingStress );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_drivingStress( drivingStress );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 9, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType dpYield;
    floatVector jac( 9, 0 );
    floatVector answer( 9, 0 );
    floatVector answer2( 9, 0 );

    TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::druckerPragerSurface( R.drivingStress, R.flowParameters[ 1 ], R.flowParameters[ 0 ], dpYield, jac, answer ) );

    TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::druckerPragerSurface( R.previousDrivingStress, R.flowParameters[ 1 ], R.flowParameters[ 0 ], dpYield, jac, answer2 ) );

    BOOST_CHECK( R.get_flowDirection( )->size( ) == 9 );

    BOOST_TEST( answer == *R.get_flowDirection( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer2 == *R.get_previousFlowDirection( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_get_flowDirection_jacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the flow direction stress
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
            stressMock elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R_ngrad( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    R.get_dFlowDirectiondCauchyStress( );

    BOOST_TEST( *R.get_flowDirection( ) == *R_ngrad.get_flowDirection( ), CHECK_PER_ELEMENT );

    // Test the jacobians
    floatType eps = 1e-6;

    floatMatrix dFlowDirectiondCauchyStress( 9, floatVector( 9, 0 ) );

    floatMatrix dFlowDirectiondF( 9, floatVector( 9, 0 ) );

    floatMatrix dFlowDirectiondSubFs( 9, floatVector( 18, 0 ) );

    floatMatrix dPreviousFlowDirectiondPreviousCauchyStress( 9, floatVector( 9, 0 ) );

    floatMatrix dPreviousFlowDirectiondPreviousF( 9, floatVector( 9, 0 ) );

    floatMatrix dPreviousFlowDirectiondPreviousSubFs( 9, floatVector( 18, 0 ) );


    // Jacobians w.r.t. the Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector + delta );

        floatVector dp = *R.get_flowDirection( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatVector dm = *R.get_flowDirection( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dFlowDirectiondCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFlowDirectiondCauchyStress ), *R.get_dFlowDirectiondCauchyStress( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_flowDirection( );

        floatVector dm = *Rm.get_flowDirection( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dFlowDirectiondF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFlowDirectiondF ), *R.get_dFlowDirectiondF( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_flowDirection( );

        floatVector dm = *Rm.get_flowDirection( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dFlowDirectiondSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dFlowDirectiondSubFs ), *R.get_dFlowDirectiondSubFs( ), 2e-5, 1e-5 ) );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousFlowDirection( );

        floatVector dm = *Rm.get_previousFlowDirection( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousFlowDirectiondPreviousCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPreviousFlowDirectiondPreviousCauchyStress ), *R.get_dPreviousFlowDirectiondPreviousCauchyStress( ) ) );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousFlowDirection( );

        floatVector dm = *Rm.get_previousFlowDirection( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousFlowDirectiondPreviousF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousFlowDirectiondPreviousF ) == *R.get_dPreviousFlowDirectiondPreviousF( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousFlowDirection( );

        floatVector dm = *Rm.get_previousFlowDirection( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousFlowDirectiondPreviousSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousFlowDirectiondPreviousSubFs ) == *R.get_dPreviousFlowDirectiondPreviousSubFs( ), CHECK_PER_ELEMENT );
}

BOOST_AUTO_TEST_CASE( test_residual_get_yieldFunction, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the yield function
     */

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

            floatVector drivingStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousDrivingStress = { 1, 0.2, 0, 0.4, 7, -1, 0.5, 3, .2 };

            floatVector yieldParameters = { 1.23, 0.46 };

        private:

            using tardigradeHydra::perzynaViscoplasticity::residual::setDrivingStress;

            virtual void setDrivingStress( const bool isPrevious ) override{

                set_yieldParameters( yieldParameters );

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_previousDrivingStress( previousDrivingStress );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_drivingStress( drivingStress );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 9, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType answer1;
    TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::druckerPragerSurface( R.drivingStress, R.yieldParameters[1], R.yieldParameters[0], answer1 ) );

    BOOST_TEST( answer1 == *R.get_yieldFunction( ) );

    floatType answer2;
    TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeStressTools::druckerPragerSurface( R.previousDrivingStress, R.yieldParameters[1], R.yieldParameters[0], answer2 ) );

    BOOST_TEST( answer2 == *R.get_previousYieldFunction( ) );

}

BOOST_AUTO_TEST_CASE( test_residual_get_yieldFunction_jacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the yield function's Jacobian
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
            stressMock elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R_ngrad( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    R.get_dYieldFunctiondCauchyStress( );

    BOOST_TEST( *R.get_yieldFunction( ) == *R_ngrad.get_yieldFunction( ) );

    // Test the jacobians
    floatType eps = 1e-6;

    floatVector dYieldFunctiondCauchyStress( 9, 0 );

    floatVector dYieldFunctiondF( 9, 0 );

    floatVector dYieldFunctiondSubFs( 18, 0 );

    floatVector dPreviousYieldFunctiondPreviousCauchyStress( 9, 0 );

    floatVector dPreviousYieldFunctiondPreviousF( 9, 0 );

    floatVector dPreviousYieldFunctiondPreviousSubFs( 18, 0 );

    // Jacobians w.r.t. the Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector + delta );

        floatType dp = *R.get_yieldFunction( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatType dm = *R.get_yieldFunction( );

        dYieldFunctiondCauchyStress[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( tolerantCheck( dYieldFunctiondCauchyStress, *R.get_dYieldFunctiondCauchyStress( ), 5e-5, 1e-5 ) );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.get_yieldFunction( );

        floatType dm = *Rm.get_yieldFunction( );

        dYieldFunctiondF[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( dYieldFunctiondF == *R.get_dYieldFunctiondF( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.get_yieldFunction( );

        floatType dm = *Rm.get_yieldFunction( );

        dYieldFunctiondSubFs[ i ] = ( dp - dm ) / ( 2 * delta[ i + 9 ] );

    }

    BOOST_TEST( tolerantCheck( dYieldFunctiondSubFs, *R.get_dYieldFunctiondSubFs( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.get_previousYieldFunction( );

        floatType dm = *Rm.get_previousYieldFunction( );

        dPreviousYieldFunctiondPreviousCauchyStress[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( dPreviousYieldFunctiondPreviousCauchyStress == *R.get_dPreviousYieldFunctiondPreviousCauchyStress( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.get_previousYieldFunction( );

        floatType dm = *Rm.get_previousYieldFunction( );

        dPreviousYieldFunctiondPreviousF[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( dPreviousYieldFunctiondPreviousF == *R.get_dPreviousYieldFunctiondPreviousF( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.get_previousYieldFunction( );

        floatType dm = *Rm.get_previousYieldFunction( );

        dPreviousYieldFunctiondPreviousSubFs[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( dPreviousYieldFunctiondPreviousSubFs == *R.get_dPreviousYieldFunctiondPreviousSubFs( ), CHECK_PER_ELEMENT );
}

BOOST_AUTO_TEST_CASE( test_residual_get_plasticThermalMultiplier, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the plastic thermal multiplier
     */

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 9, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock Rjac( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    Rjac.get_dPlasticThermalMultiplierdT( );

    Rjac.get_dPreviousPlasticThermalMultiplierdPreviousT( );

    floatType exp = ( -hydra.viscoPlasticParameters[ 3 ] * ( temperature - hydra.viscoPlasticParameters[ 5 ] ) / ( hydra.viscoPlasticParameters[ 4 ] + temperature - hydra.viscoPlasticParameters[ 5 ] ) );

    floatType answer = std::pow( 10, exp );

    exp = ( -hydra.viscoPlasticParameters[ 3 ] * ( previousTemperature - hydra.viscoPlasticParameters[ 5 ] ) / ( hydra.viscoPlasticParameters[ 4 ] + previousTemperature - hydra.viscoPlasticParameters[ 5 ] ) );

    floatType answer2 = std::pow( 10, exp );

    BOOST_TEST( answer  == *R.get_plasticThermalMultiplier( ) );

    BOOST_TEST( answer2 == *R.get_previousPlasticThermalMultiplier( ) );

    BOOST_TEST( answer  == *Rjac.get_plasticThermalMultiplier( ) );

    BOOST_TEST( answer2 == *Rjac.get_previousPlasticThermalMultiplier( ) );

    // Check the Jacobians w.r.t. the temperature
    floatType eps = 1e-6;

    floatType dPlasticThermalMultiplierdT = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + delta, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature - delta, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        residualMock Rm( &hydram, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        floatType dp = *Rp.get_plasticThermalMultiplier( );

        floatType dm = *Rm.get_plasticThermalMultiplier( );

        dPlasticThermalMultiplierdT = ( dp - dm ) / ( 2 * delta );

    }
    BOOST_TEST( dPlasticThermalMultiplierdT == *Rjac.get_dPlasticThermalMultiplierdT( ) );

    floatType dPreviousPlasticThermalMultiplierdPreviousT = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( previousTemperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature + delta, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature - delta, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        residualMock Rm( &hydram, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        floatType dp = *Rp.get_previousPlasticThermalMultiplier( );

        floatType dm = *Rm.get_previousPlasticThermalMultiplier( );

        dPreviousPlasticThermalMultiplierdPreviousT += ( dp - dm ) / ( 2 * delta );

    }
    BOOST_TEST( dPreviousPlasticThermalMultiplierdPreviousT == *Rjac.get_dPreviousPlasticThermalMultiplierdPreviousT( ) );

}

BOOST_AUTO_TEST_CASE( test_residual_get_dragStress, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the drag stress
     */

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

            floatVector stateVariables = { 0.6 };

            floatVector previousStateVariables = { 0.9 };

        private:

            using tardigradeHydra::perzynaViscoplasticity::residual::setStateVariables;

            virtual void setStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_previousStateVariables( previousStateVariables );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_stateVariables( stateVariables );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 9, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType answer = hydra.viscoPlasticParameters[ 1 ] + R.stateVariables[ 0 ] * hydra.viscoPlasticParameters[ 2 ];

    floatType answer2 = hydra.viscoPlasticParameters[ 1 ] + R.previousStateVariables[ 0 ] * hydra.viscoPlasticParameters[ 2 ];

    BOOST_TEST( answer  == *R.get_dragStress( ) );

    BOOST_TEST( answer2 == *R.get_previousDragStress( ) );

}

BOOST_AUTO_TEST_CASE( test_residual_get_dragStress_jacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the drag stress
     */

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

            floatVector stateVariables = { 0.6 };

            floatVector previousStateVariables = { 0.9 };

        private:

            using tardigradeHydra::perzynaViscoplasticity::residual::setStateVariables;

            virtual void setStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_previousStateVariables( previousStateVariables );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_stateVariables( stateVariables );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};

            floatVector _local_stateVariableDelta = { 0 };                                                   

            floatVector _local_previousStateVariableDelta = { 0 };                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 9, 1, stateVariableIndices, viscoPlasticParameters );

                viscoPlasticity.stateVariables += _local_stateVariableDelta;

                viscoPlasticity.previousStateVariables += _local_previousStateVariableDelta;

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock Rjac( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType answer = hydra.viscoPlasticParameters[ 1 ] + R.stateVariables[ 0 ] * hydra.viscoPlasticParameters[ 2 ];

    floatType answer2 = hydra.viscoPlasticParameters[ 1 ] + R.previousStateVariables[ 0 ] * hydra.viscoPlasticParameters[ 2 ];

    Rjac.get_dDragStressdStateVariables( );

    Rjac.get_dPreviousDragStressdPreviousStateVariables( );

    BOOST_TEST( answer  == *R.get_dragStress( ) );

    BOOST_TEST( answer2 == *R.get_previousDragStress( ) );

    BOOST_TEST( answer  == *Rjac.get_dragStress( ) );

    BOOST_TEST( answer2 == *Rjac.get_previousDragStress( ) );

    floatType eps = 1e-6;

    floatVector dDragStressdStateVariables( R.stateVariables.size( ), 0 );

    for ( unsigned int i = 0; i < R.stateVariables.size( ); i++ ){

        floatVector delta = floatVector( R.stateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( R.stateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        residualMock Rm( &hydram, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        Rp.stateVariables += delta;

        Rm.stateVariables -= delta;

        floatType dp = *Rp.get_dragStress( );

        floatType dm = *Rm.get_dragStress( );

        dDragStressdStateVariables[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }
    BOOST_TEST( dDragStressdStateVariables == *Rjac.get_dDragStressdStateVariables( ), CHECK_PER_ELEMENT );

    floatVector dPreviousDragStressdPreviousStateVariables( R.previousStateVariables.size( ), 0 );

    for ( unsigned int i = 0; i < R.previousStateVariables.size( ); i++ ){

        floatVector delta = floatVector( R.previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( R.previousStateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        residualMock Rm( &hydram, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        Rp.previousStateVariables += delta;

        Rm.previousStateVariables -= delta;

        floatType dp = *Rp.get_previousDragStress( );

        floatType dm = *Rm.get_previousDragStress( );

        dPreviousDragStressdPreviousStateVariables[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }
    BOOST_TEST( dPreviousDragStressdPreviousStateVariables == *Rjac.get_dPreviousDragStressdPreviousStateVariables( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_get_hardeningFunction, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the hardening function
     */

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

            floatVector stateVariables = { 0.6 };

            floatVector previousStateVariables = { 0.9 };

        private:

            using tardigradeHydra::perzynaViscoplasticity::residual::setStateVariables;

            virtual void setStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_previousStateVariables( previousStateVariables );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_stateVariables( stateVariables );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 9, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock Rjac( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatVector answer = { hydra.viscoPlasticParameters[ 9 ] + R.stateVariables[ 0 ] * hydra.viscoPlasticParameters[ 10 ] };

    floatVector answer2 = { hydra.viscoPlasticParameters[ 9 ] + R.previousStateVariables[ 0 ] * hydra.viscoPlasticParameters[ 10 ] };

    Rjac.get_dHardeningFunctiondStateVariables( );

    Rjac.get_dPreviousHardeningFunctiondPreviousStateVariables( );

    BOOST_TEST( answer  == *R.get_hardeningFunction( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer2 == *R.get_previousHardeningFunction( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer  == *Rjac.get_hardeningFunction( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer2 == *Rjac.get_previousHardeningFunction( ), CHECK_PER_ELEMENT );

    floatType eps = 1e-6;

    floatVector dHardeningFunctiondStateVariables( R.stateVariables.size( ) * R.stateVariables.size( ), 0 );

    for ( unsigned int i = 0; i < R.stateVariables.size( ); i++ ){

        floatVector delta = floatVector( R.stateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( R.stateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        residualMock Rm( &hydram, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        Rp.stateVariables += delta;

        Rm.stateVariables -= delta;

        floatVector dp = *Rp.get_hardeningFunction( );

        floatVector dm = *Rm.get_hardeningFunction( );

        for ( unsigned int j = 0; j < R.stateVariables.size( ); j++ ){
            dHardeningFunctiondStateVariables[ R.stateVariables.size( ) * j + i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );
        }

    }
    BOOST_TEST( dHardeningFunctiondStateVariables == *Rjac.get_dHardeningFunctiondStateVariables( ), CHECK_PER_ELEMENT );

    floatVector dPreviousHardeningFunctiondPreviousStateVariables( R.previousStateVariables.size( ) * R.previousStateVariables.size( ), 0 );

    for ( unsigned int i = 0; i < R.previousStateVariables.size( ); i++ ){

        floatVector delta = floatVector( R.previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( R.previousStateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              { }, { },
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        residualMock Rm( &hydram, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        Rp.previousStateVariables += delta;

        Rm.previousStateVariables -= delta;

        floatVector dp = *Rp.get_previousHardeningFunction( );

        floatVector dm = *Rm.get_previousHardeningFunction( );

        for ( unsigned int j = 0; j < R.stateVariables.size( ); j++ ){
            dPreviousHardeningFunctiondPreviousStateVariables[ R.stateVariables.size( ) * j + i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );
        }

    }
    BOOST_TEST( dPreviousHardeningFunctiondPreviousStateVariables == *Rjac.get_dPreviousHardeningFunctiondPreviousStateVariables( ), CHECK_PER_ELEMENT );
}

BOOST_AUTO_TEST_CASE( test_residual_decomposeParameters, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of decomposing the parameter vector
     */

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

            floatVector stateVariables = { 0.6 };

        private:

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 9, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatVector answer_perzynaParameters = floatVector( hydra.viscoPlasticParameters.begin( ), hydra.viscoPlasticParameters.begin( ) + 1 );

    floatVector answer_dragStressParameters = floatVector( hydra.viscoPlasticParameters.begin( ) + 1, hydra.viscoPlasticParameters.begin( ) + 3 );

    floatVector answer_thermalParameters = floatVector( hydra.viscoPlasticParameters.begin( ) + 3, hydra.viscoPlasticParameters.begin( ) + 6 );

    floatVector answer_yieldParameters = floatVector( hydra.viscoPlasticParameters.begin( ) + 6, hydra.viscoPlasticParameters.begin( ) + 8 );

    floatVector answer_flowParameters = { 0, hydra.viscoPlasticParameters[ 8 ] };

    floatVector answer_hardeningParameters = { hydra.viscoPlasticParameters[ 9 ], hydra.viscoPlasticParameters[ 10 ] };

    BOOST_TEST( answer_perzynaParameters == *R.get_perzynaParameters( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer_dragStressParameters == *R.get_dragStressParameters( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer_thermalParameters == *R.get_thermalParameters( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer_yieldParameters == *R.get_yieldParameters( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer_flowParameters == *R.get_flowParameters( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer_hardeningParameters == *R.get_hardeningParameters( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_getStateVariableIndices, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test getting the state variable indices
     */

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

            floatVector stateVariables = { 0.6 };

        private:

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 0, 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 9, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    std::vector< unsigned int > answer = { 0, 2 };

    BOOST_TEST( answer == *R.getStateVariableIndices( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_get_stateVariables, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of decomposing the parameter vector
     */

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

        private:

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 0, 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 11, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 3 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatVector answer = { 4, 6 };

    BOOST_TEST( answer == *R.get_stateVariables( ), CHECK_PER_ELEMENT );

    floatVector answer2 = { 0.01, 0.03 };

    BOOST_TEST( answer2 == *R.get_previousStateVariables( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_get_plasticMuliplier, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the plastic multiplier
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

            floatType n = 3.2;

            floatType f = 2.34;

            floatType q = 3.67;

            floatType A = 0.45;

            floatType previousf = 1.34;

            floatType previousq = 7.67;

            floatType previousA = 0.85;

            floatVector dfdCauchy = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector dfdF = { 10, 11, 12, 13, 14, 15, 16, 17, 18 };

            floatVector dfdSubFs = { 19, 20, 21, 22, 23, 24, 25, 26, 27,
                                     28, 29, 30, 31, 32, 33, 34, 35, 36 };

            floatVector dfdStateVariables = { 19 };

            floatVector dqdXi = { 27 };

            floatType   dAdT  = 28;

            floatVector dPreviousfdPreviousCauchy = { 29, 30, 31, 32, 33, 34, 35, 36, 37 };

            floatVector dPreviousfdPreviousF = { 38, 39, 40, 41, 42, 43, 44, 45, 46 };

            floatVector dPreviousfdPreviousSubFs = { 47, 48, 49, 50, 51, 52, 53, 54, 55,
                                                     56, 57, 58, 59, 60, 61, 62, 63, 64 };

            floatVector dPreviousfdPreviousStateVariables = { 47 };

            floatVector dPreviousqdPreviousXi = { 65 };

            floatType   dPreviousAdPreviousT  = 66;
            

        private:


            using tardigradeHydra::perzynaViscoplasticity::residual::setDragStress;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdDragStressdStateVariables;
            using tardigradeHydra::perzynaViscoplasticity::residual::setPlasticThermalMultiplier;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdPlasticThermalMultiplierdT;
            using tardigradeHydra::perzynaViscoplasticity::residual::setYieldFunction;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdYieldFunctiondStateVariables;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdYieldFunctiondSubFs;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdYieldFunctiondF;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdYieldFunctiondCauchyStress;

            virtual void setYieldFunction( const bool isPrevious ) override{

                set_perzynaParameters( { n } );

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_previousYieldFunction( previousf );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_yieldFunction( f );

                }

            }

            virtual void setdYieldFunctiondCauchyStress( const bool isPrevious ) override{

                set_perzynaParameters( { n } );

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousYieldFunctiondPreviousCauchyStress( dPreviousfdPreviousCauchy );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dYieldFunctiondCauchyStress( dfdCauchy );

                }

            }

            virtual void setdYieldFunctiondF( const bool isPrevious ) override{

                set_perzynaParameters( { n } );

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousYieldFunctiondPreviousF( dPreviousfdPreviousF );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dYieldFunctiondF( dfdF );

                }

            }

            virtual void setdYieldFunctiondSubFs( const bool isPrevious ) override{

                set_perzynaParameters( { n } );

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousYieldFunctiondPreviousSubFs( dPreviousfdPreviousSubFs );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dYieldFunctiondSubFs( dfdSubFs );

                }

            }

            virtual void setdYieldFunctiondStateVariables( const bool isPrevious ) override{

                set_perzynaParameters( { n } );

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousYieldFunctiondPreviousStateVariables( dPreviousfdPreviousStateVariables );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dYieldFunctiondStateVariables( dfdStateVariables );

                }

            }

            virtual void setDragStress( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_previousDragStress( previousq );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dragStress( q );

                }

            }

            virtual void setdDragStressdStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousDragStressdPreviousStateVariables( dPreviousqdPreviousXi );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dDragStressdStateVariables( dqdXi );

                }

            }

            virtual void setPlasticThermalMultiplier( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_previousPlasticThermalMultiplier( previousA );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_plasticThermalMultiplier( A );

                }

            }

            virtual void setdPlasticThermalMultiplierdT( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousPlasticThermalMultiplierdPreviousT( dPreviousAdPreviousT );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPlasticThermalMultiplierdT( dAdT );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            //tardigradeHydra::linearElasticity::residual elasticity;
            stressMock elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = stressMock( this, 9 );//tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock Rjac( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType answer = R.A * std::pow( R.f / R.q, R.n );

    floatType answer2 = R.previousA * std::pow( R.previousf / R.previousq, R.n );

    BOOST_TEST( answer == *R.get_plasticMultiplier( ) );

    BOOST_TEST( answer2 == *R.get_previousPlasticMultiplier( ) );

    Rjac.get_dPlasticMultiplierdCauchyStress( );

    Rjac.get_dPreviousPlasticMultiplierdPreviousCauchyStress( );

    BOOST_TEST( answer == *Rjac.get_plasticMultiplier( ) );

    BOOST_TEST( answer2 == *Rjac.get_previousPlasticMultiplier( ) );

}

BOOST_AUTO_TEST_CASE( test_residual_get_plasticMultiplier_jacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the plastic multiplier's Jacobian
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
            stressMock elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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
                                           0.01, 0.02, 0.005, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 2, 1, 1, 7, 1, 6, 1, 8, 3,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 0.01, 7, 8 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R_ngrad( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    R.get_dYieldFunctiondCauchyStress( );

    BOOST_TEST( *R.get_yieldFunction( ) == *R_ngrad.get_yieldFunction( ) );

    // Test the jacobians
    floatType eps = 1e-6;

    floatVector dPlasticMultiplierdCauchyStress( 9, 0 );

    floatVector dPlasticMultiplierdF( 9, 0 );

    floatVector dPlasticMultiplierdSubFs( 18, 0 );

    floatType   dPlasticMultiplierdT = 0;

    floatVector dPlasticMultiplierdStateVariables( 1, 0 );

    floatVector dPreviousPlasticMultiplierdPreviousCauchyStress( 9, 0 );

    floatVector dPreviousPlasticMultiplierdPreviousF( 9, 0 );

    floatVector dPreviousPlasticMultiplierdPreviousSubFs( 18, 0 );

    floatType   dPreviousPlasticMultiplierdPreviousT = 0;

    floatVector dPreviousPlasticMultiplierdPreviousStateVariables( 1, 0 );

    // Jacobians w.r.t. the Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector + delta );

        floatType dp = *R.get_plasticMultiplier( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatType dm = *R.get_plasticMultiplier( );

        dPlasticMultiplierdCauchyStress[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( tolerantCheck( dPlasticMultiplierdCauchyStress, *R.get_dPlasticMultiplierdCauchyStress( ), 3e-5, 1e-5 ) );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.get_plasticMultiplier( );

        floatType dm = *Rm.get_plasticMultiplier( );

        dPlasticMultiplierdF[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( tolerantCheck( dPlasticMultiplierdF, *R.get_dPlasticMultiplierdF( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.get_plasticMultiplier( );

        floatType dm = *Rm.get_plasticMultiplier( );

        dPlasticMultiplierdSubFs[ i ] = ( dp - dm ) / ( 2 * delta[ i + 9 ] );

    }

    BOOST_TEST( tolerantCheck( dPlasticMultiplierdSubFs, *R.get_dPlasticMultiplierdSubFs( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature + delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature - delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.get_plasticMultiplier( );

        floatType dm = *Rm.get_plasticMultiplier( );

        dPlasticMultiplierdT = ( dp - dm ) / ( 2 * delta );

    }

    BOOST_TEST( dPlasticMultiplierdT == *R.get_dPlasticMultiplierdT( ) );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ 27 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( unknownVector[ 27 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.get_plasticMultiplier( );

        floatType dm = *Rm.get_plasticMultiplier( );

        dPlasticMultiplierdStateVariables[ i ] = ( dp - dm ) / ( 2 * delta[ 27 + hydra.stateVariableIndices[ i ] ] );

    }

    BOOST_TEST( dPlasticMultiplierdStateVariables == *R.get_dPlasticMultiplierdStateVariables( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.get_previousPlasticMultiplier( );

        floatType dm = *Rm.get_previousPlasticMultiplier( );

        dPreviousPlasticMultiplierdPreviousCauchyStress[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( dPreviousPlasticMultiplierdPreviousCauchyStress == *R.get_dPreviousPlasticMultiplierdPreviousCauchyStress( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.get_previousPlasticMultiplier( );

        floatType dm = *Rm.get_previousPlasticMultiplier( );

        dPreviousPlasticMultiplierdPreviousF[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( dPreviousPlasticMultiplierdPreviousF == *R.get_dPreviousPlasticMultiplierdPreviousF( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.get_previousPlasticMultiplier( );

        floatType dm = *Rm.get_previousPlasticMultiplier( );

        dPreviousPlasticMultiplierdPreviousSubFs[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( dPreviousPlasticMultiplierdPreviousSubFs == *R.get_dPreviousPlasticMultiplierdPreviousSubFs( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( previousTemperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature + delta, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature - delta, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.get_previousPlasticMultiplier( );

        floatType dm = *Rm.get_previousPlasticMultiplier( );

        dPreviousPlasticMultiplierdPreviousT = ( dp - dm ) / ( 2 * delta );

    }

    BOOST_TEST( dPreviousPlasticMultiplierdPreviousT == *R.get_dPreviousPlasticMultiplierdPreviousT( ) );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ 18 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( previousStateVariables[ 18 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.get_previousPlasticMultiplier( );

        floatType dm = *Rm.get_previousPlasticMultiplier( );

        dPreviousPlasticMultiplierdPreviousStateVariables[ i ] = ( dp - dm ) / ( 2 * delta[ 18 + hydra.stateVariableIndices[ i ] ] );

    }

    BOOST_TEST( dPreviousPlasticMultiplierdPreviousStateVariables == *R.get_dPreviousPlasticMultiplierdPreviousStateVariables( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_get_velocityGradient, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the velocity gradient
     */

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

            floatType gamma = 2.4;

            floatVector dGammadCauchy = initializeVector( 9, -0.56 );

            floatVector dGammadF = initializeVector( 9, 1 );

            floatVector dGammadSubFs = initializeVector( 18, 10 );

            floatType   dGammadT = 2.4;

            floatVector dGammadXi = { 1.1 };

            floatVector nhat = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector dNhatdCauchy = initializeVector( 81, 0.1 );

            floatVector dNhatdF = initializeVector( 81, 0.3 );

            floatVector dNhatdSubFs = initializeVector( 162, -0.54 );

            floatType previousGamma = 5.6;

            floatVector dPreviousGammadPreviousCauchy = initializeVector( 9, -3 );

            floatVector dPreviousGammadPreviousF = initializeVector( 9, 0.123 );

            floatVector dPreviousGammadPreviousSubFs = initializeVector( 18, -0.224 );

            floatType   dPreviousGammadPreviousT = 0.872;

            floatVector dPreviousGammadPreviousXi = { 2.5 };

            floatVector previousNhat = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            floatVector dPreviousNhatdPreviousCauchy = initializeVector( 81, -0.37 );

            floatVector dPreviousNhatdPreviousF = initializeVector( 81, 0.373 );

            floatVector dPreviousNhatdPreviousSubFs = initializeVector( 162, -0.82 );

        private:

            floatVector initializeVector( int nvals, floatType val_0 ){

                floatVector value( nvals, 0 );

                for ( int i = 0; i < nvals; i++ ){

                    value[ i ] = i + val_0;

                }

                return value;

            }

            floatMatrix initializeMatrix( int nrows, int ncols, floatType val_0 ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                for ( int i = 0; i < nrows; i++ ){

                    for ( int j = 0; j < ncols; j++ ){

                        value[ i ][ j ] = val_0 += nrows * i + j + val_0;

                    }

                }

                return value;

            }

            using tardigradeHydra::perzynaViscoplasticity::residual::setPlasticMultiplier;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdPlasticMultiplierdCauchyStress;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdPlasticMultiplierdF;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdPlasticMultiplierdSubFs;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdPlasticMultiplierdT;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdPlasticMultiplierdStateVariables;

            using tardigradeHydra::perzynaViscoplasticity::residual::setFlowDirection;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdFlowDirectiondCauchyStress;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdFlowDirectiondF;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdFlowDirectiondSubFs;

            virtual void setPlasticMultiplier( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_previousPlasticMultiplier( previousGamma );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_plasticMultiplier( gamma );

                }

            }

            virtual void setdPlasticMultiplierdCauchyStress( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousPlasticMultiplierdPreviousCauchyStress( dPreviousGammadPreviousCauchy );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPlasticMultiplierdCauchyStress( dGammadCauchy );

                }

            }

            virtual void setdPlasticMultiplierdF( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousPlasticMultiplierdPreviousF( dPreviousGammadPreviousF );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPlasticMultiplierdF( dGammadF );

                }

            }

            virtual void setdPlasticMultiplierdSubFs( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousPlasticMultiplierdPreviousSubFs( dPreviousGammadPreviousSubFs );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPlasticMultiplierdSubFs( dGammadSubFs );

                }

            }

            virtual void setdPlasticMultiplierdT( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousPlasticMultiplierdPreviousT( dPreviousGammadPreviousT );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPlasticMultiplierdT( dGammadT );

                }

            }

            virtual void setdPlasticMultiplierdStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousPlasticMultiplierdPreviousStateVariables( dPreviousGammadPreviousXi );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPlasticMultiplierdStateVariables( dGammadXi );

                }

            }

            virtual void setFlowDirection( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_previousFlowDirection( previousNhat );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_flowDirection( nhat );

                }

            }

            virtual void setdFlowDirectiondCauchyStress( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousFlowDirectiondPreviousCauchyStress( dPreviousNhatdPreviousCauchy );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dFlowDirectiondCauchyStress( dNhatdCauchy );

                }

            }

            virtual void setdFlowDirectiondF( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousFlowDirectiondPreviousF( dPreviousNhatdPreviousF );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dFlowDirectiondF( dNhatdF );

                }

            }

            virtual void setdFlowDirectiondSubFs( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousFlowDirectiondPreviousSubFs( dPreviousNhatdPreviousSubFs );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dFlowDirectiondSubFs( dNhatdSubFs );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 9, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock Rjac( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatVector answer = R.gamma * R.nhat;

    BOOST_TEST( answer == *R.get_velocityGradient( ), CHECK_PER_ELEMENT );

    floatVector answer2 = R.previousGamma * R.previousNhat;

    BOOST_TEST( answer2 == *R.get_previousVelocityGradient( ), CHECK_PER_ELEMENT );

    Rjac.get_dVelocityGradientdCauchyStress( );

    Rjac.get_dPreviousVelocityGradientdPreviousCauchyStress( );

    BOOST_TEST( answer == *Rjac.get_velocityGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer2 == *Rjac.get_previousVelocityGradient( ), CHECK_PER_ELEMENT );


}

BOOST_AUTO_TEST_CASE( test_residual_get_velocityGradient_jacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the velocity gradient's Jacobian
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
            stressMock elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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
                                           0.01, 0.02, 0.005, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 2, 1, 1, 7, 1, 6, 1, 8, 3,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 0.01, 7, 8 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R_ngrad( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    R.get_dYieldFunctiondCauchyStress( );

    BOOST_TEST( *R.get_yieldFunction( ) == *R_ngrad.get_yieldFunction( ) );

    // Test the jacobians
    floatType eps = 1e-6;

    floatMatrix dVelocityGradientdCauchyStress( 9, floatVector( 9, 0 ) );

    floatMatrix dVelocityGradientdF( 9, floatVector( 9, 0 ) );

    floatMatrix dVelocityGradientdSubFs( 9, floatVector( 18, 0 ) );

    floatVector dVelocityGradientdT( 9, 0 );

    floatMatrix dVelocityGradientdStateVariables( 9, floatVector( 1, 0 ) );

    floatMatrix dPreviousVelocityGradientdPreviousCauchyStress( 9, floatVector( 9, 0 ) );

    floatMatrix dPreviousVelocityGradientdPreviousF( 9, floatVector( 9, 0 ) );

    floatMatrix dPreviousVelocityGradientdPreviousSubFs( 9, floatVector( 18, 0 ) );

    floatVector dPreviousVelocityGradientdPreviousT( 9, 0 );

    floatMatrix dPreviousVelocityGradientdPreviousStateVariables( 9, floatVector( 1, 0 ) );

    // Jacobians w.r.t. the Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector + delta );

        floatVector dp = *R.get_velocityGradient( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatVector dm = *R.get_velocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dVelocityGradientdCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dVelocityGradientdCauchyStress ), *R.get_dVelocityGradientdCauchyStress( ), 1e-4, 1e-5 ) );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_velocityGradient( );

        floatVector dm = *Rm.get_velocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dVelocityGradientdF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dVelocityGradientdF ), *R.get_dVelocityGradientdF( ), 5e-5, 1e-5 ) );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_velocityGradient( );

        floatVector dm = *Rm.get_velocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dVelocityGradientdSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dVelocityGradientdSubFs ), *R.get_dVelocityGradientdSubFs( ), 2e-5, 1e-5 ) );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature + delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature - delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_velocityGradient( );

        floatVector dm = *Rm.get_velocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dVelocityGradientdT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( tolerantCheck( dVelocityGradientdT, *R.get_dVelocityGradientdT( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ 27 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( unknownVector[ 27 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_velocityGradient( );

        floatVector dm = *Rm.get_velocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dVelocityGradientdStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 27 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dVelocityGradientdStateVariables ), *R.get_dVelocityGradientdStateVariables( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousVelocityGradient( );

        floatVector dm = *Rm.get_previousVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousVelocityGradientdPreviousCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousVelocityGradientdPreviousCauchyStress ) == *R.get_dPreviousVelocityGradientdPreviousCauchyStress( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousVelocityGradient( );

        floatVector dm = *Rm.get_previousVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousVelocityGradientdPreviousF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousVelocityGradientdPreviousF ) == *R.get_dPreviousVelocityGradientdPreviousF( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousVelocityGradient( );

        floatVector dm = *Rm.get_previousVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousVelocityGradientdPreviousSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousVelocityGradientdPreviousSubFs ) == *R.get_dPreviousVelocityGradientdPreviousSubFs( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( previousTemperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature + delta, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature - delta, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousVelocityGradient( );

        floatVector dm = *Rm.get_previousVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousVelocityGradientdPreviousT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dPreviousVelocityGradientdPreviousT == *R.get_dPreviousVelocityGradientdPreviousT( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ 18 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( previousStateVariables[ 18 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousVelocityGradient( );

        floatVector dm = *Rm.get_previousVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousVelocityGradientdPreviousStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 18 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousVelocityGradientdPreviousStateVariables ) == *R.get_dPreviousVelocityGradientdPreviousStateVariables( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_get_stateVariableEvolutionRate, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the velocity gradient
     */

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

            floatType gamma = 2.4;

            floatVector dGammadCauchy = initializeVector( 9, -0.56 );

            floatVector dGammadF = initializeVector( 9, 1 );

            floatVector dGammadSubFs = initializeVector( 18, 10 );

            floatType   dGammadT = 2.4;

            floatVector dGammadXi = { 1.1 };

            floatVector hardeningFunction = { 7.8 };

            floatVector dHardeningFunctiondXi = initializeVector( 1, -0.72 );

            floatType previousGamma = 5.6;

            floatVector dPreviousGammadPreviousCauchy = initializeVector( 9, -3 );

            floatVector dPreviousGammadPreviousF = initializeVector( 9, 0.123 );

            floatVector dPreviousGammadPreviousSubFs = initializeVector( 18, -0.224 );

            floatType   dPreviousGammadPreviousT = 0.872;

            floatVector dPreviousGammadPreviousXi = { 2.5 };

            floatVector previousHardeningFunction = { .56 };

            floatVector dPreviousHardeningFunctiondPreviousXi = initializeVector( 1, 0.22 );

        private:

            floatVector initializeVector( int nvals, floatType val_0 ){

                floatVector value( nvals, 0 );

                for ( int i = 0; i < nvals; i++ ){

                    value[ i ] = i + val_0;

                }

                return value;

            }

            floatMatrix initializeMatrix( int nrows, int ncols, floatType val_0 ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                for ( int i = 0; i < nrows; i++ ){

                    for ( int j = 0; j < ncols; j++ ){

                        value[ i ][ j ] = val_0 += nrows * i + j + val_0;

                    }

                }

                return value;

            }
        private:

            using tardigradeHydra::perzynaViscoplasticity::residual::setPlasticMultiplier;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdPlasticMultiplierdCauchyStress;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdPlasticMultiplierdF;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdPlasticMultiplierdSubFs;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdPlasticMultiplierdT;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdPlasticMultiplierdStateVariables;

            using tardigradeHydra::perzynaViscoplasticity::residual::setHardeningFunction;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdHardeningFunctiondStateVariables;

            virtual void setPlasticMultiplier( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_previousPlasticMultiplier( previousGamma );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_plasticMultiplier( gamma );

                }

            }

            virtual void setdPlasticMultiplierdCauchyStress( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousPlasticMultiplierdPreviousCauchyStress( dPreviousGammadPreviousCauchy );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPlasticMultiplierdCauchyStress( dGammadCauchy );

                }

            }

            virtual void setdPlasticMultiplierdF( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousPlasticMultiplierdPreviousF( dPreviousGammadPreviousF );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPlasticMultiplierdF( dGammadF );

                }

            }

            virtual void setdPlasticMultiplierdSubFs( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousPlasticMultiplierdPreviousSubFs( dPreviousGammadPreviousSubFs );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPlasticMultiplierdSubFs( dGammadSubFs );

                }

            }

            virtual void setdPlasticMultiplierdT( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousPlasticMultiplierdPreviousT( dPreviousGammadPreviousT );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPlasticMultiplierdT( dGammadT );

                }

            }

            virtual void setdPlasticMultiplierdStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousPlasticMultiplierdPreviousStateVariables( dPreviousGammadPreviousXi );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPlasticMultiplierdStateVariables( dGammadXi );

                }

            }

            virtual void setHardeningFunction( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_previousHardeningFunction( previousHardeningFunction );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_hardeningFunction( hardeningFunction );

                }

            }

            virtual void setdHardeningFunctiondStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousHardeningFunctiondPreviousStateVariables( dPreviousHardeningFunctiondPreviousXi );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dHardeningFunctiondStateVariables( dHardeningFunctiondXi );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 9, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock Rjac( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatVector answer = { R.gamma * R.hardeningFunction };

    BOOST_TEST( answer == *R.get_stateVariableEvolutionRates( ), CHECK_PER_ELEMENT );

    floatVector answer2 = { R.previousGamma * R.previousHardeningFunction };

    BOOST_TEST( answer2 == *R.get_previousStateVariableEvolutionRates( ), CHECK_PER_ELEMENT );

    Rjac.get_dStateVariableEvolutionRatesdCauchyStress( );

    Rjac.get_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( );

    BOOST_TEST( answer == *Rjac.get_stateVariableEvolutionRates( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer2 == *Rjac.get_previousStateVariableEvolutionRates( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_get_stateVariableEvolutionRates_jacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the velocity gradient's Jacobian
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
            stressMock elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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
                                           0.01, 0.02, 0.005, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 2, 1, 1, 7, 1, 6, 1, 8, 3,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 0.01, 7, 8 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R_ngrad( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    R.get_dYieldFunctiondCauchyStress( );

    BOOST_TEST( *R.get_yieldFunction( ) == *R_ngrad.get_yieldFunction( ) );

    // Test the jacobians
    floatType eps = 1e-6;

    unsigned int nvar = 1;

    floatMatrix dStateVariableEvolutionRatesdCauchyStress( nvar, floatVector( 9, 0 ) );

    floatMatrix dStateVariableEvolutionRatesdF( nvar, floatVector( 9, 0 ) );

    floatMatrix dStateVariableEvolutionRatesdSubFs( nvar, floatVector( 18, 0 ) );

    floatVector dStateVariableEvolutionRatesdT( nvar, 0 );

    floatMatrix dStateVariableEvolutionRatesdStateVariables( nvar, floatVector( 1, 0 ) );

    floatMatrix dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( nvar, floatVector( 9, 0 ) );

    floatMatrix dPreviousStateVariableEvolutionRatesdPreviousF( nvar, floatVector( 9, 0 ) );

    floatMatrix dPreviousStateVariableEvolutionRatesdPreviousSubFs( nvar, floatVector( 18, 0 ) );

    floatVector dPreviousStateVariableEvolutionRatesdPreviousT( nvar, 0 );

    floatMatrix dPreviousStateVariableEvolutionRatesdPreviousStateVariables( nvar, floatVector( 1, 0 ) );

    // Jacobians w.r.t. the Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector + delta );

        floatVector dp = *R.get_stateVariableEvolutionRates( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatVector dm = *R.get_stateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dStateVariableEvolutionRatesdCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dStateVariableEvolutionRatesdCauchyStress ), *R.get_dStateVariableEvolutionRatesdCauchyStress( ), 3e-5, 1e-5 ) );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_stateVariableEvolutionRates( );

        floatVector dm = *Rm.get_stateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dStateVariableEvolutionRatesdF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dStateVariableEvolutionRatesdF ), *R.get_dStateVariableEvolutionRatesdF( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_stateVariableEvolutionRates( );

        floatVector dm = *Rm.get_stateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dStateVariableEvolutionRatesdSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dStateVariableEvolutionRatesdSubFs ), *R.get_dStateVariableEvolutionRatesdSubFs( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature + delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature - delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_stateVariableEvolutionRates( );

        floatVector dm = *Rm.get_stateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dStateVariableEvolutionRatesdT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dStateVariableEvolutionRatesdT == *R.get_dStateVariableEvolutionRatesdT( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ 27 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( unknownVector[ 27 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_stateVariableEvolutionRates( );

        floatVector dm = *Rm.get_stateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dStateVariableEvolutionRatesdStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 27 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dStateVariableEvolutionRatesdStateVariables ) == *R.get_dStateVariableEvolutionRatesdStateVariables( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousStateVariableEvolutionRates( );

        floatVector dm = *Rm.get_previousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousStateVariableEvolutionRatesdPreviousCauchyStress ) == *R.get_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousStateVariableEvolutionRates( );

        floatVector dm = *Rm.get_previousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousStateVariableEvolutionRatesdPreviousF ) == *R.get_dPreviousStateVariableEvolutionRatesdPreviousF( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousStateVariableEvolutionRates( );

        floatVector dm = *Rm.get_previousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousStateVariableEvolutionRatesdPreviousSubFs ) == *R.get_dPreviousStateVariableEvolutionRatesdPreviousSubFs( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( previousTemperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature + delta, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature - delta, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousStateVariableEvolutionRates( );

        floatVector dm = *Rm.get_previousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dPreviousStateVariableEvolutionRatesdPreviousT == *R.get_dPreviousStateVariableEvolutionRatesdPreviousT( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ 18 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( previousStateVariables[ 18 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_previousStateVariableEvolutionRates( );

        floatVector dm = *Rm.get_previousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 18 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPreviousStateVariableEvolutionRatesdPreviousStateVariables ) == *R.get_dPreviousStateVariableEvolutionRatesdPreviousStateVariables( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_get_plasticDeformationGradient, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the plastic deformation gradient
     */

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

            floatVector velocityGradient = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousVelocityGradient = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            floatVector dLdCauchy = initializeVector( 81, 0.1 );

            floatVector dLdF      = initializeVector( 81, -0.34 );

            floatVector dLdSubFs  = initializeVector( 162, 0.23 );

            floatVector dLdT      = initializeVector( 9, 5.67 );

            floatVector dLdXi     = initializeVector( 9, -0.26 );

            floatVector dPreviousLdPreviousCauchy = initializeVector( 81, 0.27 );

            floatVector dPreviousLdPreviousF      = initializeVector( 81, -0.44 );

            floatVector dPreviousLdPreviousSubFs  = initializeVector( 162, 0.29 );

            floatVector dPreviousLdPreviousT      = initializeVector( 9, 2.11 );

            floatVector dPreviousLdPreviousXi     = initializeVector( 9, -0.19 );

        private:

            using tardigradeHydra::perzynaViscoplasticity::residual::setVelocityGradient;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdVelocityGradientdCauchyStress;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdVelocityGradientdF;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdVelocityGradientdSubFs;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdVelocityGradientdT;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdVelocityGradientdStateVariables;

            floatVector initializeVector( int nvals, floatType val_0 ){

                floatVector value( nvals, 0 );

                for ( int i = 0; i < nvals; i++ ){

                    value[ i ] = i + val_0;

                }

                return value;

            }

            floatMatrix initializeMatrix( int nrows, int ncols, floatType val_0 ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                for ( int i = 0; i < nrows; i++ ){

                    for ( int j = 0; j < ncols; j++ ){

                        value[ i ][ j ] = val_0 += nrows * i + j + val_0;

                    }

                }

                return value;

            }

            virtual void setVelocityGradient( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_previousVelocityGradient( previousVelocityGradient );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_velocityGradient( velocityGradient );

                }

            }

            virtual void setdVelocityGradientdCauchyStress( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousVelocityGradientdPreviousCauchyStress( dPreviousLdPreviousCauchy );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dVelocityGradientdCauchyStress( dLdCauchy );

                }

            }

            virtual void setdVelocityGradientdF( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousVelocityGradientdPreviousF( dPreviousLdPreviousF );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dVelocityGradientdF( dLdF );

                }

            }

            virtual void setdVelocityGradientdSubFs( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousVelocityGradientdPreviousSubFs( dPreviousLdPreviousSubFs );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dVelocityGradientdSubFs( dLdSubFs );

                }

            }

            virtual void setdVelocityGradientdT( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousVelocityGradientdPreviousT( dPreviousLdPreviousT );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dVelocityGradientdT( dLdT );

                }

            }

            virtual void setdVelocityGradientdStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousVelocityGradientdPreviousStateVariables( dPreviousLdPreviousXi );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dVelocityGradientdStateVariables( dLdXi );

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    floatVector previousStateVariables = { 0, -0.1, 0, 0, 0, 0.3, 0, 0, 0,
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

    residualMock R( &hydra, 10, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters, 0.6 );

    residualMock Rjac1( &hydra, 10, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters, 0.6 );

    residualMock Rjac2( &hydra, 10, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters, 0.6 );

    floatVector answer = { 0.41617760, -0.32005614, -0.05658586,
                          -0.33457515,  0.67081398, -0.19950498,
                          -0.08532789, -0.43831590,  0.05757591 };

    BOOST_TEST( answer == *R.get_plasticDeformationGradient( ), CHECK_PER_ELEMENT );

    Rjac1.get_dPlasticDeformationGradientdCauchyStress( );

    Rjac2.get_dPlasticDeformationGradientdPreviousCauchyStress( );

    BOOST_TEST( answer == *Rjac1.get_plasticDeformationGradient( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *Rjac2.get_plasticDeformationGradient( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_get_plasticDeformationGradient_jacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the plastic deformation gradient
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
            stressMock elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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
                                           0.01, 0.02, 0.005, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 2, 1, 1, 7, 1, 6, 1, 8, 3,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 0.01, 7, 8 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R_ngrad( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    R.get_dYieldFunctiondCauchyStress( );

    BOOST_TEST( *R.get_yieldFunction( ) == *R_ngrad.get_yieldFunction( ) );

    // Test the jacobians
    floatType eps = 1e-6;

    unsigned int nvar = 9;

    floatMatrix dPlasticDeformationGradientdCauchyStress( nvar, floatVector( 9, 0 ) );

    floatMatrix dPlasticDeformationGradientdF( nvar, floatVector( 9, 0 ) );

    floatMatrix dPlasticDeformationGradientdSubFs( nvar, floatVector( 18, 0 ) );

    floatVector dPlasticDeformationGradientdT( nvar, 0 );

    floatMatrix dPlasticDeformationGradientdStateVariables( nvar, floatVector( 1, 0 ) );

    floatMatrix dPlasticDeformationGradientdPreviousCauchyStress( nvar, floatVector( 9, 0 ) );

    floatMatrix dPlasticDeformationGradientdPreviousF( nvar, floatVector( 9, 0 ) );

    floatMatrix dPlasticDeformationGradientdPreviousSubFs( nvar, floatVector( 18, 0 ) );

    floatVector dPlasticDeformationGradientdPreviousT( nvar, 0 );

    floatMatrix dPlasticDeformationGradientdPreviousStateVariables( nvar, floatVector( 1, 0 ) );

    // Jacobians w.r.t. the Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector + delta );

        floatVector dp = *R.get_plasticDeformationGradient( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatVector dm = *R.get_plasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPlasticDeformationGradientdCauchyStress ), *R.get_dPlasticDeformationGradientdCauchyStress( ), 3e-5, 1e-5 ) );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticDeformationGradient( );

        floatVector dm = *Rm.get_plasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPlasticDeformationGradientdF ), *R.get_dPlasticDeformationGradientdF( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticDeformationGradient( );

        floatVector dm = *Rm.get_plasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPlasticDeformationGradientdSubFs ), *R.get_dPlasticDeformationGradientdSubFs( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature + delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature - delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticDeformationGradient( );

        floatVector dm = *Rm.get_plasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dPlasticDeformationGradientdT == *R.get_dPlasticDeformationGradientdT( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ 27 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( unknownVector[ 27 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticDeformationGradient( );

        floatVector dm = *Rm.get_plasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 27 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPlasticDeformationGradientdStateVariables ) == *R.get_dPlasticDeformationGradientdStateVariables( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticDeformationGradient( );

        floatVector dm = *Rm.get_plasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdPreviousCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPlasticDeformationGradientdPreviousCauchyStress ), *R.get_dPlasticDeformationGradientdPreviousCauchyStress( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticDeformationGradient( );

        floatVector dm = *Rm.get_plasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdPreviousF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPlasticDeformationGradientdPreviousF ), *R.get_dPlasticDeformationGradientdPreviousF( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticDeformationGradient( );

        floatVector dm = *Rm.get_plasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdPreviousSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPlasticDeformationGradientdPreviousSubFs ), *R.get_dPlasticDeformationGradientdPreviousSubFs( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( previousTemperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature + delta, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature - delta, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticDeformationGradient( );

        floatVector dm = *Rm.get_plasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdPreviousT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dPlasticDeformationGradientdPreviousT == *R.get_dPlasticDeformationGradientdPreviousT( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ 18 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( previousStateVariables[ 18 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticDeformationGradient( );

        floatVector dm = *Rm.get_plasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdPreviousStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 18 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPlasticDeformationGradientdPreviousStateVariables ) == *R.get_dPlasticDeformationGradientdPreviousStateVariables( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_get_plasticStateVariables, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the plastic state variables
     */

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

            floatVector stateVariableEvolutionRates = { 2.3, 1.78 };

            floatVector previousStateVariableEvolutionRates = { 4.7, 0.3 };

            floatVector previousStateVariables = { 0.3, -0.4 };

            floatVector dXidotdC = initializeVector( 18, 0.1 );

            floatVector dXidotdF = initializeVector( 18, 0.2 );

            floatVector dXidotdSubFs = initializeVector( 36, 0.3 );

            floatVector dXidotdT = initializeVector( 2, 0.4 );

            floatVector dXidotdXi = initializeVector( 4, 0.5 );

            floatVector dPreviousXidotdPreviousC = initializeVector( 18, 0.6 );

            floatVector dPreviousXidotdPreviousF = initializeVector( 18, 0.7 );

            floatVector dPreviousXidotdPreviousSubFs = initializeVector( 36, 0.8 );

            floatVector dPreviousXidotdPreviousT = initializeVector( 2, 0.9 );

            floatVector dPreviousXidotdPreviousXi = initializeVector( 4, 1.0 );

        private:

            floatVector initializeVector( int nvals, floatType val_0 ){

                floatVector value( nvals, 0 );

                for ( int i = 0; i < nvals; i++ ){

                    value[ i ] = i + val_0;

                }

                return value;

            }

            floatMatrix initializeMatrix( int nrows, int ncols, floatType val_0 ){

                floatMatrix value( nrows, floatVector( ncols, 0 ) );

                for ( int i = 0; i < nrows; i++ ){

                    for ( int j = 0; j < ncols; j++ ){

                        value[ i ][ j ] = val_0 += nrows * i + j + val_0;

                    }

                }

                return value;

            }

            using tardigradeHydra::perzynaViscoplasticity::residual::setStateVariableEvolutionRates;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdStateVariableEvolutionRatesdCauchyStress;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdStateVariableEvolutionRatesdF;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdStateVariableEvolutionRatesdSubFs;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdStateVariableEvolutionRatesdT;
            using tardigradeHydra::perzynaViscoplasticity::residual::setdStateVariableEvolutionRatesdStateVariables;

            virtual void setStateVariableEvolutionRates( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_previousStateVariableEvolutionRates( previousStateVariableEvolutionRates );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_stateVariableEvolutionRates( stateVariableEvolutionRates );

                }

            }

            virtual void setdStateVariableEvolutionRatesdCauchyStress( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousStateVariableEvolutionRatesdPreviousCauchyStress( dPreviousXidotdPreviousC );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dStateVariableEvolutionRatesdCauchyStress( dXidotdC );

                }

            }

            virtual void setdStateVariableEvolutionRatesdF( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousStateVariableEvolutionRatesdPreviousF( dPreviousXidotdPreviousF );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dStateVariableEvolutionRatesdF( dXidotdF );

                }

            }

            virtual void setdStateVariableEvolutionRatesdSubFs( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousStateVariableEvolutionRatesdPreviousSubFs( dPreviousXidotdPreviousSubFs );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dStateVariableEvolutionRatesdSubFs( dXidotdSubFs );

                }

            }

            virtual void setdStateVariableEvolutionRatesdT( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousStateVariableEvolutionRatesdPreviousT( dPreviousXidotdPreviousT );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dStateVariableEvolutionRatesdT( dXidotdT );

                }

            }

            virtual void setdStateVariableEvolutionRatesdStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dPreviousStateVariableEvolutionRatesdPreviousStateVariables( dPreviousXidotdPreviousXi );

                }
                else{

                    tardigradeHydra::perzynaViscoplasticity::residual::set_dStateVariableEvolutionRatesdStateVariables( dXidotdXi );

                }

            }

            virtual void setPreviousStateVariables( ) override{

                tardigradeHydra::perzynaViscoplasticity::residual::set_previousStateVariables( previousStateVariables );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    floatVector previousStateVariables = { 0, -0.1, 0, 0, 0, 0.3, 0, 0, 0,
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

    residualMock R( &hydra, 10, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters, 0.6 );

    residualMock R1( &hydra, 10, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters, 0.6 );

    residualMock R2( &hydra, 10, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters, 0.6 );

    floatVector answer = { 7.472, 2.2136 };

    BOOST_TEST( answer == *R.get_plasticStateVariables( ), CHECK_PER_ELEMENT );

    R1.get_dPlasticStateVariablesdCauchyStress( );

    R2.get_dPlasticStateVariablesdPreviousCauchyStress( );

    BOOST_TEST( answer == *R1.get_plasticStateVariables( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answer == *R2.get_plasticStateVariables( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_get_plasticStateVariables_jacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the plastic deformation gradient
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
            stressMock elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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
                                           0.01, 0.02, 0.005, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 2, 1, 1, 7, 1, 6, 1, 8, 3,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 0.01, 7, 8 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R_ngrad( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    R.get_dYieldFunctiondCauchyStress( );

    BOOST_TEST( *R.get_yieldFunction( ) == *R_ngrad.get_yieldFunction( ) );

    // Test the jacobians
    floatType eps = 1e-6;

    unsigned int nvar = 1;

    floatMatrix dPlasticStateVariablesdCauchyStress( nvar, floatVector( 9, 0 ) );

    floatMatrix dPlasticStateVariablesdF( nvar, floatVector( 9, 0 ) );

    floatMatrix dPlasticStateVariablesdSubFs( nvar, floatVector( 18, 0 ) );

    floatVector dPlasticStateVariablesdT( nvar, 0 );

    floatMatrix dPlasticStateVariablesdStateVariables( nvar, floatVector( 1, 0 ) );

    floatMatrix dPlasticStateVariablesdPreviousCauchyStress( nvar, floatVector( 9, 0 ) );

    floatMatrix dPlasticStateVariablesdPreviousF( nvar, floatVector( 9, 0 ) );

    floatMatrix dPlasticStateVariablesdPreviousSubFs( nvar, floatVector( 18, 0 ) );

    floatVector dPlasticStateVariablesdPreviousT( nvar, 0 );

    floatMatrix dPlasticStateVariablesdPreviousStateVariables( nvar, floatVector( 1, 0 ) );

    // Jacobians w.r.t. the Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector + delta );

        floatVector dp = *R.get_plasticStateVariables( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatVector dm = *R.get_plasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPlasticStateVariablesdCauchyStress ), *R.get_dPlasticStateVariablesdCauchyStress( ), 3e-5, 1e-5 ) );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticStateVariables( );

        floatVector dm = *Rm.get_plasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPlasticStateVariablesdF ), *R.get_dPlasticStateVariablesdF( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticStateVariables( );

        floatVector dm = *Rm.get_plasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_TEST( tolerantCheck( tardigradeVectorTools::appendVectors( dPlasticStateVariablesdSubFs ), *R.get_dPlasticStateVariablesdSubFs( ), 1e-5, 1e-5 ) );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature + delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature - delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticStateVariables( );

        floatVector dm = *Rm.get_plasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dPlasticStateVariablesdT == *R.get_dPlasticStateVariablesdT( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ 27 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( unknownVector[ 27 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticStateVariables( );

        floatVector dm = *Rm.get_plasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 27 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPlasticStateVariablesdStateVariables ) == *R.get_dPlasticStateVariablesdStateVariables( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticStateVariables( );

        floatVector dm = *Rm.get_plasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdPreviousCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPlasticStateVariablesdPreviousCauchyStress ) == *R.get_dPlasticStateVariablesdPreviousCauchyStress( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticStateVariables( );

        floatVector dm = *Rm.get_plasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdPreviousF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPlasticStateVariablesdPreviousF ) == *R.get_dPlasticStateVariablesdPreviousF( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticStateVariables( );

        floatVector dm = *Rm.get_plasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdPreviousSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPlasticStateVariablesdPreviousSubFs ) == *R.get_dPlasticStateVariablesdPreviousSubFs( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( previousTemperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature + delta, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature - delta, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticStateVariables( );

        floatVector dm = *Rm.get_plasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdPreviousT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_TEST( dPlasticStateVariablesdPreviousT == *R.get_dPlasticStateVariablesdPreviousT( ), CHECK_PER_ELEMENT );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ 18 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( previousStateVariables[ 18 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.get_plasticStateVariables( );

        floatVector dm = *Rm.get_plasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdPreviousStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 18 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dPlasticStateVariablesdPreviousStateVariables ) == *R.get_dPlasticStateVariablesdPreviousStateVariables( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_getResidual, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the plastic state variables
     */

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

            floatVector plasticDeformationGradient = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector plasticStateVariables = { 0.9 };

        private:

            virtual void setPlasticDeformationGradient( ) override{

                tardigradeHydra::perzynaViscoplasticity::residual::set_plasticDeformationGradient( plasticDeformationGradient );

            }

            virtual void setPlasticStateVariables( ) override{

                tardigradeHydra::perzynaViscoplasticity::residual::set_plasticStateVariables( plasticStateVariables );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };
    
            tardigradeHydra::linearElasticity::residual elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = tardigradeHydra::linearElasticity::residual( this, 9, elasticityParameters );

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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

    floatVector previousStateVariables = { 0, -0.1, 0, 0, 0, 0.3, 0, 0, 0,
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

    residualMock R( &hydra, 10, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters, 0.6 );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatVector answer = { -0.2, 2, 3, 4, 4, 6, 7, 8, 8, -5.1 };

    BOOST_TEST( answer == *R.getResidual( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_getJacobian, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the jacobian
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
            stressMock elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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
                                           0.01, 0.02, 0.005, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 2, 1, 1, 7, 1, 6, 1, 8, 3,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 0.01, 7, 8 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 10, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    unsigned int nvar = 10;

    floatMatrix jacobian( nvar, floatVector( unknownVector.size( ), 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < unknownVector.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        residualMock _R( &hydra, 10, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector + delta );

        floatVector dp = *_R.getResidual( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatVector dm = *_R.getResidual( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            jacobian[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( jacobian ) == *R.getJacobian( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_getdRdT, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the jacobian
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
            stressMock elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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
                                           0.01, 0.02, 0.005, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 2, 1, 1, 7, 1, 6, 1, 8, 3,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 0.01, 7, 8 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 10, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    unsigned int nvar = 10;

    floatVector dRdT( nvar, 0 );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatVector delta( 1, 0 );

        delta[ i ] = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature + delta[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature - delta[ i ], previousTemperature, deformationGradient, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 10, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 10, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getResidual( );

        floatVector dm = *Rm.getResidual( );

        dRdT = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_TEST( dRdT == *R.getdRdT( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_getdRdF, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the jacobian
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
            stressMock elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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
                                           0.01, 0.02, 0.005, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 2, 1, 1, 7, 1, 6, 1, 8, 3,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 0.01, 7, 8 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 10, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    unsigned int nvar = 10;

    floatMatrix dRdF( nvar, floatVector( 9, 0 ) );

    floatType eps = 1e-6;

    for ( unsigned int i = 0; i < deformationGradient.size( ); i++ ){

        floatVector delta( deformationGradient.size( ), 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               { }, { },
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 10, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 10, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getResidual( );

        floatVector dm = *Rm.getResidual( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dRdF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_TEST( tardigradeVectorTools::appendVectors( dRdF ) == *R.getdRdF( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_residual_getParameterizationInfo, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test of computing the jacobian
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousStress( ){

                tardigradeHydra::residualBase::setPreviousStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::perzynaViscoplasticity::residual {

        public:

            using tardigradeHydra::perzynaViscoplasticity::residual::residual;

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;
    
            floatVector elasticityParameters = { 123.4, 56.7 };

            floatVector viscoPlasticParameters = { 10.0,
                                                   1e1, 1e2,
                                                   10, 200, 293.15,
                                                   5, 0.34,
                                                   0.12,
                                                   13., 14.};
                                                   

            std::vector< unsigned int > stateVariableIndices = { 2 };

            floatVector thermalParameters = { 293.15, 1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5, 10e-5, 11e-5, 12e-5 };

            floatVector _local_deltaPreviousCauchyStress = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    
            stressMock elasticity;
    
            residualMock viscoPlasticity;

            tardigradeHydra::thermalExpansion::residual thermalExpansion;
    
            tardigradeHydra::residualBase remainder;
    
            void setResidualClasses( std::vector< tardigradeHydra::residualBase* > &residuals ){
    
                tardigradeHydra::hydraBase::setResidualClasses( residuals );
    
            }

        private:

            virtual void setResidualClasses( ){

                std::vector< tardigradeHydra::residualBase* > residuals( 4 );

                elasticity = stressMock( this, 9 );

                elasticity.previousCauchyStress += _local_deltaPreviousCauchyStress;

                viscoPlasticity = residualMock( this, 10, 1, stateVariableIndices, viscoPlasticParameters );

                thermalExpansion = tardigradeHydra::thermalExpansion::residual( this, 9, 2, thermalParameters );

                remainder = tardigradeHydra::residualBase( this, 4 );

                residuals[ 0 ] = &elasticity;

                residuals[ 1 ] = &viscoPlasticity;

                residuals[ 2 ] = &thermalExpansion;

                residuals[ 3 ] = &remainder;

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
                                           0.01, 0.02, 0.005, 0.04, 0.05 };

    floatVector parameters = { 123.4, 56.7, 293.15, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };

    floatVector unknownVector = { 2, 1, 1, 7, 1, 6, 1, 8, 3,
                                  1.2, 0.0, 0.0,
                                  0.0, 1.0, 0.0,
                                  0.0, 0.0, 1.0,
                                  1.01, 0.0, 0.0,
                                  0.0,  1.2, 0.0,
                                  0.0,  0.0, 0.9,
                                  4, 5, 0.01, 7, 8 };

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 10, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    std::string output;
    R.addParameterizationInfo( output );

}
