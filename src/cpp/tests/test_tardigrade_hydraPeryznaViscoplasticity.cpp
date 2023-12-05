/**
  * \file test_tardigrade-hydraPeryznaViscoplasticity.cpp
  *
  * Tests for tardigrade-hydraPeryznaViscoplasticity
  */

#include<tardigrade_hydraPeryznaViscoplasticity.h>
#include<tardigrade_hydraThermalExpansion.h>
#include<tardigrade_hydraLinearElasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_stress_tools.h>

#define BOOST_TEST_MODULE test_tardigrade-hydraPeryznaViscoplasticity
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

typedef tardigradeErrorTools::Node errorNode; //!< Redefinition for the error node
typedef errorNode* errorOut; //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::peryznaViscoplasticity::floatType floatType; //!< Redefinition of the floating point type
typedef tardigradeHydra::peryznaViscoplasticity::floatVector floatVector; //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::peryznaViscoplasticity::floatMatrix floatMatrix; //!< Redefinition of the matrix of floating points type

namespace tardigradeHydra{

    namespace unit_test{

        class hydraBaseTester{

            public:

                static void updateUnknownVector( tardigradeHydra::hydraBase &hydra, const floatVector &value ){

                    BOOST_CHECK_NO_THROW( hydra.updateUnknownVector( value ) );

                }

        };

    }

    namespace peryznaViscoplasticity{

        namespace unit_test{

            class residualTester{

                public:

                    static void runBasicGetTests( tardigradeHydra::peryznaViscoplasticity::residual &R ){

                        try{
    
                            BOOST_CHECK( &R._drivingStress.second == R.getDrivingStress( ) );

                            BOOST_CHECK( &R._previousDrivingStress.second == R.getPreviousDrivingStress( ) );
        
                            BOOST_CHECK( &R._flowDirection.second == R.getFlowDirection( ) );

                            BOOST_CHECK( &R._previousFlowDirection.second == R.getPreviousFlowDirection( ) );
    
                            BOOST_CHECK( &R._yieldFunction.second == R.getYieldFunction( ) );

                            BOOST_CHECK( &R._previousYieldFunction.second == R.getPreviousYieldFunction( ) );
    
                            BOOST_CHECK( &R._plasticThermalMultiplier.second == R.getPlasticThermalMultiplier( ) );

                            BOOST_CHECK( &R._previousPlasticThermalMultiplier.second == R.getPreviousPlasticThermalMultiplier( ) );

                            BOOST_CHECK( &R._dragStress.second == R.getDragStress( ) );

                            BOOST_CHECK( &R._previousDragStress.second == R.getPreviousDragStress( ) );

                            BOOST_CHECK( &R._hardeningFunction.second == R.getHardeningFunction( ) );

                            BOOST_CHECK( &R._previousHardeningFunction.second == R.getPreviousHardeningFunction( ) );
        
                            BOOST_CHECK( &R._plasticMultiplier.second == R.getPlasticMultiplier( ) );

                            BOOST_CHECK( &R._previousPlasticMultiplier.second == R.getPreviousPlasticMultiplier( ) );
        
                            BOOST_CHECK( &R._velocityGradient.second == R.getVelocityGradient( ) );

                            BOOST_CHECK( &R._previousVelocityGradient.second == R.getPreviousVelocityGradient( ) );

                            BOOST_CHECK( &R._stateVariableEvolutionRates.second == R.getStateVariableEvolutionRates( ) );

                            BOOST_CHECK( &R._previousStateVariableEvolutionRates.second == R.getPreviousStateVariableEvolutionRates( ) );

                            BOOST_CHECK( &R._plasticDeformationGradient.second == R.getPlasticDeformationGradient( ) );

                            BOOST_CHECK( &R._plasticStateVariables.second == R.getPlasticStateVariables( ) );
        
                            BOOST_CHECK( &R._stateVariables.second == R.getStateVariables( ) );
        
                            BOOST_CHECK( &R._peryznaParameters.second == R.getPeryznaParameters( ) );
        
                            BOOST_CHECK( &R._dragStressParameters.second == R.getDragStressParameters( ) );
        
                            BOOST_CHECK( &R._thermalParameters.second == R.getThermalParameters( ) );
        
                            BOOST_CHECK( &R._yieldParameters.second == R.getYieldParameters( ) );
        
                            BOOST_CHECK( &R._flowParameters.second == R.getFlowParameters( ) );
        
                            BOOST_CHECK( &R._hardeningParameters.second == R.getHardeningParameters( ) );

                            BOOST_CHECK( &R._integrationParameter == R.getIntegrationParameter( ) );

                            BOOST_CHECK( &R._dDrivingStressdCauchyStress.second == R.getdDrivingStressdCauchyStress( ) );

                            BOOST_CHECK( &R._dDrivingStressdF.second == R.getdDrivingStressdF( ) );

                            BOOST_CHECK( &R._dDrivingStressdSubFs.second == R.getdDrivingStressdSubFs( ) );

                            BOOST_CHECK( &R._dPreviousDrivingStressdPreviousCauchyStress.second == R.getdPreviousDrivingStressdPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPreviousDrivingStressdPreviousF.second == R.getdPreviousDrivingStressdPreviousF( ) );

                            BOOST_CHECK( &R._dPreviousDrivingStressdPreviousSubFs.second == R.getdPreviousDrivingStressdPreviousSubFs( ) );

                            BOOST_CHECK( &R._dFlowDirectiondCauchyStress.second == R.getdFlowDirectiondCauchyStress( ) );

                            BOOST_CHECK( &R._dFlowDirectiondF.second == R.getdFlowDirectiondF( ) );

                            BOOST_CHECK( &R._dFlowDirectiondSubFs.second == R.getdFlowDirectiondSubFs( ) );

                            BOOST_CHECK( &R._dPreviousFlowDirectiondPreviousCauchyStress.second == R.getdPreviousFlowDirectiondPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPreviousFlowDirectiondPreviousF.second == R.getdPreviousFlowDirectiondPreviousF( ) );

                            BOOST_CHECK( &R._dPreviousFlowDirectiondPreviousSubFs.second == R.getdPreviousFlowDirectiondPreviousSubFs( ) );

                            BOOST_CHECK( &R._dYieldFunctiondCauchyStress.second == R.getdYieldFunctiondCauchyStress( ) );

                            BOOST_CHECK( &R._dYieldFunctiondF.second == R.getdYieldFunctiondF( ) );

                            BOOST_CHECK( &R._dYieldFunctiondSubFs.second == R.getdYieldFunctiondSubFs( ) );

                            BOOST_CHECK( &R._dPreviousYieldFunctiondPreviousCauchyStress.second == R.getdPreviousYieldFunctiondPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPreviousYieldFunctiondPreviousF.second == R.getdPreviousYieldFunctiondPreviousF( ) );

                            BOOST_CHECK( &R._dPreviousYieldFunctiondPreviousSubFs.second == R.getdPreviousYieldFunctiondPreviousSubFs( ) );

                            BOOST_CHECK( &R._dPlasticThermalMultiplierdT.second == R.getdPlasticThermalMultiplierdT( ) );

                            BOOST_CHECK( &R._dPreviousPlasticThermalMultiplierdPreviousT.second == R.getdPreviousPlasticThermalMultiplierdPreviousT( ) );

                            BOOST_CHECK( &R._dDragStressdStateVariables.second == R.getdDragStressdStateVariables( ) );

                            BOOST_CHECK( &R._dPreviousDragStressdPreviousStateVariables.second == R.getdPreviousDragStressdPreviousStateVariables( ) );

                            BOOST_CHECK( &R._dPlasticMultiplierdCauchyStress.second == R.getdPlasticMultiplierdCauchyStress( ) );

                            BOOST_CHECK( &R._dPlasticMultiplierdF.second == R.getdPlasticMultiplierdF( ) );

                            BOOST_CHECK( &R._dPlasticMultiplierdSubFs.second == R.getdPlasticMultiplierdSubFs( ) );

                            BOOST_CHECK( &R._dPlasticMultiplierdT.second == R.getdPlasticMultiplierdT( ) );

                            BOOST_CHECK( &R._dPlasticMultiplierdStateVariables.second == R.getdPlasticMultiplierdStateVariables( ) );

                            BOOST_CHECK( &R._dPreviousPlasticMultiplierdPreviousCauchyStress.second == R.getdPreviousPlasticMultiplierdPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPreviousPlasticMultiplierdPreviousF.second == R.getdPreviousPlasticMultiplierdPreviousF( ) );

                            BOOST_CHECK( &R._dPreviousPlasticMultiplierdPreviousSubFs.second == R.getdPreviousPlasticMultiplierdPreviousSubFs( ) );

                            BOOST_CHECK( &R._dPreviousPlasticMultiplierdPreviousT.second == R.getdPreviousPlasticMultiplierdPreviousT( ) );

                            BOOST_CHECK( &R._dPreviousPlasticMultiplierdPreviousStateVariables.second == R.getdPreviousPlasticMultiplierdPreviousStateVariables( ) );

                            BOOST_CHECK( &R._dVelocityGradientdCauchyStress.second   == R.getdVelocityGradientdCauchyStress( ) );

                            BOOST_CHECK( &R._dVelocityGradientdF.second              == R.getdVelocityGradientdF( ) );

                            BOOST_CHECK( &R._dVelocityGradientdSubFs.second          == R.getdVelocityGradientdSubFs( ) );

                            BOOST_CHECK( &R._dVelocityGradientdT.second              == R.getdVelocityGradientdT( ) );

                            BOOST_CHECK( &R._dVelocityGradientdStateVariables.second == R.getdVelocityGradientdStateVariables( ) );

                            BOOST_CHECK( &R._dPreviousVelocityGradientdPreviousCauchyStress.second   == R.getdPreviousVelocityGradientdPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPreviousVelocityGradientdPreviousF.second              == R.getdPreviousVelocityGradientdPreviousF( ) );

                            BOOST_CHECK( &R._dPreviousVelocityGradientdPreviousSubFs.second          == R.getdPreviousVelocityGradientdPreviousSubFs( ) );

                            BOOST_CHECK( &R._dPreviousVelocityGradientdPreviousT.second              == R.getdPreviousVelocityGradientdPreviousT( ) );

                            BOOST_CHECK( &R._dPreviousVelocityGradientdPreviousStateVariables.second == R.getdPreviousVelocityGradientdPreviousStateVariables( ) );

                            BOOST_CHECK( &R._dStateVariableEvolutionRatesdCauchyStress.second   == R.getdStateVariableEvolutionRatesdCauchyStress( ) );

                            BOOST_CHECK( &R._dStateVariableEvolutionRatesdF.second              == R.getdStateVariableEvolutionRatesdF( ) );

                            BOOST_CHECK( &R._dStateVariableEvolutionRatesdSubFs.second          == R.getdStateVariableEvolutionRatesdSubFs( ) );

                            BOOST_CHECK( &R._dStateVariableEvolutionRatesdT.second              == R.getdStateVariableEvolutionRatesdT( ) );

                            BOOST_CHECK( &R._dStateVariableEvolutionRatesdStateVariables.second == R.getdStateVariableEvolutionRatesdStateVariables( ) );

                            BOOST_CHECK( &R._dPreviousStateVariableEvolutionRatesdPreviousCauchyStress.second   == R.getdPreviousStateVariableEvolutionRatesdPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPreviousStateVariableEvolutionRatesdPreviousF.second              == R.getdPreviousStateVariableEvolutionRatesdPreviousF( ) );

                            BOOST_CHECK( &R._dPreviousStateVariableEvolutionRatesdPreviousSubFs.second          == R.getdPreviousStateVariableEvolutionRatesdPreviousSubFs( ) );

                            BOOST_CHECK( &R._dPreviousStateVariableEvolutionRatesdPreviousT.second              == R.getdPreviousStateVariableEvolutionRatesdPreviousT( ) );

                            BOOST_CHECK( &R._dPreviousStateVariableEvolutionRatesdPreviousStateVariables.second == R.getdPreviousStateVariableEvolutionRatesdPreviousStateVariables( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdCauchyStress.second           == R.getdPlasticDeformationGradientdCauchyStress( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdF.second                      == R.getdPlasticDeformationGradientdF( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdSubFs.second                  == R.getdPlasticDeformationGradientdSubFs( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdT.second                      == R.getdPlasticDeformationGradientdT( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdStateVariables.second         == R.getdPlasticDeformationGradientdStateVariables( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdPreviousCauchyStress.second   == R.getdPlasticDeformationGradientdPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdPreviousF.second              == R.getdPlasticDeformationGradientdPreviousF( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdPreviousSubFs.second          == R.getdPlasticDeformationGradientdPreviousSubFs( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdPreviousT.second              == R.getdPlasticDeformationGradientdPreviousT( ) );

                            BOOST_CHECK( &R._dPlasticDeformationGradientdPreviousStateVariables.second == R.getdPlasticDeformationGradientdPreviousStateVariables( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdCauchyStress.second           == R.getdPlasticStateVariablesdCauchyStress( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdF.second                      == R.getdPlasticStateVariablesdF( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdSubFs.second                  == R.getdPlasticStateVariablesdSubFs( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdT.second                      == R.getdPlasticStateVariablesdT( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdStateVariables.second         == R.getdPlasticStateVariablesdStateVariables( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdPreviousCauchyStress.second   == R.getdPlasticStateVariablesdPreviousCauchyStress( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdPreviousF.second              == R.getdPlasticStateVariablesdPreviousF( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdPreviousSubFs.second          == R.getdPlasticStateVariablesdPreviousSubFs( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdPreviousT.second              == R.getdPlasticStateVariablesdPreviousT( ) );

                            BOOST_CHECK( &R._dPlasticStateVariablesdPreviousStateVariables.second == R.getdPlasticStateVariablesdPreviousStateVariables( ) );

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

BOOST_AUTO_TEST_CASE( test_residual_basicGetTests ){

    class stressResidualMock : public tardigradeHydra::residualBase {

        public:

            floatVector currentCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousCauchyStress = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

        private:

            using tardigradeHydra::residualBase::residualBase;

            virtual void setCauchyStress( ){

                tardigradeHydra::residualBase::setCauchyStress( currentCauchyStress );

            }

            virtual void setPreviousCauchyStress( ){

                tardigradeHydra::residualBase::setPreviousCauchyStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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

    try{
        hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );
    
        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );
    
        tardigradeHydra::peryznaViscoplasticity::unit_test::residualTester::runBasicGetTests( R );
    }
    catch( std::exception &e ){
        tardigradeErrorTools::printNestedExceptions( e );
        throw e;
    }

}

BOOST_AUTO_TEST_CASE( test_residual_getDrivingStress ){
    /*!
     * Test of computing the driving stress
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousCauchyStress( ){

                tardigradeHydra::residualBase::setPreviousCauchyStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock Rjac( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatVector answer = { 0.82636364, 1.0       , 0.75,
                           1.        , 1.21012101, 0.90759076,
                           0.75      , 0.90759076, 0.68069307 };

    Rjac.getdDrivingStressdCauchyStress( );

    Rjac.getdPreviousDrivingStressdPreviousCauchyStress( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getDrivingStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *Rjac.getDrivingStress( ) ) );

    answer = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getPreviousDrivingStress( ) ) ); 

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *Rjac.getPreviousDrivingStress( ) ) );

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

        floatVector dp = *R.getDrivingStress( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatVector dm = *R.getDrivingStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dDrivingStressdCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDrivingStressdCauchyStress, *Rjac.getdDrivingStressdCauchyStress( ) ) );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getDrivingStress( );

        floatVector dm = *Rm.getDrivingStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dDrivingStressdF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDrivingStressdF, *Rjac.getdDrivingStressdF( ) ) );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getDrivingStress( );

        floatVector dm = *Rm.getDrivingStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dDrivingStressdSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDrivingStressdSubFs, *Rjac.getdDrivingStressdSubFs( ) ) );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousDrivingStress( );

        floatVector dm = *Rm.getPreviousDrivingStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousDrivingStressdPreviousCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousDrivingStressdPreviousCauchyStress, *Rjac.getdPreviousDrivingStressdPreviousCauchyStress( ) ) );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousDrivingStress( );

        floatVector dm = *Rm.getPreviousDrivingStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousDrivingStressdPreviousF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousDrivingStressdPreviousF, *Rjac.getdPreviousDrivingStressdPreviousF( ) ) );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousDrivingStress( );

        floatVector dm = *Rm.getPreviousDrivingStress( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousDrivingStressdPreviousSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousDrivingStressdPreviousSubFs, *Rjac.getdPreviousDrivingStressdPreviousSubFs( ) ) );
}

BOOST_AUTO_TEST_CASE( test_residual_getFlowDirection ){
    /*!
     * Test of computing the driving stress
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

            floatVector drivingStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousDrivingStress = { 1, 0.2, 0, 0.4, 7, -1, 0.5, 3, .2 };

            floatVector flowParameters = { 1.23, 0.46 };

        private:

            virtual void setDrivingStress( const bool isPrevious ) override{

                setFlowParameters( flowParameters );

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setPreviousDrivingStress( previousDrivingStress );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setDrivingStress( drivingStress );

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType dpYield;
    floatVector jac( 9, 0 );
    floatVector answer( 9, 0 );
    floatVector answer2( 9, 0 );

    TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::druckerPragerSurface( R.drivingStress, R.flowParameters[ 1 ], R.flowParameters[ 0 ], dpYield, jac, answer ) );

    TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::druckerPragerSurface( R.previousDrivingStress, R.flowParameters[ 1 ], R.flowParameters[ 0 ], dpYield, jac, answer2 ) );

    BOOST_CHECK( R.getFlowDirection( )->size( ) == 9 );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getFlowDirection( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *R.getPreviousFlowDirection( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getFlowDirection_jacobian ){
    /*!
     * Test of computing the flow direction stress
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousCauchyStress( ){

                tardigradeHydra::residualBase::setPreviousCauchyStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R_ngrad( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    R.getdFlowDirectiondCauchyStress( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *R.getFlowDirection( ), *R_ngrad.getFlowDirection( ) ) );

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

        floatVector dp = *R.getFlowDirection( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatVector dm = *R.getFlowDirection( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dFlowDirectiondCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFlowDirectiondCauchyStress, *R.getdFlowDirectiondCauchyStress( ) ) );
    std::cout << "dFlowDirectiondCauchyStress:\n"; tardigradeVectorTools::print( dFlowDirectiondCauchyStress );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getFlowDirection( );

        floatVector dm = *Rm.getFlowDirection( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dFlowDirectiondF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFlowDirectiondF, *R.getdFlowDirectiondF( ) ) );
    std::cout << "dFlowDirectiondF:\n"; tardigradeVectorTools::print( dFlowDirectiondF );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getFlowDirection( );

        floatVector dm = *Rm.getFlowDirection( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dFlowDirectiondSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dFlowDirectiondSubFs, *R.getdFlowDirectiondSubFs( ) ) );
    std::cout << "dFlowDirectiondSubFs:\n"; tardigradeVectorTools::print( dFlowDirectiondSubFs );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousFlowDirection( );

        floatVector dm = *Rm.getPreviousFlowDirection( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousFlowDirectiondPreviousCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousFlowDirectiondPreviousCauchyStress, *R.getdPreviousFlowDirectiondPreviousCauchyStress( ) ) );
    std::cout << "dPreviousFlowDirectiondPreviousCauchyStress:\n"; tardigradeVectorTools::print( dPreviousFlowDirectiondPreviousCauchyStress );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousFlowDirection( );

        floatVector dm = *Rm.getPreviousFlowDirection( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousFlowDirectiondPreviousF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousFlowDirectiondPreviousF, *R.getdPreviousFlowDirectiondPreviousF( ) ) );
    std::cout << "dPreviousFlowDirectiondPreviousF:\n"; tardigradeVectorTools::print( dPreviousFlowDirectiondPreviousF );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousFlowDirection( );

        floatVector dm = *Rm.getPreviousFlowDirection( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousFlowDirectiondPreviousSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousFlowDirectiondPreviousSubFs, *R.getdPreviousFlowDirectiondPreviousSubFs( ) ) );
    std::cout << "dPreviousFlowDirectiondPreviousSubFs:\n"; tardigradeVectorTools::print( dPreviousFlowDirectiondPreviousSubFs );
}

BOOST_AUTO_TEST_CASE( test_residual_getYieldFunction ){
    /*!
     * Test of computing the yield function
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

            floatVector drivingStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousDrivingStress = { 1, 0.2, 0, 0.4, 7, -1, 0.5, 3, .2 };

            floatVector yieldParameters = { 1.23, 0.46 };

        private:

            virtual void setDrivingStress( const bool isPrevious ) override{

                setYieldParameters( yieldParameters );

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setPreviousDrivingStress( previousDrivingStress );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setDrivingStress( drivingStress );

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType answer;
    TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::druckerPragerSurface( R.drivingStress, R.yieldParameters[1], R.yieldParameters[0], answer ) );

    floatType answer2;
    TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeStressTools::druckerPragerSurface( R.previousDrivingStress, R.yieldParameters[1], R.yieldParameters[0], answer2 ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *R.getPreviousYieldFunction( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getYieldFunction_jacobian ){
    /*!
     * Test of computing the yield function's Jacobian
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousCauchyStress( ){

                tardigradeHydra::residualBase::setPreviousCauchyStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R_ngrad( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    R.getdYieldFunctiondCauchyStress( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *R.getYieldFunction( ), *R_ngrad.getYieldFunction( ) ) );

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

        floatType dp = *R.getYieldFunction( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatType dm = *R.getYieldFunction( );

        dYieldFunctiondCauchyStress[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dYieldFunctiondCauchyStress, *R.getdYieldFunctiondCauchyStress( ) ) );
    std::cout << "dYieldFunctiondCauchyStress:\n"; tardigradeVectorTools::print( dYieldFunctiondCauchyStress );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.getYieldFunction( );

        floatType dm = *Rm.getYieldFunction( );

        dYieldFunctiondF[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dYieldFunctiondF, *R.getdYieldFunctiondF( ) ) );
    std::cout << "dYieldFunctiondF:\n"; tardigradeVectorTools::print( dYieldFunctiondF );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.getYieldFunction( );

        floatType dm = *Rm.getYieldFunction( );

        dYieldFunctiondSubFs[ i ] = ( dp - dm ) / ( 2 * delta[ i + 9 ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dYieldFunctiondSubFs, *R.getdYieldFunctiondSubFs( ) ) );
    std::cout << "dYieldFunctiondSubFs:\n"; tardigradeVectorTools::print( dYieldFunctiondSubFs );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.getPreviousYieldFunction( );

        floatType dm = *Rm.getPreviousYieldFunction( );

        dPreviousYieldFunctiondPreviousCauchyStress[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousYieldFunctiondPreviousCauchyStress, *R.getdPreviousYieldFunctiondPreviousCauchyStress( ) ) );
    std::cout << "dPreviousYieldFunctiondPreviousCauchyStress:\n"; tardigradeVectorTools::print( dPreviousYieldFunctiondPreviousCauchyStress );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.getPreviousYieldFunction( );

        floatType dm = *Rm.getPreviousYieldFunction( );

        dPreviousYieldFunctiondPreviousF[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousYieldFunctiondPreviousF, *R.getdPreviousYieldFunctiondPreviousF( ) ) );
    std::cout << "dPreviousYieldFunctiondPreviousF:\n"; tardigradeVectorTools::print( dPreviousYieldFunctiondPreviousF );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.getPreviousYieldFunction( );

        floatType dm = *Rm.getPreviousYieldFunction( );

        dPreviousYieldFunctiondPreviousSubFs[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousYieldFunctiondPreviousSubFs, *R.getdPreviousYieldFunctiondPreviousSubFs( ) ) );
    std::cout << "dPreviousYieldFunctiondPreviousSubFs:\n"; tardigradeVectorTools::print( dPreviousYieldFunctiondPreviousSubFs );
}

BOOST_AUTO_TEST_CASE( test_residual_getPlasticThermalMultiplier ){
    /*!
     * Test of computing the plastic thermal multiplier
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock Rjac( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    Rjac.getdPlasticThermalMultiplierdT( );

    Rjac.getdPreviousPlasticThermalMultiplierdPreviousT( );

    floatType exp = ( -hydra.viscoPlasticParameters[ 3 ] * ( temperature - hydra.viscoPlasticParameters[ 5 ] ) / ( hydra.viscoPlasticParameters[ 4 ] + temperature - hydra.viscoPlasticParameters[ 5 ] ) );

    floatType answer = std::pow( 10, exp );

    exp = ( -hydra.viscoPlasticParameters[ 3 ] * ( previousTemperature - hydra.viscoPlasticParameters[ 5 ] ) / ( hydra.viscoPlasticParameters[ 4 ] + previousTemperature - hydra.viscoPlasticParameters[ 5 ] ) );

    floatType answer2 = std::pow( 10, exp );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getPlasticThermalMultiplier( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *R.getPreviousPlasticThermalMultiplier( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *Rjac.getPlasticThermalMultiplier( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *Rjac.getPreviousPlasticThermalMultiplier( ) ) );

    // Check the Jacobians w.r.t. the temperature
    floatType eps = 1e-6;

    floatType dPlasticThermalMultiplierdT = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature + delta, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature - delta, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        residualMock Rm( &hydram, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        floatType dp = *Rp.getPlasticThermalMultiplier( );

        floatType dm = *Rm.getPlasticThermalMultiplier( );

        dPlasticThermalMultiplierdT = ( dp - dm ) / ( 2 * delta );

    }
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticThermalMultiplierdT, *Rjac.getdPlasticThermalMultiplierdT( ) ) );
    std::cout << "dPlasticThermalMultiplierdT: " << dPlasticThermalMultiplierdT << "\n";

    floatType dPreviousPlasticThermalMultiplierdPreviousT = 0;

    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( previousTemperature ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature + delta, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature - delta, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        residualMock Rm( &hydram, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        floatType dp = *Rp.getPreviousPlasticThermalMultiplier( );

        floatType dm = *Rm.getPreviousPlasticThermalMultiplier( );

        dPreviousPlasticThermalMultiplierdPreviousT += ( dp - dm ) / ( 2 * delta );

    }
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousPlasticThermalMultiplierdPreviousT, *Rjac.getdPreviousPlasticThermalMultiplierdPreviousT( ) ) );
    std::cout << "dPreviousPlasticThermalMultiplierdPreviousT: " << dPreviousPlasticThermalMultiplierdPreviousT << "\n";

}

BOOST_AUTO_TEST_CASE( test_residual_getDragStress ){
    /*!
     * Test of computing the drag stress
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

            floatVector stateVariables = { 0.6 };

            floatVector previousStateVariables = { 0.9 };

        private:

            virtual void setStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setPreviousStateVariables( previousStateVariables );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setStateVariables( stateVariables );

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType answer = hydra.viscoPlasticParameters[ 1 ] + R.stateVariables[ 0 ] * hydra.viscoPlasticParameters[ 2 ];

    floatType answer2 = hydra.viscoPlasticParameters[ 1 ] + R.previousStateVariables[ 0 ] * hydra.viscoPlasticParameters[ 2 ];

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getDragStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *R.getPreviousDragStress( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getDragStress_jacobian ){
    /*!
     * Test of computing the drag stress
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

            floatVector stateVariables = { 0.6 };

            floatVector previousStateVariables = { 0.9 };

        private:

            virtual void setStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setPreviousStateVariables( previousStateVariables );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setStateVariables( stateVariables );

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock Rjac( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType answer = hydra.viscoPlasticParameters[ 1 ] + R.stateVariables[ 0 ] * hydra.viscoPlasticParameters[ 2 ];

    floatType answer2 = hydra.viscoPlasticParameters[ 1 ] + R.previousStateVariables[ 0 ] * hydra.viscoPlasticParameters[ 2 ];

    Rjac.getdDragStressdStateVariables( );

    Rjac.getdPreviousDragStressdPreviousStateVariables( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getDragStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *R.getPreviousDragStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *Rjac.getDragStress( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *Rjac.getPreviousDragStress( ) ) );

    floatType eps = 1e-6;

    floatVector dDragStressdStateVariables( R.stateVariables.size( ), 0 );

    for ( unsigned int i = 0; i < R.stateVariables.size( ); i++ ){

        floatVector delta = floatVector( R.stateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( R.stateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        residualMock Rm( &hydram, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        Rp.stateVariables += delta;

        Rm.stateVariables -= delta;

        floatType dp = *Rp.getDragStress( );

        floatType dm = *Rm.getDragStress( );

        dDragStressdStateVariables[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dDragStressdStateVariables, *Rjac.getdDragStressdStateVariables( ) ) );
    std::cout << "dDragStressdStateVariables: "; tardigradeVectorTools::print( dDragStressdStateVariables );

    floatVector dPreviousDragStressdPreviousStateVariables( R.previousStateVariables.size( ), 0 );

    for ( unsigned int i = 0; i < R.previousStateVariables.size( ); i++ ){

        floatVector delta = floatVector( R.previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( R.previousStateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        residualMock Rm( &hydram, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        Rp.previousStateVariables += delta;

        Rm.previousStateVariables -= delta;

        floatType dp = *Rp.getPreviousDragStress( );

        floatType dm = *Rm.getPreviousDragStress( );

        dPreviousDragStressdPreviousStateVariables[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousDragStressdPreviousStateVariables, *Rjac.getdPreviousDragStressdPreviousStateVariables( ) ) );
    std::cout << "dPreviousDragStressdPreviousStateVariables: "; tardigradeVectorTools::print( dPreviousDragStressdPreviousStateVariables );

}

BOOST_AUTO_TEST_CASE( test_residual_getHardeningFunction ){
    /*!
     * Test of computing the hardening function
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

            floatVector stateVariables = { 0.6 };

            floatVector previousStateVariables = { 0.9 };

        private:

            virtual void setStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setPreviousStateVariables( previousStateVariables );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setStateVariables( stateVariables );

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock Rjac( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType answer = hydra.viscoPlasticParameters[ 9 ] + R.stateVariables[ 0 ] * hydra.viscoPlasticParameters[ 10 ];

    floatType answer2 = hydra.viscoPlasticParameters[ 9 ] + R.previousStateVariables[ 0 ] * hydra.viscoPlasticParameters[ 10 ];

    Rjac.getdHardeningFunctiondStateVariables( );

    Rjac.getdPreviousHardeningFunctiondPreviousStateVariables( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getHardeningFunction( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *R.getPreviousHardeningFunction( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *Rjac.getHardeningFunction( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *Rjac.getPreviousHardeningFunction( ) ) );

    floatType eps = 1e-6;

    floatVector dHardeningFunctiondStateVariables( R.stateVariables.size( ), 0 );

    for ( unsigned int i = 0; i < R.stateVariables.size( ); i++ ){

        floatVector delta = floatVector( R.stateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( R.stateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        residualMock Rm( &hydram, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        Rp.stateVariables += delta;

        Rm.stateVariables -= delta;

        floatType dp = *Rp.getHardeningFunction( );

        floatType dm = *Rm.getHardeningFunction( );

        dHardeningFunctiondStateVariables[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dHardeningFunctiondStateVariables, *Rjac.getdHardeningFunctiondStateVariables( ) ) );
    std::cout << "dHardeningFunctiondStateVariables: "; tardigradeVectorTools::print( dHardeningFunctiondStateVariables );

    floatVector dPreviousHardeningFunctiondPreviousStateVariables( R.previousStateVariables.size( ), 0 );

    for ( unsigned int i = 0; i < R.previousStateVariables.size( ); i++ ){

        floatVector delta = floatVector( R.previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( R.previousStateVariables[ i ] ) + eps;

        hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                              previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        residualMock Rp( &hydrap, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        residualMock Rm( &hydram, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

        Rp.previousStateVariables += delta;

        Rm.previousStateVariables -= delta;

        floatType dp = *Rp.getPreviousHardeningFunction( );

        floatType dm = *Rm.getPreviousHardeningFunction( );

        dPreviousHardeningFunctiondPreviousStateVariables[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }
    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousHardeningFunctiondPreviousStateVariables, *Rjac.getdPreviousHardeningFunctiondPreviousStateVariables( ) ) );
    std::cout << "dPreviousHardeningFunctiondPreviousStateVariables: "; tardigradeVectorTools::print( dPreviousHardeningFunctiondPreviousStateVariables );
}

BOOST_AUTO_TEST_CASE( test_residual_decomposeParameters ){
    /*!
     * Test of decomposing the parameter vector
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatVector answer_peryznaParameters = floatVector( hydra.viscoPlasticParameters.begin( ), hydra.viscoPlasticParameters.begin( ) + 1 );

    floatVector answer_dragStressParameters = floatVector( hydra.viscoPlasticParameters.begin( ) + 1, hydra.viscoPlasticParameters.begin( ) + 3 );

    floatVector answer_thermalParameters = floatVector( hydra.viscoPlasticParameters.begin( ) + 3, hydra.viscoPlasticParameters.begin( ) + 6 );

    floatVector answer_yieldParameters = floatVector( hydra.viscoPlasticParameters.begin( ) + 6, hydra.viscoPlasticParameters.begin( ) + 8 );

    floatVector answer_flowParameters = { 0, hydra.viscoPlasticParameters[ 8 ] };

    floatVector answer_hardeningParameters = { hydra.viscoPlasticParameters[ 9 ], hydra.viscoPlasticParameters[ 10 ] };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer_peryznaParameters, *R.getPeryznaParameters( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer_dragStressParameters, *R.getDragStressParameters( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer_thermalParameters, *R.getThermalParameters( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer_yieldParameters, *R.getYieldParameters( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer_flowParameters, *R.getFlowParameters( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer_hardeningParameters, *R.getHardeningParameters( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getStateVariableIndices ){
    /*!
     * Test getting the state variable indices
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    std::vector< unsigned int > answer = { 0, 2 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getStateVariableIndices( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getStateVariables ){
    /*!
     * Test of decomposing the parameter vector
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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

    unsigned int numNonLinearSolveStateVariables = 3;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatVector answer = { 4, 6 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getStateVariables( ) ) );

    floatVector answer2 = { 0.01, 0.03 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *R.getPreviousStateVariables( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getPlasticMuliplier ){
    /*!
     * Test of computing the plastic multiplier
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousCauchyStress( ){

                tardigradeHydra::residualBase::setPreviousCauchyStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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

            floatVector dqdXi = { 27 };

            floatType   dAdT  = 28;

            floatVector dPreviousfdPreviousCauchy = { 29, 30, 31, 32, 33, 34, 35, 36, 37 };

            floatVector dPreviousfdPreviousF = { 38, 39, 40, 41, 42, 43, 44, 45, 46 };

            floatVector dPreviousfdPreviousSubFs = { 47, 48, 49, 50, 51, 52, 53, 54, 55,
                                                     56, 57, 58, 59, 60, 61, 62, 63, 64 };

            floatVector dPreviousqdPreviousXi = { 65 };

            floatType   dPreviousAdPreviousT  = 66;
            

        private:

            virtual void setYieldFunction( const bool isPrevious ) override{

                setPeryznaParameters( { n } );

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setPreviousYieldFunction( previousf );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setYieldFunction( f );

                }

            }

            virtual void setdYieldFunctiondCauchyStress( const bool isPrevious ) override{

                setPeryznaParameters( { n } );

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousYieldFunctiondPreviousCauchyStress( dPreviousfdPreviousCauchy );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdYieldFunctiondCauchyStress( dfdCauchy );

                }

            }

            virtual void setdYieldFunctiondF( const bool isPrevious ) override{

                setPeryznaParameters( { n } );

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousYieldFunctiondPreviousF( dPreviousfdPreviousF );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdYieldFunctiondF( dfdF );

                }

            }

            virtual void setdYieldFunctiondSubFs( const bool isPrevious ) override{

                setPeryznaParameters( { n } );

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousYieldFunctiondPreviousSubFs( dPreviousfdPreviousSubFs );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdYieldFunctiondSubFs( dfdSubFs );

                }

            }

            virtual void setDragStress( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setPreviousDragStress( previousq );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setDragStress( q );

                }

            }

            virtual void setdDragStressdStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousDragStressdPreviousStateVariables( dPreviousqdPreviousXi );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdDragStressdStateVariables( dqdXi );

                }

            }

            virtual void setPlasticThermalMultiplier( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setPreviousPlasticThermalMultiplier( previousA );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setPlasticThermalMultiplier( A );

                }

            }

            virtual void setdPlasticThermalMultiplierdT( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousPlasticThermalMultiplierdPreviousT( dPreviousAdPreviousT );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPlasticThermalMultiplierdT( dAdT );

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock Rjac( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatType answer = R.A * std::pow( R.f / R.q, R.n );

    floatType answer2 = R.previousA * std::pow( R.previousf / R.previousq, R.n );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getPlasticMultiplier( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *R.getPreviousPlasticMultiplier( ) ) );

    Rjac.getdPlasticMultiplierdCauchyStress( );

    Rjac.getdPreviousPlasticMultiplierdPreviousCauchyStress( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *Rjac.getPlasticMultiplier( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *Rjac.getPreviousPlasticMultiplier( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getPlasticMultiplier_jacobian ){
    /*!
     * Test of computing the plastic multiplier's Jacobian
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousCauchyStress( ){

                tardigradeHydra::residualBase::setPreviousCauchyStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R_ngrad( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    R.getdYieldFunctiondCauchyStress( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *R.getYieldFunction( ), *R_ngrad.getYieldFunction( ) ) );

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

        floatType dp = *R.getPlasticMultiplier( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatType dm = *R.getPlasticMultiplier( );

        dPlasticMultiplierdCauchyStress[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticMultiplierdCauchyStress, *R.getdPlasticMultiplierdCauchyStress( ) ) );
    std::cout << "dPlasticMultiplierdCauchyStress:\n"; tardigradeVectorTools::print( dPlasticMultiplierdCauchyStress );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.getPlasticMultiplier( );

        floatType dm = *Rm.getPlasticMultiplier( );

        dPlasticMultiplierdF[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticMultiplierdF, *R.getdPlasticMultiplierdF( ) ) );
    std::cout << "dPlasticMultiplierdF:\n"; tardigradeVectorTools::print( dPlasticMultiplierdF );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.getPlasticMultiplier( );

        floatType dm = *Rm.getPlasticMultiplier( );

        dPlasticMultiplierdSubFs[ i ] = ( dp - dm ) / ( 2 * delta[ i + 9 ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticMultiplierdSubFs, *R.getdPlasticMultiplierdSubFs( ) ) );
    std::cout << "dPlasticMultiplierdSubFs:\n"; tardigradeVectorTools::print( dPlasticMultiplierdSubFs );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature + delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature - delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.getPlasticMultiplier( );

        floatType dm = *Rm.getPlasticMultiplier( );

        dPlasticMultiplierdT = ( dp - dm ) / ( 2 * delta );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticMultiplierdT, *R.getdPlasticMultiplierdT( ) ) );
    std::cout << "dPlasticMultiplierdT:\n" << dPlasticMultiplierdT << "\n";

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ 27 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( unknownVector[ 27 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.getPlasticMultiplier( );

        floatType dm = *Rm.getPlasticMultiplier( );

        dPlasticMultiplierdStateVariables[ i ] = ( dp - dm ) / ( 2 * delta[ 27 + hydra.stateVariableIndices[ i ] ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticMultiplierdStateVariables, *R.getdPlasticMultiplierdStateVariables( ) ) );
    std::cout << "dPlasticMultiplierdStateVariables:\n"; tardigradeVectorTools::print( dPlasticMultiplierdStateVariables );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.getPreviousPlasticMultiplier( );

        floatType dm = *Rm.getPreviousPlasticMultiplier( );

        dPreviousPlasticMultiplierdPreviousCauchyStress[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousPlasticMultiplierdPreviousCauchyStress, *R.getdPreviousPlasticMultiplierdPreviousCauchyStress( ) ) );
    std::cout << "dPreviousPlasticMultiplierdPreviousCauchyStress:\n"; tardigradeVectorTools::print( dPreviousPlasticMultiplierdPreviousCauchyStress );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.getPreviousPlasticMultiplier( );

        floatType dm = *Rm.getPreviousPlasticMultiplier( );

        dPreviousPlasticMultiplierdPreviousF[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousPlasticMultiplierdPreviousF, *R.getdPreviousPlasticMultiplierdPreviousF( ) ) );
    std::cout << "dPreviousPlasticMultiplierdPreviousF:\n"; tardigradeVectorTools::print( dPreviousPlasticMultiplierdPreviousF );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.getPreviousPlasticMultiplier( );

        floatType dm = *Rm.getPreviousPlasticMultiplier( );

        dPreviousPlasticMultiplierdPreviousSubFs[ i ] = ( dp - dm ) / ( 2 * delta[ i ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousPlasticMultiplierdPreviousSubFs, *R.getdPreviousPlasticMultiplierdPreviousSubFs( ) ) );
    std::cout << "dPreviousPlasticMultiplierdPreviousSubFs:\n"; tardigradeVectorTools::print( dPreviousPlasticMultiplierdPreviousSubFs );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( previousTemperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature + delta, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature - delta, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.getPreviousPlasticMultiplier( );

        floatType dm = *Rm.getPreviousPlasticMultiplier( );

        dPreviousPlasticMultiplierdPreviousT = ( dp - dm ) / ( 2 * delta );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousPlasticMultiplierdPreviousT, *R.getdPreviousPlasticMultiplierdPreviousT( ) ) );
    std::cout << "dPreviousPlasticMultiplierdPreviousT:\n" << dPreviousPlasticMultiplierdPreviousT << "\n";

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ 18 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( previousStateVariables[ 18 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatType dp = *Rp.getPreviousPlasticMultiplier( );

        floatType dm = *Rm.getPreviousPlasticMultiplier( );

        dPreviousPlasticMultiplierdPreviousStateVariables[ i ] = ( dp - dm ) / ( 2 * delta[ 18 + hydra.stateVariableIndices[ i ] ] );

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousPlasticMultiplierdPreviousStateVariables, *R.getdPreviousPlasticMultiplierdPreviousStateVariables( ) ) );
    std::cout << "dPreviousPlasticMultiplierdPreviousStateVariables:\n"; tardigradeVectorTools::print( dPreviousPlasticMultiplierdPreviousStateVariables );
    std::cout << "R.getdPreviousPlasticMultiplierdPreviousStateVariables:\n"; tardigradeVectorTools::print( *R.getdPreviousPlasticMultiplierdPreviousStateVariables( ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getVelocityGradient ){
    /*!
     * Test of computing the velocity gradient
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

            floatType gamma = 2.4;

            floatVector dGammadCauchy = initializeVector( 9, -0.56 );

            floatVector dGammadF = initializeVector( 9, 1 );

            floatVector dGammadSubFs = initializeVector( 18, 10 );

            floatType   dGammadT = 2.4;

            floatVector dGammadXi = { 1.1 };

            floatVector nhat = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatMatrix dNhatdCauchy = initializeMatrix( 9, 9, 0.1 );

            floatMatrix dNhatdF = initializeMatrix( 9, 9, 0.3 );

            floatMatrix dNhatdSubFs = initializeMatrix( 9, 18, -0.54 );

            floatType previousGamma = 5.6;

            floatVector dPreviousGammadPreviousCauchy = initializeVector( 9, -3 );

            floatVector dPreviousGammadPreviousF = initializeVector( 9, 0.123 );

            floatVector dPreviousGammadPreviousSubFs = initializeVector( 18, -0.224 );

            floatType   dPreviousGammadPreviousT = 0.872;

            floatVector dPreviousGammadPreviousXi = { 2.5 };

            floatVector previousNhat = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            floatMatrix dPreviousNhatdPreviousCauchy = initializeMatrix( 9, 9, -0.37 );

            floatMatrix dPreviousNhatdPreviousF = initializeMatrix( 9, 9, 0.373 );

            floatMatrix dPreviousNhatdPreviousSubFs = initializeMatrix( 9, 18, -0.82 );

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

            virtual void setPlasticMultiplier( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setPreviousPlasticMultiplier( previousGamma );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setPlasticMultiplier( gamma );

                }

            }

            virtual void setdPlasticMultiplierdCauchyStress( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousPlasticMultiplierdPreviousCauchyStress( dPreviousGammadPreviousCauchy );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPlasticMultiplierdCauchyStress( dGammadCauchy );

                }

            }

            virtual void setdPlasticMultiplierdF( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousPlasticMultiplierdPreviousF( dPreviousGammadPreviousF );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPlasticMultiplierdF( dGammadF );

                }

            }

            virtual void setdPlasticMultiplierdSubFs( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousPlasticMultiplierdPreviousSubFs( dPreviousGammadPreviousSubFs );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPlasticMultiplierdSubFs( dGammadSubFs );

                }

            }

            virtual void setdPlasticMultiplierdT( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousPlasticMultiplierdPreviousT( dPreviousGammadPreviousT );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPlasticMultiplierdT( dGammadT );

                }

            }

            virtual void setdPlasticMultiplierdStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousPlasticMultiplierdPreviousStateVariables( dPreviousGammadPreviousXi );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPlasticMultiplierdStateVariables( dGammadXi );

                }

            }

            virtual void setFlowDirection( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setPreviousFlowDirection( previousNhat );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setFlowDirection( nhat );

                }

            }

            virtual void setdFlowDirectiondCauchyStress( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousFlowDirectiondPreviousCauchyStress( dPreviousNhatdPreviousCauchy );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdFlowDirectiondCauchyStress( dNhatdCauchy );

                }

            }

            virtual void setdFlowDirectiondF( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousFlowDirectiondPreviousF( dPreviousNhatdPreviousF );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdFlowDirectiondF( dNhatdF );

                }

            }

            virtual void setdFlowDirectiondSubFs( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousFlowDirectiondPreviousSubFs( dPreviousNhatdPreviousSubFs );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdFlowDirectiondSubFs( dNhatdSubFs );

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock Rjac( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatVector answer = R.gamma * R.nhat;

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getVelocityGradient( ) ) );

    floatVector answer2 = R.previousGamma * R.previousNhat;

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *R.getPreviousVelocityGradient( ) ) );

    try{
    Rjac.getdVelocityGradientdCauchyStress( );

    Rjac.getdPreviousVelocityGradientdPreviousCauchyStress( );
    }
    catch(std::exception &e){tardigradeErrorTools::printNestedExceptions(e);}

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *Rjac.getVelocityGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *Rjac.getPreviousVelocityGradient( ) ) );


}

BOOST_AUTO_TEST_CASE( test_residual_getVelocityGradient_jacobian ){
    /*!
     * Test of computing the velocity gradient's Jacobian
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousCauchyStress( ){

                tardigradeHydra::residualBase::setPreviousCauchyStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R_ngrad( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    R.getdYieldFunctiondCauchyStress( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *R.getYieldFunction( ), *R_ngrad.getYieldFunction( ) ) );

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

        floatVector dp = *R.getVelocityGradient( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatVector dm = *R.getVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dVelocityGradientdCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dVelocityGradientdCauchyStress, *R.getdVelocityGradientdCauchyStress( ) ) );
    std::cout << "dVelocityGradientdCauchyStress:\n"; tardigradeVectorTools::print( dVelocityGradientdCauchyStress );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getVelocityGradient( );

        floatVector dm = *Rm.getVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dVelocityGradientdF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dVelocityGradientdF, *R.getdVelocityGradientdF( ) ) );
    std::cout << "dVelocityGradientdF:\n"; tardigradeVectorTools::print( dVelocityGradientdF );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getVelocityGradient( );

        floatVector dm = *Rm.getVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dVelocityGradientdSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dVelocityGradientdSubFs, *R.getdVelocityGradientdSubFs( ) ) );
    std::cout << "dVelocityGradientdSubFs:\n"; tardigradeVectorTools::print( dVelocityGradientdSubFs );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature + delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature - delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getVelocityGradient( );

        floatVector dm = *Rm.getVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dVelocityGradientdT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dVelocityGradientdT, *R.getdVelocityGradientdT( ) ) );
    std::cout << "dVelocityGradientdT:\n"; tardigradeVectorTools::print( dVelocityGradientdT );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ 27 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( unknownVector[ 27 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getVelocityGradient( );

        floatVector dm = *Rm.getVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dVelocityGradientdStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 27 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dVelocityGradientdStateVariables, *R.getdVelocityGradientdStateVariables( ), 1e-5, 1e-5 ) );
    std::cout << "dVelocityGradientdStateVariables:\n"; tardigradeVectorTools::print( dVelocityGradientdStateVariables );
    std::cout << "R.getdVelocityGradientdStateVariables:\n"; tardigradeVectorTools::print( *R.getdVelocityGradientdStateVariables( ) );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousVelocityGradient( );

        floatVector dm = *Rm.getPreviousVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousVelocityGradientdPreviousCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousVelocityGradientdPreviousCauchyStress, *R.getdPreviousVelocityGradientdPreviousCauchyStress( ) ) );
    std::cout << "dPreviousVelocityGradientdPreviousCauchyStress:\n"; tardigradeVectorTools::print( dPreviousVelocityGradientdPreviousCauchyStress );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousVelocityGradient( );

        floatVector dm = *Rm.getPreviousVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousVelocityGradientdPreviousF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousVelocityGradientdPreviousF, *R.getdPreviousVelocityGradientdPreviousF( ) ) );
    std::cout << "dPreviousVelocityGradientdPreviousF:\n"; tardigradeVectorTools::print( dPreviousVelocityGradientdPreviousF );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousVelocityGradient( );

        floatVector dm = *Rm.getPreviousVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousVelocityGradientdPreviousSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousVelocityGradientdPreviousSubFs, *R.getdPreviousVelocityGradientdPreviousSubFs( ) ) );
    std::cout << "dPreviousVelocityGradientdPreviousSubFs:\n"; tardigradeVectorTools::print( dPreviousVelocityGradientdPreviousSubFs );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( previousTemperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature + delta, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature - delta, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousVelocityGradient( );

        floatVector dm = *Rm.getPreviousVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousVelocityGradientdPreviousT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousVelocityGradientdPreviousT, *R.getdPreviousVelocityGradientdPreviousT( ) ) );
    std::cout << "dPreviousVelocityGradientdPreviousT:\n"; tardigradeVectorTools::print( dPreviousVelocityGradientdPreviousT );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ 18 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( previousStateVariables[ 18 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousVelocityGradient( );

        floatVector dm = *Rm.getPreviousVelocityGradient( );

        for ( unsigned int j = 0; j < 9; j++ ){

            dPreviousVelocityGradientdPreviousStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 18 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousVelocityGradientdPreviousStateVariables, *R.getdPreviousVelocityGradientdPreviousStateVariables( ) ) );
    std::cout << "dPreviousVelocityGradientdPreviousStateVariables:\n"; tardigradeVectorTools::print( dPreviousVelocityGradientdPreviousStateVariables );
    std::cout << "R.getdPreviousVelocityGradientdPreviousStateVariables:\n"; tardigradeVectorTools::print( *R.getdPreviousVelocityGradientdPreviousStateVariables( ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getStateVariableEvolutionRate ){
    /*!
     * Test of computing the velocity gradient
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

            floatType gamma = 2.4;

            floatVector dGammadCauchy = initializeVector( 9, -0.56 );

            floatVector dGammadF = initializeVector( 9, 1 );

            floatVector dGammadSubFs = initializeVector( 18, 10 );

            floatType   dGammadT = 2.4;

            floatVector dGammadXi = { 1.1 };

            floatType hardeningFunction = 7.8;

            floatVector dHardeningFunctiondXi = initializeVector( 1, -0.72 );

            floatType previousGamma = 5.6;

            floatVector dPreviousGammadPreviousCauchy = initializeVector( 9, -3 );

            floatVector dPreviousGammadPreviousF = initializeVector( 9, 0.123 );

            floatVector dPreviousGammadPreviousSubFs = initializeVector( 18, -0.224 );

            floatType   dPreviousGammadPreviousT = 0.872;

            floatVector dPreviousGammadPreviousXi = { 2.5 };

            floatType previousHardeningFunction = .56;

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

            virtual void setPlasticMultiplier( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setPreviousPlasticMultiplier( previousGamma );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setPlasticMultiplier( gamma );

                }

            }

            virtual void setdPlasticMultiplierdCauchyStress( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousPlasticMultiplierdPreviousCauchyStress( dPreviousGammadPreviousCauchy );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPlasticMultiplierdCauchyStress( dGammadCauchy );

                }

            }

            virtual void setdPlasticMultiplierdF( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousPlasticMultiplierdPreviousF( dPreviousGammadPreviousF );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPlasticMultiplierdF( dGammadF );

                }

            }

            virtual void setdPlasticMultiplierdSubFs( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousPlasticMultiplierdPreviousSubFs( dPreviousGammadPreviousSubFs );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPlasticMultiplierdSubFs( dGammadSubFs );

                }

            }

            virtual void setdPlasticMultiplierdT( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousPlasticMultiplierdPreviousT( dPreviousGammadPreviousT );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPlasticMultiplierdT( dGammadT );

                }

            }

            virtual void setdPlasticMultiplierdStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousPlasticMultiplierdPreviousStateVariables( dPreviousGammadPreviousXi );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPlasticMultiplierdStateVariables( dGammadXi );

                }

            }

            virtual void setHardeningFunction( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setPreviousHardeningFunction( previousHardeningFunction );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setHardeningFunction( hardeningFunction );

                }

            }

            virtual void setdHardeningFunctiondStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousHardeningFunctiondPreviousStateVariables( dPreviousHardeningFunctiondPreviousXi );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdHardeningFunctiondStateVariables( dHardeningFunctiondXi );

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock Rjac( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    floatVector answer = { R.gamma * R.hardeningFunction };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getStateVariableEvolutionRates( ) ) );

    floatVector answer2 = { R.previousGamma * R.previousHardeningFunction };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *R.getPreviousStateVariableEvolutionRates( ) ) );

    Rjac.getdStateVariableEvolutionRatesdCauchyStress( );

    Rjac.getdPreviousStateVariableEvolutionRatesdPreviousCauchyStress( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *Rjac.getStateVariableEvolutionRates( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer2, *Rjac.getPreviousStateVariableEvolutionRates( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getStateVariableEvolutionRates_jacobian ){
    /*!
     * Test of computing the velocity gradient's Jacobian
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousCauchyStress( ){

                tardigradeHydra::residualBase::setPreviousCauchyStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R_ngrad( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    R.getdYieldFunctiondCauchyStress( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *R.getYieldFunction( ), *R_ngrad.getYieldFunction( ) ) );

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

        floatVector dp = *R.getStateVariableEvolutionRates( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatVector dm = *R.getStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dStateVariableEvolutionRatesdCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dStateVariableEvolutionRatesdCauchyStress, *R.getdStateVariableEvolutionRatesdCauchyStress( ) ) );
    std::cout << "dStateVariableEvolutionRatesdCauchyStress:\n"; tardigradeVectorTools::print( dStateVariableEvolutionRatesdCauchyStress );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getStateVariableEvolutionRates( );

        floatVector dm = *Rm.getStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dStateVariableEvolutionRatesdF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dStateVariableEvolutionRatesdF, *R.getdStateVariableEvolutionRatesdF( ), 1e-5, 1e-5 ) );
    std::cout << "dStateVariableEvolutionRatesdF:\n"; tardigradeVectorTools::print( dStateVariableEvolutionRatesdF );
    std::cout << "R.getdStateVariableEvolutionRatesdF:\n"; tardigradeVectorTools::print( *R.getdStateVariableEvolutionRatesdF( ) );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getStateVariableEvolutionRates( );

        floatVector dm = *Rm.getStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dStateVariableEvolutionRatesdSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dStateVariableEvolutionRatesdSubFs, *R.getdStateVariableEvolutionRatesdSubFs( ), 1e-5, 1e-5 ) );
    std::cout << "dStateVariableEvolutionRatesdSubFs:\n"; tardigradeVectorTools::print( dStateVariableEvolutionRatesdSubFs );
    std::cout << "R.getdStateVariableEvolutionRatesdSubFs:\n"; tardigradeVectorTools::print( *R.getdStateVariableEvolutionRatesdSubFs( ) );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature + delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature - delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getStateVariableEvolutionRates( );

        floatVector dm = *Rm.getStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dStateVariableEvolutionRatesdT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dStateVariableEvolutionRatesdT, *R.getdStateVariableEvolutionRatesdT( ) ) );
    std::cout << "dStateVariableEvolutionRatesdT:\n"; tardigradeVectorTools::print( dStateVariableEvolutionRatesdT );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ 27 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( unknownVector[ 27 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getStateVariableEvolutionRates( );

        floatVector dm = *Rm.getStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dStateVariableEvolutionRatesdStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 27 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dStateVariableEvolutionRatesdStateVariables, *R.getdStateVariableEvolutionRatesdStateVariables( ), 1e-5, 1e-5 ) );
    std::cout << "dStateVariableEvolutionRatesdStateVariables:\n"; tardigradeVectorTools::print( dStateVariableEvolutionRatesdStateVariables );
    std::cout << "R.getdStateVariableEvolutionRatesdStateVariables:\n"; tardigradeVectorTools::print( *R.getdStateVariableEvolutionRatesdStateVariables( ) );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousStateVariableEvolutionRates( );

        floatVector dm = *Rm.getPreviousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousStateVariableEvolutionRatesdPreviousCauchyStress, *R.getdPreviousStateVariableEvolutionRatesdPreviousCauchyStress( ) ) );
    std::cout << "dPreviousStateVariableEvolutionRatesdPreviousCauchyStress:\n"; tardigradeVectorTools::print( dPreviousStateVariableEvolutionRatesdPreviousCauchyStress );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousStateVariableEvolutionRates( );

        floatVector dm = *Rm.getPreviousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousStateVariableEvolutionRatesdPreviousF, *R.getdPreviousStateVariableEvolutionRatesdPreviousF( ) ) );
    std::cout << "dPreviousStateVariableEvolutionRatesdPreviousF:\n"; tardigradeVectorTools::print( dPreviousStateVariableEvolutionRatesdPreviousF );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousStateVariableEvolutionRates( );

        floatVector dm = *Rm.getPreviousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousStateVariableEvolutionRatesdPreviousSubFs, *R.getdPreviousStateVariableEvolutionRatesdPreviousSubFs( ) ) );
    std::cout << "dPreviousStateVariableEvolutionRatesdPreviousSubFs:\n"; tardigradeVectorTools::print( dPreviousStateVariableEvolutionRatesdPreviousSubFs );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( previousTemperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature + delta, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature - delta, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousStateVariableEvolutionRates( );

        floatVector dm = *Rm.getPreviousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousStateVariableEvolutionRatesdPreviousT, *R.getdPreviousStateVariableEvolutionRatesdPreviousT( ) ) );
    std::cout << "dPreviousStateVariableEvolutionRatesdPreviousT:\n"; tardigradeVectorTools::print( dPreviousStateVariableEvolutionRatesdPreviousT );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ 18 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( previousStateVariables[ 18 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPreviousStateVariableEvolutionRates( );

        floatVector dm = *Rm.getPreviousStateVariableEvolutionRates( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPreviousStateVariableEvolutionRatesdPreviousStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 18 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPreviousStateVariableEvolutionRatesdPreviousStateVariables, *R.getdPreviousStateVariableEvolutionRatesdPreviousStateVariables( ) ) );
    std::cout << "dPreviousStateVariableEvolutionRatesdPreviousStateVariables:\n"; tardigradeVectorTools::print( dPreviousStateVariableEvolutionRatesdPreviousStateVariables );
    std::cout << "R.getdPreviousStateVariableEvolutionRatesdPreviousStateVariables:\n"; tardigradeVectorTools::print( *R.getdPreviousStateVariableEvolutionRatesdPreviousStateVariables( ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getPlasticDeformationGradient ){
    /*!
     * Test of computing the plastic deformation gradient
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

            floatVector velocityGradient = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector previousVelocityGradient = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9 };

            floatMatrix dLdCauchy = initializeMatrix( 9, 9, 0.1 );

            floatMatrix dLdF      = initializeMatrix( 9, 9, -0.34 );

            floatMatrix dLdSubFs  = initializeMatrix( 9, 18, 0.23 );

            floatVector dLdT      = initializeVector( 9, 5.67 );

            floatMatrix dLdXi     = initializeMatrix( 9, 1, -0.26 );

            floatMatrix dPreviousLdPreviousCauchy = initializeMatrix( 9, 9, 0.27 );

            floatMatrix dPreviousLdPreviousF      = initializeMatrix( 9, 9, -0.44 );

            floatMatrix dPreviousLdPreviousSubFs  = initializeMatrix( 9, 18, 0.29 );

            floatVector dPreviousLdPreviousT      = initializeVector( 9, 2.11 );

            floatMatrix dPreviousLdPreviousXi     = initializeMatrix( 9, 1, -0.19 );

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

            virtual void setVelocityGradient( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setPreviousVelocityGradient( previousVelocityGradient );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setVelocityGradient( velocityGradient );

                }

            }

            virtual void setdVelocityGradientdCauchyStress( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousVelocityGradientdPreviousCauchyStress( dPreviousLdPreviousCauchy );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdVelocityGradientdCauchyStress( dLdCauchy );

                }

            }

            virtual void setdVelocityGradientdF( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousVelocityGradientdPreviousF( dPreviousLdPreviousF );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdVelocityGradientdF( dLdF );

                }

            }

            virtual void setdVelocityGradientdSubFs( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousVelocityGradientdPreviousSubFs( dPreviousLdPreviousSubFs );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdVelocityGradientdSubFs( dLdSubFs );

                }

            }

            virtual void setdVelocityGradientdT( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousVelocityGradientdPreviousT( dPreviousLdPreviousT );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdVelocityGradientdT( dLdT );

                }

            }

            virtual void setdVelocityGradientdStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousVelocityGradientdPreviousStateVariables( dPreviousLdPreviousXi );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdVelocityGradientdStateVariables( dLdXi );

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters, 0.6 );

    residualMock Rjac1( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters, 0.6 );

    residualMock Rjac2( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters, 0.6 );

    floatVector answer = { 0.41617760, -0.32005614, -0.05658586,
                          -0.33457515,  0.67081398, -0.19950498,
                          -0.08532789, -0.43831590,  0.05757591 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getPlasticDeformationGradient( ) ) );

    try{
    Rjac1.getdPlasticDeformationGradientdCauchyStress( );

    Rjac2.getdPlasticDeformationGradientdPreviousCauchyStress( );
    }
    catch(std::exception &e){tardigradeErrorTools::printNestedExceptions(e);}

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *Rjac1.getPlasticDeformationGradient( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *Rjac2.getPlasticDeformationGradient( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getPlasticDeformationGradient_jacobian ){
    /*!
     * Test of computing the plastic deformation gradient
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousCauchyStress( ){

                tardigradeHydra::residualBase::setPreviousCauchyStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R_ngrad( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    R.getdYieldFunctiondCauchyStress( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *R.getYieldFunction( ), *R_ngrad.getYieldFunction( ) ) );

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

        floatVector dp = *R.getPlasticDeformationGradient( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatVector dm = *R.getPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticDeformationGradientdCauchyStress, *R.getdPlasticDeformationGradientdCauchyStress( ) ) );
    std::cout << "dPlasticDeformationGradientdCauchyStress:\n"; tardigradeVectorTools::print( dPlasticDeformationGradientdCauchyStress );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticDeformationGradient( );

        floatVector dm = *Rm.getPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticDeformationGradientdF, *R.getdPlasticDeformationGradientdF( ), 1e-5, 1e-5 ) );
    std::cout << "dPlasticDeformationGradientdF:\n"; tardigradeVectorTools::print( dPlasticDeformationGradientdF );
    std::cout << "R.getdPlasticDeformationGradientdF:\n"; tardigradeVectorTools::print( *R.getdPlasticDeformationGradientdF( ) );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticDeformationGradient( );

        floatVector dm = *Rm.getPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticDeformationGradientdSubFs, *R.getdPlasticDeformationGradientdSubFs( ), 1e-5, 1e-5 ) );
    std::cout << "dPlasticDeformationGradientdSubFs:\n"; tardigradeVectorTools::print( dPlasticDeformationGradientdSubFs );
    std::cout << "R.getdPlasticDeformationGradientdSubFs:\n"; tardigradeVectorTools::print( *R.getdPlasticDeformationGradientdSubFs( ) );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature + delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature - delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticDeformationGradient( );

        floatVector dm = *Rm.getPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticDeformationGradientdT, *R.getdPlasticDeformationGradientdT( ) ) );
    std::cout << "dPlasticDeformationGradientdT:\n"; tardigradeVectorTools::print( dPlasticDeformationGradientdT );
    std::cout << "R.getdPlasticDeformationGradientdT:\n"; tardigradeVectorTools::print( *R.getdPlasticDeformationGradientdT( ) );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ 27 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( unknownVector[ 27 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticDeformationGradient( );

        floatVector dm = *Rm.getPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 27 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticDeformationGradientdStateVariables, *R.getdPlasticDeformationGradientdStateVariables( ), 1e-5, 1e-5 ) );
    std::cout << "dPlasticDeformationGradientdStateVariables:\n"; tardigradeVectorTools::print( dPlasticDeformationGradientdStateVariables );
    std::cout << "R.getdPlasticDeformationGradientdStateVariables:\n"; tardigradeVectorTools::print( *R.getdPlasticDeformationGradientdStateVariables( ) );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticDeformationGradient( );

        floatVector dm = *Rm.getPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdPreviousCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticDeformationGradientdPreviousCauchyStress, *R.getdPlasticDeformationGradientdPreviousCauchyStress( ) ) );
    std::cout << "dPlasticDeformationGradientdPreviousCauchyStress:\n"; tardigradeVectorTools::print( dPlasticDeformationGradientdPreviousCauchyStress );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticDeformationGradient( );

        floatVector dm = *Rm.getPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdPreviousF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticDeformationGradientdPreviousF, *R.getdPlasticDeformationGradientdPreviousF( ) ) );
    std::cout << "dPlasticDeformationGradientdPreviousF:\n"; tardigradeVectorTools::print( dPlasticDeformationGradientdPreviousF );
    std::cout << "R.getdPlasticDeformationGradientdPreviousF:\n"; tardigradeVectorTools::print( *R.getdPlasticDeformationGradientdPreviousF( ) );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticDeformationGradient( );

        floatVector dm = *Rm.getPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdPreviousSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticDeformationGradientdPreviousSubFs, *R.getdPlasticDeformationGradientdPreviousSubFs( ) ) );
    std::cout << "dPlasticDeformationGradientdPreviousSubFs:\n"; tardigradeVectorTools::print( dPlasticDeformationGradientdPreviousSubFs );
    std::cout << "R.getdPlasticDeformationGradientdPreviousSubFs:\n"; tardigradeVectorTools::print( *R.getdPlasticDeformationGradientdPreviousSubFs( ) );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( previousTemperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature + delta, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature - delta, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticDeformationGradient( );

        floatVector dm = *Rm.getPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdPreviousT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticDeformationGradientdPreviousT, *R.getdPlasticDeformationGradientdPreviousT( ) ) );
    std::cout << "dPlasticDeformationGradientdPreviousT:\n"; tardigradeVectorTools::print( dPlasticDeformationGradientdPreviousT );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ 18 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( previousStateVariables[ 18 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticDeformationGradient( );

        floatVector dm = *Rm.getPlasticDeformationGradient( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticDeformationGradientdPreviousStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 18 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticDeformationGradientdPreviousStateVariables, *R.getdPlasticDeformationGradientdPreviousStateVariables( ) ) );
    std::cout << "dPlasticDeformationGradientdPreviousStateVariables:\n"; tardigradeVectorTools::print( dPlasticDeformationGradientdPreviousStateVariables );
    std::cout << "R.getdPlasticDeformationGradientdPreviousStateVariables:\n"; tardigradeVectorTools::print( *R.getdPlasticDeformationGradientdPreviousStateVariables( ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getPlasticStateVariables ){
    /*!
     * Test of computing the plastic state variables
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

            floatVector stateVariableEvolutionRates = { 2.3, 1.78 };

            floatVector previousStateVariableEvolutionRates = { 4.7, 0.3 };

            floatVector previousStateVariables = { 0.3, -0.4 };

            floatMatrix dXidotdC = initializeMatrix( 2, 9, 0.1 );

            floatMatrix dXidotdF = initializeMatrix( 2, 9, 0.2 );

            floatMatrix dXidotdSubFs = initializeMatrix( 2, 18, 0.3 );

            floatVector dXidotdT = initializeVector( 2, 0.4 );

            floatMatrix dXidotdXi = initializeMatrix( 2, 2, 0.5 );

            floatMatrix dPreviousXidotdPreviousC = initializeMatrix( 2, 9, 0.6 );

            floatMatrix dPreviousXidotdPreviousF = initializeMatrix( 2, 9, 0.7 );

            floatMatrix dPreviousXidotdPreviousSubFs = initializeMatrix( 2, 18, 0.8 );

            floatVector dPreviousXidotdPreviousT = initializeVector( 2, 0.9 );

            floatMatrix dPreviousXidotdPreviousXi = initializeMatrix( 2, 2, 1.0 );

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

            virtual void setStateVariableEvolutionRates( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setPreviousStateVariableEvolutionRates( previousStateVariableEvolutionRates );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setStateVariableEvolutionRates( stateVariableEvolutionRates );

                }

            }

            virtual void setdStateVariableEvolutionRatesdCauchyStress( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousStateVariableEvolutionRatesdPreviousCauchyStress( dPreviousXidotdPreviousC );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdStateVariableEvolutionRatesdCauchyStress( dXidotdC );

                }

            }

            virtual void setdStateVariableEvolutionRatesdF( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousStateVariableEvolutionRatesdPreviousF( dPreviousXidotdPreviousF );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdStateVariableEvolutionRatesdF( dXidotdF );

                }

            }

            virtual void setdStateVariableEvolutionRatesdSubFs( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousStateVariableEvolutionRatesdPreviousSubFs( dPreviousXidotdPreviousSubFs );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdStateVariableEvolutionRatesdSubFs( dXidotdSubFs );

                }

            }

            virtual void setdStateVariableEvolutionRatesdT( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousStateVariableEvolutionRatesdPreviousT( dPreviousXidotdPreviousT );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdStateVariableEvolutionRatesdT( dXidotdT );

                }

            }

            virtual void setdStateVariableEvolutionRatesdStateVariables( const bool isPrevious ) override{

                if ( isPrevious ){

                    tardigradeHydra::peryznaViscoplasticity::residual::setdPreviousStateVariableEvolutionRatesdPreviousStateVariables( dPreviousXidotdPreviousXi );

                }
                else{

                    tardigradeHydra::peryznaViscoplasticity::residual::setdStateVariableEvolutionRatesdStateVariables( dXidotdXi );

                }

            }

            virtual void setPreviousStateVariables( ) override{

                tardigradeHydra::peryznaViscoplasticity::residual::setPreviousStateVariables( previousStateVariables );

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters, 0.6 );

    residualMock R1( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters, 0.6 );

    residualMock R2( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters, 0.6 );

    floatVector answer = { 7.472, 2.2136 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getPlasticStateVariables( ) ) );

    try{
    R1.getdPlasticStateVariablesdCauchyStress( );

    R2.getdPlasticStateVariablesdPreviousCauchyStress( );
    }
    catch(std::exception &e){tardigradeErrorTools::printNestedExceptions( e );}

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R1.getPlasticStateVariables( ) ) );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R2.getPlasticStateVariables( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getPlasticStateVariables_jacobian ){
    /*!
     * Test of computing the plastic deformation gradient
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousCauchyStress( ){

                tardigradeHydra::residualBase::setPreviousCauchyStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R_ngrad( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    residualMock R( &hydra, 9, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    R.getdYieldFunctiondCauchyStress( );

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( *R.getYieldFunction( ), *R_ngrad.getYieldFunction( ) ) );

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

        floatVector dp = *R.getPlasticStateVariables( );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector - delta );

        floatVector dm = *R.getPlasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticStateVariablesdCauchyStress, *R.getdPlasticStateVariablesdCauchyStress( ) ) );
    std::cout << "dPlasticStateVariablesdCauchyStress:\n"; tardigradeVectorTools::print( dPlasticStateVariablesdCauchyStress );

    // Jacobians w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( deformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient + delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient - delta, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticStateVariables( );

        floatVector dm = *Rm.getPlasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticStateVariablesdF, *R.getdPlasticStateVariablesdF( ), 1e-5, 1e-5 ) );
    std::cout << "dPlasticStateVariablesdF:\n"; tardigradeVectorTools::print( dPlasticStateVariablesdF );
    std::cout << "R.getdPlasticStateVariablesdF:\n"; tardigradeVectorTools::print( *R.getdPlasticStateVariablesdF( ) );

    // Jacobians w.r.t. the sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ i + 9 ] = eps * std::fabs( unknownVector[ i + 9 ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticStateVariables( );

        floatVector dm = *Rm.getPlasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i + 9 ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticStateVariablesdSubFs, *R.getdPlasticStateVariablesdSubFs( ), 1e-5, 1e-5 ) );
    std::cout << "dPlasticStateVariablesdSubFs:\n"; tardigradeVectorTools::print( dPlasticStateVariablesdSubFs );
    std::cout << "R.getdPlasticStateVariablesdSubFs:\n"; tardigradeVectorTools::print( *R.getdPlasticStateVariablesdSubFs( ) );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( temperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature + delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature - delta, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticStateVariables( );

        floatVector dm = *Rm.getPlasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticStateVariablesdT, *R.getdPlasticStateVariablesdT( ) ) );
    std::cout << "dPlasticStateVariablesdT:\n"; tardigradeVectorTools::print( dPlasticStateVariablesdT );
    std::cout << "R.getdPlasticStateVariablesdT:\n"; tardigradeVectorTools::print( *R.getdPlasticStateVariablesdT( ) );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( unknownVector.size( ), 0 );

        delta[ 27 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( unknownVector[ 27 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector + delta );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector - delta );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticStateVariables( );

        floatVector dm = *Rm.getPlasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 27 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticStateVariablesdStateVariables, *R.getdPlasticStateVariablesdStateVariables( ), 1e-5, 1e-5 ) );
    std::cout << "dPlasticStateVariablesdStateVariables:\n"; tardigradeVectorTools::print( dPlasticStateVariablesdStateVariables );
    std::cout << "R.getdPlasticStateVariablesdStateVariables:\n"; tardigradeVectorTools::print( *R.getdPlasticStateVariablesdStateVariables( ) );

    // Jacobians w.r.t. the Previous Cauchy stress
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( hydra.elasticity.previousCauchyStress[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydra_p._local_deltaPreviousCauchyStress = delta;

        hydra_m._local_deltaPreviousCauchyStress = -delta;

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticStateVariables( );

        floatVector dm = *Rm.getPlasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdPreviousCauchyStress[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticStateVariablesdPreviousCauchyStress, *R.getdPlasticStateVariablesdPreviousCauchyStress( ) ) );
    std::cout << "dPlasticStateVariablesdPreviousCauchyStress:\n"; tardigradeVectorTools::print( dPlasticStateVariablesdPreviousCauchyStress );

    // Jacobians w.r.t. the previous deformation gradient
    for ( unsigned int i = 0; i < 9; i++ ){

        floatVector delta( 9, 0 );

        delta[ i ] = eps * std::fabs( previousDeformationGradient[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient + delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient - delta,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticStateVariables( );

        floatVector dm = *Rm.getPlasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdPreviousF[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticStateVariablesdPreviousF, *R.getdPlasticStateVariablesdPreviousF( ) ) );
    std::cout << "dPlasticStateVariablesdPreviousF:\n"; tardigradeVectorTools::print( dPlasticStateVariablesdPreviousF );
    std::cout << "R.getdPlasticStateVariablesdPreviousF:\n"; tardigradeVectorTools::print( *R.getdPlasticStateVariablesdPreviousF( ) );

    // Jacobians w.r.t. the previous sub-deformation gradients
    for ( unsigned int i = 0; i < 18; i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ i ] = eps * std::fabs( unknownVector[ i ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticStateVariables( );

        floatVector dm = *Rm.getPlasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdPreviousSubFs[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ i ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticStateVariablesdPreviousSubFs, *R.getdPlasticStateVariablesdPreviousSubFs( ) ) );
    std::cout << "dPlasticStateVariablesdPreviousSubFs:\n"; tardigradeVectorTools::print( dPlasticStateVariablesdPreviousSubFs );
    std::cout << "R.getdPlasticStateVariablesdPreviousSubFs:\n"; tardigradeVectorTools::print( *R.getdPlasticStateVariablesdPreviousSubFs( ) );

    // Jacobians w.r.t. the temperature
    for ( unsigned int i = 0; i < 1; i++ ){

        floatType delta = eps * std::fabs( previousTemperature ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature + delta, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature - delta, deformationGradient, previousDeformationGradient,
                               previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticStateVariables( );

        floatVector dm = *Rm.getPlasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdPreviousT[ j ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticStateVariablesdPreviousT, *R.getdPlasticStateVariablesdPreviousT( ) ) );
    std::cout << "dPlasticStateVariablesdPreviousT:\n"; tardigradeVectorTools::print( dPlasticStateVariablesdPreviousT );

    // Jacobians w.r.t. the state variables
    for ( unsigned int i = 0; i < hydra.stateVariableIndices.size( ); i++ ){

        floatVector delta( previousStateVariables.size( ), 0 );

        delta[ 18 + hydra.stateVariableIndices[ i ] ] += eps * std::fabs( previousStateVariables[ 18 + hydra.stateVariableIndices[ i ] ] ) + eps;

        hydraBaseMock hydra_p( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        hydraBaseMock hydra_m( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                               previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_p, unknownVector );

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra_m, unknownVector );

        residualMock Rp( &hydra_p, 9, 1, hydra_p.stateVariableIndices, hydra_p.viscoPlasticParameters );

        residualMock Rm( &hydra_m, 9, 1, hydra_m.stateVariableIndices, hydra_m.viscoPlasticParameters );

        floatVector dp = *Rp.getPlasticStateVariables( );

        floatVector dm = *Rm.getPlasticStateVariables( );

        for ( unsigned int j = 0; j < nvar; j++ ){

            dPlasticStateVariablesdPreviousStateVariables[ j ][ i ] = ( dp[ j ] - dm[ j ] ) / ( 2 * delta[ 18 + hydra.stateVariableIndices[ i ] ] );

        }

    }

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( dPlasticStateVariablesdPreviousStateVariables, *R.getdPlasticStateVariablesdPreviousStateVariables( ) ) );
    std::cout << "dPlasticStateVariablesdPreviousStateVariables:\n"; tardigradeVectorTools::print( dPlasticStateVariablesdPreviousStateVariables );
    std::cout << "R.getdPlasticStateVariablesdPreviousStateVariables:\n"; tardigradeVectorTools::print( *R.getdPlasticStateVariablesdPreviousStateVariables( ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getResidual ){
    /*!
     * Test of computing the plastic state variables
     */

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

            floatVector plasticDeformationGradient = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

            floatVector plasticStateVariables = { 0.9 };

        private:

            virtual void setPlasticDeformationGradient( ) override{

                tardigradeHydra::peryznaViscoplasticity::residual::setPlasticDeformationGradient( plasticDeformationGradient );

            }

            virtual void setPlasticStateVariables( ) override{

                tardigradeHydra::peryznaViscoplasticity::residual::setPlasticStateVariables( plasticStateVariables );

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
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    residualMock R( &hydra, 10, 1, hydra.stateVariableIndices, hydra.viscoPlasticParameters, 0.6 );

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector( hydra, unknownVector );

    floatVector answer = { -0.2, 2, 3, 4, 4, 6, 7, 8, 8, -5.1 };

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( answer, *R.getResidual( ) ) );

}

BOOST_AUTO_TEST_CASE( test_residual_getJacobian ){
    /*!
     * Test of computing the jacobian
     */

    class stressMock : public tardigradeHydra::residualBase {

        public:

            using tardigradeHydra::residualBase::residualBase;

            floatVector previousCauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        private:

            virtual void setPreviousCauchyStress( ){

                tardigradeHydra::residualBase::setPreviousCauchyStress( previousCauchyStress );

            }

    };

    class residualMock : public tardigradeHydra::peryznaViscoplasticity::residual {

        public:

            using tardigradeHydra::peryznaViscoplasticity::residual::residual;

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

    BOOST_CHECK( tardigradeVectorTools::fuzzyEquals( jacobian, *R.getJacobian( ) ) );
    std::cout << "jacobian:\n"; tardigradeVectorTools::print( jacobian );
    std::cout << "R.getJacobian:\n"; tardigradeVectorTools::print( *R.getJacobian( ) );

}
