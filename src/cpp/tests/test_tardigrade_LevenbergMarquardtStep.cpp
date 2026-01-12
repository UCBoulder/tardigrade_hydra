/**
  * \file test_tardigrade_LevenbergMarquardtStep.cpp
  *
  * Tests for tardigrade_LevenbergMarquardtStep
  */

#include"tardigrade_LevenbergMarquardtStep.h"
#include"tardigrade_ResidualBase.h"
#include"tardigrade_hydra.h"

#define BOOST_TEST_MODULE test_tardigrade_LevenbergMarquardtStep
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element( )

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

namespace tardigradeHydra{

    namespace unit_test{

        class SolverStepBaseTester{

            public:

                static void setMuk( SolverStepBase &step, const floatType &value ){

                    step.setMuk( value );

                }

        };

        class hydraBaseTester{

            public:

                static void set_residual( hydraBase &hydra, const floatVector &value ){

                    hydra._residual.second = value;
                    hydra._residual.first = true;

                    hydra.addIterationData( &hydra._residual );

                }

                static void set_unknownVector( hydraBase &hydra, const floatVector &value ){

                    hydra._X.second = value;
                    hydra._X.first = true;

                }

                static void set_flatJacobian( hydraBase &hydra, const floatVector &value ){

                    hydra._jacobian.second = value;
                    hydra._jacobian.first = true;

                    hydra.addIterationData( &hydra._jacobian );
                }

        };

    }

}

BOOST_AUTO_TEST_CASE( test_LevenbergMarquardtStep_getNonlinearLMTerms, * boost::unit_test::tolerance( 1e-5 ) ){
    /*!
     * Test checking the gradient convergence
     */

    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                        0.10262954,  0.43893794, -0.15378708,
                                        0.9615284 ,  0.36965948, -0.0381362 };

    tardigradeHydra::floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081, -0.63501654, -0.64909649 };

    tardigradeHydra::floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                           0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                           0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                           0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                           0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                           0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                           0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                           0.54506801, 0.34276383, 0.30412079 }; 

    tardigradeHydra::floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::floatVector residual = { 1, 2, 3, 4, 5 };

            tardigradeHydra::floatVector jacobian = { -0.15378708,  0.9615284 ,  0.36965948, -0.0381362 , -0.21576496,
                                     -0.31364397,  0.45809941, -0.12285551, -0.88064421, -0.20391149,
                                      0.47599081, -0.63501654, -0.64909649,  0.06310275,  0.06365517,
                                      0.26880192,  0.69886359,  0.44891065,  0.22204702,  0.44488677,
                                     -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225 };

            tardigradeHydra::floatType mu_k = 1.34;

            using tardigradeHydra::hydraBase::hydraBase;

            virtual void formNonLinearResidual( ) override{

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

            }

            virtual void formNonLinearDerivatives( ) override{

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, jacobian );

                tardigradeHydra::unit_test::SolverStepBaseTester::setMuk( *(solver->step), mu_k );

            }

            virtual void decomposeUnknownVector( ) override{ return; }
            virtual const unsigned int getNumUnknowns( ) override{ return 5; }

            tardigradeHydra::SolverBase *getSolver( ){ return solver; }

    };

    class LevenbergMarquardtStepMock : public tardigradeHydra::LevenbergMarquardtStep{

        public:

            using tardigradeHydra::LevenbergMarquardtStep::LevenbergMarquardtStep;

    };

    tardigradeHydra::floatVector answerRHS = { -0.04830576,  1.38601851, -2.74506611, -2.78478784,  2.6566859 };

    tardigradeHydra::floatVector answerLHS = {  1.88621891, -0.30808058, -0.01417759,  0.51788095,  0.15427055,
                              -0.30808058,  3.44245776,  1.17530079, -0.21093811, -0.10279234,
                              -0.01417759,  1.17530079,  2.40995212,  0.377036  , -0.03867597,
                               0.51788095, -0.21093811,  0.377036  ,  2.34049101,  0.18253039,
                               0.15427055, -0.10279234, -0.03867597,  0.18253039,  1.69872961 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    LevenbergMarquardtStepMock step;

    hydra.getSolver( )->step = &step;
    step.setSolver( hydra.getSolver( ) );

    tardigradeHydra::floatVector unknownVector = { 0.39293837, -0.42772133, -0.54629709,  0.10262954,  0.43893794 };

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra, unknownVector );

    BOOST_TEST( answerRHS == *hydra.getSolver( )->step->getNonlinearRHS( ), CHECK_PER_ELEMENT );

    BOOST_TEST( answerLHS == *hydra.getSolver( )->step->getFlatNonlinearLHS( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_LevenbergMarquardtStep_setBaseQuantities, * boost::unit_test::tolerance( 1e-5 ) ){
    /*!
     * Test checking the gradient convergence
     */

    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = { 0.39293837, -0.42772133, -0.54629709,
                                                         0.10262954,  0.43893794, -0.15378708,
                                                         0.9615284 ,  0.36965948, -0.0381362 };

    tardigradeHydra::floatVector previousDeformationGradient = { -0.21576496, -0.31364397,  0.45809941,
                                                                 -0.12285551, -0.88064421, -0.20391149,
                                                                  0.47599081, -0.63501654, -0.64909649 };

    tardigradeHydra::floatVector previousStateVariables = { 0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532,
                                                            0.61102351, 0.72244338, 0.32295891, 0.36178866, 0.22826323,
                                                            0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276,
                                                            0.4936851 , 0.42583029, 0.31226122, 0.42635131, 0.89338916,
                                                            0.94416002, 0.50183668, 0.62395295, 0.1156184 , 0.31728548,
                                                            0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979,
                                                            0.51948512, 0.61289453, 0.12062867, 0.8263408 , 0.60306013,
                                                            0.54506801, 0.34276383, 0.30412079 }; 

    tardigradeHydra::floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    class LevenbergMarquardtStepMock : public tardigradeHydra::LevenbergMarquardtStep{

        public:

            tardigradeHydra::floatType rnorm = 10.3;

            tardigradeHydra::floatVector dRNormdX = { 1, 2, 3 };

            virtual void setResidualNorm( ) override{

                solver->step->set_residualNorm( rnorm );

            }

            virtual void setdResidualNormdX( ) override{

                solver->step->set_dResidualNormdX( dRNormdX );;

            }

            void runSetBaseQuantities( ){ setBaseQuantities( ); }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            tardigradeHydra::floatVector residual = { 1, 2, 3, 4, 5 };

            tardigradeHydra::floatVector jacobian = { -0.15378708,  0.9615284 ,  0.36965948, -0.0381362 , -0.21576496,
                                                      -0.31364397,  0.45809941, -0.12285551, -0.88064421, -0.20391149,
                                                       0.47599081, -0.63501654, -0.64909649,  0.06310275,  0.06365517,
                                                       0.26880192,  0.69886359,  0.44891065,  0.22204702,  0.44488677,
                                                      -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225 };

            tardigradeHydra::floatType mu_k = 1.34;

            using tardigradeHydra::hydraBase::hydraBase;

            void set_solver( tardigradeHydra::SolverBase *_solver ){ solver = _solver; }

            virtual void formNonLinearResidual( ) override{

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

            }

            virtual void formNonLinearDerivatives( ) override{

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, jacobian );

                tardigradeHydra::unit_test::SolverStepBaseTester::setMuk( *(solver->step), mu_k );

            }

            virtual void decomposeUnknownVector( ) override{ return; }
            virtual const unsigned int getNumUnknowns( ) override{ return 5; }

    };

    tardigradeHydra::floatType answer1 = 0.5 * 1e-8 * 10.3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    LevenbergMarquardtStepMock step;
    tardigradeHydra::SolverBase solver;
    hydra.set_solver( &solver );
    solver.hydra = &hydra;
    solver.step = &step;
    step.setSolver( &solver );

    step.runSetBaseQuantities( );

    BOOST_TEST( answer1 == step.getMuk( ) );

    BOOST_TEST( step.rnorm == *step.get_baseResidualNorm( ) );

    BOOST_TEST( step.dRNormdX == *step.get_basedResidualNormdX( ) );

    step.rnorm = 1e-9;

    step.setResidualNorm( );

    step.runSetBaseQuantities( );

    BOOST_TEST( step.rnorm == step.getMuk( ) );

}

