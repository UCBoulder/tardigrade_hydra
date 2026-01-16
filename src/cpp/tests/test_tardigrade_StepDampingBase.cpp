/**
  * \file test_tardigrade_StepDampingBase.cpp
  *
  * Tests for tardigrade_StepDampingBase
  */

#include"tardigrade_StepDampingBase.h"
#include"tardigrade_hydra.h"
#include"tardigrade_SolverBase.h"
#include"tardigrade_SolverStepBase.h"

#define BOOST_TEST_MODULE test_tardigrade_StepDampingBase
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

        class hydraBaseTester{

            public:

                static void set_residual( hydraBase &hydra, const tardigradeHydra::floatVector &value ){

                    hydra._residual.second = value;
                    hydra._residual.first = true;

                    hydra.addIterationData( &hydra._residual );

                }

                static void set_unknownVector( hydraBase &hydra, const tardigradeHydra::floatVector &value ){

                    hydra._X.second = value;
                    hydra._X.first = true;

                }

                static void set_flatJacobian( hydraBase &hydra, const tardigradeHydra::floatVector &value ){

                    hydra._jacobian.second = value;
                    hydra._jacobian.first = true;

                    hydra.addIterationData( &hydra._jacobian );
                }

                static void initializeUnknownVector( hydraBase &hydra ){

                    BOOST_CHECK_NO_THROW( hydra.initializeUnknownVector( ) );

                }

        };

        class StepDampingBaseTester{

            public:

                static void checkLSAlpha( StepDampingBase &damping ){

                    BOOST_CHECK( damping._lsAlpha == damping.getLSAlpha( ) );

                }

                static void checkGradientRho( StepDampingBase &damping ){

                    BOOST_CHECK( damping._gradientRho == damping.getGradientRho( ) );

                }

                static void checkGradientP( StepDampingBase &damping ){

                    BOOST_CHECK( damping._gradientP == damping.getGradientP( ) );

                }

                static void checkGradientBeta( StepDampingBase &damping ){

                    BOOST_CHECK( damping._gradientBeta == damping.getGradientBeta( ) );

                }

                static void checkGradientSigma( StepDampingBase &damping ){

                    BOOST_CHECK( damping._gradientSigma == damping.getGradientSigma( ) );

                }

                static void checkMaxGradientIterations( StepDampingBase &damping ){

                    BOOST_CHECK( damping._maxGradientIterations == damping.getMaxGradientIterations( ) );

                }

                static void checkMuk( StepDampingBase &damping ){

                    BOOST_CHECK( damping._mu_k == damping.getMuk( ) );

                }

//                static void checkLMMu( SolverStepBase &step ){
//
//                    BOOST_CHECK( step._lm_mu == step.getLMMu( ) );
//
//                }
//
//                static void checkUseLevenbergMarquardt( SolverStepBase &step ){
//
//                    BOOST_CHECK( step._use_LM_step == step.getUseLevenbergMarquardt( ) );
//
//                }

                static void checkUseGradientDescent( StepDampingBase &damping ){

                    BOOST_CHECK( damping._use_gradient_descent == damping.getUseGradientDescent( ) );

                }

//                static void setMuk( SolverStepBase &step, const tardigradeHydra::floatType &value ){
//
//                    step.setMuk( value );
//
//                }

        };

    }

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_getLSAlpha, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::StepDampingBase damping;
    damping.setLSAlpha( 1.23 );

    tardigradeHydra::unit_test::StepDampingBaseTester::checkLSAlpha( damping );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_getLSResidualNorm, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            void setSolver( tardigradeHydra::SolverBase *_solver ){ solver = _solver; }
            tardigradeHydra::SolverBase *getSolver( ){ return solver; }

    };

    hydraBaseMock hydra;
    tardigradeHydra::SolverBase solver;
    tardigradeHydra::SolverStepBase step;
    tardigradeHydra::PreconditionerBase preconditioner;
    tardigradeHydra::StepDampingBase damping;

    hydra.setSolver( &solver );

    solver.hydra = &hydra;
    solver.step  = &step;
    solver.preconditioner = &preconditioner;

    step.setSolver( &solver );
    preconditioner.setSolver( &solver );

    step.damping = &damping;
    damping.step = &step;

    tardigradeHydra::floatVector residual = { 1, 2, 3 };

    tardigradeHydra::floatType lsResidualNormAnswer = tardigradeVectorTools::l2norm( residual );
      tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    BOOST_TEST( *damping.getLSResidualNorm( ) == lsResidualNormAnswer );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_checkLSConvergence, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            tardigradeHydra::SolverBase *getSolver( ){ return solver; }

    };

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

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::SolverStepBase step;
    tardigradeHydra::StepDampingBase damping;

    step.damping = &damping;
    damping.step = &step;

    damping.setLSAlpha( 1e-4 );
    damping.setMaxLSIterations( 5 );

    hydra.getSolver( )->step = &step;
    step.setSolver( hydra.getSolver( ) );

    tardigradeHydra::floatVector residual = { 1, 2, -3, 0 };

    tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    BOOST_CHECK( !damping.checkLSConvergence( ) );

    residual = { 0.1, 2, -3, 0 };

    tardigradeHydra::unit_test::hydraBaseTester::set_residual( hydra, residual );

    BOOST_CHECK( damping.checkLSConvergence( ) );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_getUseGradientDescent, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::StepDampingBase damping;

    tardigradeHydra::unit_test::StepDampingBaseTester::checkUseGradientDescent( damping );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_setUseGradientDescent, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::StepDampingBase damping;

    damping.setUseGradientDescent( true );

    BOOST_TEST( true == damping.getUseGradientDescent( ) );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_getGradientRho, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::StepDampingBase damping;

    tardigradeHydra::unit_test::StepDampingBaseTester::checkGradientRho( damping );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_getGradientP, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::StepDampingBase damping;

    tardigradeHydra::unit_test::StepDampingBaseTester::checkGradientP( damping );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_getGradientBeta, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::StepDampingBase damping;

    tardigradeHydra::unit_test::StepDampingBaseTester::checkGradientBeta( damping );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_getGradientSigma, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::StepDampingBase damping;

    tardigradeHydra::unit_test::StepDampingBaseTester::checkGradientSigma( damping );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_setGradientRho, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::StepDampingBase damping;

    damping.setGradientRho( 123.4 );

    BOOST_TEST( 123.4 == damping.getGradientRho( ) );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_setGradientP, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::StepDampingBase damping;

    damping.setGradientP( 123.4 );

    BOOST_TEST( 123.4 == damping.getGradientP( ) );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_setGradientBeta, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::StepDampingBase damping;

    damping.setGradientBeta( 123.4 );

    BOOST_TEST( 123.4 == damping.getGradientBeta( ) );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_setGradientSigma, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::StepDampingBase damping;

    damping.setGradientSigma( 123.4 );

    BOOST_TEST( 123.4 == damping.getGradientSigma( ) );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_getMaxGradientIterations, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::StepDampingBase damping;

    tardigradeHydra::unit_test::StepDampingBaseTester::checkMaxGradientIterations( damping );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_setMaxGradientIterations, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::StepDampingBase damping;

    damping.setMaxGradientIterations( 123 );

    BOOST_TEST( 123 == damping.getMaxGradientIterations( ) );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_getMuk, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::StepDampingBase step;

    tardigradeHydra::unit_test::StepDampingBaseTester::checkMuk( step );

}

BOOST_AUTO_TEST_CASE( test_StepDampingBase_setMuk, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class StepDampingBaseMock : public tardigradeHydra::StepDampingBase{

        using tardigradeHydra::StepDampingBase::StepDampingBase;

        public:

            void public_setMuk( const tardigradeHydra::floatType &v ){

                setMuk( v );

            }

    };

    StepDampingBaseMock damping;

    damping.public_setMuk( 123.4 );

    BOOST_TEST( 123.4 == damping.getMuk( ) );

}

BOOST_AUTO_TEST_CASE( test_setResidualNorm, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
    /*!
     * Test setting the norm of the residual
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

            tardigradeHydra::floatVector jacobian = { -0.15378708,  0.9615284 ,  0.36965948, -0.0381362 , -0.21576496,
                                     -0.31364397,  0.45809941, -0.12285551, -0.88064421, -0.20391149,
                                      0.47599081, -0.63501654, -0.64909649,  0.06310275,  0.06365517,
                                      0.26880192,  0.69886359,  0.44891065,  0.22204702,  0.44488677,
                                     -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225 };

            using tardigradeHydra::hydraBase::hydraBase;

            void set_solver( tardigradeHydra::SolverBase *_solver ){ solver = _solver; }

            virtual void formNonLinearResidual( ) override{

                tardigradeHydra::floatVector residual( 5, 0 );

                auto *X = getUnknownVector( );

                for ( unsigned int i = 0; i < 5; i++ ){
                    for ( unsigned int j = 0; j < 5; j++ ){
                        residual[ i ] += jacobian[5*i+j] * ( *X )[j];
                    }
                }

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

            }

            virtual void formNonLinearDerivatives( ) override{

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, jacobian );

            }

            virtual const unsigned int getNumUnknowns( ) override{ return 5; }

    };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::floatVector unknownVector = { 0.39293837, -0.42772133, -0.54629709,  0.10262954,  0.43893794 };

    tardigradeHydra::floatType answer = 1.6716509825117496;

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra, unknownVector );

    tardigradeHydra::SolverBase solver;
    tardigradeHydra::SolverStepBase step;
    tardigradeHydra::StepDampingBase damping;

    hydra.set_solver( &solver );
    solver.hydra = &hydra;
    solver.step = &step;
    step.setSolver( &solver );
    step.damping = &damping;
    damping.step = &step;

    BOOST_TEST( answer == *damping.get_residualNorm( ) );

    tardigradeHydra::floatVector dResidualNormdX( 5, 0 );

    {

        constexpr unsigned int NUM_INPUTS = 5;
        tardigradeHydra::floatVector x = unknownVector;

        tardigradeHydra::floatType eps = 1e-6;

        for ( unsigned int i = 0; i < NUM_INPUTS; ++i ){

            tardigradeHydra::floatType delta = eps * std::fabs( x[ i ] ) + eps;

            tardigradeHydra::floatVector xp = x;

            tardigradeHydra::floatVector xm = x;

            xp[ i ] += delta;

            xm[ i ] -= delta;

            hydraBaseMock hydrap( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                  { }, { },
                                  previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

            hydraBaseMock hydram( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                                  { }, { },
                                  previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

            tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydrap, xp );

            tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydram, xm );

            tardigradeHydra::SolverBase solverp;
            tardigradeHydra::SolverBase solverm;

            tardigradeHydra::SolverStepBase stepp;
            tardigradeHydra::SolverStepBase stepm;

            tardigradeHydra::StepDampingBase dampingp;
            tardigradeHydra::StepDampingBase dampingm;

            hydrap.set_solver( &solverp );
            solverp.hydra = &hydrap;
            solverp.step  = &stepp;
            stepp.setSolver( &solverp );
            stepp.damping = &dampingp;
            dampingp.step = &stepp;

            hydram.set_solver( &solverm );
            solverm.hydra = &hydram;
            solverm.step  = &stepm;
            stepm.setSolver( &solverm );
            stepm.damping = &dampingm;
            dampingm.step = &stepm;

            dResidualNormdX[ i ] = ( *dampingp.get_residualNorm( ) - *dampingm.get_residualNorm( ) ) / ( 2 * delta );

        }

        BOOST_TEST( dResidualNormdX == *damping.get_dResidualNormdX( ), CHECK_PER_ELEMENT );

    }

}
