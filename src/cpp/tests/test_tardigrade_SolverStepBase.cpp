
/**
  * \file test_tardigrade_SolverStepBase.cpp
  *
  * Tests for tardigrade_SolverStepBase
  */

#include"tardigrade_SolverStepBase.h"
#include"tardigrade_hydra.h"
#include"tardigrade_SolverBase.h"

#define BOOST_TEST_MODULE test_tardigrade_SolverStepBase
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

                static void setMuk( StepDampingBase &damping, const tardigradeHydra::floatType &value ){

                    damping.setMuk( value );

                }

        };

        class SolverStepBaseTester{

            public:

                static void checkUseLevenbergMarquardt( SolverStepBase &step ){

                    BOOST_CHECK( step._use_LM_step == step.getUseLevenbergMarquardt( ) );

                }

                static void solveConstrainedQP( SolverStepBase &step, tardigradeHydra::floatVector &dx, const unsigned int kmax=100 ){

                    step.solveConstrainedQP( dx, kmax );

                }

                static void solveNewtonUpdate( SolverStepBase &step, tardigradeHydra::floatVector &deltaX ){

                    step.solveNewtonUpdate( deltaX );

                }

        };

    }

}

BOOST_AUTO_TEST_CASE( test_SolverStepBase_getUseLevenbergMarquardt, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::SolverStepBase step;

    tardigradeHydra::unit_test::SolverStepBaseTester::checkUseLevenbergMarquardt( step );

}

BOOST_AUTO_TEST_CASE( test_SolverStepBase_setUseLevenbergMarquardt, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::SolverStepBase step;

    step.setUseLevenbergMarquardt( true );

    BOOST_TEST( true == step.getUseLevenbergMarquardt( ) );

    BOOST_TEST( true == step.damping->getUseGradientDescent( ) );
}

BOOST_AUTO_TEST_CASE( test_SolverStepBase_solveConstrainedQP, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase{

        public:

            using tardigradeHydra::SolverStepBase::SolverStepBase;

            tardigradeHydra::floatVector initialUnknownVector = { 2, 0 };

        protected:

            virtual void initializeActiveConstraints( std::vector< bool > &active_constraints ) override{

                active_constraints = { false, false, true, false, true };

            }

            virtual void assembleKKTRHSVector( const tardigradeHydra::floatVector &dx, tardigradeHydra::floatVector &RHS, const std::vector< bool > &active_constraints ) override{

                RHS = tardigradeHydra::floatVector( 7, 0 );

                RHS[ 0 ] = 2 * ( initialUnknownVector[ 0 ] + dx[ 0 ] - 1.0 );
                RHS[ 1 ] = 2 * ( initialUnknownVector[ 1 ] + dx[ 1 ] - 2.5 );

                RHS[ 2 + 0 ] =  ( initialUnknownVector[ 0 ] + dx[ 0 ] ) - 2 * ( initialUnknownVector[ 1 ] + dx[ 1 ] ) + 2;
                RHS[ 2 + 1 ] = -( initialUnknownVector[ 0 ] + dx[ 0 ] ) - 2 * ( initialUnknownVector[ 1 ] + dx[ 1 ] ) + 6;
                RHS[ 2 + 2 ] = -( initialUnknownVector[ 0 ] + dx[ 0 ] ) + 2 * ( initialUnknownVector[ 1 ] + dx[ 1 ] ) + 2;
                RHS[ 2 + 3 ] =  initialUnknownVector[ 0 ] + dx[ 0 ];
                RHS[ 2 + 4 ] =  initialUnknownVector[ 1 ] + dx[ 1 ];

                for ( unsigned int i = 0; i < active_constraints.size( ); i++ ){

                    if ( !active_constraints[ i ] ){

                        RHS[ 2 + i ] = 0;

                    }

                }

            }

            virtual void assembleKKTMatrix( tardigradeHydra::floatVector &K, const std::vector< bool > &active_constraints ) override{

                const unsigned int numConstraints = solver->hydra->getNumConstraints( );

                K = tardigradeHydra::floatVector( ( 2 + 5 ) * ( 2 + 5 ), 0 );

                K[ 7 * 0 + 0 ] = 2;
                K[ 7 * 1 + 1 ] = 2;

                for ( unsigned int i = 0; i < numConstraints; i++ ){

                    if ( active_constraints[ i ] ){

                        K[ 7 * 0 + i + 2 ] = ( *solver->hydra->getConstraintJacobians( ) )[ 2 * i + 0 ];
                        K[ 7 * 1 + i + 2 ] = ( *solver->hydra->getConstraintJacobians( ) )[ 2 * i + 1 ];

                        K[ 7 * ( i + 2 ) + 0 ] = ( *solver->hydra->getConstraintJacobians( ) )[ 2 * i + 0 ];
                        K[ 7 * ( i + 2 ) + 1 ] = ( *solver->hydra->getConstraintJacobians( ) )[ 2 * i + 1 ];

                    }
                    else{

                        K[ 7 * ( i + 2 ) + i + 2 ] = 1;

                    }

                }

            }

            virtual void updateKKTMatrix( tardigradeHydra::floatVector &K, const std::vector< bool > &active_constraints ) override{

                const unsigned int numConstraints = solver->hydra->getNumConstraints( );

                for ( unsigned int i = 0; i < numConstraints; i++ ){

                    if ( active_constraints[ i ] ){

                        K[ 7 * 0 + i + 2 ] = ( *solver->hydra->getConstraintJacobians( ) )[ 2 * i + 0 ];
                        K[ 7 * 1 + i + 2 ] = ( *solver->hydra->getConstraintJacobians( ) )[ 2 * i + 1 ];

                        K[ 7 * ( i + 2 ) + 0 ] = ( *solver->hydra->getConstraintJacobians( ) )[ 2 * i + 0 ];
                        K[ 7 * ( i + 2 ) + 1 ] = ( *solver->hydra->getConstraintJacobians( ) )[ 2 * i + 1 ];

                        K[ 7 * ( i + 2 ) + i + 2 ] = 0;

                    }
                    else{

                        K[ 7 * 0 + i + 2 ] = 0;
                        K[ 7 * 1 + i + 2 ] = 0;

                        K[ 7 * ( i + 2 ) + 0 ] = 0;
                        K[ 7 * ( i + 2 ) + 1 ] = 0;

                        K[ 7 * ( i + 2 ) + i + 2 ] = 1;

                    }

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            tardigradeHydra::floatVector initialUnknownVector = { 2, 0 };

            void public_solveConstrainedQP( tardigradeHydra::floatVector &dx ){

                tardigradeHydra::unit_test::SolverStepBaseTester::solveConstrainedQP( *(solver->step), dx );

            }

            void set_solver( tardigradeHydra::SolverBase *_solver ){ solver = _solver; }

        protected:

            virtual void setConstraints( ) override{

                auto constraints = get_SetDataStorage_constraints( );

                *constraints.value = { 2, 6, 2, 0, 0 };

                for ( unsigned int i = 0; i < 5; i++ ){

                    for ( unsigned int j = 0; j < 2; j++ ){

                        ( *constraints.value )[ i ] += ( *getConstraintJacobians( ) )[ 2 * i + j ] * initialUnknownVector[ j ];

                    }

                }

            }

            virtual void setConstraintJacobians( ) override{

                auto constraintJacobians = get_SetDataStorage_constraintJacobians( );

                *constraintJacobians.value = { 1, -2,
                                              -1, -2,
                                              -1,  2,
                                               1,  0,
                                               0,  1 };

            }

            virtual const unsigned int getNumUnknowns( ) override{ return initialUnknownVector.size( ); }

            virtual const unsigned int getNumConstraints( ) override{ return 5; }

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
                                           0.54506801, 0.34276383, 0.30412079, 0.0,        0.1 }; 

    tardigradeHydra::floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::SolverBase solver;
    SolverStepBaseMock step;

    solver.hydra = &hydra;
    solver.step  = &step;
    step.setSolver( &solver );

    hydra.set_solver( &solver );

    tardigradeHydra::floatVector result = { 0, 0 };

    tardigradeHydra::floatVector answer = { 1.4, 1.7 };

    hydra.public_solveConstrainedQP( result );

    BOOST_TEST( ( result + hydra.initialUnknownVector ) == answer, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_SolverStepBase_initializeActiveConstraints, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase{

        public:

            using tardigradeHydra::SolverStepBase::SolverStepBase;

            void public_initializeActiveConstraints( std::vector< bool > &active_constraints ){

                initializeActiveConstraints( active_constraints );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            auto set_solver( tardigradeHydra::SolverBase *_solver ){ solver = _solver; }

        protected:

            virtual void setConstraints( ) override{

                auto constraints = get_SetDataStorage_constraints( );

                *constraints.value = { 2, 6, -2, 0.1, 2 };

            }

            virtual const unsigned int getNumConstraints( ) override{ return 5; }

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
                                           0.54506801, 0.34276383, 0.30412079, 0.0,        0.1 }; 

    tardigradeHydra::floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::SolverBase solver;
    SolverStepBaseMock step;

    hydra.set_solver( &solver );

    solver.hydra = &hydra;
    solver.step = &step;

    step.setSolver( &solver );

    std::vector< bool > answer = { false, false, true, false, false };

    std::vector< bool > result;

    step.public_initializeActiveConstraints( result );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_SolverStepBase_assembleKKTRHSVector, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase{

        public:

            using tardigradeHydra::SolverStepBase::SolverStepBase;

            virtual void public_assembleKKTRHSVector( const tardigradeHydra::floatVector &dx, tardigradeHydra::floatVector &RHS, const std::vector< bool > &active_constraints ){

                assembleKKTRHSVector( dx, RHS, active_constraints );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            tardigradeHydra::floatVector initialUnknownVector = { 2, 1 };

            auto public_setMuk( const tardigradeHydra::floatType &value ){ tardigradeHydra::unit_test::StepDampingBaseTester::setMuk( *(solver->step->damping), value ); }

            auto set_solver( tardigradeHydra::SolverBase *_solver ){ solver = _solver; }

        protected:

            virtual void setConstraints( ) override{

                auto constraints = get_SetDataStorage_constraints( );

                *constraints.value = { 2, 6, 2, 0, 0 };

                for ( unsigned int i = 0; i < 5; i++ ){

                    for ( unsigned int j = 0; j < 2; j++ ){

                        ( *constraints.value )[ i ] += ( *getConstraintJacobians( ) )[ 2 * i + j ] * initialUnknownVector[ j ];

                    }

                }

            }

            virtual void setConstraintJacobians( ) override{

                auto constraintJacobians = get_SetDataStorage_constraintJacobians( );

                *constraintJacobians.value = { 1, -2,
                                              -1, -2,
                                              -1,  2,
                                               1,  0,
                                               0,  1 };

            }

            virtual void formNonLinearResidual( ) override{

                tardigradeHydra::floatVector residual = { 1., 2. };

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

            }

            virtual void formNonLinearDerivatives( ) override{

                tardigradeHydra::floatVector jacobian = { std::pow( 2, 0.5 ), 0.4, -0.1, std::pow( 2, 0.5 ) };

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, jacobian );

            }

            virtual const unsigned int getNumUnknowns( ) override{ return initialUnknownVector.size( ); }

            virtual const unsigned int getNumConstraints( ) override{ return 5; }

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
                                           0.54506801, 0.34276383, 0.30412079, 0.0,        0.1 }; 

    tardigradeHydra::floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::SolverBase solver;
    SolverStepBaseMock step;

    hydra.set_solver( &solver );
    solver.hydra = &hydra;
    solver.step = &step;
    step.setSolver( &solver );

    tardigradeHydra::floatVector dx = { -0.2, 1.4 };

    hydra.public_setMuk( 0.1 );

    tardigradeHydra::floatVector result_KKTRHSVector;
    std::vector< bool > active_constraints( 5, false );

    tardigradeHydra::floatVector answer1_KKTRHSVector = { 1.38618326,  6.30757431,  0.0       ,  0.0       ,  0.0       ,  0.0       ,  0.0       };

    tardigradeHydra::floatVector answer2_KKTRHSVector = { 1.38618326,  6.30757431,  0.0       ,  0.0       ,  5.        ,  0.0       ,  2.4       };

    step.public_assembleKKTRHSVector( dx, result_KKTRHSVector, active_constraints );

    BOOST_TEST( answer1_KKTRHSVector == result_KKTRHSVector, CHECK_PER_ELEMENT );

    active_constraints[ 2 ] = true;
    active_constraints[ 4 ] = true;

    result_KKTRHSVector.clear( );

    step.public_assembleKKTRHSVector( dx, result_KKTRHSVector, active_constraints );

    BOOST_TEST( answer2_KKTRHSVector == result_KKTRHSVector, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_SolverStepBase_assembleKKTMatrix, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase{

        public:

            using tardigradeHydra::SolverStepBase::SolverStepBase;

            virtual void public_assembleKKTMatrix( tardigradeHydra::floatVector &K, const std::vector< bool > &active_constraints ){

                assembleKKTMatrix( K, active_constraints );

            }

            virtual void public_updateKKTMatrix( tardigradeHydra::floatVector &K, const std::vector< bool > &active_constraints ){

                updateKKTMatrix( K, active_constraints );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            tardigradeHydra::floatVector initialUnknownVector = { 2, 1 };

            auto public_setMuk( const tardigradeHydra::floatType &value ){ tardigradeHydra::unit_test::StepDampingBaseTester::setMuk( *(solver->step->damping), value ); }

            auto set_solver( tardigradeHydra::SolverBase *_solver ){ solver = _solver; }

        protected:

            virtual void setConstraints( ) override{

                auto constraints = get_SetDataStorage_constraints( );

                *constraints.value = { 2, 6, 2, 0, 0 };

                for ( unsigned int i = 0; i < 5; i++ ){

                    for ( unsigned int j = 0; j < 2; j++ ){

                        ( *constraints.value )[ i ] += ( *getConstraintJacobians( ) )[ 2 * i + j ] * initialUnknownVector[ j ];

                    }

                }

            }

            virtual void setConstraintJacobians( ) override{

                auto constraintJacobians = get_SetDataStorage_constraintJacobians( );

                *constraintJacobians.value = { 1, -2,
                                              -1, -2,
                                              -1,  2,
                                               1,  0,
                                               0,  1 };

            }

            virtual void formNonLinearDerivatives( ) override{

                tardigradeHydra::floatVector jacobian = { std::pow( 2, 0.5 ), 0.4, -0.1, std::pow( 2, 0.5 ) };

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, jacobian );

            }

            virtual const unsigned int getNumUnknowns( ) override{ return initialUnknownVector.size( ); }

            virtual const unsigned int getNumConstraints( ) override{ return 5; }

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
                                           0.54506801, 0.34276383, 0.30412079, 0.0,        0.1 }; 

    tardigradeHydra::floatVector parameters = { 1, 2, 3, 4, 5 };

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::SolverBase solver;
    SolverStepBaseMock step;

    hydra.set_solver( &solver );

    solver.hydra = &hydra;
    solver.step = &step;

    step.setSolver( &solver );

    hydra.public_setMuk( 0.1 );

    tardigradeHydra::floatVector result_KKT;
    std::vector< bool > active_constraints( 5, false );

    tardigradeHydra::floatVector answer1_KKTMatrix = { 2.11      , 0.42426407, 0., 0., 0., 0., 0.,
                                      0.42426407, 2.26      , 0., 0., 0., 0., 0.,
                                      0.        , 0.        , 1., 0., 0., 0., 0.,
                                      0.        , 0.        , 0., 1., 0., 0., 0.,
                                      0.        , 0.        , 0., 0., 1., 0., 0.,
                                      0.        , 0.        , 0., 0., 0., 1., 0.,
                                      0.        , 0.        , 0., 0., 0., 0., 1. };

    tardigradeHydra::floatVector answer2_KKTMatrix = { 2.11      ,  0.42426407,  0.,  0., -1.,  0.,  0.,
                                      0.42426407,  2.26      ,  0.,  0.,  2.,  0.,  1.,
                                      0.        ,  0.        ,  1.,  0.,  0.,  0.,  0.,
                                      0.        ,  0.        ,  0.,  1.,  0.,  0.,  0.,
                                     -1.        ,  2.        ,  0.,  0.,  0.,  0.,  0.,
                                      0.        ,  0.        ,  0.,  0.,  0.,  1.,  0.,
                                      0.        ,  1.        ,  0.,  0.,  0.,  0.,  0. };

    step.public_assembleKKTMatrix( result_KKT, active_constraints );

    BOOST_TEST( answer1_KKTMatrix == result_KKT, CHECK_PER_ELEMENT );

    active_constraints[ 2 ] = true;
    active_constraints[ 4 ] = true;

    step.public_updateKKTMatrix( result_KKT, active_constraints );

    BOOST_TEST( answer2_KKTMatrix == result_KKT, CHECK_PER_ELEMENT );

    result_KKT.clear( );

    step.public_assembleKKTMatrix( result_KKT, active_constraints );

    BOOST_TEST( answer2_KKTMatrix == result_KKT, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_SolverStepBase_solveNewtonUpdate, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            tardigradeHydra::floatVector flatJacobian = { 1, 2, 3, 4, 5, 6, 7, 8, 2 };

            tardigradeHydra::floatVector residual = { 1, 2, 3 };

            virtual const unsigned int getNumUnknowns( ) override{ return residual.size( ); }

            auto set_solver( tardigradeHydra::SolverBase *_solver ){ solver = _solver; }

        protected:

            virtual void initializeUnknownVector( ){

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, flatJacobian );

            }


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

    tardigradeHydra::SolverBase solver;
    tardigradeHydra::SolverStepBase step;
    tardigradeHydra::PreconditionerBase preconditioner; //TODO: Replace with null preconditioner
    preconditioner._use_preconditioner = false;

    hydra.set_solver( &solver );

    solver.hydra = &hydra;
    solver.step  = &step;
    solver.preconditioner = &preconditioner;

    step.setSolver( &solver );
    preconditioner.setSolver( &solver );

    tardigradeHydra::floatVector answer = { 1./3, -2./3, 0 };

    tardigradeHydra::floatVector result( 3, 0 );

    tardigradeHydra::unit_test::hydraBaseTester::initializeUnknownVector( hydra );
    tardigradeHydra::unit_test::SolverStepBaseTester::solveNewtonUpdate( step, result );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

    result = { 0, 0, 0 };

    hydraBaseMock hydra_pre( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                             { }, { },
                             previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension,
                             9, 1e-9, 1e-9, 20, 5, 1e-4, true, 0 );

    tardigradeHydra::SolverBase solver_pre;
    tardigradeHydra::SolverStepBase step_pre;
    tardigradeHydra::PreconditionerBase preconditioner_pre;
    preconditioner_pre._use_preconditioner = true;

    hydra_pre.set_solver( &solver_pre );

    solver_pre.hydra = &hydra_pre;
    solver_pre.step  = &step_pre;
    solver_pre.preconditioner = &preconditioner_pre;

    step_pre.setSolver( &solver_pre );
    preconditioner_pre.setSolver( &solver_pre );

    tardigradeHydra::unit_test::hydraBaseTester::initializeUnknownVector( hydra_pre );
    tardigradeHydra::unit_test::SolverStepBaseTester::solveNewtonUpdate( step_pre, result );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_SolverStepBase_performPreconditionedSolve, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class hydraBaseMock : public tardigradeHydra::hydraBase {

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            tardigradeHydra::floatVector flatJacobian = { 1, 2, 3, 4, 5, 6, 7, 8, 2 };

            tardigradeHydra::floatVector residual = { 1, 2, 3 };

            virtual const unsigned int getNumUnknowns( ) override{ return residual.size( ); }

            void setSolver( tardigradeHydra::SolverBase *_solver ){ solver = _solver; }
            tardigradeHydra::SolverBase *getSolver( ){ return solver; }

        protected:

            virtual void initializeUnknownVector( ){

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, flatJacobian );

            }


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

    tardigradeHydra::SolverBase solver;

    tardigradeHydra::SolverStepBase step;

    tardigradeHydra::PreconditionerBase preconditioner;

    hydra.setSolver( &solver );

    solver.hydra = &hydra;
    solver.step  = &step;
    solver.preconditioner = &preconditioner;

    step.setSolver( &solver );
    preconditioner.setSolver( &solver );

    tardigradeHydra::floatVector answer = { 1./3, -2./3, 0 };

    tardigradeHydra::floatVector result( 3, 0 );

    tardigradeHydra::unit_test::hydraBaseTester::initializeUnknownVector( hydra );

    step.performPreconditionedSolve( result );

    BOOST_TEST( result == answer, CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_SolverStepBase_checkDescentDirection, * boost::unit_test::tolerance( 1e-5 ) ){
    /*!
     * Test checking the descent direction
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

    class StepDampingMock : public tardigradeHydra::StepDampingBase{

        public:

            using tardigradeHydra::StepDampingBase::StepDampingBase;

            tardigradeHydra::floatType residualNorm = 0.2408779076031648;

            tardigradeHydra::floatVector dResidualNormdX = { -1.17899799,  0.07843952, -0.01708813, -0.01779959, -0.06410942 };

            void mockInitialize( ){

                set_residualNorm( residualNorm );

                set_basedResidualNormdX( dResidualNormdX );

            }


    };

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase{

        public:

            using tardigradeHydra::SolverStepBase::SolverStepBase;

            virtual bool runCheckDescentDirection( const tardigradeHydra::floatVector &dx ){

                return checkDescentDirection( dx );

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            void set_solver( tardigradeHydra::SolverBase *_solver ){ solver = _solver; }

            virtual const unsigned int getNumUnknowns( ) override{ return 5; }

    };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    tardigradeHydra::SolverBase solver;
    SolverStepBaseMock step;
    StepDampingMock damping;
    hydra.set_solver( &solver );
    solver.hydra = &hydra;
    solver.step = &step;
    step.setSolver( &solver );
    step.damping = &damping;
    damping.step = &step;

    tardigradeHydra::floatVector dx = damping.dResidualNormdX;

    damping.mockInitialize( );

    BOOST_TEST( !step.runCheckDescentDirection( dx ) );

    dx = -damping.dResidualNormdX;

    BOOST_TEST( step.runCheckDescentDirection( dx ) );

}

BOOST_AUTO_TEST_CASE( test_SolverStepBase_checkGradientConvergence, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){
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

            tardigradeHydra::floatVector X0 = { -0.81579012, -0.13259765, -0.13827447, -0.0126298 , -0.14833942 };

            tardigradeHydra::floatVector jacobian = { -0.15378708,  0.9615284 ,  0.36965948, -0.0381362 , -0.21576496,
                                     -0.31364397,  0.45809941, -0.12285551, -0.88064421, -0.20391149,
                                      0.47599081, -0.63501654, -0.64909649,  0.06310275,  0.06365517,
                                      0.26880192,  0.69886359,  0.44891065,  0.22204702,  0.44488677,
                                     -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225 };

            using tardigradeHydra::hydraBase::hydraBase;

            virtual void formNonLinearResidual( ) override{

                tardigradeHydra::floatVector residual( 5, 0 );

                const tardigradeHydra::floatVector *X = getUnknownVector( );

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

            tardigradeHydra::SolverBase *getSolver( ){ return solver; }

    };

    class StepDampingBaseMock : public tardigradeHydra::StepDampingBase{

        public:

            using tardigradeHydra::StepDampingBase::StepDampingBase;

            tardigradeHydra::floatType baseResidualNorm = 0.4459139462561169;

            tardigradeHydra::floatVector basedResidualNormdX = { -0.86442794, -0.34410741, -0.58249594, -0.97271835, -0.32478706 };

            virtual void mockInitialize( ){

                set_baseResidualNorm( baseResidualNorm );

                set_basedResidualNormdX( basedResidualNormdX );

            }

    };

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase{

        public:

            using tardigradeHydra::SolverStepBase::SolverStepBase;

    };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    SolverStepBaseMock step;
    StepDampingBaseMock damping;
    damping.step = &step;
    step.damping = &damping;

    hydra.getSolver( )->step = &step;
    step.setSolver( hydra.getSolver( ) );

    tardigradeHydra::floatVector unknownVector = { 0.39293837, -0.42772133, -0.54629709,  0.10262954,  0.43893794 };

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra, unknownVector );

    damping.mockInitialize( );

    BOOST_TEST( !step.checkGradientConvergence( hydra.X0 ) );

    hydraBaseMock hydra2( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    SolverStepBaseMock step2;
    StepDampingBaseMock damping2;

    hydra2.getSolver( )->step = &step2;
    step2.setSolver( hydra2.getSolver( ) );
    damping2.step = &step2;
    step2.damping = &damping2;

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra2, 0.5 * hydra.X0 );

    damping2.mockInitialize( );

    BOOST_TEST( step2.checkGradientConvergence( hydra2.X0 ) );

}

BOOST_AUTO_TEST_CASE( test_SolverStepBase_performGradientStep, * boost::unit_test::tolerance( 1e-5 ) ){
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

            tardigradeHydra::floatVector A = { -0.15378708,  0.9615284 ,  0.36965948, -0.0381362 , -0.21576496,
                              -0.31364397,  0.45809941, -0.12285551, -0.88064421, -0.20391149,
                               0.47599081, -0.63501654, -0.64909649,  0.06310275,  0.06365517,
                               0.26880192,  0.69886359,  0.44891065,  0.22204702,  0.44488677,
                              -0.35408217, -0.27642269, -0.54347354, -0.41257191,  0.26195225 };

            using tardigradeHydra::hydraBase::hydraBase;

            virtual void formNonLinearResidual( ) override{

                tardigradeHydra::floatVector residual( 5, 0 );

                const tardigradeHydra::floatVector *X = getUnknownVector( );

                for ( unsigned int i = 0; i < 5; i++ ){
                    for ( unsigned int j = 0; j < 5; j++ ){
                        residual[ i ] += A[ 5 * i + j ] * ( ( *X )[ j ] * ( *X )[ j ] );
                    }
                }

                tardigradeHydra::unit_test::hydraBaseTester::set_residual( *this, residual );

            }

            virtual void formNonLinearDerivatives( ) override{

                tardigradeHydra::floatVector jacobian( 25, 0 );

                const tardigradeHydra::floatVector *X = getUnknownVector( );

                for ( unsigned int i = 0; i < 5; i++ ){
                    for ( unsigned int j = 0; j < 5; j++ ){
                        jacobian[ 5 * i + j ] += 2 * A[ 5 * i + j ] * ( *X )[ j ];
                    }
                }

                tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian( *this, jacobian );

            }

            virtual void decomposeUnknownVector( ) override{ return; }
            virtual const unsigned int getNumUnknowns( ) override{ return 5; }
            tardigradeHydra::SolverBase *getSolver( ){ return solver; }

    };

    class StepDampingBaseMock : public tardigradeHydra::StepDampingBase {

        public:

            using tardigradeHydra::StepDampingBase::StepDampingBase;

            tardigradeHydra::floatType baseResidualNorm = 0.2408779076031648;

            tardigradeHydra::floatVector basedResidualNormdX = { -1.17899799,  0.07843952, -0.01708813, -0.01779959, -0.06410942 };

            virtual void mockInitialize( ){

                set_baseResidualNorm( baseResidualNorm );

                set_basedResidualNormdX( basedResidualNormdX );

            }
    };

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase {

        public:

            using tardigradeHydra::SolverStepBase::SolverStepBase;

    };

    tardigradeHydra::floatVector X0 = { -0.81579012, -0.13259765, -0.13827447, -0.0126298 , -0.14833942 };

    tardigradeHydra::floatVector answer = { 0.36320787, -0.21103717, -0.12118634,  0.00516978, -0.08423 };

    hydraBaseMock hydra( time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient,
                         { }, { },
                         previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables, dimension );

    SolverStepBaseMock step;
    StepDampingBaseMock damping;

    step.damping = &damping;
    damping.step = &step;

    hydra.getSolver( )->step = &step;
    step.setSolver( hydra.getSolver( ) );

    tardigradeHydra::floatVector unknownVector = { 0.39293837, -0.42772133, -0.54629709,  0.10262954,  0.43893794 };

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra, unknownVector );

    damping.mockInitialize( );

    step.performGradientStep( X0 );

    BOOST_TEST( answer == *hydra.getUnknownVector( ), CHECK_PER_ELEMENT );

}

BOOST_AUTO_TEST_CASE( test_SolverStepBase_getNonlinearTerms, * boost::unit_test::tolerance( 1e-5 ) ){
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

                tardigradeHydra::unit_test::StepDampingBaseTester::setMuk( *(solver->step->damping), mu_k );

            }

            virtual void decomposeUnknownVector( ) override{ return; }
            virtual const unsigned int getNumUnknowns( ) override{ return 5; }

            tardigradeHydra::SolverBase *getSolver( ){ return solver; }

    };

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase{

        public:

            using tardigradeHydra::SolverStepBase::SolverStepBase;

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

    SolverStepBaseMock step;

    hydra.getSolver( )->step = &step;
    step.setSolver( hydra.getSolver( ) );

    tardigradeHydra::floatVector unknownVector = { 0.39293837, -0.42772133, -0.54629709,  0.10262954,  0.43893794 };

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector( hydra, unknownVector );

    hydra.getSolver( )->step->setUseLevenbergMarquardt( false );

    BOOST_TEST( hydra.residual == *hydra.getSolver( )->step->getNonlinearRHS( ), CHECK_PER_ELEMENT );

    BOOST_TEST( hydra.jacobian == *hydra.getSolver( )->step->getFlatNonlinearLHS( ), CHECK_PER_ELEMENT );

}
