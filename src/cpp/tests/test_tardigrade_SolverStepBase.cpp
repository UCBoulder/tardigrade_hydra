
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

        };

        class SolverStepBaseTester{

            public:

                static void checkMuk( SolverStepBase &step ){

                    BOOST_CHECK( step._mu_k == step.getMuk( ) );

                }

                static void checkLMMu( SolverStepBase &step ){

                    BOOST_CHECK( step._lm_mu == step.getLMMu( ) );

                }

                static void setMuk( SolverStepBase &step, const tardigradeHydra::floatType &value ){

                    step.setMuk( value );

                }

                static void solveConstrainedQP( SolverStepBase &step, tardigradeHydra::floatVector &dx, const unsigned int kmax=100 ){

                    step.solveConstrainedQP( dx, kmax );

                }

        };

    }


}

BOOST_AUTO_TEST_CASE( test_hydraBase_getMuk, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::SolverStepBase step;

    tardigradeHydra::unit_test::SolverStepBaseTester::checkMuk( step );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_getLMMu, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    tardigradeHydra::SolverStepBase step;

    tardigradeHydra::unit_test::SolverStepBaseTester::checkLMMu( step );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setMuk, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase{

        using tardigradeHydra::SolverStepBase::SolverStepBase;

        public:

            void public_setMuk( const tardigradeHydra::floatType &v ){

                setMuk( v );

            }

    };

    SolverStepBaseMock step;

    step.public_setMuk( 123.4 );

    BOOST_TEST( 123.4 == step.getMuk( ) );

}

BOOST_AUTO_TEST_CASE( test_hydraBase_setLMMu, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase{

        using tardigradeHydra::SolverStepBase::SolverStepBase;

        public:

            void public_setLMMu( const tardigradeHydra::floatType &v ){

                setLMMu( v );

            }

    };

    SolverStepBaseMock step;

    step.public_setLMMu( 123.4 );

    BOOST_TEST( 123.4 == step.getLMMu( ) );

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

    hydra.set_solver( &solver );
    solver.hydra = &hydra;
    solver.step = &step;
    step.setSolver( &solver );

    BOOST_TEST( answer == *step.get_residualNorm( ) );

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

            hydrap.set_solver( &solverp );
            solverp.hydra = &hydrap;
            solverp.step  = &stepp;
            stepp.setSolver( &solverp );

            hydram.set_solver( &solverm );
            solverm.hydra = &hydram;
            solverm.step  = &stepm;
            stepm.setSolver( &solverm );

            dResidualNormdX[ i ] = ( *stepp.get_residualNorm( ) - *stepm.get_residualNorm( ) ) / ( 2 * delta );

        }

        BOOST_TEST( dResidualNormdX == *step.get_dResidualNormdX( ), CHECK_PER_ELEMENT );

    }

}

BOOST_AUTO_TEST_CASE( test_setBaseQuantities, * boost::unit_test::tolerance( 1e-5 ) ){
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

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase{

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

            const tardigradeHydra::floatType *getBaseResidualNorm( ){ return get_baseResidualNorm( ); }

            const tardigradeHydra::floatVector *getBasedResidualNormdX( ){ return get_basedResidualNormdX( ); }

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

    SolverStepBaseMock step;
    tardigradeHydra::SolverBase solver;
    hydra.set_solver( &solver );
    solver.hydra = &hydra;
    solver.step = &step;
    step.setSolver( &solver );

    step.runSetBaseQuantities( );

    BOOST_TEST( answer1 == step.getMuk( ) );

    BOOST_TEST( step.rnorm == *hydra.getBaseResidualNorm( ) );

    BOOST_TEST( step.dRNormdX == *hydra.getBasedResidualNormdX( ) );

    step.rnorm = 1e-9;

    step.setResidualNorm( );

    step.runSetBaseQuantities( );

    BOOST_TEST( step.rnorm == step.getMuk( ) );

}

BOOST_AUTO_TEST_CASE( test_SolverStepBase_solveConstrainedQP, * boost::unit_test::tolerance( DEFAULT_TEST_TOLERANCE ) ){

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase{

        public:

            using tardigradeHydra::SolverStepBase::SolverStepBase;

            tardigradeHydra::floatVector initialUnknownVector = { 2, 0 };

        protected:

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

            virtual void updateKKTMatrix( tardigradeHydra::floatVector &K, const std::vector< bool > &active_constraints ) override{

                const unsigned int numConstraints = getNumConstraints( );

                for ( unsigned int i = 0; i < numConstraints; i++ ){

                    if ( active_constraints[ i ] ){

                        K[ 7 * 0 + i + 2 ] = ( *getConstraintJacobians( ) )[ 2 * i + 0 ];
                        K[ 7 * 1 + i + 2 ] = ( *getConstraintJacobians( ) )[ 2 * i + 1 ];

                        K[ 7 * ( i + 2 ) + 0 ] = ( *getConstraintJacobians( ) )[ 2 * i + 0 ];
                        K[ 7 * ( i + 2 ) + 1 ] = ( *getConstraintJacobians( ) )[ 2 * i + 1 ];

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

            virtual void initializeActiveConstraints( std::vector< bool > &active_constraints ) override{

                active_constraints = { false, false, true, false, true };

            }

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

            auto public_setMuk( const tardigradeHydra::floatType &value ){ tardigradeHydra::unit_test::SolverStepBaseTester::setMuk( *(solver->step), value ); }

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

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase{

        public:

            using tardigradeHydra::hydraBase::hydraBase;

            using tardigradeHydra::hydraBase::setResidualClasses;

            tardigradeHydra::floatVector initialUnknownVector = { 2, 1 };

            virtual void public_updateKKTMatrix( tardigradeHydra::floatVector &K, const std::vector< bool > &active_constraints ){

                updateKKTMatrix( K, active_constraints );

            }

            auto public_setMuk( const tardigradeHydra::floatType &value ){ tardigradeHydra::unit_test::SolverStepBaseTester::setMuk( *(solver->step), value ); }

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

    hydra.public_updateKKTMatrix( result_KKT, active_constraints );

    BOOST_TEST( answer2_KKTMatrix == result_KKT, CHECK_PER_ELEMENT );

    result_KKT.clear( );

    step.public_assembleKKTMatrix( result_KKT, active_constraints );

    BOOST_TEST( answer2_KKTMatrix == result_KKT, CHECK_PER_ELEMENT );

}

