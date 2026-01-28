/**
 * \file test_tardigrade_IterativeSolverBase.cpp
 *
 * Tests for tardigrade_IterativeSolverBase
 */

#include "tardigrade_ArmijoGradientDamping.h"
#include "tardigrade_IterativeSolverBase.h"
#include "tardigrade_ResidualBase.h"
#include "tardigrade_NewtonStep.h"
#include "tardigrade_hydra.h"

#define BOOST_TEST_MODULE test_tardigrade_IterativeSolverBase
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element()

bool tolerantCheck(const std::vector<double> &v1, const std::vector<double> &v2, double eps = 1e-6, double tol = 1e-9) {
    if (v1.size() != v2.size()) {
        return false;
    }

    BOOST_CHECK(v1.size() == v2.size());

    const unsigned int len = v1.size();

    for (unsigned int i = 0; i < len; i++) {
        if (std::fabs(v1[i]) < tol) {
            if (std::fabs(v1[i] - v2[i]) > eps) {
                return false;
            }

            BOOST_CHECK(std::fabs(v1[i] - v2[i]) <= eps);

        } else {
            if ((std::fabs(v1[i] - v2[i]) / std::fabs(v1[i]) > eps) ||
                (std::fabs(v1[i] - v2[i]) / std::fabs(v2[i]) > eps)) {
                return false;
            }

            BOOST_TEST(v1[i] == v2[i]);
        }
    }

    return true;
}

namespace tardigradeHydra {

    namespace unit_test {

        class hydraBaseTester {
           public:
            static void set_residual(hydraBase &hydra, const tardigradeHydra::floatVector &value) {
                hydra._residual.second = value;
                hydra._residual.first  = true;

                hydra.addIterationData(&hydra._residual);
            }

            static void set_unknownVector(hydraBase &hydra, const tardigradeHydra::floatVector &value) {
                hydra._X.second = value;
                hydra._X.first  = true;
            }

            static void set_flatJacobian(hydraBase &hydra, const floatVector &value) {
                hydra._jacobian.second = value;
                hydra._jacobian.first  = true;

                hydra.addIterationData(&hydra._jacobian);
            }

            static void resetIterationData(hydraBase &hydra) { hydra.resetIterationData(); }
        };

    }  // namespace unit_test

}  // namespace tardigradeHydra

unsigned int test_SolverBase_solve_in_gradient_convergence = 0;

BOOST_AUTO_TEST_CASE(test_IterativeSolverBase_solve, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        tardigradeHydra::floatVector initialUnknownVector = {1, 2, 3};

        std::vector<bool> isGradient = {0, 1, 0, 0};

        std::vector<unsigned int> numLSIterations = {1, 2, 1, 3};

        std::vector<std::vector<tardigradeHydra::floatVector> > residual = {
            {{1, 2, 3}},
            {{4, 5, 6}, {7, 8, 9}},
            {{10, 9, 8}},
            {{7, 6, 5}, {4, 3, 2}, {1, 1, 1}},
        };

        std::vector<tardigradeHydra::floatType> expectedBaseResidualNorms = {14, 77, 245, 110};

        std::vector<tardigradeHydra::floatVector> expectedBasedResidualNormdXs = {
            {2.97757248,   -1.95288068, -3.03183756 },
            {-16.10008822, -6.7048463,  14.58211116 },
            {-33.99029958, 28.301765,   -20.74369594},
            {-4.96839724,  18.87890456, 2.15863168  }
        };

        std::vector<tardigradeHydra::floatType> expectedMuk = {7e-8, 7e-8, 7e-8, 7e-8};

        std::vector<std::vector<tardigradeHydra::floatVector> > flatJacobian = {
            {{0.99951474, -0.18854971, -0.59377647, 0.73798389, 0.50649978, -0.24789123, -0.32889876, -0.60029673,
              -0.14211995}},
            {{0.39531964, -0.92103617, 0.15529267, -0.88403283, -0.29325643, 0.33220876, -0.86852642, 0.29966728,
              0.83480685},
             {-0.52407245, 0.35366894, 0.08551378, 0.70688925, -0.14848006, 0.24923792, 0.67051005, -0.08228108,
              0.65747464}},
            {{-0.81643355, 0.85279545, -0.18966435, -0.42616781, 0.53676336, -0.92714295, -0.624413, 0.09900722,
              -0.01636474}},
            {{-0.2448239, 0.74837591, 0.68140055, -0.14194502, 0.67937811, -0.28307126, 0.01624776, 0.02491045,
              -0.39841209},
             {-0.65594132, -0.44378971, -0.43630527, 0.13856471, 0.95333808, -0.89484459, -0.26311552, -0.94030312,
              0.69678361},
             {-0.82400549, -0.05626004, -0.17170742, 0.34131108, 0.49817018, -0.78483747, 0.95727384, 0.23324588,
              -0.22667228}},
        };

        std::vector<std::vector<tardigradeHydra::floatVector> > expectedXVectors = {
            {{17.15005309, -9.55454694, 35.53888217}},
            {{18.95259726, -5.57445709, 28.79822658},
             {33.25014131, -2.84970064, 20.95677101},
             {31.64013249, -3.52018527, 22.41498213}},
            {{44.6038325, -1.77731706, 27.17238891}},
            {{420.89746176, 81.81218142, 60.29432705},
             {232.75064713, 40.01743218, 43.73335798},
             {138.67723982, 19.12005756, 35.45287345}},
        };

        virtual const unsigned int getNumUnknowns() override { return initialUnknownVector.size(); }

        unsigned int num_derivative_calls = 0;

        unsigned int num_residual_calls = 0;

        void                         setSolver(tardigradeHydra::SolverBase *_solver) { solver = _solver; }
        tardigradeHydra::SolverBase *getSolver() { return solver; }

       private:
        using tardigradeHydra::hydraBase::getResidual;

        virtual void initializeUnknownVector() override {
            tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector(*this, initialUnknownVector);

            tardigradeHydra::unit_test::hydraBaseTester::set_residual(*this, residual[0][0]);

            tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian(*this, flatJacobian[0][0]);
        }

        virtual void formNonLinearResidual() override {
            auto local_solver = dynamic_cast<tardigradeHydra::IterativeSolverBase*>(solver);
            unsigned int iteration = local_solver->getIteration();

            auto local_damping = dynamic_cast<tardigradeHydra::ArmijoGradientDamping *>(solver->step->damping);

            unsigned int LSIteration = local_damping->getLSIteration();

            unsigned int gradIteration = local_damping->getGradientIteration();

            unsigned int iterationOffset = 0;

            unsigned int LSoffset = 0;

            if (!isGradient[iteration] && (LSIteration < residual[iteration].size() - 1)) {
                LSoffset += 1;

            } else if (isGradient[iteration] && (gradIteration < residual[iteration].size() - 2)) {
                LSoffset += 1;

            } else if (iteration < residual.size() - 1) {
                iterationOffset += 1;
                LSIteration   = 0;
                gradIteration = 0;
            }

            tardigradeHydra::unit_test::hydraBaseTester::set_residual(
                *this, residual[iteration + iterationOffset][LSIteration + LSoffset + gradIteration]);

            num_residual_calls++;
        }

        virtual void formNonLinearDerivatives() override {
            auto local_solver = dynamic_cast<tardigradeHydra::IterativeSolverBase*>(solver);
            unsigned int iteration = local_solver->getIteration();

            tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian(*this, flatJacobian[iteration][0]);

            num_derivative_calls++;
        }

        virtual void updateUnknownVector(const tardigradeHydra::floatVector &newUnknownVector) override {
            tardigradeHydra::unit_test::hydraBaseTester::resetIterationData(*this);

            auto local_solver = dynamic_cast<tardigradeHydra::IterativeSolverBase*>(solver);
            unsigned int iteration = local_solver->getIteration();

            auto local_damping = dynamic_cast<tardigradeHydra::ArmijoGradientDamping *>(solver->step->damping);

            unsigned int LSIteration = local_damping->getLSIteration();

            unsigned int gradIteration = local_damping->getGradientIteration();

            unsigned int subIteration = LSIteration + gradIteration + test_SolverBase_solve_in_gradient_convergence;

            BOOST_TEST(expectedXVectors[iteration][subIteration] == newUnknownVector, CHECK_PER_ELEMENT);

            BOOST_TEST(expectedBaseResidualNorms[iteration] == *solver->step->damping->get_baseResidualNorm());

            BOOST_TEST(expectedBasedResidualNormdXs[iteration] == *solver->step->damping->get_basedResidualNormdX(),
                       CHECK_PER_ELEMENT);

            BOOST_TEST(expectedMuk[iteration] == solver->step->damping->getMuk());

            tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector(*this, newUnknownVector);

            test_SolverBase_solve_in_gradient_convergence = 0;
        }
    };

    class IterativeSolverBaseMock : public tardigradeHydra::IterativeSolverBase {
       public:
        std::vector<std::vector<tardigradeHydra::floatVector> > residual = {
            {{1, 2, 3}},
            {{4, 5, 6}, {7, 8, 9}},
            {{10, 9, 8}},
            {{7, 6, 5}, {4, 3, 2}, {1, 1, 1}},
        };

        unsigned int num_pre_iterativesolve_calls = 0;

        unsigned int num_successful_iterativestep_calls = 0;

        unsigned int num_post_iterativesolve_calls = 0;

        using tardigradeHydra::IterativeSolverBase::IterativeSolverBase;

        virtual void callResidualPreIterativeSolve() override { num_pre_iterativesolve_calls++; }

        virtual void callResidualSuccessfulIterativeStep() override {
            tardigradeHydra::IterativeSolverBase::callResidualSuccessfulIterativeStep();

            num_successful_iterativestep_calls++;
        }

        virtual void callResidualPostIterativeSolve() override { num_post_iterativesolve_calls++; }

        virtual bool checkConvergence() override {
            getResidual();

            unsigned int iteration = getIteration();

            if (iteration < residual.size()) {
                return false;
            }

            return true;
        }
    };

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase {
       public:
        using tardigradeHydra::SolverStepBase::SolverStepBase;

        std::vector<std::vector<tardigradeHydra::floatVector> > residual = {
            {{1, 2, 3}},
            {{4, 5, 6}, {7, 8, 9}},
            {{10, 9, 8}},
            {{7, 6, 5}, {4, 3, 2}, {1, 1, 1}},
        };

        tardigradeHydra::SolverBase* getSolver(){ return solver; }

        unsigned int getIteration(){
            auto local_solver = dynamic_cast<IterativeSolverBaseMock*>(solver);
            return local_solver->getIteration();
        }
    };

    class NewtonStepMock : public tardigradeHydra::NewtonStep {
       public:
        using tardigradeHydra::NewtonStep::NewtonStep;

        unsigned int getIteration(){
            auto local_step   = dynamic_cast<SolverStepBaseMock*>(step);
            return local_step->getIteration();
        }

    };

    class ArmijoGradientDampingMock : public tardigradeHydra::ArmijoGradientDamping {
       public:
        using tardigradeHydra::ArmijoGradientDamping::ArmijoGradientDamping;

        std::vector<std::vector<tardigradeHydra::floatVector> > residual = {
            {{1, 2, 3}},
            {{4, 5, 6}, {7, 8, 9}},
            {{10, 9, 8}},
            {{7, 6, 5}, {4, 3, 2}, {1, 1, 1}},
        };

        virtual bool checkLSConvergence() override {
            auto local_step = dynamic_cast<SolverStepBaseMock*>(step);
            unsigned int iteration = local_step->getIteration();

            unsigned int LSIteration = getLSIteration();

            getResidual();

            if (LSIteration < residual[iteration].size() - 1) {
                return false;
            }

            return true;
        }

        virtual bool checkDescentDirection(const tardigradeHydra::floatVector &dx) override {
            auto local_step = dynamic_cast<SolverStepBaseMock*>(step);
            unsigned int iteration = local_step->getIteration();

            if (iteration == 1) {
                return false;
            }

            return true;
        }

        virtual bool checkGradientConvergence(const tardigradeHydra::floatVector &X0) override {
            auto local_step = dynamic_cast<SolverStepBaseMock*>(step);
            unsigned int iteration = local_step->getIteration();

            unsigned int gradientIteration = getGradientIteration();

            test_SolverBase_solve_in_gradient_convergence = 1;

            getResidual();

            if (gradientIteration < residual[iteration].size() - 1) {
                return false;
            }

            test_SolverBase_solve_in_gradient_convergence = 0;

            return true;
        }

        virtual void performGradientStep(const tardigradeHydra::floatVector &X0) override {
            test_SolverBase_solve_in_gradient_convergence = 1;

            tardigradeHydra::StepDampingBase::performGradientStep(X0);
        }
    };

    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    test_SolverBase_solve_in_gradient_convergence = 0;
    hydraBaseMock hydra(time, deltaTime, temperature, previousTemperature, deformationGradient,
                        previousDeformationGradient, {}, {}, previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables, dimension);

    IterativeSolverBaseMock   solver;
    SolverStepBaseMock        step;
    ArmijoGradientDampingMock damping;
    NewtonStepMock     trial_step;

    step.trial_step = &trial_step;
    step.damping    = &damping;

    trial_step.step = &step;
    damping.step    = &step;

    damping.setMaxLSIterations(5);
    damping.setLSAlpha(1e-4);
    tardigradeHydra::PreconditionerBase preconditioner;

    hydra.setSolver(&solver);

    solver.hydra              = &hydra;
    solver.step               = &step;
    trial_step.preconditioner = &preconditioner;

    step.setSolver(&solver);
    preconditioner.trial_step = &trial_step;

    hydra.getSolver()->step->damping->setUseGradientDescent(true);

    solver.solve();

    BOOST_TEST(step.getNumUndamped() == 2);

    BOOST_TEST(solver.num_pre_iterativesolve_calls == 1);

    BOOST_TEST(solver.num_post_iterativesolve_calls == 1);

    BOOST_TEST(solver.num_successful_iterativestep_calls == 4);

    BOOST_TEST(damping.getNumLS() == 1);

    BOOST_TEST(damping.getNumGrad() == 1);

    BOOST_TEST(hydra.num_residual_calls == 9);  // 9 because we initialize the residual

    BOOST_TEST(hydra.num_derivative_calls == 3);  // 3 because we initialize the jacobian

    test_SolverBase_solve_in_gradient_convergence = 0;
    hydraBaseMock hydra_pre(time, deltaTime, temperature, previousTemperature, deformationGradient,
                            previousDeformationGradient, {}, {}, previousStateVariables, parameters, numConfigurations,
                            numNonLinearSolveStateVariables, dimension, 9, 1e-9, 1e-9);

    IterativeSolverBaseMock   solver_pre;
    SolverStepBaseMock        step_pre;
    ArmijoGradientDampingMock damping_pre;
    NewtonStepMock     trial_step_pre;

    step_pre.trial_step = &trial_step_pre;
    step_pre.damping    = &damping_pre;

    trial_step_pre.step = &step_pre;
    damping_pre.step    = &step_pre;

    damping_pre.setMaxLSIterations(5);
    damping_pre.setLSAlpha(1e-4);
    tardigradeHydra::MaxRowPreconditioner preconditioner_pre;

    hydra_pre.setSolver(&solver_pre);

    solver_pre.hydra              = &hydra_pre;
    solver_pre.step               = &step_pre;
    trial_step_pre.preconditioner = &preconditioner_pre;

    step_pre.setSolver(&solver_pre);
    preconditioner_pre.trial_step = &trial_step_pre;

    hydra_pre.getSolver()->step->damping->setUseGradientDescent(true);

    solver_pre.solve();

    BOOST_TEST(step_pre.getNumUndamped() == 2);

    BOOST_TEST(solver_pre.num_pre_iterativesolve_calls == 1);

    BOOST_TEST(solver_pre.num_post_iterativesolve_calls == 1);

    BOOST_TEST(solver_pre.num_successful_iterativestep_calls == 4);

    BOOST_TEST(damping.getNumLS() == 1);

    BOOST_TEST(damping_pre.getNumGrad() == 1);

    BOOST_TEST(hydra.num_residual_calls == 9);  // 9 because we initialize the residual

    BOOST_TEST(hydra.num_derivative_calls == 3);  // 3 because we initialize the jacobian
}

BOOST_AUTO_TEST_CASE(test_IterativeSolverBase_callResidualSuccessfulIterativeStep,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> {
       public:
        using tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>::ResidualBase;

        unsigned int numSuccessfulIterativeStepCalls = 0;

        virtual void successfulIterativeStep() override {
            BOOST_TEST(hydra->getMutableResidual());

            numSuccessfulIterativeStepCalls++;
        }

       protected:
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        residualMock r1;

        residualMock r2;

        residualMock r3;

        unsigned int s1 = 36;

        unsigned int s2 = 2;

        unsigned int s3 = 3;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            r1 = residualMock(this, s1);

            r2 = residualMock(this, s2);

            r3 = residualMock(this, s3);

            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            residuals[0] = &r1;

            residuals[1] = &r2;

            residuals[2] = &r3;

            setResidualClasses(residuals);
        }

        void setSolver(tardigradeHydra::SolverBase *_solver) { solver = solver; }
    };

    class IterativeSolverBaseMock : public tardigradeHydra::IterativeSolverBase {
       public:
        using tardigradeHydra::IterativeSolverBase::IterativeSolverBase;

        virtual void public_callResidualSuccessfulIterativeStep() { callResidualSuccessfulIterativeStep(); }
    };

    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079, 0.0,        0.1};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra(time, deltaTime, temperature, previousTemperature, deformationGradient,
                        previousDeformationGradient, {}, {}, previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables, dimension);

    IterativeSolverBaseMock solver;

    hydra.setSolver(&solver);
    solver.hydra = &hydra;

    solver.public_callResidualSuccessfulIterativeStep();

    BOOST_TEST(hydra.r1.numSuccessfulIterativeStepCalls == 1);

    BOOST_TEST(hydra.r2.numSuccessfulIterativeStepCalls == 1);

    BOOST_TEST(hydra.r3.numSuccessfulIterativeStepCalls == 1);
}

BOOST_AUTO_TEST_CASE(test_IterativeSolverBase_callResidualPreIterativeSolve,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> {
       public:
        using tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>::ResidualBase;

        unsigned int numPreIterativeSolveCalls = 0;

        virtual void preIterativeSolve() override { numPreIterativeSolveCalls++; }

       protected:
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        residualMock r1;

        residualMock r2;

        residualMock r3;

        unsigned int s1 = 36;

        unsigned int s2 = 2;

        unsigned int s3 = 3;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            r1 = residualMock(this, s1);

            r2 = residualMock(this, s2);

            r3 = residualMock(this, s3);

            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            residuals[0] = &r1;

            residuals[1] = &r2;

            residuals[2] = &r3;

            setResidualClasses(residuals);
        }

        void                         setSolver(tardigradeHydra::SolverBase *_solver) { solver = _solver; }
        tardigradeHydra::SolverBase *getSolver() { return solver; }
    };

    class IterativeSolverBaseMock : public tardigradeHydra::IterativeSolverBase {
       public:
        using tardigradeHydra::IterativeSolverBase::IterativeSolverBase;

        void public_callResidualPreIterativeSolve() { callResidualPreIterativeSolve(); }
    };

    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079, 0.0,        0.1};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra(time, deltaTime, temperature, previousTemperature, deformationGradient,
                        previousDeformationGradient, {}, {}, previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables, dimension);

    IterativeSolverBaseMock solver;
    solver.hydra = &hydra;
    hydra.setSolver(&solver);

    solver.public_callResidualPreIterativeSolve();

    BOOST_TEST(hydra.r1.numPreIterativeSolveCalls == 1);

    BOOST_TEST(hydra.r2.numPreIterativeSolveCalls == 1);

    BOOST_TEST(hydra.r3.numPreIterativeSolveCalls == 1);
}

BOOST_AUTO_TEST_CASE(test_IterativeSolverBase_callResidualPostIterativeSolve,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> {
       public:
        using tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>::ResidualBase;

        unsigned int numPostIterativeSolveCalls = 0;

        virtual void postIterativeSolve() override { numPostIterativeSolveCalls++; }

       protected:
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        residualMock r1;

        residualMock r2;

        residualMock r3;

        unsigned int s1 = 36;

        unsigned int s2 = 2;

        unsigned int s3 = 3;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            r1 = residualMock(this, s1);

            r2 = residualMock(this, s2);

            r3 = residualMock(this, s3);

            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            residuals[0] = &r1;

            residuals[1] = &r2;

            residuals[2] = &r3;

            setResidualClasses(residuals);
        }

        void setSolver(tardigradeHydra::SolverBase *_solver) { solver = _solver; }
    };

    class IterativeSolverBaseMock : public tardigradeHydra::IterativeSolverBase {
       public:
        using tardigradeHydra::IterativeSolverBase::IterativeSolverBase;

        virtual void public_callResidualPostIterativeSolve() { callResidualPostIterativeSolve(); }
    };

    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079, 0.0,        0.1};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra(time, deltaTime, temperature, previousTemperature, deformationGradient,
                        previousDeformationGradient, {}, {}, previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables, dimension);

    IterativeSolverBaseMock solver;
    hydra.setSolver(&solver);
    solver.hydra = &hydra;

    solver.public_callResidualPostIterativeSolve();

    BOOST_TEST(hydra.r1.numPostIterativeSolveCalls == 1);

    BOOST_TEST(hydra.r2.numPostIterativeSolveCalls == 1);

    BOOST_TEST(hydra.r3.numPostIterativeSolveCalls == 1);
}

BOOST_AUTO_TEST_CASE(test_IterativeSolverBase_checkConvergence, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        void setSolver(tardigradeHydra::SolverBase *_solver) { solver = solver; }
    };

    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra(time, deltaTime, temperature, previousTemperature, deformationGradient,
                        previousDeformationGradient, {}, {}, previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables, dimension);

    tardigradeHydra::IterativeSolverBase solver;

    hydra.setSolver(&solver);
    solver.hydra = &hydra;

    tardigradeHydra::floatVector residual = {1, 2, -3, 0};

    tardigradeHydra::floatVector unknownVector = {-2, 5, 10, 0.3};

    tardigradeHydra::unit_test::hydraBaseTester::set_residual(hydra, residual);

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector(hydra, unknownVector);

    BOOST_CHECK(!solver.checkConvergence());

    residual = {0, 2, -3, 0};

    tardigradeHydra::unit_test::hydraBaseTester::set_residual(hydra, residual);

    BOOST_CHECK(!solver.checkConvergence());

    residual = {0, 0, 0, 0};

    tardigradeHydra::unit_test::hydraBaseTester::set_residual(hydra, residual);

    BOOST_CHECK(solver.checkConvergence());
}

BOOST_AUTO_TEST_CASE(test_IterativeSolverBase_setTolerance, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        void setSolver(tardigradeHydra::SolverBase *_solver) { solver = _solver; }
    };

    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra(time, deltaTime, temperature, previousTemperature, deformationGradient,
                        previousDeformationGradient, {}, {}, previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables, dimension);

    tardigradeHydra::IterativeSolverBase solver;

    hydra.setSolver(&solver);
    solver.hydra = &hydra;

    tardigradeHydra::floatVector residual = {1, 2, -3, 0};

    tardigradeHydra::floatVector unknownVector = {-2, 5, 10, 0.3};

    tardigradeHydra::floatVector toleranceAnswer =
        1e-9 * (tardigradeVectorTools::abs(residual) + tardigradeVectorTools::abs(unknownVector)) + 1e-9;

    tardigradeHydra::unit_test::hydraBaseTester::set_residual(hydra, residual);

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector(hydra, unknownVector);

    BOOST_TEST(*solver.getTolerance() == toleranceAnswer, CHECK_PER_ELEMENT);
}
