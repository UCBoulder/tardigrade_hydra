
/**
 * \file test_tardigrade_SolverStepBase.cpp
 *
 * Tests for tardigrade_SolverStepBase
 */

#include "tardigrade_SolverBase.h"
#include "tardigrade_SolverStepBase.h"
#include "tardigrade_hydra.h"

#define BOOST_TEST_MODULE test_tardigrade_SolverStepBase
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

            static void set_flatJacobian(hydraBase &hydra, const tardigradeHydra::floatVector &value) {
                hydra._jacobian.second = value;
                hydra._jacobian.first  = true;

                hydra.addIterationData(&hydra._jacobian);
            }

            static void initializeUnknownVector(hydraBase &hydra) {
                BOOST_CHECK_NO_THROW(hydra.initializeUnknownVector());
            }
        };

        class StepDampingBaseTester {
           public:
            static void setMuk(StepDampingBase &damping, const tardigradeHydra::floatType &value) {
                damping.setMuk(value);
            }
        };

        class SolverStepBaseTester {
           public:
            static void checkUseLevenbergMarquardt(SolverStepBase &step) {
                BOOST_CHECK(step._use_LM_step == step.getUseLevenbergMarquardt());
            }
        };

        class TrialStepBaseTester {
           public:
            static void solveConstrainedQP(TrialStepBase &trial_step, tardigradeHydra::floatVector &dx,
                                           const unsigned int kmax = 100) {
                trial_step.solveConstrainedQP(dx, kmax);
            }

            static void solveNewtonUpdate(TrialStepBase &trial_step, tardigradeHydra::floatVector &deltaX) {
                trial_step.solveNewtonUpdate(deltaX);
            }
        };

    }  // namespace unit_test

}  // namespace tardigradeHydra

BOOST_AUTO_TEST_CASE(test_SolverStepBase_getUseLevenbergMarquardt,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::SolverStepBase step;

    tardigradeHydra::unit_test::SolverStepBaseTester::checkUseLevenbergMarquardt(step);
}

BOOST_AUTO_TEST_CASE(test_SolverStepBase_setUseLevenbergMarquardt,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::SolverStepBase step;

    step.setUseLevenbergMarquardt(true);

    BOOST_TEST(true == step.getUseLevenbergMarquardt());

    BOOST_TEST(true == step.damping->getUseGradientDescent());
}

BOOST_AUTO_TEST_CASE(test_SolverStepBase_solveNewtonUpdate, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        tardigradeHydra::floatVector flatJacobian = {1, 2, 3, 4, 5, 6, 7, 8, 2};

        tardigradeHydra::floatVector residual = {1, 2, 3};

        virtual const unsigned int getNumUnknowns() override { return residual.size(); }

        auto set_solver(tardigradeHydra::SolverBase *_solver) { solver = _solver; }

       protected:
        virtual void initializeUnknownVector() override {
            tardigradeHydra::unit_test::hydraBaseTester::set_residual(*this, residual);

            tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian(*this, flatJacobian);
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

    hydraBaseMock hydra(time, deltaTime, temperature, previousTemperature, deformationGradient,
                        previousDeformationGradient, {}, {}, previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables, dimension);

    tardigradeHydra::SolverBase         solver;
    tardigradeHydra::SolverStepBase     step;
    tardigradeHydra::PreconditionerBase preconditioner;  // TODO: Replace with null preconditioner
    preconditioner._use_preconditioner = false;

    hydra.set_solver(&solver);

    solver.hydra                    = &hydra;
    solver.step                     = &step;
    step.trial_step->preconditioner = &preconditioner;

    step.setSolver(&solver);
    preconditioner.trial_step = step.trial_step;

    tardigradeHydra::floatVector answer = {1. / 3, -2. / 3, 0};

    tardigradeHydra::floatVector result(3, 0);

    tardigradeHydra::unit_test::hydraBaseTester::initializeUnknownVector(hydra);
    tardigradeHydra::unit_test::TrialStepBaseTester::solveNewtonUpdate(*step.trial_step, result);

    BOOST_TEST(result == answer, CHECK_PER_ELEMENT);

    result = {0, 0, 0};

    hydraBaseMock hydra_pre(time, deltaTime, temperature, previousTemperature, deformationGradient,
                            previousDeformationGradient, {}, {}, previousStateVariables, parameters, numConfigurations,
                            numNonLinearSolveStateVariables, dimension, 9, 1e-9, 1e-9, 20, 5, 1e-4, true, 0);

    tardigradeHydra::SolverBase         solver_pre;
    tardigradeHydra::SolverStepBase     step_pre;
    tardigradeHydra::PreconditionerBase preconditioner_pre;
    preconditioner_pre._use_preconditioner = true;

    hydra_pre.set_solver(&solver_pre);

    solver_pre.hydra                    = &hydra_pre;
    solver_pre.step                     = &step_pre;
    step_pre.trial_step->preconditioner = &preconditioner_pre;

    step_pre.setSolver(&solver_pre);
    preconditioner_pre.trial_step = step_pre.trial_step;

    tardigradeHydra::unit_test::hydraBaseTester::initializeUnknownVector(hydra_pre);
    tardigradeHydra::unit_test::TrialStepBaseTester::solveNewtonUpdate(*step_pre.trial_step, result);

    BOOST_TEST(result == answer, CHECK_PER_ELEMENT);
}
