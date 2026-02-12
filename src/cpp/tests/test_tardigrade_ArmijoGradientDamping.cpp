/**
 * \file test_tardigrade_ArmijoGradientDamping.cpp
 *
 * Tests for tardigrade_ArmijoGradientDamping
 */

#include "tardigrade_ArmijoGradientDamping.h"
#include "tardigrade_SolverBase.h"
#include "tardigrade_SolverStepBase.h"
#include "tardigrade_hydra.h"

#define BOOST_TEST_MODULE test_tardigrade_ArmijoGradientDamping
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

        class ArmijoGradientDampingTester {
           public:
            static void checkLSAlpha(ArmijoGradientDamping &damping) {
                BOOST_CHECK(damping._lsAlpha == damping.getLSAlpha());
            }
        };

        class StepDampingBaseTester {
           public:
            static void checkGradientRho(StepDampingBase &damping) {
                BOOST_CHECK(damping._gradientRho == damping.getGradientRho());
            }

            static void checkGradientP(StepDampingBase &damping) {
                BOOST_CHECK(damping._gradientP == damping.getGradientP());
            }

            static void checkGradientBeta(StepDampingBase &damping) {
                BOOST_CHECK(damping._gradientBeta == damping.getGradientBeta());
            }

            static void checkGradientSigma(StepDampingBase &damping) {
                BOOST_CHECK(damping._gradientSigma == damping.getGradientSigma());
            }

            static void checkMaxGradientIterations(StepDampingBase &damping) {
                BOOST_CHECK(damping._maxGradientIterations == damping.getMaxGradientIterations());
            }

            static void checkMuk(StepDampingBase &damping) { BOOST_CHECK(damping._mu_k == damping.getMuk()); }

            static void checkLMMu(StepDampingBase &damping) { BOOST_CHECK(damping._lm_mu == damping.getLMMu()); }

            static void checkUseGradientDescent(StepDampingBase &damping) {
                BOOST_CHECK(damping._use_gradient_descent == damping.getUseGradientDescent());
            }

            //                static void setMuk( SolverStepBase &step, const tardigradeHydra::floatType &value ){
            //
            //                    step.setMuk( value );
            //
            //                }
        };

    }  // namespace unit_test

}  // namespace tardigradeHydra

BOOST_AUTO_TEST_CASE(test_ArmijoGradientDamping_getLSAlpha, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::ArmijoGradientDamping damping;
    damping.setLSAlpha(1.23);

    tardigradeHydra::unit_test::ArmijoGradientDampingTester::checkLSAlpha(damping);
}

BOOST_AUTO_TEST_CASE(test_ArmijoGradientDamping_getLSResidualNorm,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;
    };

    hydraBaseMock                   hydra;
    tardigradeHydra::SolverBase     solver;
    tardigradeHydra::SolverStepBase step;
    //    tardigradeHydra::PreconditionerBase    preconditioner;
    tardigradeHydra::ArmijoGradientDamping damping;

    hydra.solver = &solver;

    solver.hydra = &hydra;
    solver.step  = &step;
    //    step.trial_step->preconditioner = &preconditioner;

    step.solver = &solver;
    //    preconditioner.trial_step = step.trial_step;

    step.damping = &damping;
    damping.step = &step;

    tardigradeHydra::floatVector residual = {1, 2, 3};

    tardigradeHydra::floatType lsResidualNormAnswer = tardigradeVectorTools::l2norm(residual);
    tardigradeHydra::unit_test::hydraBaseTester::set_residual(hydra, residual);

    BOOST_TEST(*damping.getLSResidualNorm() == lsResidualNormAnswer);
}

BOOST_AUTO_TEST_CASE(test_ArmijoGradientDamping_checkLSConvergence,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;
    };

    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature);

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

    hydraBaseMock hydra(dof, deformationGradient, previousDeformationGradient, {}, {},
                        previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::SolverStepBase        step;
    tardigradeHydra::ArmijoGradientDamping damping;

    step.damping = &damping;
    damping.step = &step;

    damping.setLSAlpha(1e-4);
    damping.setMaxLSIterations(5);

    hydra.solver->step = &step;
    step.solver        = hydra.solver;

    tardigradeHydra::floatVector residual = {1, 2, -3, 0};

    tardigradeHydra::unit_test::hydraBaseTester::set_residual(hydra, residual);

    BOOST_CHECK(!damping.checkLSConvergence());

    residual = {0.1, 2, -3, 0};

    tardigradeHydra::unit_test::hydraBaseTester::set_residual(hydra, residual);

    BOOST_CHECK(damping.checkLSConvergence());
}
