/**
 * \file test_tardigrade_NonlinearStepBase.cpp
 *
 * Tests for tardigrade_NonlinearStepBase
 */

#include "tardigrade_NonlinearStepBase.h"
#include "tardigrade_hydra.h"
#include "tardigrade_MaxRowPreconditioner.h"

#define BOOST_TEST_MODULE test_tardigrade_NonlinearStepBase
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

    }  // namespace unit_test

}  // namespace tardigradeHydra

BOOST_AUTO_TEST_CASE(test_NonlinearStepBase_getNonlinearTerms, *boost::unit_test::tolerance(1e-5)) {
    /*!
     * Test checking the gradient convergence
     */

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

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        tardigradeHydra::floatVector residual = {1, 2, 3, 4, 5};

        tardigradeHydra::floatVector jacobian = {-0.15378708, 0.9615284,   0.36965948,  -0.0381362,  -0.21576496,
                                                 -0.31364397, 0.45809941,  -0.12285551, -0.88064421, -0.20391149,
                                                 0.47599081,  -0.63501654, -0.64909649, 0.06310275,  0.06365517,
                                                 0.26880192,  0.69886359,  0.44891065,  0.22204702,  0.44488677,
                                                 -0.35408217, -0.27642269, -0.54347354, -0.41257191, 0.26195225};

        tardigradeHydra::floatType mu_k = 1.34;

        using tardigradeHydra::hydraBase::hydraBase;

        virtual void formNonLinearResidual() override {
            tardigradeHydra::unit_test::hydraBaseTester::set_residual(*this, residual);
        }

        virtual void formNonLinearDerivatives() override {
            tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian(*this, jacobian);

            tardigradeHydra::unit_test::StepDampingBaseTester::setMuk(*(solver->step->damping), mu_k);
        }

        virtual void               decomposeUnknownVector() override { return; }
        virtual const unsigned int getNumUnknowns() override { return 5; }

        tardigradeHydra::SolverBase *getSolver() { return solver; }
    };

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase {
       public:
        using tardigradeHydra::SolverStepBase::SolverStepBase;
    };

    class NonlinearStepBaseMock : public tardigradeHydra::NonlinearStepBase {
       public:
        using tardigradeHydra::NonlinearStepBase::NonlinearStepBase;
    };

    hydraBaseMock hydra(time, deltaTime, temperature, previousTemperature, deformationGradient,
                        previousDeformationGradient, {}, {}, previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables, dimension);

    SolverStepBaseMock step;
    NonlinearStepBaseMock  trial_step;

    step.trial_step = &trial_step;
    trial_step.step = &step;

    hydra.getSolver()->step = &step;
    step.setSolver(hydra.getSolver());

    tardigradeHydra::floatVector unknownVector = {0.39293837, -0.42772133, -0.54629709, 0.10262954, 0.43893794};

    tardigradeHydra::unit_test::hydraBaseTester::set_unknownVector(hydra, unknownVector);

    hydra.getSolver()->step->setUseLevenbergMarquardt(false);

    BOOST_TEST(hydra.residual == *trial_step.getNonlinearRHS(), CHECK_PER_ELEMENT);

    BOOST_TEST(hydra.jacobian == *trial_step.getFlatNonlinearLHS(), CHECK_PER_ELEMENT);
}
