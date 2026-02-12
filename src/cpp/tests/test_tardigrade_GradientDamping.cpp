/**
 * \file test_tardigrade_GradientDamping.cpp
 *
 * Tests for tardigrade_GradientDamping
 */

#include "tardigrade_GradientDamping.h"
#include "tardigrade_ResidualBase.h"
#include "tardigrade_hydra.h"

#define BOOST_TEST_MODULE test_tardigrade_GradientDamping
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

        class GradientDampingTester {
           public:
            static void setMuk(GradientDamping &damping, const floatType &value) { damping.setMuk(value); }
        };

        class hydraBaseTester {
           public:
            static void set_residual(hydraBase &hydra, const floatVector &value) {
                hydra._residual.second = value;
                hydra._residual.first  = true;

                hydra.addIterationData(&hydra._residual);
            }

            static void set_unknownVector(hydraBase &hydra, const floatVector &value) {
                hydra._X.second = value;
                hydra._X.first  = true;
            }

            static void set_flatJacobian(hydraBase &hydra, const floatVector &value) {
                hydra._jacobian.second = value;
                hydra._jacobian.first  = true;

                hydra.addIterationData(&hydra._jacobian);
            }
        };

    }  // namespace unit_test

}  // namespace tardigradeHydra

BOOST_AUTO_TEST_CASE(test_GradientDamping_setBaseQuantities, *boost::unit_test::tolerance(1e-5)) {
    /*!
     * Test checking the gradient convergence
     */

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

    class GradientDampingMock : public tardigradeHydra::GradientDamping {
       public:
        tardigradeHydra::floatType rnorm = 10.3;

        tardigradeHydra::floatVector dRNormdX = {1, 2, 3};

        virtual void setResidualNorm() override { set_residualNorm(rnorm); }

        virtual void setdResidualNormdX() override {
            set_dResidualNormdX(dRNormdX);
            ;
        }

        void public_setMuk(const tardigradeHydra::floatType &value) { setMuk(value); }

        void runSetBaseQuantities() { setBaseQuantities(); }
    };

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

            auto local_damping = dynamic_cast<GradientDampingMock *>(solver->step->damping);

            if (local_damping == nullptr) {
                throw std::runtime_error("dynamic cast failed");
            }

            local_damping->public_setMuk(mu_k);
        }

        virtual void               decomposeUnknownVector() override { return; }
        virtual const unsigned int getNumUnknowns() override { return 5; }
    };

    tardigradeHydra::floatType answer1 = 0.5 * 1e-8 * 10.3;

    hydraBaseMock hydra(dof, deformationGradient, previousDeformationGradient, {}, {},
                        previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    GradientDampingMock             damping;
    tardigradeHydra::SolverStepBase step;
    tardigradeHydra::SolverBase     solver;
    hydra.solver = &solver;
    solver.hydra = &hydra;
    solver.step  = &step;
    step.solver  = &solver;
    step.damping = &damping;
    damping.step = &step;

    damping.runSetBaseQuantities();

    BOOST_TEST(answer1 == damping.getMuk());

    BOOST_TEST(damping.rnorm == *damping.get_baseResidualNorm());

    BOOST_TEST(damping.dRNormdX == *damping.get_basedResidualNormdX());

    damping.rnorm = 1e-9;

    damping.setResidualNorm();

    damping.runSetBaseQuantities();

    BOOST_TEST(damping.rnorm == damping.getMuk());
}
