/**
 * \file test_tardigrade_PreconditionerBase.cpp
 *
 * Tests for tardigrade_PreconditionerBase
 */

#include <tardigrade_PreconditionerBase.h>
#include <tardigrade_hydra.h>

#define BOOST_TEST_MODULE test_tardigrade_PreconditionerBase
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

        class PreconditionerBaseTester {
           public:
            static void checkGetFlatPreconditioner(PreconditionerBase &preconditioner) {
                BOOST_CHECK(&preconditioner._preconditioner.second == preconditioner.getFlatPreconditioner());
            }

            static void set_preconditioner_nohydra(PreconditionerBase &preconditioner, floatVector &value) {
                preconditioner._preconditioner.second = value;
                preconditioner._preconditioner.first  = true;
            }

            static void set_preconditioner(PreconditionerBase &preconditioner, floatVector &value) {
                set_preconditioner_nohydra(preconditioner, value);

                preconditioner.addIterationData(&preconditioner._preconditioner);
            }

            static void set_flatPreconditioner(PreconditionerBase &preconditioner, const floatVector &value) {
                preconditioner._preconditioner.second = value;
                preconditioner._preconditioner.first  = true;

                preconditioner.addIterationData(&preconditioner._preconditioner);
            }
        };

        class hydraBaseTester {
           public:
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

        class SolverBaseTester {
           public:
            static hydraBase *get_hydra(SolverBase &solver) { return solver.hydra; }
        };

        class SolverStepBaseTester {
           public:
            static SolverBase *get_solver(SolverStepBase &step) { return step.solver; }
        };

    }  // namespace unit_test

}  // namespace tardigradeHydra

BOOST_AUTO_TEST_CASE(test_PreconditionerBase_getFlatPreconditioner,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class PreconditionerMock : public tardigradeHydra::PreconditionerBase {
       public:
        using tardigradeHydra::PreconditionerBase::PreconditionerBase;

        tardigradeHydra::floatVector expected_value = {1, 2, 3};

        virtual void formPreconditioner() override {
            tardigradeHydra::unit_test::PreconditionerBaseTester::set_preconditioner_nohydra(*this, expected_value);
        }
    };

    PreconditionerMock preconditioner;

    tardigradeHydra::unit_test::PreconditionerBaseTester::checkGetFlatPreconditioner(preconditioner);

    std::vector<double> answer = {1, 2, 3};

    BOOST_TEST(answer == *preconditioner.getFlatPreconditioner(), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_PreconditionerBase_formPreconditioner, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    /*!
     * Boost test of the get-configuration command
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

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient, previousDeformationGradient);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    class PreconditionerBaseMock : public tardigradeHydra::PreconditionerBase {
       public:
        using tardigradeHydra::PreconditionerBase::PreconditionerBase;

        tardigradeHydra::floatVector expected_preconditioner = {1, 2, 3, 4, 5, 6, 7};

        virtual void formPreconditioner() override {
            tardigradeHydra::unit_test::PreconditionerBaseTester::set_flatPreconditioner(*this,
                                                                                         expected_preconditioner);
        }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;
    };

    hydraBaseMock hydra(dof,  {}, {}, previousStateVariables,
                        parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::IterativeSolverBase solver(&hydra);
    tardigradeHydra::SolverStepBase      step(&solver);
    solver.step = &step;
    tardigradeHydra::NonlinearStepBase trial_step(&step);
    PreconditionerBaseMock             preconditioner(&trial_step);
    tardigradeHydra::StepDampingBase   damping(&step);

    hydra.solver = &solver;

    BOOST_TEST(preconditioner.expected_preconditioner == *preconditioner.getFlatPreconditioner(), CHECK_PER_ELEMENT);
}
