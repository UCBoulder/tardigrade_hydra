/**
 * \file test_tardigrade_MaxRowPreconditioner.cpp
 *
 * Tests for tardigrade_MaxRowPreconditioner
 */

#include <tardigrade_MaxRowPreconditioner.h>
#include <tardigrade_hydra.h>

#define BOOST_TEST_MODULE test_tardigrade_MaxRowPreconditioner
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

            static void set_preconditionerType(PreconditionerBase &preconditioner, const unsigned int &value) {
                preconditioner._preconditioner_type = value;
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


BOOST_AUTO_TEST_CASE(test_placeholder,*boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {}
