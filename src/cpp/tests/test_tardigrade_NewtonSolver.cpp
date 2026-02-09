/**
 * \file test_tardigrade_NewtonSolver.cpp
 *
 * Tests for tardigrade_NewtonSolver
 */

#include "tardigrade_NewtonSolver.h"

#define BOOST_TEST_MODULE test_tardigrade_NewtonSolver
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

BOOST_AUTO_TEST_CASE(test_configuration, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {

    tardigradeHydra::NewtonSolver solver;

    auto step    = dynamic_cast<tardigradeHydra::NewtonStep*>(solver.step->trial_step);
    auto damping = dynamic_cast<tardigradeHydra::ArmijoGradientDamping*>(solver.step->damping);

    BOOST_TEST(step != nullptr);

    BOOST_TEST(damping != nullptr);

}
