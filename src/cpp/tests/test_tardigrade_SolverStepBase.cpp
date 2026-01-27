
/**
 * \file test_tardigrade_SolverStepBase.cpp
 *
 * Tests for tardigrade_SolverStepBase
 */

#include "tardigrade_SolverBase.h"
#include "tardigrade_SolverStepBase.h"
#include "tardigrade_NonlinearStepBase.h"
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

    }  // namespace unit_test

}  // namespace tardigradeHydra

BOOST_AUTO_TEST_CASE(test_SolverStepBase_placeholder, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) { }
