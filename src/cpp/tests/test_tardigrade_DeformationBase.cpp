/**
 * \file test_tardigrade_DeformationBase.cpp
 *
 * Tests for tardigrade_DeformationBase
 */

#include <tardigrade_DeformationBase.h>

#define BOOST_TEST_MODULE test_tardigrade_DeformationBase
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element()

typedef tardigradeErrorTools::Node errorNode;  //!< Redefinition for the error node
typedef errorNode                 *errorOut;   //!< Redefinition for a pointer to the error node

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

        class DeformationBaseTester {
           public:

            static void checkConfigurations(DeformationBase &deformation) {
                BOOST_CHECK(&deformation._configurations.second == deformation.get_configurations());
            }

            static void checkPreviousConfigurations(DeformationBase &deformation) {
                BOOST_CHECK(&deformation._previousConfigurations.second == deformation.get_previousConfigurations());
            }

            static void checkInverseConfigurations(DeformationBase &deformation) {
                BOOST_CHECK(&deformation._inverseConfigurations.second == deformation.get_inverseConfigurations());
            }

            static void checkPreviousInverseConfigurations(DeformationBase &deformation) {
                BOOST_CHECK(&deformation._previousInverseConfigurations.second == deformation.get_previousInverseConfigurations());
            }

        };

    }

}

BOOST_AUTO_TEST_CASE(test_DeformationBase_get_configurations, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::DeformationBase deformation;

    tardigradeHydra::unit_test::DeformationBaseTester::checkConfigurations(deformation);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_get_previousConfigurations, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::DeformationBase deformation;

    tardigradeHydra::unit_test::DeformationBaseTester::checkPreviousConfigurations(deformation);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_get_inverseConfigurations, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::DeformationBase deformation;

    tardigradeHydra::unit_test::DeformationBaseTester::checkInverseConfigurations(deformation);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_get_previousInverseConfigurations,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::DeformationBase deformation;

    tardigradeHydra::unit_test::DeformationBaseTester::checkPreviousInverseConfigurations(deformation);
}

