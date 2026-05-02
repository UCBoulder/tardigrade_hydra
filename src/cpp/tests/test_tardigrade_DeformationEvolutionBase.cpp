/**
 * \file test_tardigrade_DeformationEvolutionBase.cpp
 *
 * Tests for tardigrade_DeformationEvolutionBase
 */

#include "tardigrade_DeformationEvolutionBase.h"
#include "tardigrade_hydra.h"

#define BOOST_TEST_MODULE test_tardigrade_DeformationEvolutionBase
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

        class hydraBaseTester {
           public:
            static auto getIterationData(hydraBase &hydra) { return hydra._iterationData; }
        };

        class DeformationEvolutionBaseTester {
           public:
        };

    }  // namespace unit_test

}  // namespace tardigradeHydra

BOOST_AUTO_TEST_CASE(test_DeformationEvolutionBase, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {

    class DeformationEvolution : public tardigradeHydra::DeformationEvolutionBase<tardigradeHydra::hydraBase, 3> { };

    tardigradeHydra::floatType dt = 1.45;

    tardigradeHydra::floatType alpha = 0.56;
    
    std::array<tardigradeHydra::floatType, 3 * 3> Ft = {1.1, 0.2, 0.3, 0.4, 1.5, 0.6, 0.7, 0.8, 1.9};

    std::array<tardigradeHydra::floatType, 3 * 3> Lt = {0.01, -0.02, -0.02, -0.04, -0.02, -0.  ,  0.02, -0.01, -0.04};

    std::array<tardigradeHydra::floatType, 3 * 3> Ltp1 = {0.04,  0.03,  0.05, -0.  , -0.01, -0.04,  0.01, -0.05,  0.02};

    std::array<tardigradeHydra::floatType, 3 * 3> answer = {1.1676671 , 0.24451068, 0.37037461, 0.34129644, 1.44041789, 0.51993393, 0.70062294, 0.72789142, 1.86367623};

    std::array<tardigradeHydra::floatType, 3 * 3> result;

    DeformationEvolution de;

    de.integration_parameter = alpha;

    de.computeDeformation(dt, std::begin(Lt), std::end(Lt), std::begin(Ltp1), std::end(Ltp1), std::begin(Ft), std::end(Ft), std::begin(result), std::end(result));

    BOOST_TEST(result == answer, CHECK_PER_ELEMENT);

}
