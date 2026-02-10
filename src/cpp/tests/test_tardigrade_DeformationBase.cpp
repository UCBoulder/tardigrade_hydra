/**
 * \file test_tardigrade_DeformationBase.cpp
 *
 * Tests for tardigrade_DeformationBase
 */

#include "tardigrade_DeformationBase.h"
#include "tardigrade_hydra.h"

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

        class hydraBaseTester {

            public:

                static auto getIterationData(hydraBase &hydra) {
                    return hydra._iterationData;
                }

        };

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

            static void checkdF1dF(hydraBase &hydra, DeformationBase &deformation) {
                auto iterationData = tardigradeHydra::unit_test::hydraBaseTester::getIterationData(hydra);
                BOOST_CHECK(std::find(iterationData.begin(), iterationData.end(), &deformation._dF1dF) !=
                            iterationData.end());
            }

            static void checkdF1dFn(hydraBase &hydra, DeformationBase &deformation) {
                auto iterationData = tardigradeHydra::unit_test::hydraBaseTester::getIterationData(hydra);
                BOOST_CHECK(std::find(iterationData.begin(), iterationData.end(), &deformation._dF1dFn) !=
                            iterationData.end());
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

BOOST_AUTO_TEST_CASE(test_DeformationBase_setFirstConfigurationGradients,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                       -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,  -0.12285551, -0.88064421,
                                               -0.20391149, 0.47599081,  -0.63501654, -0.64909649};

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

    tardigradeHydra::hydraBase hydra(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                     previousDeformationGradient, {}, {}, previousStateVariables, parameters,
                                     numConfigurations, numNonLinearSolveStateVariables, dimension);

    hydra.initialize();

    tardigradeHydra::floatType eps = 1e-6;

    tardigradeHydra::floatMatrix dF1dF_answer(deformationGradient.size(), tardigradeHydra::floatVector(deformationGradient.size(), 0));

    tardigradeHydra::floatMatrix dF1dFn_answer(deformationGradient.size(),
                              tardigradeHydra::floatVector(deformationGradient.size() * (numConfigurations - 1), 0));

    for (unsigned int i = 0; i < deformationGradient.size(); i++) {
        tardigradeHydra::floatVector delta(deformationGradient.size(), 0);

        delta[i] = eps * std::fabs(deformationGradient[i]) + eps;

        tardigradeHydra::hydraBase hydra_p(time, deltaTime, temperature, previousTemperature,
                                           deformationGradient + delta, previousDeformationGradient, {}, {},
                                           previousStateVariables, parameters, numConfigurations,
                                           numNonLinearSolveStateVariables, dimension);

        tardigradeHydra::hydraBase hydra_m(time, deltaTime, temperature, previousTemperature,
                                           deformationGradient - delta, previousDeformationGradient, {}, {},
                                           previousStateVariables, parameters, numConfigurations,
                                           numNonLinearSolveStateVariables, dimension);

        hydra_p.initialize();

        hydra_m.initialize();

        tardigradeHydra::floatVector F1_p, F1_m;

        F1_p = tardigradeHydra::floatVector(hydra_p.deformation->get_configurations()->begin(), hydra_p.deformation->get_configurations()->begin() + 9);

        F1_m = tardigradeHydra::floatVector(hydra_m.deformation->get_configurations()->begin(), hydra_m.deformation->get_configurations()->begin() + 9);

        for (unsigned int j = 0; j < F1_p.size(); j++) {
            dF1dF_answer[j][i] = (F1_p[j] - F1_m[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(dF1dF_answer) == *hydra.deformation->get_dF1dF(), CHECK_PER_ELEMENT);

    tardigradeHydra::unit_test::DeformationBaseTester::checkdF1dF(hydra, *hydra.deformation);

    for (unsigned int i = 0; i < (numConfigurations - 1) * deformationGradient.size(); i++) {
        tardigradeHydra::floatVector delta(previousStateVariables.size(), 0);

        delta[i] += eps * std::fabs(previousStateVariables[i]) + eps;

        tardigradeHydra::hydraBase hydra_p(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                           previousDeformationGradient, {}, {}, previousStateVariables + delta,
                                           parameters, numConfigurations, numNonLinearSolveStateVariables, dimension);

        tardigradeHydra::hydraBase hydra_m(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                           previousDeformationGradient, {}, {}, previousStateVariables - delta,
                                           parameters, numConfigurations, numNonLinearSolveStateVariables, dimension);

        hydra_p.initialize();

        hydra_m.initialize();

        tardigradeHydra::floatVector F1_p, F1_m;

        F1_p = tardigradeHydra::floatVector(hydra_p.deformation->get_configurations()->begin(), hydra_p.deformation->get_configurations()->begin() + 9);

        F1_m = tardigradeHydra::floatVector(hydra_m.deformation->get_configurations()->begin(), hydra_m.deformation->get_configurations()->begin() + 9);

        for (unsigned int j = 0; j < F1_p.size(); j++) {
            dF1dFn_answer[j][i] = (F1_p[j] - F1_m[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(dF1dFn_answer) == *hydra.deformation->get_dF1dFn(), CHECK_PER_ELEMENT);

    tardigradeHydra::unit_test::DeformationBaseTester::checkdF1dFn(hydra, *hydra.deformation);
}
