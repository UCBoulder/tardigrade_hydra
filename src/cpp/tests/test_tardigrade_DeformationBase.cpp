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
            static auto getIterationData(hydraBase &hydra) { return hydra._iterationData; }
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
                BOOST_CHECK(&deformation._previousInverseConfigurations.second ==
                            deformation.get_previousInverseConfigurations());
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

    }  // namespace unit_test

}  // namespace tardigradeHydra

BOOST_AUTO_TEST_CASE(test_DeformationBase_get_configurations, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::DeformationBase deformation;

    tardigradeHydra::unit_test::DeformationBaseTester::checkConfigurations(deformation);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_get_previousConfigurations,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::DeformationBase deformation;

    tardigradeHydra::unit_test::DeformationBaseTester::checkPreviousConfigurations(deformation);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_get_inverseConfigurations,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::DeformationBase deformation;

    tardigradeHydra::unit_test::DeformationBaseTester::checkInverseConfigurations(deformation);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_get_previousInverseConfigurations,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::DeformationBase deformation;

    tardigradeHydra::unit_test::DeformationBaseTester::checkPreviousInverseConfigurations(deformation);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_getSubConfiguration, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    auto result = hydra.deformation->getSubConfiguration(0, 4);

    BOOST_TEST(result == *hydra.getDeformationGradient(), CHECK_PER_ELEMENT);

    tardigradeHydra::floatVector answer1 = {2.24332648, 1.48246714, 2.02801682, 1.50380989, 2.98203598,
                                            2.08079721, 1.58939152, 1.2551092,  2.38201794};

    auto result1 = hydra.deformation->getSubConfiguration(1, 3);

    BOOST_TEST(result1 == answer1, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_getPrecedingConfiguration,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    BOOST_TEST(hydra.deformation->getPrecedingConfiguration(4) == *hydra.getDeformationGradient(), CHECK_PER_ELEMENT);

    tardigradeHydra::floatVector answer1 = {0.73947165, -0.24328161, -0.68904986, 0.04049558, 0.25403825,
                                            -0.1748363, 0.97015752,  -0.04452644, -0.77301275};

    BOOST_TEST(hydra.deformation->getPrecedingConfiguration(2) == answer1, CHECK_PER_ELEMENT);

    tardigradeHydra::floatVector answer2 = {0.54568475,  -0.42501821, -0.54244544, -0.01317666, 0.30165847,
                                            -0.09442353, 0.80588282,  -0.10806097, -0.42143322};

    BOOST_TEST(hydra.deformation->getPrecedingConfiguration(3) == answer2, CHECK_PER_ELEMENT);

    tardigradeHydra::floatVector answer3 = {1, 0, 0, 0, 1, 0, 0, 0, 1};

    BOOST_TEST(hydra.deformation->getPrecedingConfiguration(0) == answer3, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_getFollowingConfiguration,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::floatVector answer1 = {2.09953091, 1.83604029, 2.3712323,  0.98756433, 2.58928197,
                                            1.05684715, 1.33422708, 1.67694162, 2.96443669};

    BOOST_TEST(hydra.deformation->getFollowingConfiguration(1) == answer1, CHECK_PER_ELEMENT);

    tardigradeHydra::floatVector answer2 = {1.42635131, 0.89338916, 0.94416002, 0.50183668, 1.62395295,
                                            0.1156184,  0.31728548, 0.41482621, 1.86630916};

    BOOST_TEST(hydra.deformation->getFollowingConfiguration(2) == answer2, CHECK_PER_ELEMENT);

    tardigradeHydra::floatVector answer3 = {1, 0, 0, 0, 1, 0, 0, 0, 1};

    BOOST_TEST(hydra.deformation->getFollowingConfiguration(3) == answer3, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_getPreviousSubConfiguration,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    BOOST_TEST(hydra.deformation->getPreviousSubConfiguration(0, 4) == *hydra.getPreviousDeformationGradient(),
               CHECK_PER_ELEMENT);

    tardigradeHydra::floatVector answer1 = {2.24332648, 1.48246714, 2.02801682, 1.50380989, 2.98203598,
                                            2.08079721, 1.58939152, 1.2551092,  2.38201794};

    BOOST_TEST(hydra.deformation->getPreviousSubConfiguration(1, 3) == answer1, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_getPreviousPrecedingConfiguration,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    BOOST_TEST(hydra.deformation->getPreviousPrecedingConfiguration(4) == *hydra.getPreviousDeformationGradient(),
               CHECK_PER_ELEMENT);

    tardigradeHydra::floatVector answer1 = {-0.30350143, -0.21223491, 0.47296395,  0.18299993, -0.42974886,
                                            -0.06195712, 0.92470041,  -0.36418391, -0.8287879};

    BOOST_TEST(hydra.deformation->getPreviousPrecedingConfiguration(2) == answer1, CHECK_PER_ELEMENT);

    tardigradeHydra::floatVector answer2 = {-0.15883228, -0.1920217, 0.33770597,  0.15460279, -0.58876502,
                                            -0.15099813, 0.69307214, -0.60345639, -0.66103563};

    BOOST_TEST(hydra.deformation->getPreviousPrecedingConfiguration(3) == answer2, CHECK_PER_ELEMENT);

    tardigradeHydra::floatVector answer3 = {1, 0, 0, 0, 1, 0, 0, 0, 1};

    BOOST_TEST(hydra.deformation->getPreviousPrecedingConfiguration(0) == answer3, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_getPreviousFollowingConfiguration,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::floatVector answer1 = {2.09953091, 1.83604029, 2.3712323,  0.98756433, 2.58928197,
                                            1.05684715, 1.33422708, 1.67694162, 2.96443669};

    BOOST_TEST(hydra.deformation->getPreviousFollowingConfiguration(1) == answer1, CHECK_PER_ELEMENT);

    tardigradeHydra::floatVector answer2 = {1.42635131, 0.89338916, 0.94416002, 0.50183668, 1.62395295,
                                            0.1156184,  0.31728548, 0.41482621, 1.86630916};

    BOOST_TEST(hydra.deformation->getPreviousFollowingConfiguration(2) == answer2, CHECK_PER_ELEMENT);

    tardigradeHydra::floatVector answer3 = {1, 0, 0, 0, 1, 0, 0, 0, 1};

    BOOST_TEST(hydra.deformation->getPreviousFollowingConfiguration(3) == answer3, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_getSubConfigurationJacobian,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::floatVector configurations = *hydra.deformation->get_configurations();

    tardigradeHydra::floatVector x = configurations;

    tardigradeHydra::floatType eps = 1e-6;

    const unsigned int dimension = 3;

    tardigradeHydra::floatMatrix gradient(dimension * dimension, tardigradeHydra::floatVector(x.size(), 0));

    unsigned int lower = 0;

    unsigned int upper = 4;

    for (unsigned int i = 0; i < x.size(); i++) {
        tardigradeHydra::floatVector delta(numConfigurations * dimension * dimension, 0);

        delta[i] = std::fabs(eps * configurations[i]) + eps;

        tardigradeHydra::floatVector Fscp;

        tardigradeHydra::floatVector Fscm;

        Fscp = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations + delta, lower, upper);

        Fscm = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations - delta, lower, upper);

        for (unsigned int j = 0; j < (dimension * dimension); j++) {
            gradient[j][i] = (Fscp[j] - Fscm[j]) / (2 * delta[i]);
        }
    }

    auto result = hydra.deformation->getSubConfigurationJacobian<3, 3, 3>(configurations, lower, upper);

    BOOST_TEST(tardigradeVectorTools::appendVectors(gradient) == result, CHECK_PER_ELEMENT);

    lower = 1;
    upper = 3;

    for (unsigned int i = 0; i < x.size(); i++) {
        tardigradeHydra::floatVector delta(numConfigurations * dimension * dimension, 0);

        delta[i] = std::fabs(eps * configurations[i]) + eps;

        tardigradeHydra::floatVector Fscp;

        tardigradeHydra::floatVector Fscm;

        Fscp = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations + delta, lower, upper);

        Fscm = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations - delta, lower, upper);

        for (unsigned int j = 0; j < (dimension * dimension); j++) {
            gradient[j][i] = (Fscp[j] - Fscm[j]) / (2 * delta[i]);
        }
    }

    result = hydra.deformation->getSubConfigurationJacobian<3, 3, 3>(configurations, lower, upper);

    BOOST_TEST(tardigradeVectorTools::appendVectors(gradient) == result, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_getSubConfigurationJacobian2,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::floatVector configurations = *hydra.deformation->get_configurations();

    tardigradeHydra::floatVector x = configurations;

    tardigradeHydra::floatType eps = 1e-6;

    const unsigned int dimension = 3;

    tardigradeHydra::floatMatrix gradient(dimension * dimension, tardigradeHydra::floatVector(x.size(), 0));

    unsigned int lower = 0;

    unsigned int upper = 4;

    for (unsigned int i = 0; i < x.size(); i++) {
        tardigradeHydra::floatVector delta(numConfigurations * dimension * dimension, 0);

        delta[i] = std::fabs(eps * configurations[i]) + eps;

        tardigradeHydra::floatVector Fscp;

        tardigradeHydra::floatVector Fscm;

        Fscp = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations + delta, lower, upper);

        Fscm = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations - delta, lower, upper);

        for (unsigned int j = 0; j < (dimension * dimension); j++) {
            gradient[j][i] = (Fscp[j] - Fscm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(gradient) ==
                   hydra.deformation->getSubConfigurationJacobian(lower, upper),
               CHECK_PER_ELEMENT);

    lower = 1;
    upper = 3;

    for (unsigned int i = 0; i < x.size(); i++) {
        tardigradeHydra::floatVector delta(numConfigurations * dimension * dimension, 0);

        delta[i] = std::fabs(eps * configurations[i]) + eps;

        tardigradeHydra::floatVector Fscp;

        tardigradeHydra::floatVector Fscm;

        Fscp = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations + delta, lower, upper);

        Fscm = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations - delta, lower, upper);

        for (unsigned int j = 0; j < (dimension * dimension); j++) {
            gradient[j][i] = (Fscp[j] - Fscm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(gradient) ==
                   hydra.deformation->getSubConfigurationJacobian(lower, upper),
               CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_getPrecedingConfigurationJacobian,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::floatVector configurations = *hydra.deformation->get_configurations();

    tardigradeHydra::floatVector x = configurations;

    tardigradeHydra::floatType eps = 1e-6;

    const unsigned int dimension = 3;

    tardigradeHydra::floatMatrix gradient(dimension * dimension, tardigradeHydra::floatVector(x.size(), 0));

    unsigned int lower = 0;

    unsigned int upper = 4;

    for (unsigned int i = 0; i < x.size(); i++) {
        tardigradeHydra::floatVector delta(numConfigurations * dimension * dimension, 0);

        delta[i] = std::fabs(eps * configurations[i]) + eps;

        tardigradeHydra::floatVector Fscp;

        tardigradeHydra::floatVector Fscm;

        Fscp = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations + delta, lower, upper);

        Fscm = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations - delta, lower, upper);

        for (unsigned int j = 0; j < (dimension * dimension); j++) {
            gradient[j][i] = (Fscp[j] - Fscm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(gradient) ==
                   hydra.deformation->getPrecedingConfigurationJacobian(upper),
               CHECK_PER_ELEMENT);

    lower = 0;
    upper = 3;

    for (unsigned int i = 0; i < x.size(); i++) {
        tardigradeHydra::floatVector delta(numConfigurations * dimension * dimension, 0);

        delta[i] = std::fabs(eps * configurations[i]) + eps;

        tardigradeHydra::floatVector Fscp;

        tardigradeHydra::floatVector Fscm;

        Fscp = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations + delta, lower, upper);

        Fscm = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations - delta, lower, upper);

        for (unsigned int j = 0; j < (dimension * dimension); j++) {
            gradient[j][i] = (Fscp[j] - Fscm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(gradient) ==
                   hydra.deformation->getPrecedingConfigurationJacobian(upper),
               CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_getFollowingConfigurationJacobian,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::floatVector configurations = *hydra.deformation->get_configurations();

    tardigradeHydra::floatVector x = configurations;

    tardigradeHydra::floatType eps = 1e-6;

    const unsigned int dimension = 3;

    tardigradeHydra::floatMatrix gradient(dimension * dimension, tardigradeHydra::floatVector(x.size(), 0));

    unsigned int lower = 1;

    unsigned int upper = 4;

    for (unsigned int i = 0; i < x.size(); i++) {
        tardigradeHydra::floatVector delta(numConfigurations * dimension * dimension, 0);

        delta[i] = std::fabs(eps * configurations[i]) + eps;

        tardigradeHydra::floatVector Fscp;

        tardigradeHydra::floatVector Fscm;

        Fscp = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations + delta, lower + 1, upper);

        Fscm = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations - delta, lower + 1, upper);

        for (unsigned int j = 0; j < (dimension * dimension); j++) {
            gradient[j][i] = (Fscp[j] - Fscm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(gradient) ==
                   hydra.deformation->getFollowingConfigurationJacobian(lower),
               CHECK_PER_ELEMENT);

    lower = 2;
    upper = 4;

    for (unsigned int i = 0; i < x.size(); i++) {
        tardigradeHydra::floatVector delta(numConfigurations * dimension * dimension, 0);

        delta[i] = std::fabs(eps * configurations[i]) + eps;

        tardigradeHydra::floatVector Fscp;

        tardigradeHydra::floatVector Fscm;

        Fscp = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations + delta, lower + 1, upper);

        Fscm = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations - delta, lower + 1, upper);

        for (unsigned int j = 0; j < (dimension * dimension); j++) {
            gradient[j][i] = (Fscp[j] - Fscm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(gradient) ==
                   hydra.deformation->getFollowingConfigurationJacobian(lower),
               CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_getPreviousSubConfigurationJacobian,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::floatVector configurations = *hydra.deformation->get_previousConfigurations();

    tardigradeHydra::floatVector x = configurations;

    tardigradeHydra::floatType eps = 1e-6;

    const unsigned int dimension = 3;

    tardigradeHydra::floatMatrix gradient(dimension * dimension, tardigradeHydra::floatVector(x.size(), 0));

    unsigned int lower = 0;

    unsigned int upper = 4;

    for (unsigned int i = 0; i < x.size(); i++) {
        tardigradeHydra::floatVector delta(numConfigurations * dimension * dimension, 0);

        delta[i] = std::fabs(eps * configurations[i]) + eps;

        tardigradeHydra::floatVector Fscp;

        tardigradeHydra::floatVector Fscm;

        Fscp = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations + delta, lower, upper);

        Fscm = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations - delta, lower, upper);

        for (unsigned int j = 0; j < (dimension * dimension); j++) {
            gradient[j][i] = (Fscp[j] - Fscm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(gradient) ==
                   hydra.deformation->getPreviousSubConfigurationJacobian(lower, upper),
               CHECK_PER_ELEMENT);

    lower = 1;
    upper = 3;

    for (unsigned int i = 0; i < x.size(); i++) {
        tardigradeHydra::floatVector delta(numConfigurations * dimension * dimension, 0);

        delta[i] = std::fabs(eps * configurations[i]) + eps;

        tardigradeHydra::floatVector Fscp;

        tardigradeHydra::floatVector Fscm;

        Fscp = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations + delta, lower, upper);

        Fscm = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations - delta, lower, upper);

        for (unsigned int j = 0; j < (dimension * dimension); j++) {
            gradient[j][i] = (Fscp[j] - Fscm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(gradient) ==
                   hydra.deformation->getPreviousSubConfigurationJacobian(lower, upper),
               CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_getPreviousPrecedingConfigurationJacobian,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::floatVector configurations = *hydra.deformation->get_previousConfigurations();

    tardigradeHydra::floatVector x = configurations;

    tardigradeHydra::floatType eps = 1e-6;

    const unsigned int dimension = 3;

    tardigradeHydra::floatMatrix gradient(dimension * dimension, tardigradeHydra::floatVector(x.size(), 0));

    unsigned int lower = 0;

    unsigned int upper = 4;

    for (unsigned int i = 0; i < x.size(); i++) {
        tardigradeHydra::floatVector delta(numConfigurations * dimension * dimension, 0);

        delta[i] = std::fabs(eps * configurations[i]) + eps;

        tardigradeHydra::floatVector Fscp;

        tardigradeHydra::floatVector Fscm;

        Fscp = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations + delta, lower, upper);

        Fscm = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations - delta, lower, upper);

        for (unsigned int j = 0; j < (dimension * dimension); j++) {
            gradient[j][i] = (Fscp[j] - Fscm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(gradient) ==
                   hydra.deformation->getPreviousPrecedingConfigurationJacobian(upper),
               CHECK_PER_ELEMENT);

    lower = 0;
    upper = 3;

    for (unsigned int i = 0; i < x.size(); i++) {
        tardigradeHydra::floatVector delta(numConfigurations * dimension * dimension, 0);

        delta[i] = std::fabs(eps * configurations[i]) + eps;

        tardigradeHydra::floatVector Fscp;

        tardigradeHydra::floatVector Fscm;

        Fscp = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations + delta, lower, upper);

        Fscm = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations - delta, lower, upper);

        for (unsigned int j = 0; j < (dimension * dimension); j++) {
            gradient[j][i] = (Fscp[j] - Fscm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(gradient) ==
                   hydra.deformation->getPreviousPrecedingConfigurationJacobian(upper),
               CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_getPreviousFollowingConfigurationJacobian,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::floatVector configurations = *hydra.deformation->get_previousConfigurations();

    tardigradeHydra::floatVector x = configurations;

    tardigradeHydra::floatType eps = 1e-6;

    const unsigned int dimension = 3;

    tardigradeHydra::floatMatrix gradient(dimension * dimension, tardigradeHydra::floatVector(x.size(), 0));

    unsigned int lower = 1;

    unsigned int upper = 4;

    for (unsigned int i = 0; i < x.size(); i++) {
        tardigradeHydra::floatVector delta(numConfigurations * dimension * dimension);

        delta[i] = std::fabs(eps * configurations[i]) + eps;

        tardigradeHydra::floatVector Fscp;

        tardigradeHydra::floatVector Fscm;

        Fscp = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations + delta, lower + 1, upper);

        Fscm = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations - delta, lower + 1, upper);

        for (unsigned int j = 0; j < (dimension * dimension); j++) {
            gradient[j][i] = (Fscp[j] - Fscm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(gradient) ==
                   hydra.deformation->getPreviousFollowingConfigurationJacobian(lower),
               CHECK_PER_ELEMENT);

    lower = 2;
    upper = 4;

    for (unsigned int i = 0; i < x.size(); i++) {
        tardigradeHydra::floatVector delta(numConfigurations * dimension * dimension);

        delta[i] = std::fabs(eps * configurations[i]) + eps;

        tardigradeHydra::floatVector Fscp;

        tardigradeHydra::floatVector Fscm;

        Fscp = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations + delta, lower + 1, upper);

        Fscm = hydra.deformation->getSubConfiguration<3, 3, 3>(configurations - delta, lower + 1, upper);

        for (unsigned int j = 0; j < (dimension * dimension); j++) {
            gradient[j][i] = (Fscp[j] - Fscm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(gradient) ==
                   hydra.deformation->getPreviousFollowingConfigurationJacobian(lower),
               CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_setPreviousFirstConfigurationGradients,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::floatType eps = 1e-6;

    tardigradeHydra::floatMatrix previousdF1dF_answer(deformationGradient.size(),
                                                      tardigradeHydra::floatVector(deformationGradient.size(), 0));

    tardigradeHydra::floatMatrix previousdF1dFn_answer(
        deformationGradient.size(),
        tardigradeHydra::floatVector(deformationGradient.size() * (numConfigurations - 1), 0));

    for (unsigned int i = 0; i < previousDeformationGradient.size(); i++) {
        tardigradeHydra::floatVector delta(previousDeformationGradient.size(), 0);

        delta[i] = eps * std::fabs(previousDeformationGradient[i]) + eps;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient + delta, additionalDOF, previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient - delta, additionalDOF, previousAdditionalDOF);

        tardigradeHydra::hydraBase hydra_p(dofp, model_configuration);

        tardigradeHydra::hydraBase hydra_m(dofm, model_configuration);

        hydra_p.initialize();

        hydra_m.initialize();

        tardigradeHydra::floatVector F1_p, F1_m;

        F1_p = tardigradeHydra::floatVector(hydra_p.deformation->get_previousConfigurations()->begin(),
                                            hydra_p.deformation->get_previousConfigurations()->begin() + 9);

        F1_m = tardigradeHydra::floatVector(hydra_m.deformation->get_previousConfigurations()->begin(),
                                            hydra_m.deformation->get_previousConfigurations()->begin() + 9);

        for (unsigned int j = 0; j < F1_p.size(); j++) {
            previousdF1dF_answer[j][i] = (F1_p[j] - F1_m[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(previousdF1dF_answer) == *hydra.deformation->get_previousdF1dF(),
               CHECK_PER_ELEMENT);

    for (unsigned int i = 0; i < (numConfigurations - 1) * deformationGradient.size(); i++) {
        tardigradeHydra::floatVector delta(previousStateVariables.size(), 0);

        delta[i] += eps * std::fabs(previousStateVariables[i]) + eps;

        tardigradeHydra::ModelConfigurationBase model_configurationp(previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables);

        tardigradeHydra::ModelConfigurationBase model_configurationm(previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables);

        tardigradeHydra::hydraBase hydra_p(dof, model_configurationp);

        tardigradeHydra::hydraBase hydra_m(dof, model_configurationm);

        hydra_p.initialize();

        hydra_m.initialize();

        tardigradeHydra::floatVector F1_p, F1_m;

        F1_p = tardigradeHydra::floatVector(hydra_p.deformation->get_previousConfigurations()->begin(),
                                            hydra_p.deformation->get_previousConfigurations()->begin() + 9);

        F1_m = tardigradeHydra::floatVector(hydra_m.deformation->get_previousConfigurations()->begin(),
                                            hydra_m.deformation->get_previousConfigurations()->begin() + 9);

        for (unsigned int j = 0; j < F1_p.size(); j++) {
            previousdF1dFn_answer[j][i] = (F1_p[j] - F1_m[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(previousdF1dFn_answer) == *hydra.deformation->get_previousdF1dFn(),
               CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_DeformationBase_setFirstConfigurationGradients,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::hydraBase hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::floatType eps = 1e-6;

    tardigradeHydra::floatMatrix dF1dF_answer(deformationGradient.size(),
                                              tardigradeHydra::floatVector(deformationGradient.size(), 0));

    tardigradeHydra::floatMatrix dF1dFn_answer(
        deformationGradient.size(),
        tardigradeHydra::floatVector(deformationGradient.size() * (numConfigurations - 1), 0));

    for (unsigned int i = 0; i < deformationGradient.size(); i++) {
        tardigradeHydra::floatVector delta(deformationGradient.size(), 0);

        delta[i] = eps * std::fabs(deformationGradient[i]) + eps;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature,
                                             deformationGradient + delta, previousDeformationGradient, additionalDOF,
                                             previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature,
                                             deformationGradient - delta, previousDeformationGradient, additionalDOF,
                                             previousAdditionalDOF);

        tardigradeHydra::hydraBase hydra_p(dofp, model_configuration);

        tardigradeHydra::hydraBase hydra_m(dofm, model_configuration);

        hydra_p.initialize();

        hydra_m.initialize();

        tardigradeHydra::floatVector F1_p, F1_m;

        F1_p = tardigradeHydra::floatVector(hydra_p.deformation->get_configurations()->begin(),
                                            hydra_p.deformation->get_configurations()->begin() + 9);

        F1_m = tardigradeHydra::floatVector(hydra_m.deformation->get_configurations()->begin(),
                                            hydra_m.deformation->get_configurations()->begin() + 9);

        for (unsigned int j = 0; j < F1_p.size(); j++) {
            dF1dF_answer[j][i] = (F1_p[j] - F1_m[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(dF1dF_answer) == *hydra.deformation->get_dF1dF(),
               CHECK_PER_ELEMENT);

    tardigradeHydra::unit_test::DeformationBaseTester::checkdF1dF(hydra, *hydra.deformation);

    for (unsigned int i = 0; i < (numConfigurations - 1) * deformationGradient.size(); i++) {
        tardigradeHydra::floatVector delta(previousStateVariables.size(), 0);

        delta[i] += eps * std::fabs(previousStateVariables[i]) + eps;

        tardigradeHydra::ModelConfigurationBase model_configurationp(previousStateVariables + delta, parameters, numConfigurations, numNonLinearSolveStateVariables);

        tardigradeHydra::ModelConfigurationBase model_configurationm(previousStateVariables - delta, parameters, numConfigurations, numNonLinearSolveStateVariables);

        tardigradeHydra::hydraBase hydra_p(dof, model_configurationp);

        tardigradeHydra::hydraBase hydra_m(dof, model_configurationm);

        hydra_p.initialize();

        hydra_m.initialize();

        tardigradeHydra::floatVector F1_p, F1_m;

        F1_p = tardigradeHydra::floatVector(hydra_p.deformation->get_configurations()->begin(),
                                            hydra_p.deformation->get_configurations()->begin() + 9);

        F1_m = tardigradeHydra::floatVector(hydra_m.deformation->get_configurations()->begin(),
                                            hydra_m.deformation->get_configurations()->begin() + 9);

        for (unsigned int j = 0; j < F1_p.size(); j++) {
            dF1dFn_answer[j][i] = (F1_p[j] - F1_m[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(dF1dFn_answer) == *hydra.deformation->get_dF1dFn(),
               CHECK_PER_ELEMENT);

    tardigradeHydra::unit_test::DeformationBaseTester::checkdF1dFn(hydra, *hydra.deformation);
}
