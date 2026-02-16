/**
 * \file test_tardigrade_hydraDOFVelocityGradientDeformation.cpp
 *
 * Tests for tardigrade_hydraDOFVelocityGradientDeformation
 */

#include <tardigrade_constitutive_tools.h>
#include <tardigrade_hydraDOFVelocityGradientDeformation.h>
#include <tardigrade_hydraLinearElasticity.h>

#define BOOST_TEST_MODULE test_tardigrade_hydraDOFVelocityGradientDeformation
#include <boost/test/included/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>

#define DEFAULT_TEST_TOLERANCE 1e-6
#define CHECK_PER_ELEMENT boost::test_tools::per_element()

typedef tardigradeErrorTools::Node errorNode;  //!< Redefinition for the error node
typedef errorNode                 *errorOut;   //!< Redefinition for a pointer to the error node
typedef tardigradeHydra::dofVelocityGradientDeformation::floatType
    floatType;  //!< Redefinition of the floating point type
typedef tardigradeHydra::dofVelocityGradientDeformation::floatVector
    floatVector;  //!< Redefinition of the vector of floating points type
typedef tardigradeHydra::dofVelocityGradientDeformation::floatMatrix
    floatMatrix;  //!< Redefinition of the matrix of floating points type

namespace tardigradeHydra {

    namespace unit_test {

        class hydraBaseTester {
           public:
            static void updateUnknownVector(tardigradeHydra::hydraBase &hydra, const floatVector &value) {
                BOOST_CHECK_NO_THROW(hydra.updateUnknownVector(value));
            }
        };

    }  // namespace unit_test

    namespace dofVelocityGradientDeformation {

        namespace unit_test {

            class residualTester {
               public:
                static void runBasicGetTests(tardigradeHydra::dofVelocityGradientDeformation::residual &R) {
                    BOOST_CHECK(R._dofConfigurationIndex == R.getDOFConfigurationIndex());

                    BOOST_CHECK(R._dofVelocityGradientIndex == R.getDOFVelocityGradientIndex());

                    BOOST_CHECK(R._integrationParameter == R.getIntegrationParameter());

                    BOOST_CHECK(&R._dofVelocityGradient.second == R.get_dofVelocityGradient());

                    BOOST_CHECK(&R._previousDOFVelocityGradient.second == R.get_previousDOFVelocityGradient());
                }
            };

        }  // namespace unit_test

    }  // namespace dofVelocityGradientDeformation

}  // namespace tardigradeHydra

void adaptive_tolerance_test(const floatVector &result, const floatVector &answer, floatType atol = 1e-9) {
    /*!
     * Test which allows for very small values to be compared absolutely while larger values are compared relatively
     */

    BOOST_TEST(result.size() == answer.size());
    for (unsigned int i = 0; i < result.size(); i++) {
        if (std::abs(answer[i]) < atol) {
            BOOST_REQUIRE_SMALL(std::abs(answer[i] - result[i]), atol);
        } else {
            BOOST_TEST(answer[i] == result[i]);
        }
    }
}

BOOST_AUTO_TEST_CASE(test_residual_basicGetTests, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            dofDeformation = residualMock(this, 11, 1, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &dofDeformation;

            residuals[2] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousStateVariables = {1, 2, 3, 4, 5, 6, 7, 8, 9, -1, -2, -3, -4, -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 5, 6, 7};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 5;

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    residualMock R(&hydra, 11, 1, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    tardigradeHydra::dofVelocityGradientDeformation::unit_test::residualTester::runBasicGetTests(R);

    floatVector dofVelocityGradientAnswer = {0.44, 0.55, 0.66, 0.77, 0.88, 0.99, 1.11, 1.22, 1.33};

    floatVector previousDOFVelocityGradientAnswer = {-0.44, -0.55, -0.66, -0.77, -0.88, -0.99, -1.11, -1.22, -1.33};

    floatType densityAnswer = 0.22;

    floatType internalEnergyAnswer = 1.44;

    BOOST_TEST(dofVelocityGradientAnswer == *R.get_dofVelocityGradient(), CHECK_PER_ELEMENT);

    BOOST_TEST(R.getDensityIndex() == 1);

    BOOST_TEST(densityAnswer == *R.get_density());

    BOOST_TEST(R.getInternalEnergyIndex() == 12);

    BOOST_TEST(R.getInternalEnergyScaledByDensity());

    BOOST_TEST(internalEnergyAnswer == *R.get_internalEnergy());

    BOOST_TEST(previousDOFVelocityGradientAnswer == *R.get_previousDOFVelocityGradient(), CHECK_PER_ELEMENT);

    BOOST_TEST(hydra.dofDeformationParameters[0] == *R.get_massChangeRateFactor());

    BOOST_TEST(hydra.dofDeformationParameters[1] == *R.get_internalHeatGenerationRateFactor());

    BOOST_TEST(hydra.stateVariableIndices == *R.getStateVariableIndices(), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_residual_dofPrecedingDeformationGradient_1,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            dofDeformation = residualMock(this, 11, 1, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &dofDeformation;

            residuals[2] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1, -2, -3, -4, -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,   1,    1,    1,    1,   1, 1, 1, 1.2, 0.21, 0.29,
                                 0.23, 0.9, 0.11, 0.30, 0.25, 1.1, 3, 4, 5, 6,   7};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    floatVector answer = {0.99004464,  -0.04875021, -0.25613675, -0.21008905, 1.17745363,
                          -0.06235825, -0.20292692, -0.2263232,  0.98522215};

    floatVector previousAnswer = {1.02155172,  -0.06034483, -0.14224138, -0.14655172, 0.81034483,
                                  -0.23275862, -0.31465517, -0.31896552, 0.67672414};

    residualMock R(&hydra, 11, 1, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

    residualMock Rgrad(&hydra, 11, 1, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

    Rgrad.get_dPrecedingDeformationGradientdDeformationGradient();

    Rgrad.get_dPreviousPrecedingDeformationGradientdPreviousDeformationGradient();

    BOOST_TEST(answer == *R.get_precedingDeformationGradient(), CHECK_PER_ELEMENT);

    BOOST_TEST(previousAnswer == *R.get_previousPrecedingDeformationGradient(), CHECK_PER_ELEMENT);

    BOOST_TEST(answer == *Rgrad.get_precedingDeformationGradient(), CHECK_PER_ELEMENT);

    BOOST_TEST(previousAnswer == *Rgrad.get_previousPrecedingDeformationGradient(), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_residual_dofPrecedingDeformationGradient_2,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }

       protected:
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> plasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(4);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            plasticity = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 9);

            dofDeformation = residualMock(this, 11, 2, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &plasticity;

            residuals[2] = &dofDeformation;

            residuals[3] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3,
                                          0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1,  -2,  -3,  -4,  -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,    1,    1,    1,    1,    1,   1,   1,    1.1,  0.12,
                                 0.13, 0.14, 1.3,  0.15, 0.16, 0.17, 1.4, 1.2, 0.21, 0.29, 0.23,
                                 0.9,  0.11, 0.30, 0.25, 1.1,  3,    4,   5,   6,    7};

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    residualMock R(&hydra, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

    floatType eps = 1e-6;

    floatVector dpFdF(81, 0);

    floatVector dpFdFn(81 * 2, 0);

    floatVector previousdpFdF(81, 0);

    floatVector previousdpFdFn(81 * 2, 0);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(deformationGradient[i]) + eps;

        floatVector Fp = deformationGradient;

        floatVector Fm = deformationGradient;

        Fp[i] += delta;

        Fm[i] -= delta;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, Fp,
                                             previousDeformationGradient, additionalDOF, previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, Fm,
                                             previousDeformationGradient, additionalDOF, previousAdditionalDOF);

        hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        floatVector vp = *Rp.get_precedingDeformationGradient();

        floatVector vm = *Rm.get_precedingDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            dpFdF[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(dpFdF == *R.get_dPrecedingDeformationGradientdDeformationGradient(), CHECK_PER_ELEMENT);

    unsigned int offset = 9;

    for (unsigned int i = 0; i < 18; i++) {
        floatType delta = eps * std::fabs(unknownVector[i + offset]) + eps;

        floatVector xp = unknownVector;

        floatVector xm = unknownVector;

        xp[i + offset] += delta;

        xm[i + offset] -= delta;

        hydraBaseMock hydrap(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, xp);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, xm);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        floatVector vp = *Rp.get_precedingDeformationGradient();

        floatVector vm = *Rm.get_precedingDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            dpFdFn[18 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    adaptive_tolerance_test(*R.get_dPrecedingDeformationGradientdSubDeformationGradients(), dpFdFn);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(deformationGradient[i]) + eps;

        floatVector Fp = previousDeformationGradient;

        floatVector Fm = previousDeformationGradient;

        Fp[i] += delta;

        Fm[i] -= delta;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient, Fp,
                                             additionalDOF, previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient, Fm,
                                             additionalDOF, previousAdditionalDOF);

        hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        floatVector vp = *Rp.get_previousPrecedingDeformationGradient();

        floatVector vm = *Rm.get_previousPrecedingDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            previousdpFdF[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(previousdpFdF == *R.get_dPreviousPrecedingDeformationGradientdPreviousDeformationGradient(),
               CHECK_PER_ELEMENT);

    offset = 0;

    for (unsigned int i = 0; i < 18; i++) {
        floatType delta = eps * std::fabs(previousStateVariables[i + offset]) + eps;

        floatVector xp = previousStateVariables;

        floatVector xm = previousStateVariables;

        xp[i + offset] += delta;

        xm[i + offset] -= delta;

        tardigradeHydra::ModelConfigurationBase model_configurationp(xp, parameters);

        tardigradeHydra::ModelConfigurationBase model_configurationm(xm, parameters);

        hydraBaseMock hydrap(dof, model_configurationp, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dof, model_configurationm, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        floatVector vp = *Rp.get_previousPrecedingDeformationGradient();

        floatVector vm = *Rm.get_previousPrecedingDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            previousdpFdFn[18 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    adaptive_tolerance_test(*R.get_dPreviousPrecedingDeformationGradientdPreviousSubDeformationGradients(),
                            previousdpFdFn);
}

BOOST_AUTO_TEST_CASE(test_residual_dofIntermediateVelocityGradient_1,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

       protected:
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            dofDeformation = residualMock(this, 11, 1, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &dofDeformation;

            residuals[2] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00,
                                 7.00, 8.00, 9.00, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, 0.10,  0.20,  0.30,  0.40,  0.50,  0.60,
                                         0.70,  0.80,  0.90,  -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1, -2, -3, -4, -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,   1,    1,    1,    1,   1, 1, 1, 1.2, 0.21, 0.29,
                                 0.23, 0.9, 0.11, 0.30, 0.25, 1.1, 3, 4, 5, 6,   7};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    floatVector answer = {1.12104841, 4.2545126,  5.17734672, 1.85599459, 4.92506006,
                          5.28347558, 4.13186343, 9.15476705, 8.95389153};

    floatVector previousAnswer = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

    residualMock R(&hydra, 11, 1, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

    residualMock Rgrad(&hydra, 11, 1, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

    Rgrad.get_dDOFIntermediateVelocityGradientdDOFVelocityGradient();

    Rgrad.get_dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient();

    BOOST_TEST(answer == *R.get_dofIntermediateVelocityGradient(), CHECK_PER_ELEMENT);

    BOOST_TEST(previousAnswer == *R.get_previousDOFIntermediateVelocityGradient(), CHECK_PER_ELEMENT);

    BOOST_TEST(answer == *Rgrad.get_dofIntermediateVelocityGradient(), CHECK_PER_ELEMENT);

    BOOST_TEST(previousAnswer == *Rgrad.get_previousDOFIntermediateVelocityGradient(), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_residual_dofIntermediateVelocityGradient_2,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }

       protected:
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> plasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(4);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            plasticity = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 9);

            dofDeformation = residualMock(this, 11, 2, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &plasticity;

            residuals[2] = &dofDeformation;

            residuals[3] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3,
                                          0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1,  -2,  -3,  -4,  -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,    1,    1,    1,    1,    1,   1,   1,    1.1,  0.12,
                                 0.13, 0.14, 1.3,  0.15, 0.16, 0.17, 1.4, 1.2, 0.21, 0.29, 0.23,
                                 0.9,  0.11, 0.30, 0.25, 1.1,  3,    4,   5,   6,    7};

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    residualMock R(&hydra, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

    floatType eps = 1e-6;

    floatVector dILdF(81, 0);

    floatVector dILdFn(81 * 2, 0);

    floatVector dILdL(81, 0);

    floatVector previousdILdF(81, 0);

    floatVector previousdILdFn(81 * 2, 0);

    floatVector previousdILdL(81, 0);

    unsigned int offset = 0;

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(deformationGradient[i]) + eps;

        floatVector Fp = deformationGradient;

        floatVector Fm = deformationGradient;

        Fp[i] += delta;

        Fm[i] -= delta;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, Fp,
                                             previousDeformationGradient, additionalDOF, previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, Fm,
                                             previousDeformationGradient, additionalDOF, previousAdditionalDOF);

        hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        floatVector vp = *Rp.get_dofIntermediateVelocityGradient();

        floatVector vm = *Rm.get_dofIntermediateVelocityGradient();

        for (unsigned int j = 0; j < 9; j++) {
            dILdF[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(dILdF == *R.get_dDOFIntermediateVelocityGradientdDeformationGradient(), CHECK_PER_ELEMENT);

    offset = 9;

    for (unsigned int i = 0; i < 18; i++) {
        floatType delta = eps * std::fabs(unknownVector[i + offset]) + eps;

        floatVector xp = unknownVector;

        floatVector xm = unknownVector;

        xp[i + offset] += delta;

        xm[i + offset] -= delta;

        hydraBaseMock hydrap(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, xp);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, xm);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        floatVector vp = *Rp.get_dofIntermediateVelocityGradient();

        floatVector vm = *Rm.get_dofIntermediateVelocityGradient();

        for (unsigned int j = 0; j < 9; j++) {
            dILdFn[18 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    adaptive_tolerance_test(*R.get_dDOFIntermediateVelocityGradientdSubDeformationGradients(), dILdFn);

    offset = 3;

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(additionalDOF[i + offset]) + eps;

        floatVector xp = additionalDOF;

        floatVector xm = additionalDOF;

        xp[i + offset] += delta;

        xm[i + offset] -= delta;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient, xp, previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient, xm, previousAdditionalDOF);

        hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        floatVector vp = *Rp.get_dofIntermediateVelocityGradient();

        floatVector vm = *Rm.get_dofIntermediateVelocityGradient();

        for (unsigned int j = 0; j < 9; j++) {
            dILdL[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(*R.get_dDOFIntermediateVelocityGradientdDOFVelocityGradient() == dILdL, CHECK_PER_ELEMENT);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(deformationGradient[i]) + eps;

        floatVector Fp = previousDeformationGradient;

        floatVector Fm = previousDeformationGradient;

        Fp[i] += delta;

        Fm[i] -= delta;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient, Fp,
                                             additionalDOF, previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient, Fm,
                                             additionalDOF, previousAdditionalDOF);

        hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        floatVector vp = *Rp.get_previousDOFIntermediateVelocityGradient();

        floatVector vm = *Rm.get_previousDOFIntermediateVelocityGradient();

        for (unsigned int j = 0; j < 9; j++) {
            previousdILdF[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(previousdILdF == *R.get_dPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient(),
               CHECK_PER_ELEMENT);

    offset = 0;

    for (unsigned int i = 0; i < 18; i++) {
        floatType delta = eps * std::fabs(previousStateVariables[i + offset]) + eps;

        floatVector xp = previousStateVariables;

        floatVector xm = previousStateVariables;

        xp[i + offset] += delta;

        xm[i + offset] -= delta;

        tardigradeHydra::ModelConfigurationBase model_configurationp(xp, parameters);

        tardigradeHydra::ModelConfigurationBase model_configurationm(xm, parameters);

        hydraBaseMock hydrap(dof, model_configurationp, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dof, model_configurationm, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        floatVector vp = *Rp.get_previousDOFIntermediateVelocityGradient();

        floatVector vm = *Rm.get_previousDOFIntermediateVelocityGradient();

        for (unsigned int j = 0; j < 9; j++) {
            previousdILdFn[18 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    adaptive_tolerance_test(*R.get_dPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients(),
                            previousdILdFn);

    offset = 3;

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(previousAdditionalDOF[i + offset]) + eps;

        floatVector xp = previousAdditionalDOF;

        floatVector xm = previousAdditionalDOF;

        xp[i + offset] += delta;

        xm[i + offset] -= delta;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient, additionalDOF, xp);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient, additionalDOF, xm);

        hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

        floatVector vp = *Rp.get_previousDOFIntermediateVelocityGradient();

        floatVector vm = *Rm.get_previousDOFIntermediateVelocityGradient();

        for (unsigned int j = 0; j < 9; j++) {
            previousdILdL[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(*R.get_dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient() == previousdILdL,
               CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_residual_dofDeformationGradient_1, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        floatVector dofIntermediateVelocityGradient = {1, 2, 3, 4, 5, 6, 7, 8, 9};

        floatVector previousDOFIntermediateVelocityGradient = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

        floatVector dDOFIntermediateVelocityGradientdDOFVelocityGradient = initializeVector(81);

        floatVector dDOFIntermediateVelocityGradientdDeformationGradient = initializeVector(81);

        floatVector dDOFIntermediateVelocityGradientdSubDeformationGradients = initializeVector(81);

        floatVector dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient = initializeVector(81);

        floatVector dPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient = initializeVector(81);

        floatVector dPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients = initializeVector(81);

        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }

       protected:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::setDOFIntermediateVelocityGradient;

        virtual void setDOFIntermediateVelocityGradient(const bool &isPrevious) override {
            if (isPrevious) {
                set_previousDOFIntermediateVelocityGradient(previousDOFIntermediateVelocityGradient);

            } else {
                set_dofIntermediateVelocityGradient(dofIntermediateVelocityGradient);
            }
        }

        virtual void setDOFIntermediateVelocityGradientDerivatives(const bool &isPrevious) override {
            if (isPrevious) {
                set_dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient(
                    dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient);

                set_dPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient(
                    dPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient);

                set_dPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients(
                    dPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients);

            } else {
                set_dDOFIntermediateVelocityGradientdDOFVelocityGradient(
                    dDOFIntermediateVelocityGradientdDOFVelocityGradient);

                set_dDOFIntermediateVelocityGradientdDeformationGradient(
                    dDOFIntermediateVelocityGradientdDeformationGradient);

                set_dDOFIntermediateVelocityGradientdSubDeformationGradients(
                    dDOFIntermediateVelocityGradientdSubDeformationGradients);
            }
        }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            dofDeformation = residualMock(this, 11, 1, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &dofDeformation;

            residuals[2] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1, -2, -3, -4, -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,   1,    1,    1,    1,   1, 1, 1, 1.2, 0.21, 0.29,
                                 0.23, 0.9, 0.11, 0.30, 0.25, 1.1, 3, 4, 5, 6,   7};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    floatVector answer = {1.95711900e+10, 2.40473864e+10, 2.85235828e+10, 4.43210244e+10, 5.44578434e+10,
                          6.45946624e+10, 6.90708588e+10, 8.48683004e+10, 1.00665742e+11};

    residualMock R(&hydra, 11, 1, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters, 0.67);

    residualMock Rgrad(&hydra, 11, 1, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters, 0.67);

    Rgrad.get_dDOFDeformationGradientdDOFVelocityGradient();

    Rgrad.get_dDOFDeformationGradientdPreviousDOFVelocityGradient();

    BOOST_TEST(answer == *R.get_dofDeformationGradient(), CHECK_PER_ELEMENT);

    BOOST_TEST(answer == *Rgrad.get_dofDeformationGradient(), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_residual_dofDeformationGradient_2, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> plasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(4);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            plasticity = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 9);

            dofDeformation = residualMock(this, 11, 2, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &plasticity;

            residuals[2] = &dofDeformation;

            residuals[3] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3,
                                          0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1,  -2,  -3,  -4,  -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,    1,    1,    1,    1,    1,   1,   1,    1.1,  0.12,
                                 0.13, 0.14, 1.3,  0.15, 0.16, 0.17, 1.4, 1.2, 0.21, 0.29, 0.23,
                                 0.9,  0.11, 0.30, 0.25, 1.1,  3,    4,   5,   6,    7};

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    floatType alpha = 0.67;

    residualMock R(&hydra, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters, alpha);

    floatVector dFmdL(81, 0);

    floatVector dFmdF(81, 0);

    floatVector dFmdFn(81 * 2, 0);

    floatVector previousdFmdL(81, 0);

    floatVector previousdFmdF(81, 0);

    floatVector previousdFmdFn(81 * 2, 0);

    floatType eps = 1e-6;

    unsigned int offset = 0;

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(deformationGradient[i]) + eps;

        floatVector Fp = deformationGradient;

        floatVector Fm = deformationGradient;

        Fp[i] += delta;

        Fm[i] -= delta;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, Fp,
                                             previousDeformationGradient, additionalDOF, previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, Fm,
                                             previousDeformationGradient, additionalDOF, previousAdditionalDOF);

        hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        floatVector vp = *Rp.get_dofDeformationGradient();

        floatVector vm = *Rm.get_dofDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            dFmdF[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(dFmdF == *R.get_dDOFDeformationGradientdDeformationGradient(), CHECK_PER_ELEMENT);

    offset = 9;

    for (unsigned int i = 0; i < 18; i++) {
        floatType delta = eps * std::fabs(unknownVector[i + offset]) + eps;

        floatVector xp = unknownVector;

        floatVector xm = unknownVector;

        xp[i + offset] += delta;

        xm[i + offset] -= delta;

        hydraBaseMock hydrap(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, xp);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, xm);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        floatVector vp = *Rp.get_dofDeformationGradient();

        floatVector vm = *Rm.get_dofDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            dFmdFn[18 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    adaptive_tolerance_test(*R.get_dDOFDeformationGradientdSubDeformationGradients(), dFmdFn, 1e-8);

    offset = 3;

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(additionalDOF[i + offset]) + eps;

        floatVector xp = additionalDOF;

        floatVector xm = additionalDOF;

        xp[i + offset] += delta;

        xm[i + offset] -= delta;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient, xp, previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient, xm, previousAdditionalDOF);

        hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        floatVector vp = *Rp.get_dofDeformationGradient();

        floatVector vm = *Rm.get_dofDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            dFmdL[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(*R.get_dDOFDeformationGradientdDOFVelocityGradient() == dFmdL, CHECK_PER_ELEMENT);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(deformationGradient[i]) + eps;

        floatVector Fp = previousDeformationGradient;

        floatVector Fm = previousDeformationGradient;

        Fp[i] += delta;

        Fm[i] -= delta;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient, Fp,
                                             additionalDOF, previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient, Fm,
                                             additionalDOF, previousAdditionalDOF);

        hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        floatVector vp = *Rp.get_dofDeformationGradient();

        floatVector vm = *Rm.get_dofDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            previousdFmdF[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(previousdFmdF == *R.get_dDOFDeformationGradientdPreviousDeformationGradient(), CHECK_PER_ELEMENT);

    offset = 0;

    for (unsigned int i = 0; i < 18; i++) {
        floatType delta = eps * std::fabs(previousStateVariables[i + offset]) + eps;

        floatVector xp = previousStateVariables;

        floatVector xm = previousStateVariables;

        xp[i + offset] += delta;

        xm[i + offset] -= delta;

        tardigradeHydra::ModelConfigurationBase model_configurationp(xp, parameters);

        tardigradeHydra::ModelConfigurationBase model_configurationm(xm, parameters);

        hydraBaseMock hydrap(dof, model_configurationp, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dof, model_configurationm, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        floatVector vp = *Rp.get_dofDeformationGradient();

        floatVector vm = *Rm.get_dofDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            previousdFmdFn[18 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    adaptive_tolerance_test(*R.get_dDOFDeformationGradientdPreviousSubDeformationGradients(), previousdFmdFn, 1e-7);

    offset = 3;

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(previousAdditionalDOF[i + offset]) + eps;

        floatVector xp = previousAdditionalDOF;

        floatVector xm = previousAdditionalDOF;

        xp[i + offset] += delta;

        xm[i + offset] -= delta;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient, additionalDOF, xp);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient, additionalDOF, xm);

        hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        floatVector vp = *Rp.get_dofDeformationGradient();

        floatVector vm = *Rm.get_dofDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            previousdFmdL[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(*R.get_dDOFDeformationGradientdPreviousDOFVelocityGradient() == previousdFmdL, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_residual_dofDeformationGradient_3, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        floatVector dofIntermediateVelocityGradient = {1, 2, 3, 4, 5, 6, 7, 8, 9};

        floatVector previousDOFIntermediateVelocityGradient = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

        floatVector dDOFIntermediateVelocityGradientdDOFVelocityGradient = initializeVector(81);

        floatVector dDOFIntermediateVelocityGradientdDeformationGradient = initializeVector(81);

        floatVector dDOFIntermediateVelocityGradientdSubDeformationGradients = initializeVector(81);

        floatVector dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient = initializeVector(81);

        floatVector dPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient = initializeVector(81);

        floatVector dPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients = initializeVector(81);

        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }

       protected:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::setDOFIntermediateVelocityGradient;

        virtual void setDOFIntermediateVelocityGradient(const bool &isPrevious) override {
            if (isPrevious) {
                set_previousDOFIntermediateVelocityGradient(previousDOFIntermediateVelocityGradient);

            } else {
                set_dofIntermediateVelocityGradient(dofIntermediateVelocityGradient);
            }
        }

        virtual void setDOFIntermediateVelocityGradientDerivatives(const bool &isPrevious) override {
            if (isPrevious) {
                set_dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient(
                    dPreviousDOFIntermediateVelocityGradientdPreviousDOFVelocityGradient);

                set_dPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient(
                    dPreviousDOFIntermediateVelocityGradientdPreviousDeformationGradient);

                set_dPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients(
                    dPreviousDOFIntermediateVelocityGradientdPreviousSubDeformationGradients);

            } else {
                set_dDOFIntermediateVelocityGradientdDOFVelocityGradient(
                    dDOFIntermediateVelocityGradientdDOFVelocityGradient);

                set_dDOFIntermediateVelocityGradientdDeformationGradient(
                    dDOFIntermediateVelocityGradientdDeformationGradient);

                set_dDOFIntermediateVelocityGradientdSubDeformationGradients(
                    dDOFIntermediateVelocityGradientdSubDeformationGradients);
            }
        }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            dofDeformation = residualMock(this, 11, 1, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);
            dofDeformation.setUseTrapezoidalIntegration(true);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &dofDeformation;

            residuals[2] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1, -2, -3, -4, -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,   1,    1,    1,    1,   1, 1, 1, 1.2, 0.21, 0.29,
                                 0.23, 0.9, 0.11, 0.30, 0.25, 1.1, 3, 4, 5, 6,   7};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    floatVector answer = {3.612667141e-01,  -3.072550304e-01, 2.422322512e-02,  -3.719571754e-01, 5.956917520e-01,
                          -4.366593205e-01, -1.051810650e-01, -5.013614655e-01, 1.024581339e-01};

    residualMock R(&hydra, 11, 1, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters, 0.67);

    residualMock Rgrad(&hydra, 11, 1, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters, 0.67);

    R.setUseTrapezoidalIntegration(true);
    Rgrad.setUseTrapezoidalIntegration(true);

    Rgrad.get_dDOFDeformationGradientdDOFVelocityGradient();

    Rgrad.get_dDOFDeformationGradientdPreviousDOFVelocityGradient();

    BOOST_TEST(answer == *R.get_dofDeformationGradient(), CHECK_PER_ELEMENT);

    BOOST_TEST(answer == *Rgrad.get_dofDeformationGradient(), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_residual_dofDeformationGradient_4, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> plasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(4);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            plasticity = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 9);

            dofDeformation = residualMock(this, 11, 2, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &plasticity;

            residuals[2] = &dofDeformation;
            dofDeformation.setUseTrapezoidalIntegration(true);

            residuals[3] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3,
                                          0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1,  -2,  -3,  -4,  -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,    1,    1,    1,    1,    1,   1,   1,    1.1,  0.12,
                                 0.13, 0.14, 1.3,  0.15, 0.16, 0.17, 1.4, 1.2, 0.21, 0.29, 0.23,
                                 0.9,  0.11, 0.30, 0.25, 1.1,  3,    4,   5,   6,    7};

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    floatType alpha = 0.67;

    residualMock R(&hydra, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters, alpha);
    R.setUseTrapezoidalIntegration(true);

    floatVector dFmdL(81, 0);

    floatVector dFmdF(81, 0);

    floatVector dFmdFn(81 * 2, 0);

    floatVector previousdFmdL(81, 0);

    floatVector previousdFmdF(81, 0);

    floatVector previousdFmdFn(81 * 2, 0);

    floatType eps = 1e-6;

    unsigned int offset = 0;

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(deformationGradient[i]) + eps;

        floatVector Fp = deformationGradient;

        floatVector Fm = deformationGradient;

        Fp[i] += delta;

        Fm[i] -= delta;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, Fp,
                                             previousDeformationGradient, additionalDOF, previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, Fm,
                                             previousDeformationGradient, additionalDOF, previousAdditionalDOF);

        hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        Rp.setUseTrapezoidalIntegration(true);

        Rm.setUseTrapezoidalIntegration(true);

        floatVector vp = *Rp.get_dofDeformationGradient();

        floatVector vm = *Rm.get_dofDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            dFmdF[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(dFmdF == *R.get_dDOFDeformationGradientdDeformationGradient(), CHECK_PER_ELEMENT);

    offset = 9;

    for (unsigned int i = 0; i < 18; i++) {
        floatType delta = eps * std::fabs(unknownVector[i + offset]) + eps;

        floatVector xp = unknownVector;

        floatVector xm = unknownVector;

        xp[i + offset] += delta;

        xm[i + offset] -= delta;

        hydraBaseMock hydrap(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, xp);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, xm);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        Rp.setUseTrapezoidalIntegration(true);

        Rm.setUseTrapezoidalIntegration(true);

        floatVector vp = *Rp.get_dofDeformationGradient();

        floatVector vm = *Rm.get_dofDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            dFmdFn[18 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    adaptive_tolerance_test(*R.get_dDOFDeformationGradientdSubDeformationGradients(), dFmdFn, 1e-8);

    offset = 3;

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(additionalDOF[i + offset]) + eps;

        floatVector xp = additionalDOF;

        floatVector xm = additionalDOF;

        xp[i + offset] += delta;

        xm[i + offset] -= delta;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient, xp, previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient, xm, previousAdditionalDOF);

        hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        Rp.setUseTrapezoidalIntegration(true);

        Rm.setUseTrapezoidalIntegration(true);

        floatVector vp = *Rp.get_dofDeformationGradient();

        floatVector vm = *Rm.get_dofDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            dFmdL[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(*R.get_dDOFDeformationGradientdDOFVelocityGradient() == dFmdL, CHECK_PER_ELEMENT);

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(deformationGradient[i]) + eps;

        floatVector Fp = previousDeformationGradient;

        floatVector Fm = previousDeformationGradient;

        Fp[i] += delta;

        Fm[i] -= delta;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient, Fp,
                                             additionalDOF, previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient, Fm,
                                             additionalDOF, previousAdditionalDOF);

        hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        Rp.setUseTrapezoidalIntegration(true);

        Rm.setUseTrapezoidalIntegration(true);

        floatVector vp = *Rp.get_dofDeformationGradient();

        floatVector vm = *Rm.get_dofDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            previousdFmdF[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    BOOST_TEST(previousdFmdF == *R.get_dDOFDeformationGradientdPreviousDeformationGradient(), CHECK_PER_ELEMENT);

    offset = 0;

    for (unsigned int i = 0; i < 18; i++) {
        floatType delta = eps * std::fabs(previousStateVariables[i + offset]) + eps;

        floatVector xp = previousStateVariables;

        floatVector xm = previousStateVariables;

        xp[i + offset] += delta;

        xm[i + offset] -= delta;

        tardigradeHydra::ModelConfigurationBase model_configurationp(xp, parameters);

        tardigradeHydra::ModelConfigurationBase model_configurationm(xm, parameters);

        hydraBaseMock hydrap(dof, model_configurationp, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dof, model_configurationm, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        Rp.setUseTrapezoidalIntegration(true);

        Rm.setUseTrapezoidalIntegration(true);

        floatVector vp = *Rp.get_dofDeformationGradient();

        floatVector vm = *Rm.get_dofDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            previousdFmdFn[18 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    adaptive_tolerance_test(*R.get_dDOFDeformationGradientdPreviousSubDeformationGradients(), previousdFmdFn, 1e-7);

    offset = 3;

    for (unsigned int i = 0; i < 9; i++) {
        floatType delta = eps * std::fabs(previousAdditionalDOF[i + offset]) + eps;

        floatVector xp = previousAdditionalDOF;

        floatVector xm = previousAdditionalDOF;

        xp[i + offset] += delta;

        xm[i + offset] -= delta;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient, additionalDOF, xp);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient, additionalDOF, xm);

        hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

        tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

        residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                        alpha);

        Rp.setUseTrapezoidalIntegration(true);

        Rm.setUseTrapezoidalIntegration(true);

        floatVector vp = *Rp.get_dofDeformationGradient();

        floatVector vm = *Rm.get_dofDeformationGradient();

        for (unsigned int j = 0; j < 9; j++) {
            previousdFmdL[9 * j + i] = (vp[j] - vm[j]) / (2 * delta);
        }
    }

    adaptive_tolerance_test(*R.get_dDOFDeformationGradientdPreviousDOFVelocityGradient(), previousdFmdL, 1e-7);
}

BOOST_AUTO_TEST_CASE(test_residual_massChangeRate_1, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }

       protected:
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            dofDeformation = residualMock(this, 11, 1, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &dofDeformation;

            residuals[2] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1, -2, -3, -4, -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,   1,    1,    1,    1,   1, 1, 1, 1.2, 0.21, 0.29,
                                 0.23, 0.9, 0.11, 0.30, 0.25, 1.1, 3, 4, 5, 6,   7};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    floatType answer = 0.78 * 0.22 * (0.44 + 0.88 + 1.33);

    residualMock R(&hydra, 11, 1, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters, 0.67);

    residualMock Rgrad(&hydra, 11, 1, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters, 0.67);

    Rgrad.get_dMassChangeRatedDensity();

    Rgrad.get_dMassChangeRatedDOFVelocityGradient();

    BOOST_TEST(answer == *R.get_massChangeRate());

    BOOST_TEST(answer == *Rgrad.get_massChangeRate());
}

BOOST_AUTO_TEST_CASE(test_residual_massChangeRate_2, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> plasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(4);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            plasticity = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 9);

            dofDeformation = residualMock(this, 11, 2, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &plasticity;

            residuals[2] = &dofDeformation;

            residuals[3] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3,
                                          0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1,  -2,  -3,  -4,  -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,    1,    1,    1,    1,    1,   1,   1,    1.1,  0.12,
                                 0.13, 0.14, 1.3,  0.15, 0.16, 0.17, 1.4, 1.2, 0.21, 0.29, 0.23,
                                 0.9,  0.11, 0.30, 0.25, 1.1,  3,    4,   5,   6,    7};

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    floatType alpha = 0.67;

    residualMock R(&hydra, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters, alpha);

    floatVector dCdAdditionalDOF(18, 0);

    floatType eps = 1e-6;

    {
        constexpr unsigned int OUT_SIZE = 1;

        constexpr unsigned int VAR_SIZE = 18;

        std::vector<double> X(std::begin(additionalDOF), std::end(additionalDOF));

        for (unsigned int i = 0; i < VAR_SIZE; i++) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;

            floatVector xm = X;

            xp[i] += delta;

            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, xp, previousAdditionalDOF);

            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, xm, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydrap.initialize();

            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

            residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                            alpha);

            residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                            alpha);

            floatVector vp = {*Rp.get_massChangeRate()};

            floatVector vm = {*Rm.get_massChangeRate()};

            for (unsigned int j = 0; j < OUT_SIZE; j++) {
                dCdAdditionalDOF[VAR_SIZE * j + i] = (vp[j] - vm[j]) / (2 * delta);
            }
        }
    }

    floatVector answer(18, 0);
    answer[18 * 0 + R.getDensityIndex()] = *R.get_dMassChangeRatedDensity();
    std::copy(std::begin(*R.get_dMassChangeRatedDOFVelocityGradient()),
              std::end(*R.get_dMassChangeRatedDOFVelocityGradient()),
              std::begin(answer) + R.getDOFVelocityGradientIndex());

    BOOST_TEST(dCdAdditionalDOF == answer, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_residual_internalHeatGenerationRate_1, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }

       protected:
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            dofDeformation = residualMock(this, 11, 1, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &dofDeformation;

            residuals[2] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1, -2, -3, -4, -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,   1,    1,    1,    1,   1, 1, 1, 1.2, 0.21, 0.29,
                                 0.23, 0.9, 0.11, 0.30, 0.25, 1.1, 3, 4, 5, 6,   7};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    floatType answer = 0.89 * 1.44 / 0.22 * (0.44 + 0.88 + 1.33);

    residualMock R(&hydra, 11, 1, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters, 0.67);

    residualMock Rgrad(&hydra, 11, 1, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters, 0.67);

    Rgrad.get_dInternalHeatGenerationRatedInternalEnergy();

    Rgrad.get_dInternalHeatGenerationRatedDensity();

    Rgrad.get_dInternalHeatGenerationRatedDOFVelocityGradient();

    BOOST_TEST(answer == *R.get_internalHeatGenerationRate());

    BOOST_TEST(answer == *Rgrad.get_internalHeatGenerationRate());
}

BOOST_AUTO_TEST_CASE(test_residual_internalHeatGenerationRate_2, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> plasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(4);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            plasticity = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 9);

            dofDeformation = residualMock(this, 11, 2, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &plasticity;

            residuals[2] = &dofDeformation;

            residuals[3] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3,
                                          0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1,  -2,  -3,  -4,  -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,    1,    1,    1,    1,    1,   1,   1,    1.1,  0.12,
                                 0.13, 0.14, 1.3,  0.15, 0.16, 0.17, 1.4, 1.2, 0.21, 0.29, 0.23,
                                 0.9,  0.11, 0.30, 0.25, 1.1,  3,    4,   5,   6,    7};

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    floatType alpha = 0.67;

    residualMock R(&hydra, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters, alpha);

    floatVector drdAdditionalDOF(18, 0);

    floatType eps = 1e-6;

    {
        constexpr unsigned int OUT_SIZE = 1;

        constexpr unsigned int VAR_SIZE = 18;

        std::vector<double> X(std::begin(additionalDOF), std::end(additionalDOF));

        for (unsigned int i = 0; i < VAR_SIZE; i++) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;

            floatVector xm = X;

            xp[i] += delta;

            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, xp, previousAdditionalDOF);

            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, xm, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydrap.initialize();

            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

            residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                            alpha);

            residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                            alpha);

            floatVector vp = {*Rp.get_internalHeatGenerationRate()};

            floatVector vm = {*Rm.get_internalHeatGenerationRate()};

            for (unsigned int j = 0; j < OUT_SIZE; j++) {
                drdAdditionalDOF[VAR_SIZE * j + i] = (vp[j] - vm[j]) / (2 * delta);
            }
        }
    }

    floatVector answer(18, 0);
    answer[18 * 0 + R.getDensityIndex()] = *R.get_dInternalHeatGenerationRatedDensity();
    std::copy(std::begin(*R.get_dInternalHeatGenerationRatedDOFVelocityGradient()),
              std::end(*R.get_dInternalHeatGenerationRatedDOFVelocityGradient()),
              std::begin(answer) + R.getDOFVelocityGradientIndex());
    answer[18 * 0 + R.getInternalEnergyIndex()] = *R.get_dInternalHeatGenerationRatedInternalEnergy();

    BOOST_TEST(drdAdditionalDOF == answer, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_residual_internalHeatGenerationRate_3, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }

       protected:
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            dofDeformation = residualMock(this, 11, 1, 1, 12, 3, false, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &dofDeformation;

            residuals[2] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1, -2, -3, -4, -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,   1,    1,    1,    1,   1, 1, 1, 1.2, 0.21, 0.29,
                                 0.23, 0.9, 0.11, 0.30, 0.25, 1.1, 3, 4, 5, 6,   7};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    floatType answer = 0.89 * 1.44 * (0.44 + 0.88 + 1.33);

    residualMock R(&hydra, 11, 1, 1, 12, 3, false, hydra.stateVariableIndices, hydra.dofDeformationParameters, 0.67);

    residualMock Rgrad(&hydra, 11, 1, 1, 12, 3, false, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                       0.67);

    Rgrad.get_dInternalHeatGenerationRatedInternalEnergy();

    Rgrad.get_dInternalHeatGenerationRatedDensity();

    Rgrad.get_dInternalHeatGenerationRatedDOFVelocityGradient();

    BOOST_TEST(answer == *R.get_internalHeatGenerationRate());

    BOOST_TEST(answer == *Rgrad.get_internalHeatGenerationRate());
}

BOOST_AUTO_TEST_CASE(test_residual_internalHeatGenerationRate_4, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> plasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(4);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            plasticity = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 9);

            dofDeformation = residualMock(this, 11, 2, 1, 12, 3, false, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &plasticity;

            residuals[2] = &dofDeformation;

            residuals[3] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3,
                                          0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1,  -2,  -3,  -4,  -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,    1,    1,    1,    1,    1,   1,   1,    1.1,  0.12,
                                 0.13, 0.14, 1.3,  0.15, 0.16, 0.17, 1.4, 1.2, 0.21, 0.29, 0.23,
                                 0.9,  0.11, 0.30, 0.25, 1.1,  3,    4,   5,   6,    7};

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    floatType alpha = 0.67;

    residualMock R(&hydra, 11, 2, 1, 12, 3, false, hydra.stateVariableIndices, hydra.dofDeformationParameters, alpha);

    floatVector drdAdditionalDOF(18, 0);

    floatType eps = 1e-6;

    {
        constexpr unsigned int OUT_SIZE = 1;

        constexpr unsigned int VAR_SIZE = 18;

        std::vector<double> X(std::begin(additionalDOF), std::end(additionalDOF));

        for (unsigned int i = 0; i < VAR_SIZE; i++) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;

            floatVector xm = X;

            xp[i] += delta;

            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, xp, previousAdditionalDOF);

            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, xm, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydrap.initialize();

            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

            residualMock Rp(&hydrap, 11, 2, 1, 12, 3, false, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                            alpha);

            residualMock Rm(&hydram, 11, 2, 1, 12, 3, false, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                            alpha);

            floatVector vp = {*Rp.get_internalHeatGenerationRate()};

            floatVector vm = {*Rm.get_internalHeatGenerationRate()};

            for (unsigned int j = 0; j < OUT_SIZE; j++) {
                drdAdditionalDOF[VAR_SIZE * j + i] = (vp[j] - vm[j]) / (2 * delta);
            }
        }
    }

    floatVector answer(18, 0);
    answer[18 * 0 + R.getDensityIndex()] = *R.get_dInternalHeatGenerationRatedDensity();
    std::copy(std::begin(*R.get_dInternalHeatGenerationRatedDOFVelocityGradient()),
              std::end(*R.get_dInternalHeatGenerationRatedDOFVelocityGradient()),
              std::begin(answer) + R.getDOFVelocityGradientIndex());
    answer[18 * 0 + R.getInternalEnergyIndex()] = *R.get_dInternalHeatGenerationRatedInternalEnergy();

    BOOST_TEST(drdAdditionalDOF == answer, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_residual_dofResidual, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        floatVector dofDeformationGradient = {0.111, 0.222, 0.333, 0.444, 0.555, 0.666, 0.777, 0.888, 0.999};

        floatType massChangeRate = 1.23;

        floatType internalHeatGenerationRate = 2.34;

       protected:
        virtual void setDOFDeformationGradient() override { set_dofDeformationGradient(dofDeformationGradient); }

        virtual void setMassChangeRate() override {
            auto mcr = get_SetDataStorage_massChangeRate();

            *mcr.value = massChangeRate;
        }

        virtual void setInternalHeatGenerationRate() override {
            auto ihgr = get_SetDataStorage_internalHeatGenerationRate();

            *ihgr.value = internalHeatGenerationRate;
        }

        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> plasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(4);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            plasticity = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 9);

            dofDeformation = residualMock(this, 11, 2, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &plasticity;

            residuals[2] = &dofDeformation;

            residuals[3] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3,
                                          0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1,  -2,  -3,  -4,  -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,    1,    1,    1,    1,    1,   1,   1,    1.1,  0.12,
                                 0.13, 0.14, 1.3,  0.15, 0.16, 0.17, 1.4, 1.2, 0.21, 0.29, 0.23,
                                 0.9,  0.11, 0.30, 0.25, 1.1,  3,    4,   5,   6,    7};

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    floatType alpha = 0.67;

    residualMock R(&hydra, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters, alpha);

    floatVector answer = {-1.089, 0.012, 0.043, 0.214, -0.345, 0.556, 0.477, 0.638, -0.101, 1.23 - 6, 2.34 - 7};

    BOOST_TEST(answer == *R.getResidual(), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_residual_dofResidual_2, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

       protected:
        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> plasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(4);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            plasticity = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 9);

            dofDeformation = residualMock(this, 11, 2, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &plasticity;

            residuals[2] = &dofDeformation;

            residuals[3] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3,
                                          0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1,  -2,  -3,  -4,  -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,    1,    1,    1,    1,    1,   1,   1,    1.1,  0.12,
                                 0.13, 0.14, 1.3,  0.15, 0.16, 0.17, 1.4, 1.2, 0.21, 0.29, 0.23,
                                 0.9,  0.11, 0.30, 0.25, 1.1,  3,    4,   5,   6,    7};

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    floatType alpha = 0.67;

    residualMock R(&hydra, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters, alpha);

    floatVector jacobian(11 * unknownVector.size(), 0);

    floatVector dRdF(11 * 9, 0);

    floatVector dRdT(11, 0);

    floatVector dRdAdditionalDOF(11 * 18, 0);

    floatType eps = 1e-6;

    {
        constexpr unsigned int VAR_DIM = 32;

        constexpr unsigned int OUT_DIM = 11;

        floatVector X(std::begin(unknownVector), std::end(unknownVector));

        for (unsigned int i = 0; i < VAR_DIM; i++) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;

            floatVector xm = X;

            xp[i] += delta;

            xm[i] -= delta;

            hydraBaseMock hydrap(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydraBaseMock hydram(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydrap.initialize();

            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, xp);

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, xm);

            residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                            alpha);

            residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                            alpha);

            floatVector vp = *Rp.getResidual();

            floatVector vm = *Rm.getResidual();

            for (unsigned int j = 0; j < OUT_DIM; j++) {
                jacobian[VAR_DIM * j + i] += (vp[j] - vm[j]) / (2 * delta);
            }
        }
    }

    adaptive_tolerance_test(*R.getJacobian(), jacobian, 1e-8);

    {
        constexpr unsigned int VAR_DIM = 9;

        constexpr unsigned int OUT_DIM = 11;

        floatVector X(std::begin(deformationGradient), std::end(deformationGradient));

        for (unsigned int i = 0; i < VAR_DIM; i++) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;

            floatVector xm = X;

            xp[i] += delta;

            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, xp,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, xm,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydrap.initialize();

            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

            residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                            alpha);

            residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                            alpha);

            floatVector vp = *Rp.getResidual();

            floatVector vm = *Rm.getResidual();

            for (unsigned int j = 0; j < OUT_DIM; j++) {
                dRdF[VAR_DIM * j + i] += (vp[j] - vm[j]) / (2 * delta);
            }
        }
    }

    adaptive_tolerance_test(*R.getdRdF(), dRdF, 1e-8);

    {
        constexpr unsigned int VAR_DIM = 1;

        constexpr unsigned int OUT_DIM = 11;

        floatVector X = {temperature};

        for (unsigned int i = 0; i < VAR_DIM; i++) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;

            floatVector xm = X;

            xp[i] += delta;

            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, xp[0], previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, xm[0], previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydrap.initialize();

            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

            residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                            alpha);

            residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                            alpha);

            floatVector vp = *Rp.getResidual();

            floatVector vm = *Rm.getResidual();

            for (unsigned int j = 0; j < OUT_DIM; j++) {
                dRdT[VAR_DIM * j + i] += (vp[j] - vm[j]) / (2 * delta);
            }
        }
    }

    adaptive_tolerance_test(*R.getdRdT(), dRdT, 1e-8);

    eps = 1e-7;

    {
        constexpr unsigned int VAR_DIM = 18;

        constexpr unsigned int OUT_DIM = 11;

        floatVector X = additionalDOF;

        for (unsigned int i = 0; i < VAR_DIM; i++) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;

            floatVector xm = X;

            xp[i] += delta;

            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, xp, previousAdditionalDOF);

            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, xm, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydrap.initialize();

            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

            residualMock Rp(&hydrap, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                            alpha);

            residualMock Rm(&hydram, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters,
                            alpha);

            floatVector vp = *Rp.getResidual();

            floatVector vm = *Rm.getResidual();

            for (unsigned int j = 0; j < OUT_DIM; j++) {
                dRdAdditionalDOF[VAR_DIM * j + i] += (vp[j] - vm[j]) / (2 * delta);
            }
        }
    }

    adaptive_tolerance_test(*R.getdRdAdditionalDOF(), dRdAdditionalDOF, 1e-8);
}

BOOST_AUTO_TEST_CASE(test_residual_exampleModel, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {0, 1};

        tardigradeHydra::linearElasticity::residual elasticity;

        residualMock dofDeformation;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(2);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            dofDeformation = residualMock(this, 11, 1, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            residuals[0] = &elasticity;

            residuals[1] = &dofDeformation;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    floatVector parameters = {123.4, 56.7, 0.78, 0.89};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 2;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    floatVector dXdF(20 * 9, 0);

    floatVector dXdT(20, 0);

    floatVector dXdAdditionalDOF(20 * 18, 0);

    hydra.evaluate();

    floatType eps = 1e-6;

    {
        constexpr unsigned int VAR_SIZE = 9;

        constexpr unsigned int OUT_SIZE = 20;

        std::vector<double> X(std::begin(deformationGradient), std::end(deformationGradient));

        for (unsigned int i = 0; i < VAR_SIZE; i++) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;

            floatVector xm = X;

            xp[i] += delta;

            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, xp,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, xm,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydrap.evaluate();

            hydram.evaluate();

            floatVector vp = *hydrap.getUnknownVector();

            floatVector vm = *hydram.getUnknownVector();

            for (unsigned int j = 0; j < OUT_SIZE; j++) {
                dXdF[VAR_SIZE * j + i] = (vp[j] - vm[j]) / (2 * delta);
            }
        }
    }

    BOOST_TEST(*hydra.getFlatdXdF() == dXdF, CHECK_PER_ELEMENT);

    {
        constexpr unsigned int VAR_SIZE = 1;

        constexpr unsigned int OUT_SIZE = 20;

        std::vector<double> X = {temperature};

        for (unsigned int i = 0; i < VAR_SIZE; i++) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;

            floatVector xm = X;

            xp[i] += delta;

            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, xp[0], previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, xm[0], previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydrap.evaluate();

            hydram.evaluate();

            floatVector vp = *hydrap.getUnknownVector();

            floatVector vm = *hydram.getUnknownVector();

            for (unsigned int j = 0; j < OUT_SIZE; j++) {
                dXdT[VAR_SIZE * j + i] = (vp[j] - vm[j]) / (2 * delta);
            }
        }
    }

    BOOST_TEST(*hydra.getFlatdXdT() == dXdT, CHECK_PER_ELEMENT);

    {
        constexpr unsigned int VAR_SIZE = 18;

        constexpr unsigned int OUT_SIZE = 20;

        std::vector<double> X = additionalDOF;

        for (unsigned int i = 0; i < VAR_SIZE; i++) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;

            floatVector xm = X;

            xp[i] += delta;

            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, xp, previousAdditionalDOF);

            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, xm, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydrap.evaluate();

            hydram.evaluate();

            floatVector vp = *hydrap.getUnknownVector();

            floatVector vm = *hydram.getUnknownVector();

            for (unsigned int j = 0; j < OUT_SIZE; j++) {
                dXdAdditionalDOF[VAR_SIZE * j + i] = (vp[j] - vm[j]) / (2 * delta);
            }
        }
    }

    BOOST_TEST(*hydra.getFlatdXdAdditionalDOF() == dXdAdditionalDOF, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_residual_exampleModel2, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {0, 1};

        tardigradeHydra::linearElasticity::residual elasticity;

        residualMock dofDeformation;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(2);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            dofDeformation = residualMock(this, 11, 1, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);
            dofDeformation.setUseTrapezoidalIntegration(true);

            residuals[0] = &elasticity;

            residuals[1] = &dofDeformation;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    floatVector parameters = {123.4, 56.7, 0.78, 0.89};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 2;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    floatVector dXdF(20 * 9, 0);

    floatVector dXdT(20, 0);

    floatVector dXdAdditionalDOF(20 * 18, 0);

    hydra.evaluate();

    floatType eps = 1e-6;

    {
        constexpr unsigned int VAR_SIZE = 9;

        constexpr unsigned int OUT_SIZE = 20;

        std::vector<double> X(std::begin(deformationGradient), std::end(deformationGradient));

        for (unsigned int i = 0; i < VAR_SIZE; i++) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;

            floatVector xm = X;

            xp[i] += delta;

            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, xp,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, xm,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydrap.evaluate();

            hydram.evaluate();

            floatVector vp = *hydrap.getUnknownVector();

            floatVector vm = *hydram.getUnknownVector();

            for (unsigned int j = 0; j < OUT_SIZE; j++) {
                dXdF[VAR_SIZE * j + i] = (vp[j] - vm[j]) / (2 * delta);
            }
        }
    }

    BOOST_TEST(*hydra.getFlatdXdF() == dXdF, CHECK_PER_ELEMENT);

    {
        constexpr unsigned int VAR_SIZE = 1;

        constexpr unsigned int OUT_SIZE = 20;

        std::vector<double> X = {temperature};

        for (unsigned int i = 0; i < VAR_SIZE; i++) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;

            floatVector xm = X;

            xp[i] += delta;

            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, xp[0], previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, xm[0], previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydrap.evaluate();

            hydram.evaluate();

            floatVector vp = *hydrap.getUnknownVector();

            floatVector vm = *hydram.getUnknownVector();

            for (unsigned int j = 0; j < OUT_SIZE; j++) {
                dXdT[VAR_SIZE * j + i] = (vp[j] - vm[j]) / (2 * delta);
            }
        }
    }

    BOOST_TEST(*hydra.getFlatdXdT() == dXdT, CHECK_PER_ELEMENT);

    {
        constexpr unsigned int VAR_SIZE = 18;

        constexpr unsigned int OUT_SIZE = 20;

        std::vector<double> X = additionalDOF;

        for (unsigned int i = 0; i < VAR_SIZE; i++) {
            floatType delta = eps * std::fabs(X[i]) + eps;

            floatVector xp = X;

            floatVector xm = X;

            xp[i] += delta;

            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, xp, previousAdditionalDOF);

            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, xm, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydraBaseMock hydram(dofm, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

            hydrap.evaluate();

            hydram.evaluate();

            floatVector vp = *hydrap.getUnknownVector();

            floatVector vm = *hydram.getUnknownVector();

            for (unsigned int j = 0; j < OUT_SIZE; j++) {
                dXdAdditionalDOF[VAR_SIZE * j + i] = (vp[j] - vm[j]) / (2 * delta);
            }
        }
    }

    BOOST_TEST(*hydra.getFlatdXdAdditionalDOF() == dXdAdditionalDOF, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_residual_suggestInitialIterateValues, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::dofVelocityGradientDeformation::residual {
       public:
        using tardigradeHydra::dofVelocityGradientDeformation::residual::residual;

        floatVector dofDeformationGradient = {0.111, 0.222, 0.333, 0.444, 0.555, 0.666, 0.777, 0.888, 0.999};

       protected:
        virtual void setDOFDeformationGradient() override { set_dofDeformationGradient(dofDeformationGradient); }

        floatVector initializeVector(unsigned int size) { return floatVector(size, 0); }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        floatVector elasticityParameters = {123.4, 56.7};

        floatVector dofDeformationParameters = {0.78, 0.89};

        std::vector<unsigned int> stateVariableIndices = {3, 4};

        tardigradeHydra::linearElasticity::residual elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> plasticity;

        residualMock dofDeformation;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        void setResidualClasses(std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> &residuals) {
            tardigradeHydra::hydraBase::setResidualClasses(residuals);
        }

       private:
        virtual void setResidualClasses() {
            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(4);

            elasticity = tardigradeHydra::linearElasticity::residual(this, 9, elasticityParameters);

            plasticity = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 9);

            dofDeformation = residualMock(this, 11, 2, 1, 12, 3, true, stateVariableIndices, dofDeformationParameters);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 3);

            residuals[0] = &elasticity;

            residuals[1] = &plasticity;

            residuals[2] = &dofDeformation;

            residuals[3] = &remainder;

            setResidualClasses(residuals);
        }
    };

    floatType time = 1.1;

    floatType deltaTime = 2.2;

    floatType temperature = 300.0;

    floatType previousTemperature = 320.4;

    floatVector deformationGradient = {1.1, 0.1, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector previousDeformationGradient = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    floatVector additionalDOF = {0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 0.99,
                                 1.11, 1.22, 1.33, 1.44, 1.55, 1.66, 1.77, 1.88, 1.99};

    floatVector previousAdditionalDOF = {-0.11, -0.22, -0.33, -0.44, -0.55, -0.66, -0.77, -0.88, -0.99,
                                         -1.11, -1.22, -1.33, -1.44, -1.55, -1.66, -1.77, -1.88, -1.99};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    floatVector previousStateVariables = {0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.3,
                                          0.4, 0.5, 0.6, 0.7, 0.8, 0.9, -1,  -2,  -3,  -4,  -5};

    floatVector parameters = {123.4, 56.7, 0.78};

    floatVector unknownVector = {1,    1,    1,    1,    1,    1,    1,   1,   1,    1.1,  0.12,
                                 0.13, 0.14, 1.3,  0.15, 0.16, 0.17, 1.4, 1.2, 0.21, 0.29, 0.23,
                                 0.9,  0.11, 0.30, 0.25, 1.1,  3,    4,   5,   6,    7};

    const unsigned int numConfigurations = 3;

    const unsigned int numNonLinearSolveStateVariables = 5;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    residualMock R(&hydra, 11, 2, 1, 12, 3, true, hydra.stateVariableIndices, hydra.dofDeformationParameters);

    std::vector<unsigned int> indices = {18, 19, 20, 21, 22, 23, 24, 25, 26, 30, 31};

    floatVector values = {0.111, 0.222, 0.333, 0.444, 0.555, 0.666, 0.777, 0.888, 0.999, 0.45474, 15.437454545};

    std::vector<unsigned int> result_1;
    std::vector<floatType>    result_2;

    R.suggestInitialIterateValues(result_1, result_2);

    BOOST_TEST(indices == result_1, CHECK_PER_ELEMENT);

    BOOST_TEST(values == result_2, CHECK_PER_ELEMENT);
}
