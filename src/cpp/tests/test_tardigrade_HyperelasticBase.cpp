/**
 * \file test_tardigrade_HyperelasticBase.cpp
 *
 * Tests for tardigrade_HyperelasticBase
 */

#include "tardigrade_HyperelasticBase.h"
#include "tardigrade_SetDataStorage.h"

#define BOOST_TEST_MODULE test_tardigrade_HyperelasticBase
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

BOOST_AUTO_TEST_CASE(test_residual_setFe, *boost::unit_test::tolerance(1e-5)) {
    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        tardigradeHydra::HyperelasticBase elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        unsigned int elasticitySize = 9;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            elasticity = tardigradeHydra::HyperelasticBase(this, elasticitySize);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 9);

            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(2);

            residuals[0] = &elasticity;

            residuals[1] = &remainder;

            setResidualClasses(residuals);
        }
    };

    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {9.92294371e-01, -1.32912834e-02, 3.57093896e-02,  2.70156033e-02, 9.77585948e-01,
                                       4.76270457e-04, 4.42875587e-02,  -4.95172679e-02, 9.54139963e-01};

    tardigradeHydra::floatVector previousDeformationGradient = {1.02150915, -0.01813721, -0.01347062, 0.0406867, 1.0450375,
                                               0.0038319,  0.07107588,  0.0013805,   0.98514977};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {0.00315514, 0.00318276, 0.0134401,   0.03494318, 0.02244553,
                                          0.01110235, 0.02224434, -0.01770411, -0.01382113};

    tardigradeHydra::floatVector parameters = {123.4, 56.7};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    tardigradeHydra::floatVector FeAnswer = {0.98921175,  -0.0156822, 0.02290497,  -0.00614278, 0.95596779,
                            -0.01019557, 0.02379954, -0.03175083, 0.96754518};

    tardigradeHydra::floatVector previousFeAnswer = {1.01964692,  -0.02138607, -0.02731485, 0.00513148, 1.0219469,
                                    -0.00768935, 0.04807642,  0.01848297,  0.99809319};

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations,
                                                                numNonLinearSolveStateVariables);

    hydraBaseMock hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::HyperelasticBase R(&hydra, 9);

    tardigradeHydra::HyperelasticBase RJ(&hydra, 9);

    RJ.get_dFedF();

    BOOST_TEST(FeAnswer == *R.get_Fe(), CHECK_PER_ELEMENT);

    BOOST_TEST(previousFeAnswer == *R.get_previousFe(), CHECK_PER_ELEMENT);

    BOOST_TEST(FeAnswer == *RJ.get_Fe(), CHECK_PER_ELEMENT);

    BOOST_TEST(previousFeAnswer == *RJ.get_previousFe(), CHECK_PER_ELEMENT);

    // Test the Jacobians
    tardigradeHydra::floatType   eps = 1e-6;
    tardigradeHydra::floatMatrix jac(deformationGradient.size(), tardigradeHydra::floatVector(deformationGradient.size(), 0));
    for (unsigned int i = 0; i < deformationGradient.size(); i++) {
        tardigradeHydra::floatVector delta(deformationGradient.size(), 0);

        delta[i] = eps * std::fabs(deformationGradient[i]) + eps;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature,
                                             deformationGradient + delta, previousDeformationGradient, additionalDOF,
                                             previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature,
                                             deformationGradient - delta, previousDeformationGradient, additionalDOF,
                                             previousAdditionalDOF);

        hydraBaseMock hydrap(dofp, model_configuration);

        hydraBaseMock hydram(dofm, model_configuration);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::HyperelasticBase Rp(&hydrap, 9);

        tardigradeHydra::HyperelasticBase Rm(&hydram, 9);

        tardigradeHydra::floatVector vp = *Rp.get_Fe();

        tardigradeHydra::floatVector vm = *Rm.get_Fe();

        for (unsigned int j = 0; j < deformationGradient.size(); j++) {
            jac[j][i] = (vp[j] - vm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(jac) == *R.get_dFedF(), CHECK_PER_ELEMENT);

    jac = tardigradeHydra::floatMatrix(deformationGradient.size(), tardigradeHydra::floatVector(deformationGradient.size(), 0));
    for (unsigned int i = 0; i < deformationGradient.size(); i++) {
        tardigradeHydra::floatVector delta(deformationGradient.size(), 0);

        delta[i] = eps * std::fabs(deformationGradient[i]) + eps;

        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient + delta, additionalDOF, previousAdditionalDOF);

        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                             previousDeformationGradient - delta, additionalDOF, previousAdditionalDOF);

        hydraBaseMock hydrap(dofp, model_configuration);

        hydraBaseMock hydram(dofm, model_configuration);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::HyperelasticBase Rp(&hydrap, 9);

        tardigradeHydra::HyperelasticBase Rm(&hydram, 9);

        tardigradeHydra::floatVector vp = *Rp.get_previousFe();

        tardigradeHydra::floatVector vm = *Rm.get_previousFe();

        for (unsigned int j = 0; j < deformationGradient.size(); j++) {
            jac[j][i] = (vp[j] - vm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(jac) == *R.get_previousdFedF(), CHECK_PER_ELEMENT);

    jac = tardigradeHydra::floatMatrix(deformationGradient.size(), tardigradeHydra::floatVector(previousStateVariables.size(), 0));
    for (unsigned int i = 0; i < previousStateVariables.size(); i++) {
        tardigradeHydra::floatVector delta(previousStateVariables.size(), 0);

        delta[i] = eps * std::fabs(previousStateVariables[i]) + eps;

        tardigradeHydra::ModelConfigurationBase model_configurationp(previousStateVariables + delta, parameters,
                                                                     numConfigurations,
                                                                     numNonLinearSolveStateVariables);

        tardigradeHydra::ModelConfigurationBase model_configurationm(previousStateVariables - delta, parameters,
                                                                     numConfigurations,
                                                                     numNonLinearSolveStateVariables);

        hydraBaseMock hydrap(dof, model_configurationp);

        hydraBaseMock hydram(dof, model_configurationm);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::HyperelasticBase Rp(&hydrap, 9);

        tardigradeHydra::HyperelasticBase Rm(&hydram, 9);

        tardigradeHydra::floatVector vp = *Rp.get_Fe();

        tardigradeHydra::floatVector vm = *Rm.get_Fe();

        for (unsigned int j = 0; j < deformationGradient.size(); j++) {
            jac[j][i] = (vp[j] - vm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(jac) == *R.get_dFedFn(), CHECK_PER_ELEMENT);

    jac = tardigradeHydra::floatMatrix(deformationGradient.size(), tardigradeHydra::floatVector(previousStateVariables.size(), 0));
    for (unsigned int i = 0; i < previousStateVariables.size(); i++) {
        tardigradeHydra::floatVector delta(previousStateVariables.size(), 0);

        delta[i] = eps * std::fabs(previousStateVariables[i]) + eps;

        tardigradeHydra::ModelConfigurationBase model_configurationp(previousStateVariables + delta, parameters,
                                                                     numConfigurations,
                                                                     numNonLinearSolveStateVariables);

        tardigradeHydra::ModelConfigurationBase model_configurationm(previousStateVariables - delta, parameters,
                                                                     numConfigurations,
                                                                     numNonLinearSolveStateVariables);

        hydraBaseMock hydrap(dof, model_configurationp);

        hydraBaseMock hydram(dof, model_configurationm);

        hydrap.initialize();

        hydram.initialize();

        tardigradeHydra::HyperelasticBase Rp(&hydrap, 9);

        tardigradeHydra::HyperelasticBase Rm(&hydram, 9);

        tardigradeHydra::floatVector vp = *Rp.get_previousFe();

        tardigradeHydra::floatVector vm = *Rm.get_previousFe();

        for (unsigned int j = 0; j < deformationGradient.size(); j++) {
            jac[j][i] = (vp[j] - vm[j]) / (2 * delta[i]);
        }
    }

    BOOST_TEST(tardigradeVectorTools::appendVectors(jac) == *R.get_previousdFedFn(), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_residual_setCauchyStress, *boost::unit_test::tolerance(1e-6)) {

    class HyperelasticBaseMock : public tardigradeHydra::HyperelasticBase {

        public:
            using tardigradeHydra::HyperelasticBase::HyperelasticBase;

            tardigradeHydra::floatVector Fe = {2, 2, 3, 4, 6, 6, 7, 8, 10};

            tardigradeHydra::floatVector previousFe = {1.1, 0.2, 0.3, 0.4, 1.5, 0.6, 0.7, 0.8, 1.9};

            tardigradeHydra::floatVector dStrainEnergydFe = {.11, .22, .33, .44, .55, .66, .77, .88, .99};

            tardigradeHydra::floatVector dPreviousStrainEnergydPreviousFe = {-.11, -.22, -.33, -.44, -.55, -.66, -.77, -.88, -.99};

        protected:
            virtual void setFe(const bool isPrevious) override {

                if ( isPrevious ){

                    auto v = get_SetDataStorage_previousFe();
                    *v.value = previousFe;

                }
                else{

                    auto v = get_SetDataStorage_Fe();
                    *v.value = Fe;

                }

            }

            virtual void setStrainEnergyJacobians(const bool isPrevious) override {

                if ( isPrevious ){

                    auto v = get_SetDataStorage_dPreviousStrainEnergydPreviousFe();
                    *v.value = dPreviousStrainEnergydPreviousFe;

                }
                else{

                    auto v = get_SetDataStorage_dStrainEnergydFe();
                    *v.value = dStrainEnergydFe;

                }

            }

    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        HyperelasticBaseMock elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        unsigned int elasticitySize = 9;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            elasticity = HyperelasticBaseMock(this, elasticitySize);

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 9);

            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(2);

            residuals[0] = &elasticity;

            residuals[1] = &remainder;

            setResidualClasses(residuals);
        }
    };

    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {9.92294371e-01, -1.32912834e-02, 3.57093896e-02,  2.70156033e-02, 9.77585948e-01,
                                       4.76270457e-04, 4.42875587e-02,  -4.95172679e-02, 9.54139963e-01};

    tardigradeHydra::floatVector previousDeformationGradient = {1.02150915, -0.01813721, -0.01347062, 0.0406867, 1.0450375,
                                               0.0038319,  0.07107588,  0.0013805,   0.98514977};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {0.00315514, 0.00318276, 0.0134401,   0.03494318, 0.02244553,
                                          0.01110235, 0.02224434, -0.01770411, -0.01382113};

    tardigradeHydra::floatVector parameters = {123.4, 56.7};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations,
                                                                numNonLinearSolveStateVariables);

    hydraBaseMock hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::floatVector answer = { -0.825,  -1.87 ,  -2.915,
                                            -1.98 ,  -4.51 ,  -7.04 ,
                                            -3.135,  -7.15 , -11.165 };

    tardigradeHydra::floatVector previousAnswer = { -0.1137931 , -0.24655172, -0.37931034,
                                                    -0.34137931, -0.60215517, -0.86293103,
                                                    -0.56896552, -0.95775862, -1.34655172 };

    HyperelasticBaseMock R(&hydra, 9);

    BOOST_TEST(answer == *R.get_cauchyStress(), CHECK_PER_ELEMENT);

    BOOST_TEST(previousAnswer == *R.get_previousCauchyStress(), CHECK_PER_ELEMENT);

//    HyperelasticBaseMock RJ(&hydra, 9);
//
//    RJ.get_dFedF();
//
//    BOOST_TEST(FeAnswer == *R.get_Fe(), CHECK_PER_ELEMENT);
//
//    BOOST_TEST(previousFeAnswer == *R.get_previousFe(), CHECK_PER_ELEMENT);
//
//    BOOST_TEST(FeAnswer == *RJ.get_Fe(), CHECK_PER_ELEMENT);
//
//    BOOST_TEST(previousFeAnswer == *RJ.get_previousFe(), CHECK_PER_ELEMENT);
//
//    // Test the Jacobians
//    tardigradeHydra::floatType   eps = 1e-6;
//    tardigradeHydra::floatMatrix jac(deformationGradient.size(), tardigradeHydra::floatVector(deformationGradient.size(), 0));
//    for (unsigned int i = 0; i < deformationGradient.size(); i++) {
//        tardigradeHydra::floatVector delta(deformationGradient.size(), 0);
//
//        delta[i] = eps * std::fabs(deformationGradient[i]) + eps;
//
//        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature,
//                                             deformationGradient + delta, previousDeformationGradient, additionalDOF,
//                                             previousAdditionalDOF);
//
//        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature,
//                                             deformationGradient - delta, previousDeformationGradient, additionalDOF,
//                                             previousAdditionalDOF);
//
//        hydraBaseMock hydrap(dofp, model_configuration);
//
//        hydraBaseMock hydram(dofm, model_configuration);
//
//        hydrap.initialize();
//
//        hydram.initialize();
//
//        tardigradeHydra::HyperelasticBase Rp(&hydrap, 9);
//
//        tardigradeHydra::HyperelasticBase Rm(&hydram, 9);
//
//        tardigradeHydra::floatVector vp = *Rp.get_Fe();
//
//        tardigradeHydra::floatVector vm = *Rm.get_Fe();
//
//        for (unsigned int j = 0; j < deformationGradient.size(); j++) {
//            jac[j][i] = (vp[j] - vm[j]) / (2 * delta[i]);
//        }
//    }
//
//    BOOST_TEST(tardigradeVectorTools::appendVectors(jac) == *R.get_dFedF(), CHECK_PER_ELEMENT);
//
//    jac = tardigradeHydra::floatMatrix(deformationGradient.size(), tardigradeHydra::floatVector(deformationGradient.size(), 0));
//    for (unsigned int i = 0; i < deformationGradient.size(); i++) {
//        tardigradeHydra::floatVector delta(deformationGradient.size(), 0);
//
//        delta[i] = eps * std::fabs(deformationGradient[i]) + eps;
//
//        tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
//                                             previousDeformationGradient + delta, additionalDOF, previousAdditionalDOF);
//
//        tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
//                                             previousDeformationGradient - delta, additionalDOF, previousAdditionalDOF);
//
//        hydraBaseMock hydrap(dofp, model_configuration);
//
//        hydraBaseMock hydram(dofm, model_configuration);
//
//        hydrap.initialize();
//
//        hydram.initialize();
//
//        tardigradeHydra::HyperelasticBase Rp(&hydrap, 9);
//
//        tardigradeHydra::HyperelasticBase Rm(&hydram, 9);
//
//        tardigradeHydra::floatVector vp = *Rp.get_previousFe();
//
//        tardigradeHydra::floatVector vm = *Rm.get_previousFe();
//
//        for (unsigned int j = 0; j < deformationGradient.size(); j++) {
//            jac[j][i] = (vp[j] - vm[j]) / (2 * delta[i]);
//        }
//    }
//
//    BOOST_TEST(tardigradeVectorTools::appendVectors(jac) == *R.get_previousdFedF(), CHECK_PER_ELEMENT);
//
//    jac = tardigradeHydra::floatMatrix(deformationGradient.size(), tardigradeHydra::floatVector(previousStateVariables.size(), 0));
//    for (unsigned int i = 0; i < previousStateVariables.size(); i++) {
//        tardigradeHydra::floatVector delta(previousStateVariables.size(), 0);
//
//        delta[i] = eps * std::fabs(previousStateVariables[i]) + eps;
//
//        tardigradeHydra::ModelConfigurationBase model_configurationp(previousStateVariables + delta, parameters,
//                                                                     numConfigurations,
//                                                                     numNonLinearSolveStateVariables);
//
//        tardigradeHydra::ModelConfigurationBase model_configurationm(previousStateVariables - delta, parameters,
//                                                                     numConfigurations,
//                                                                     numNonLinearSolveStateVariables);
//
//        hydraBaseMock hydrap(dof, model_configurationp);
//
//        hydraBaseMock hydram(dof, model_configurationm);
//
//        hydrap.initialize();
//
//        hydram.initialize();
//
//        tardigradeHydra::HyperelasticBase Rp(&hydrap, 9);
//
//        tardigradeHydra::HyperelasticBase Rm(&hydram, 9);
//
//        tardigradeHydra::floatVector vp = *Rp.get_Fe();
//
//        tardigradeHydra::floatVector vm = *Rm.get_Fe();
//
//        for (unsigned int j = 0; j < deformationGradient.size(); j++) {
//            jac[j][i] = (vp[j] - vm[j]) / (2 * delta[i]);
//        }
//    }
//
//    BOOST_TEST(tardigradeVectorTools::appendVectors(jac) == *R.get_dFedFn(), CHECK_PER_ELEMENT);
//
//    jac = tardigradeHydra::floatMatrix(deformationGradient.size(), tardigradeHydra::floatVector(previousStateVariables.size(), 0));
//    for (unsigned int i = 0; i < previousStateVariables.size(); i++) {
//        tardigradeHydra::floatVector delta(previousStateVariables.size(), 0);
//
//        delta[i] = eps * std::fabs(previousStateVariables[i]) + eps;
//
//        tardigradeHydra::ModelConfigurationBase model_configurationp(previousStateVariables + delta, parameters,
//                                                                     numConfigurations,
//                                                                     numNonLinearSolveStateVariables);
//
//        tardigradeHydra::ModelConfigurationBase model_configurationm(previousStateVariables - delta, parameters,
//                                                                     numConfigurations,
//                                                                     numNonLinearSolveStateVariables);
//
//        hydraBaseMock hydrap(dof, model_configurationp);
//
//        hydraBaseMock hydram(dof, model_configurationm);
//
//        hydrap.initialize();
//
//        hydram.initialize();
//
//        tardigradeHydra::HyperelasticBase Rp(&hydrap, 9);
//
//        tardigradeHydra::HyperelasticBase Rm(&hydram, 9);
//
//        tardigradeHydra::floatVector vp = *Rp.get_previousFe();
//
//        tardigradeHydra::floatVector vm = *Rm.get_previousFe();
//
//        for (unsigned int j = 0; j < deformationGradient.size(); j++) {
//            jac[j][i] = (vp[j] - vm[j]) / (2 * delta[i]);
//        }
//    }
//
//    BOOST_TEST(tardigradeVectorTools::appendVectors(jac) == *R.get_previousdFedFn(), CHECK_PER_ELEMENT);
}
