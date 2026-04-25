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

namespace tardigradeHydra {

    namespace unit_test {

        class hydraBaseTester {
           public:
            static void updateUnknownVector(tardigradeHydra::hydraBase &hydra, const floatVector &value) {
                BOOST_CHECK_NO_THROW(hydra.updateUnknownVector(value));
            }
        };

    }  // namespace unit_test

}  // namespace tardigradeHydra

BOOST_AUTO_TEST_CASE(test_HyperelasticBase_setFe, *boost::unit_test::tolerance(1e-5)) {
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

    tardigradeHydra::floatVector deformationGradient = {9.92294371e-01, -1.32912834e-02, 3.57093896e-02,
                                                        2.70156033e-02, 9.77585948e-01,  4.76270457e-04,
                                                        4.42875587e-02, -4.95172679e-02, 9.54139963e-01};

    tardigradeHydra::floatVector previousDeformationGradient = {
        1.02150915, -0.01813721, -0.01347062, 0.0406867, 1.0450375, 0.0038319, 0.07107588, 0.0013805, 0.98514977};

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
    tardigradeHydra::floatMatrix jac(deformationGradient.size(),
                                     tardigradeHydra::floatVector(deformationGradient.size(), 0));
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

    jac = tardigradeHydra::floatMatrix(deformationGradient.size(),
                                       tardigradeHydra::floatVector(deformationGradient.size(), 0));
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

    jac = tardigradeHydra::floatMatrix(deformationGradient.size(),
                                       tardigradeHydra::floatVector(previousStateVariables.size(), 0));
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

    jac = tardigradeHydra::floatMatrix(deformationGradient.size(),
                                       tardigradeHydra::floatVector(previousStateVariables.size(), 0));
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

BOOST_AUTO_TEST_CASE(test_HyperelasticBase_setJe, *boost::unit_test::tolerance(1e-6)) {
    class HyperelasticBaseMock : public tardigradeHydra::HyperelasticBase {
       public:
        using tardigradeHydra::HyperelasticBase::HyperelasticBase;

        tardigradeHydra::floatVector Fe = {2, 1, 2, 3, 6, 3, 2, 1, 10};

        tardigradeHydra::floatVector previousFe = {1.1, 0.2, 0.3, 0.4, 1.5, 0.6, 0.7, 0.8, 1.9};

       protected:
        virtual void setFe(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousFe();
                *v.value = previousFe;

            } else {
                auto v   = get_SetDataStorage_Fe();
                *v.value = Fe;
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

    tardigradeHydra::floatVector deformationGradient = {9.92294371e-01, -1.32912834e-02, 3.57093896e-02,
                                                        2.70156033e-02, 9.77585948e-01,  4.76270457e-04,
                                                        4.42875587e-02, -4.95172679e-02, 9.54139963e-01};

    tardigradeHydra::floatVector previousDeformationGradient = {
        1.02150915, -0.01813721, -0.01347062, 0.0406867, 1.0450375, 0.0038319, 0.07107588, 0.0013805, 0.98514977};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {0.00315514, 0.00318276, 0.0134401,   0.03494318, 0.02244553,
                                                           0.01110235, 0.02224434, -0.01770411, -0.01382113};

    tardigradeHydra::floatVector parameters = {60., 100.0, 0.94, 50., 0.7, 1e3, 1.2, 1.4};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations,
                                                                numNonLinearSolveStateVariables);

    hydraBaseMock hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::floatType answer = 72.0;

    tardigradeHydra::floatType previousAnswer = 2.32;

    HyperelasticBaseMock R(&hydra, 9);

    BOOST_TEST(answer == *R.get_Je());

    BOOST_TEST(previousAnswer == *R.get_previousJe());
}

BOOST_AUTO_TEST_CASE(test_HyperelasticBase_setJe_derivatives, *boost::unit_test::tolerance(1e-6)) {
    class HyperelasticBaseMock : public tardigradeHydra::HyperelasticBase {
       public:
        using tardigradeHydra::HyperelasticBase::HyperelasticBase;

        tardigradeHydra::floatVector Fe = {2, 1, 2, 3, 6, 3, 2, 1, 10};

        tardigradeHydra::floatVector previousFe = {1.1, 0.2, 0.3, 0.4, 1.5, 0.6, 0.7, 0.8, 1.9};

       protected:
        virtual void setFe(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousFe();
                *v.value = previousFe;

            } else {
                auto v   = get_SetDataStorage_Fe();
                *v.value = Fe;
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

    tardigradeHydra::floatVector deformationGradient = {9.92294371e-01, -1.32912834e-02, 3.57093896e-02,
                                                        2.70156033e-02, 9.77585948e-01,  4.76270457e-04,
                                                        4.42875587e-02, -4.95172679e-02, 9.54139963e-01};

    tardigradeHydra::floatVector previousDeformationGradient = {
        1.02150915, -0.01813721, -0.01347062, 0.0406867, 1.0450375, 0.0038319, 0.07107588, 0.0013805, 0.98514977};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {0.00315514, 0.00318276, 0.0134401,   0.03494318, 0.02244553,
                                                           0.01110235, 0.02224434, -0.01770411, -0.01382113};

    tardigradeHydra::floatVector parameters = {60., 100.0, 0.94, 50., 0.7, 1e3, 1.2, 1.4};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations,
                                                                numNonLinearSolveStateVariables);

    hydraBaseMock hydra(dof, model_configuration);

    hydra.initialize();

    HyperelasticBaseMock R(&hydra, 9);

    {
        double                 eps      = 1e-6;
        constexpr unsigned int VAR_SIZE = 9;
        constexpr unsigned int OUT_SIZE = 1;
        std::vector<double>    x        = R.Fe;
        std::vector<double>    answer(VAR_SIZE * OUT_SIZE, 0);

        for (unsigned int i = 0; i < VAR_SIZE; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);
            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            tardigradeHydra::ModelConfigurationBase model_configurationp(previousStateVariables, parameters,
                                                                         numConfigurations,
                                                                         numNonLinearSolveStateVariables);
            tardigradeHydra::ModelConfigurationBase model_configurationm(previousStateVariables, parameters,
                                                                         numConfigurations,
                                                                         numNonLinearSolveStateVariables);

            hydraBaseMock hydrap(dofp, model_configurationp);
            hydraBaseMock hydram(dofm, model_configurationm);

            hydrap.initialize();
            hydram.initialize();

            HyperelasticBaseMock Rp(&hydrap, 9);
            HyperelasticBaseMock Rm(&hydram, 9);

            Rp.Fe = xp;
            Rm.Fe = xm;

            auto rp = *Rp.get_Je();
            auto rm = *Rm.get_Je();

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                answer[VAR_SIZE * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_dJedFe(), CHECK_PER_ELEMENT);
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int VAR_SIZE = 9;
        constexpr unsigned int OUT_SIZE = 9;
        std::vector<double>    x        = R.Fe;
        std::vector<double>    answer(VAR_SIZE * OUT_SIZE, 0);

        for (unsigned int i = 0; i < VAR_SIZE; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);
            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            tardigradeHydra::ModelConfigurationBase model_configurationp(previousStateVariables, parameters,
                                                                         numConfigurations,
                                                                         numNonLinearSolveStateVariables);
            tardigradeHydra::ModelConfigurationBase model_configurationm(previousStateVariables, parameters,
                                                                         numConfigurations,
                                                                         numNonLinearSolveStateVariables);

            hydraBaseMock hydrap(dofp, model_configurationp);
            hydraBaseMock hydram(dofm, model_configurationm);

            hydrap.initialize();
            hydram.initialize();

            HyperelasticBaseMock Rp(&hydrap, 9);
            HyperelasticBaseMock Rm(&hydram, 9);

            Rp.Fe = xp;
            Rm.Fe = xm;

            auto rp = *Rp.get_dJedFe();
            auto rm = *Rm.get_dJedFe();

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                answer[VAR_SIZE * j + i] = (rp[j] - rm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_d2JedFe2(), CHECK_PER_ELEMENT);
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int VAR_SIZE = 9;
        constexpr unsigned int OUT_SIZE = 1;
        std::vector<double>    x        = R.previousFe;
        std::vector<double>    answer(VAR_SIZE * OUT_SIZE, 0);

        for (unsigned int i = 0; i < VAR_SIZE; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);
            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            tardigradeHydra::ModelConfigurationBase model_configurationp(previousStateVariables, parameters,
                                                                         numConfigurations,
                                                                         numNonLinearSolveStateVariables);
            tardigradeHydra::ModelConfigurationBase model_configurationm(previousStateVariables, parameters,
                                                                         numConfigurations,
                                                                         numNonLinearSolveStateVariables);

            hydraBaseMock hydrap(dofp, model_configurationp);
            hydraBaseMock hydram(dofm, model_configurationm);

            hydrap.initialize();
            hydram.initialize();

            HyperelasticBaseMock Rp(&hydrap, 9);
            HyperelasticBaseMock Rm(&hydram, 9);

            Rp.previousFe = xp;
            Rm.previousFe = xm;

            auto rp = *Rp.get_previousJe();
            auto rm = *Rm.get_previousJe();

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                answer[VAR_SIZE * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_dPreviousJedPreviousFe(), CHECK_PER_ELEMENT);
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int VAR_SIZE = 9;
        constexpr unsigned int OUT_SIZE = 9;
        std::vector<double>    x        = R.previousFe;
        std::vector<double>    answer(VAR_SIZE * OUT_SIZE, 0);

        for (unsigned int i = 0; i < VAR_SIZE; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);
            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            tardigradeHydra::ModelConfigurationBase model_configurationp(previousStateVariables, parameters,
                                                                         numConfigurations,
                                                                         numNonLinearSolveStateVariables);
            tardigradeHydra::ModelConfigurationBase model_configurationm(previousStateVariables, parameters,
                                                                         numConfigurations,
                                                                         numNonLinearSolveStateVariables);

            hydraBaseMock hydrap(dofp, model_configurationp);
            hydraBaseMock hydram(dofm, model_configurationm);

            hydrap.initialize();
            hydram.initialize();

            HyperelasticBaseMock Rp(&hydrap, 9);
            HyperelasticBaseMock Rm(&hydram, 9);

            Rp.previousFe = xp;
            Rm.previousFe = xm;

            auto rp = *Rp.get_dPreviousJedPreviousFe();
            auto rm = *Rm.get_dPreviousJedPreviousFe();

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                answer[VAR_SIZE * j + i] = (rp[j] - rm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_d2PreviousJedPreviousFe2(), CHECK_PER_ELEMENT);
    }
}

BOOST_AUTO_TEST_CASE(test_HyperelasticBase_setCauchyStress, *boost::unit_test::tolerance(1e-6)) {
    class HyperelasticBaseMock : public tardigradeHydra::HyperelasticBase {
       public:
        using tardigradeHydra::HyperelasticBase::HyperelasticBase;

        tardigradeHydra::floatVector Fe = {2, 2, 3, 4, 6, 6, 7, 8, 10};

        tardigradeHydra::floatVector previousFe = {1.1, 0.2, 0.3, 0.4, 1.5, 0.6, 0.7, 0.8, 1.9};

        tardigradeHydra::floatVector dStrainEnergydFe = {.11, .22, .33, .44, .55, .66, .77, .88, .99};

        tardigradeHydra::floatVector dPreviousStrainEnergydPreviousFe = {-.11, -.22, -.33, -.44, -.55,
                                                                         -.66, -.77, -.88, -.99};

       protected:
        virtual void setFe(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousFe();
                *v.value = previousFe;

            } else {
                auto v   = get_SetDataStorage_Fe();
                *v.value = Fe;
            }
        }

        virtual void setStrainEnergyJacobians(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_dPreviousStrainEnergydPreviousFe();
                *v.value = dPreviousStrainEnergydPreviousFe;

            } else {
                auto v   = get_SetDataStorage_dStrainEnergydFe();
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

    tardigradeHydra::floatVector deformationGradient = {9.92294371e-01, -1.32912834e-02, 3.57093896e-02,
                                                        2.70156033e-02, 9.77585948e-01,  4.76270457e-04,
                                                        4.42875587e-02, -4.95172679e-02, 9.54139963e-01};

    tardigradeHydra::floatVector previousDeformationGradient = {
        1.02150915, -0.01813721, -0.01347062, 0.0406867, 1.0450375, 0.0038319, 0.07107588, 0.0013805, 0.98514977};

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

    tardigradeHydra::floatVector answer = {-0.825, -1.87, -2.915, -1.98, -4.51, -7.04, -3.135, -7.15, -11.165};

    tardigradeHydra::floatVector previousAnswer = {-0.1137931,  -0.24655172, -0.37931034, -0.34137931, -0.60215517,
                                                   -0.86293103, -0.56896552, -0.95775862, -1.34655172};

    HyperelasticBaseMock R(&hydra, 9);

    BOOST_TEST(answer == *R.get_cauchyStress(), CHECK_PER_ELEMENT);

    BOOST_TEST(previousAnswer == *R.get_previousCauchyStress(), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_HyperelasticBase_setCauchyStressJacobians, *boost::unit_test::tolerance(1e-6)) {
    class HyperelasticBaseMock : public tardigradeHydra::HyperelasticBase {
       public:
        using tardigradeHydra::HyperelasticBase::HyperelasticBase;

        tardigradeHydra::floatVector _d2StrainEnergydFe2 = {
            1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
            28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
            55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81};

        tardigradeHydra::floatVector _d2StrainEnergydFedT = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

       protected:
        virtual void setStrainEnergyJacobians(const bool isPrevious) override {
            const tardigradeHydra::floatVector *Fe;

            tardigradeHydra::floatType T;

            tardigradeHydra::SetDataStorageBase<tardigradeHydra::floatVector> dStrainEnergydFe;

            if (isPrevious) {
                Fe = get_previousFe();

                T = hydra->getPreviousTemperature();

                dStrainEnergydFe = get_SetDataStorage_dPreviousStrainEnergydPreviousFe();

            } else {
                Fe = get_Fe();

                T = hydra->getTemperature();

                dStrainEnergydFe = get_SetDataStorage_dStrainEnergydFe();
            }

            dStrainEnergydFe.zero(9);

            for (unsigned int I = 0; I < 9; ++I) {
                (*dStrainEnergydFe.value)[I] += _d2StrainEnergydFedT[I] * T;

                for (unsigned int J = 0; J < 9; ++J) {
                    (*dStrainEnergydFe.value)[I] += _d2StrainEnergydFe2[9 * I + J] * (*Fe)[J];
                }
            }
        }

        virtual void setStrainEnergyHessians(const bool isPrevious) override {
            tardigradeHydra::SetDataStorageBase<tardigradeHydra::floatVector> d2StrainEnergydFe2;

            tardigradeHydra::SetDataStorageBase<tardigradeHydra::floatVector> d2StrainEnergydFedT;

            if (isPrevious) {
                d2StrainEnergydFe2 = get_SetDataStorage_d2PreviousStrainEnergydPreviousFe2();

                d2StrainEnergydFedT = get_SetDataStorage_d2PreviousStrainEnergydPreviousFedPreviousT();

            } else {
                d2StrainEnergydFe2 = get_SetDataStorage_d2StrainEnergydFe2();

                d2StrainEnergydFedT = get_SetDataStorage_d2StrainEnergydFedT();
            }

            *d2StrainEnergydFe2.value = _d2StrainEnergydFe2;

            *d2StrainEnergydFedT.value = _d2StrainEnergydFedT;
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

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 18);

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

    tardigradeHydra::floatVector deformationGradient = {9.92294371e-01, -1.32912834e-02, 3.57093896e-02,
                                                        2.70156033e-02, 9.77585948e-01,  4.76270457e-04,
                                                        4.42875587e-02, -4.95172679e-02, 9.54139963e-01};

    tardigradeHydra::floatVector previousDeformationGradient = {
        1.02150915, -0.01813721, -0.01347062, 0.0406867, 1.0450375, 0.0038319, 0.07107588, 0.0013805, 0.98514977};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.00315514, 0.00318276, 0.0134401, 0.03494318, 0.02244553, 0.01110235, 0.02224434, -0.01770411, -0.01382113,
        0,          0,          0,         0,          0,          0,          0,          0,           0};

    tardigradeHydra::floatVector parameters = {123.4, 56.7};

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 0;

    tardigradeHydra::floatVector unknownVector = {1,    1,    1,   1,    1,    1,    1,     1,    1,
                                                  1.02, .01,  .00, .00,  1.00, .00,  .00,   .03,  1.00,
                                                  0.95, -0.2, .01, 0.05, 1.03, 0.02, -0.01, 0.01, 0.99};

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations,
                                                                numNonLinearSolveStateVariables);

    hydraBaseMock hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    HyperelasticBaseMock R(&hydra, 9);

    {
        double                 eps         = 1e-6;
        constexpr unsigned int NUM_INPUTS  = 9;
        constexpr unsigned int NUM_OUTPUTS = 9;
        std::vector<double>    x           = deformationGradient;

        std::vector<double> answer(NUM_INPUTS * NUM_OUTPUTS, 0);

        for (unsigned int i = 0; i < NUM_INPUTS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, xp,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);
            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, xm,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration);
            hydraBaseMock hydram(dofm, model_configuration);

            hydrap.initialize();
            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);
            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

            HyperelasticBaseMock Rp(&hydrap, 9);
            HyperelasticBaseMock Rm(&hydram, 9);

            auto vp = *Rp.get_cauchyStress();
            auto vm = *Rm.get_cauchyStress();

            for (unsigned int j = 0; j < NUM_OUTPUTS; ++j) {
                answer[NUM_INPUTS * j + i] += (vp[j] - vm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_dCauchyStressdF(), CHECK_PER_ELEMENT);
    }

    {
        double                 eps         = 1e-6;
        constexpr unsigned int NUM_INPUTS  = 1;
        constexpr unsigned int NUM_OUTPUTS = 9;
        double                 x           = temperature;

        std::vector<double> answer(NUM_INPUTS * NUM_OUTPUTS, 0);

        for (unsigned int i = 0; i < NUM_INPUTS; ++i) {
            double delta = eps * std::fabs(x) + eps;

            double xp = x;
            double xm = x;

            xp += delta;
            xm -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, xp, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);
            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, xm, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration);
            hydraBaseMock hydram(dofm, model_configuration);

            hydrap.initialize();
            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);
            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

            HyperelasticBaseMock Rp(&hydrap, 9);
            HyperelasticBaseMock Rm(&hydram, 9);

            auto vp = *Rp.get_cauchyStress();
            auto vm = *Rm.get_cauchyStress();

            for (unsigned int j = 0; j < NUM_OUTPUTS; ++j) {
                answer[NUM_INPUTS * j + i] += (vp[j] - vm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_dCauchyStressdT(), CHECK_PER_ELEMENT);
    }

    {
        double                 eps         = 1e-6;
        constexpr unsigned int NUM_INPUTS  = 18;
        constexpr unsigned int NUM_OUTPUTS = 9;
        std::vector<double>    x           = unknownVector;

        std::vector<double> answer(NUM_INPUTS * NUM_OUTPUTS, 0);

        for (unsigned int i = 0; i < NUM_INPUTS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i + 9] += delta;
            xm[i + 9] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);
            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration);
            hydraBaseMock hydram(dofm, model_configuration);

            hydrap.initialize();
            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, xp);
            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, xm);

            HyperelasticBaseMock Rp(&hydrap, 9);
            HyperelasticBaseMock Rm(&hydram, 9);

            auto vp = *Rp.get_cauchyStress();
            auto vm = *Rm.get_cauchyStress();

            for (unsigned int j = 0; j < NUM_OUTPUTS; ++j) {
                answer[NUM_INPUTS * j + i] += (vp[j] - vm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_dCauchyStressdFn(), CHECK_PER_ELEMENT);
    }

    {
        double                 eps         = 1e-6;
        constexpr unsigned int NUM_INPUTS  = 9;
        constexpr unsigned int NUM_OUTPUTS = 9;
        std::vector<double>    x           = previousDeformationGradient;

        std::vector<double> answer(NUM_INPUTS * NUM_OUTPUTS, 0);

        for (unsigned int i = 0; i < NUM_INPUTS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 xp, additionalDOF, previousAdditionalDOF);
            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 xm, additionalDOF, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration);
            hydraBaseMock hydram(dofm, model_configuration);

            hydrap.initialize();
            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);
            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

            HyperelasticBaseMock Rp(&hydrap, 9);
            HyperelasticBaseMock Rm(&hydram, 9);

            auto vp = *Rp.get_previousCauchyStress();
            auto vm = *Rm.get_previousCauchyStress();

            for (unsigned int j = 0; j < NUM_OUTPUTS; ++j) {
                answer[NUM_INPUTS * j + i] += (vp[j] - vm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_dPreviousCauchyStressdPreviousF(), CHECK_PER_ELEMENT);
    }

    {
        double                 eps         = 1e-6;
        constexpr unsigned int NUM_INPUTS  = 1;
        constexpr unsigned int NUM_OUTPUTS = 9;
        double                 x           = previousTemperature;

        std::vector<double> answer(NUM_INPUTS * NUM_OUTPUTS, 0);

        for (unsigned int i = 0; i < NUM_INPUTS; ++i) {
            double delta = eps * std::fabs(x) + eps;

            double xp = x;
            double xm = x;

            xp += delta;
            xm -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, xp, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);
            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, xm, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration);
            hydraBaseMock hydram(dofm, model_configuration);

            hydrap.initialize();
            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);
            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

            HyperelasticBaseMock Rp(&hydrap, 9);
            HyperelasticBaseMock Rm(&hydram, 9);

            auto vp = *Rp.get_previousCauchyStress();
            auto vm = *Rm.get_previousCauchyStress();

            for (unsigned int j = 0; j < NUM_OUTPUTS; ++j) {
                answer[NUM_INPUTS * j + i] += (vp[j] - vm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_dPreviousCauchyStressdPreviousT(), CHECK_PER_ELEMENT);
    }

    {
        double                 eps         = 1e-6;
        constexpr unsigned int NUM_INPUTS  = 18;
        constexpr unsigned int NUM_OUTPUTS = 9;
        std::vector<double>    x           = previousStateVariables;

        std::vector<double> answer(NUM_INPUTS * NUM_OUTPUTS, 0);

        for (unsigned int i = 0; i < NUM_INPUTS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);
            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            tardigradeHydra::ModelConfigurationBase model_configurationp(xp, parameters, numConfigurations,
                                                                         numNonLinearSolveStateVariables);

            tardigradeHydra::ModelConfigurationBase model_configurationm(xm, parameters, numConfigurations,
                                                                         numNonLinearSolveStateVariables);

            hydraBaseMock hydrap(dofp, model_configurationp);
            hydraBaseMock hydram(dofm, model_configurationm);

            hydrap.initialize();
            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);
            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

            HyperelasticBaseMock Rp(&hydrap, 9);
            HyperelasticBaseMock Rm(&hydram, 9);

            auto vp = *Rp.get_previousCauchyStress();
            auto vm = *Rm.get_previousCauchyStress();

            for (unsigned int j = 0; j < NUM_OUTPUTS; ++j) {
                answer[NUM_INPUTS * j + i] += (vp[j] - vm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_dPreviousCauchyStressdPreviousFn(), CHECK_PER_ELEMENT);
    }
}

BOOST_AUTO_TEST_CASE(test_HyperelasticBase_setStress, *boost::unit_test::tolerance(1e-6)) {
    class HyperelasticBaseMock : public tardigradeHydra::HyperelasticBase {
       public:
        using tardigradeHydra::HyperelasticBase::HyperelasticBase;

        tardigradeHydra::floatVector cauchyStress = {1, 2, 3, 4, 5, 6, 7, 8, 9};

        tardigradeHydra::floatVector previousCauchyStress = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

       protected:
        virtual void setCauchyStress(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousCauchyStress();
                *v.value = previousCauchyStress;

            } else {
                auto v   = get_SetDataStorage_cauchyStress();
                *v.value = cauchyStress;
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

    tardigradeHydra::floatVector deformationGradient = {9.92294371e-01, -1.32912834e-02, 3.57093896e-02,
                                                        2.70156033e-02, 9.77585948e-01,  4.76270457e-04,
                                                        4.42875587e-02, -4.95172679e-02, 9.54139963e-01};

    tardigradeHydra::floatVector previousDeformationGradient = {
        1.02150915, -0.01813721, -0.01347062, 0.0406867, 1.0450375, 0.0038319, 0.07107588, 0.0013805, 0.98514977};

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

    HyperelasticBaseMock R(&hydra, 9);

    BOOST_TEST(R.cauchyStress == *R.getStress(), CHECK_PER_ELEMENT);

    BOOST_TEST(R.previousCauchyStress == *R.getPreviousStress(), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_HyperelasticBase_setResidual, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class HyperelasticBaseMock : public tardigradeHydra::HyperelasticBase {
       public:
        using tardigradeHydra::HyperelasticBase::HyperelasticBase;

        void setStress(tardigradeHydra::floatVector &cauchyStress) {
            tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>::setStress(cauchyStress);
        }

       private:
        void setStress() {
            tardigradeHydra::floatVector cauchyStress = {2, 3, 4, 5, 6, 7, 8, 9, 10};

            setStress(cauchyStress);
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

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 18);

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

    tardigradeHydra::floatVector deformationGradient = {0.39293837,  -0.42772133, -0.54629709, 0.10262954, 0.43893794,
                                                        -0.15378708, 0.9615284,   0.36965948,  -0.0381362};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {0.1,  0.2, 0.3, 0.2,  0.4,  -0.2, 0.3, 0.2, -0.1,
                                                           -0.5, 0.4, 0.2, -0.2, -0.1, 0.1,  0.1, 0.3, 0.4};

    tardigradeHydra::floatVector parameters = {123.4, 56.7};

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 0;

    tardigradeHydra::floatVector unknownVector = {1,  1,  1,  1,  1,   1,   1,   1,   1,   .1,  .2,  .3,  .4, .5,
                                                  .6, .7, .8, .9, .10, .11, .12, .13, .14, .15, .16, .17, .18};

    tardigradeHydra::floatVector residualAnswer = {1, 2, 3, 4, 5, 6, 7, 8, 9};

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations,
                                                                numNonLinearSolveStateVariables);

    hydraBaseMock hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    HyperelasticBaseMock R(&hydra, 9);

    BOOST_TEST(residualAnswer == *R.getResidual(), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_HyperelasticBase_Jacobians, *boost::unit_test::tolerance(1e-6)) {
    class HyperelasticBaseMock : public tardigradeHydra::HyperelasticBase {
       public:
        using tardigradeHydra::HyperelasticBase::HyperelasticBase;

        tardigradeHydra::floatVector _d2StrainEnergydFe2 = {
            1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
            28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
            55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81};

        tardigradeHydra::floatVector _d2StrainEnergydFedT = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

       protected:
        virtual void setStrainEnergyJacobians(const bool isPrevious) override {
            const tardigradeHydra::floatVector *Fe;

            tardigradeHydra::floatType T;

            tardigradeHydra::SetDataStorageBase<tardigradeHydra::floatVector> dStrainEnergydFe;

            if (isPrevious) {
                Fe = get_previousFe();

                T = hydra->getPreviousTemperature();

                dStrainEnergydFe = get_SetDataStorage_dPreviousStrainEnergydPreviousFe();

            } else {
                Fe = get_Fe();

                T = hydra->getTemperature();

                dStrainEnergydFe = get_SetDataStorage_dStrainEnergydFe();
            }

            dStrainEnergydFe.zero(9);

            for (unsigned int I = 0; I < 9; ++I) {
                (*dStrainEnergydFe.value)[I] += _d2StrainEnergydFedT[I] * T;

                for (unsigned int J = 0; J < 9; ++J) {
                    (*dStrainEnergydFe.value)[I] += _d2StrainEnergydFe2[9 * I + J] * (*Fe)[J];
                }
            }
        }

        virtual void setStrainEnergyHessians(const bool isPrevious) override {
            tardigradeHydra::SetDataStorageBase<tardigradeHydra::floatVector> d2StrainEnergydFe2;

            tardigradeHydra::SetDataStorageBase<tardigradeHydra::floatVector> d2StrainEnergydFedT;

            if (isPrevious) {
                d2StrainEnergydFe2 = get_SetDataStorage_d2PreviousStrainEnergydPreviousFe2();

                d2StrainEnergydFedT = get_SetDataStorage_d2PreviousStrainEnergydPreviousFedPreviousT();

            } else {
                d2StrainEnergydFe2 = get_SetDataStorage_d2StrainEnergydFe2();

                d2StrainEnergydFedT = get_SetDataStorage_d2StrainEnergydFedT();
            }

            *d2StrainEnergydFe2.value = _d2StrainEnergydFe2;

            *d2StrainEnergydFedT.value = _d2StrainEnergydFedT;
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

            remainder = tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(this, 18);

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

    tardigradeHydra::floatVector deformationGradient = {9.92294371e-01, -1.32912834e-02, 3.57093896e-02,
                                                        2.70156033e-02, 9.77585948e-01,  4.76270457e-04,
                                                        4.42875587e-02, -4.95172679e-02, 9.54139963e-01};

    tardigradeHydra::floatVector previousDeformationGradient = {
        1.02150915, -0.01813721, -0.01347062, 0.0406867, 1.0450375, 0.0038319, 0.07107588, 0.0013805, 0.98514977};

    tardigradeHydra::floatVector additionalDOF = {};

    tardigradeHydra::floatVector previousAdditionalDOF = {};

    tardigradeHydra::DOFStorageBase dof(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                        previousDeformationGradient, additionalDOF, previousAdditionalDOF);

    tardigradeHydra::floatVector previousStateVariables = {
        0.00315514, 0.00318276, 0.0134401, 0.03494318, 0.02244553, 0.01110235, 0.02224434, -0.01770411, -0.01382113,
        0,          0,          0,         0,          0,          0,          0,          0,           0};

    tardigradeHydra::floatVector parameters = {123.4, 56.7};

    unsigned int numConfigurations = 3;

    unsigned int numNonLinearSolveStateVariables = 0;

    tardigradeHydra::floatVector unknownVector = {1,    1,    1,   1,    1,    1,    1,     1,    1,
                                                  1.02, .01,  .00, .00,  1.00, .00,  .00,   .03,  1.00,
                                                  0.95, -0.2, .01, 0.05, 1.03, 0.02, -0.01, 0.01, 0.99};

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters, numConfigurations,
                                                                numNonLinearSolveStateVariables);

    hydraBaseMock hydra(dof, model_configuration);

    hydra.initialize();

    tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydra, unknownVector);

    HyperelasticBaseMock R(&hydra, 9);

    {
        double                 eps         = 1e-6;
        constexpr unsigned int NUM_INPUTS  = 9;
        constexpr unsigned int NUM_OUTPUTS = 9;
        std::vector<double>    x           = deformationGradient;

        std::vector<double> answer(NUM_INPUTS * NUM_OUTPUTS, 0);

        for (unsigned int i = 0; i < NUM_INPUTS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, xp,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);
            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, xm,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration);
            hydraBaseMock hydram(dofm, model_configuration);

            hydrap.initialize();
            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);
            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

            HyperelasticBaseMock Rp(&hydrap, 9);
            HyperelasticBaseMock Rm(&hydram, 9);

            auto vp = *Rp.getResidual();
            auto vm = *Rm.getResidual();

            for (unsigned int j = 0; j < NUM_OUTPUTS; ++j) {
                answer[NUM_INPUTS * j + i] += (vp[j] - vm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.getdRdF(), CHECK_PER_ELEMENT);
    }

    {
        double                 eps         = 1e-6;
        constexpr unsigned int NUM_INPUTS  = 1;
        constexpr unsigned int NUM_OUTPUTS = 9;
        double                 x           = temperature;

        std::vector<double> answer(NUM_INPUTS * NUM_OUTPUTS, 0);

        for (unsigned int i = 0; i < NUM_INPUTS; ++i) {
            double delta = eps * std::fabs(x) + eps;

            double xp = x;
            double xm = x;

            xp += delta;
            xm -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, xp, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);
            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, xm, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration);
            hydraBaseMock hydram(dofm, model_configuration);

            hydrap.initialize();
            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, unknownVector);
            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, unknownVector);

            HyperelasticBaseMock Rp(&hydrap, 9);
            HyperelasticBaseMock Rm(&hydram, 9);

            auto vp = *Rp.getResidual();
            auto vm = *Rm.getResidual();

            for (unsigned int j = 0; j < NUM_OUTPUTS; ++j) {
                answer[NUM_INPUTS * j + i] += (vp[j] - vm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.getdRdT(), CHECK_PER_ELEMENT);
    }

    {
        double                 eps         = 1e-6;
        constexpr unsigned int NUM_INPUTS  = 27;
        constexpr unsigned int NUM_OUTPUTS = 9;
        std::vector<double>    x           = unknownVector;

        std::vector<double> answer(NUM_INPUTS * NUM_OUTPUTS, 0);

        for (unsigned int i = 0; i < NUM_INPUTS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            tardigradeHydra::DOFStorageBase dofp(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);
            tardigradeHydra::DOFStorageBase dofm(time, deltaTime, temperature, previousTemperature, deformationGradient,
                                                 previousDeformationGradient, additionalDOF, previousAdditionalDOF);

            hydraBaseMock hydrap(dofp, model_configuration);
            hydraBaseMock hydram(dofm, model_configuration);

            hydrap.initialize();
            hydram.initialize();

            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydrap, xp);
            tardigradeHydra::unit_test::hydraBaseTester::updateUnknownVector(hydram, xm);

            HyperelasticBaseMock Rp(&hydrap, 9);
            HyperelasticBaseMock Rm(&hydram, 9);

            auto vp = *Rp.getResidual();
            auto vm = *Rm.getResidual();

            for (unsigned int j = 0; j < NUM_OUTPUTS; ++j) {
                answer[NUM_INPUTS * j + i] += (vp[j] - vm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.getJacobian(), CHECK_PER_ELEMENT);
    }
}
