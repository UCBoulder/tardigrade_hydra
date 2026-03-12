/**
 * \file test_tardigrade_CHIPFoamStrainEnergy.cpp
 *
 * Tests for tardigrade_CHIPFoamStrainEnergy
 */

#include "tardigrade_CHIPFoamStrainEnergy.h"
#include "tardigrade_SetDataStorage.h"

#define BOOST_TEST_MODULE test_tardigrade_CHIPFoamStrainEnergy
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

BOOST_AUTO_TEST_CASE(test_CHIPFoamStrainEnergy_setJe, *boost::unit_test::tolerance(1e-6)) {
    class CHIPFoamStrainEnergyMock : public tardigradeHydra::CHIPFoamStrainEnergy {
       public:
        using tardigradeHydra::CHIPFoamStrainEnergy::CHIPFoamStrainEnergy;

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
        CHIPFoamStrainEnergyMock elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        unsigned int elasticitySize = 9;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            elasticity = CHIPFoamStrainEnergyMock(this, elasticitySize, *getParameters());

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

    CHIPFoamStrainEnergyMock R(&hydra, 9, parameters);

    BOOST_TEST(answer == *R.get_Je());

    BOOST_TEST(previousAnswer == *R.get_previousJe());
}

BOOST_AUTO_TEST_CASE(test_CHIPFoamStrainEnergy_setJe_derivatives, *boost::unit_test::tolerance(1e-6)) {
    class CHIPFoamStrainEnergyMock : public tardigradeHydra::CHIPFoamStrainEnergy {
       public:
        using tardigradeHydra::CHIPFoamStrainEnergy::CHIPFoamStrainEnergy;

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
        CHIPFoamStrainEnergyMock elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        unsigned int elasticitySize = 9;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            elasticity = CHIPFoamStrainEnergyMock(this, elasticitySize, *getParameters());

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

    CHIPFoamStrainEnergyMock R(&hydra, 9, parameters);

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

            CHIPFoamStrainEnergyMock Rp(&hydrap, 9, parameters);
            CHIPFoamStrainEnergyMock Rm(&hydram, 9, parameters);

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

            CHIPFoamStrainEnergyMock Rp(&hydrap, 9, parameters);
            CHIPFoamStrainEnergyMock Rm(&hydram, 9, parameters);

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

            CHIPFoamStrainEnergyMock Rp(&hydrap, 9, parameters);
            CHIPFoamStrainEnergyMock Rm(&hydram, 9, parameters);

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

            CHIPFoamStrainEnergyMock Rp(&hydrap, 9, parameters);
            CHIPFoamStrainEnergyMock Rm(&hydram, 9, parameters);

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

BOOST_AUTO_TEST_CASE(test_CHIPFoamStrainEnergy_setIbar1, *boost::unit_test::tolerance(1e-6)) {
    class CHIPFoamStrainEnergyMock : public tardigradeHydra::CHIPFoamStrainEnergy {
       public:
        using tardigradeHydra::CHIPFoamStrainEnergy::CHIPFoamStrainEnergy;

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
        CHIPFoamStrainEnergyMock elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        unsigned int elasticitySize = 9;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            elasticity = CHIPFoamStrainEnergyMock(this, elasticitySize, *getParameters());

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

    tardigradeHydra::floatType answer = 9.707057840908885;

    tardigradeHydra::floatType previousAnswer = 5.049921348278632;

    CHIPFoamStrainEnergyMock R(&hydra, 9, parameters);

    BOOST_TEST(answer == *R.get_Ibar1());

    BOOST_TEST(previousAnswer == *R.get_previousIbar1());
}

BOOST_AUTO_TEST_CASE(test_CHIPFoamStrainEnergy_setIbar1_derivatives, *boost::unit_test::tolerance(1e-6)) {
    class CHIPFoamStrainEnergyMock : public tardigradeHydra::CHIPFoamStrainEnergy {
       public:
        using tardigradeHydra::CHIPFoamStrainEnergy::CHIPFoamStrainEnergy;

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
        CHIPFoamStrainEnergyMock elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        unsigned int elasticitySize = 9;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            elasticity = CHIPFoamStrainEnergyMock(this, elasticitySize, *getParameters());

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

    CHIPFoamStrainEnergyMock R(&hydra, 9, parameters);

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

            CHIPFoamStrainEnergyMock Rp(&hydrap, 9, parameters);
            CHIPFoamStrainEnergyMock Rm(&hydram, 9, parameters);

            Rp.Fe = xp;
            Rm.Fe = xm;

            auto rp = *Rp.get_Ibar1();
            auto rm = *Rm.get_Ibar1();

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                answer[VAR_SIZE * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_dIbar1dFe(), CHECK_PER_ELEMENT);
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

            CHIPFoamStrainEnergyMock Rp(&hydrap, 9, parameters);
            CHIPFoamStrainEnergyMock Rm(&hydram, 9, parameters);

            Rp.Fe = xp;
            Rm.Fe = xm;

            auto rp = *Rp.get_dIbar1dFe();
            auto rm = *Rm.get_dIbar1dFe();

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                answer[VAR_SIZE * j + i] = (rp[j] - rm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_d2Ibar1dFe2(), CHECK_PER_ELEMENT);
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

            CHIPFoamStrainEnergyMock Rp(&hydrap, 9, parameters);
            CHIPFoamStrainEnergyMock Rm(&hydram, 9, parameters);

            Rp.previousFe = xp;
            Rm.previousFe = xm;

            auto rp = *Rp.get_previousIbar1();
            auto rm = *Rm.get_previousIbar1();

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                answer[VAR_SIZE * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_dPreviousIbar1dPreviousFe(), CHECK_PER_ELEMENT);
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

            CHIPFoamStrainEnergyMock Rp(&hydrap, 9, parameters);
            CHIPFoamStrainEnergyMock Rm(&hydram, 9, parameters);

            Rp.previousFe = xp;
            Rm.previousFe = xm;

            auto rp = *Rp.get_dPreviousIbar1dPreviousFe();
            auto rm = *Rm.get_dPreviousIbar1dPreviousFe();

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                answer[VAR_SIZE * j + i] = (rp[j] - rm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_d2PreviousIbar1dPreviousFe2(), CHECK_PER_ELEMENT);
    }
}

BOOST_AUTO_TEST_CASE(test_CHIPFoamStrainEnergy_setWLB, *boost::unit_test::tolerance(1e-6)) {
    class CHIPFoamStrainEnergyMock : public tardigradeHydra::CHIPFoamStrainEnergy {
       public:
        using tardigradeHydra::CHIPFoamStrainEnergy::CHIPFoamStrainEnergy;

        tardigradeHydra::floatType Je = 0.9;

        tardigradeHydra::floatType Ibar1 = 1.2;

        tardigradeHydra::floatType previousJe = 0.98;

        tardigradeHydra::floatType previousIbar1 = 2.2;

       protected:
        virtual void setJe(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousJe();
                *v.value = previousJe;

            } else {
                auto v   = get_SetDataStorage_Je();
                *v.value = Je;
            }
        }

        virtual void setIbar1(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousIbar1();
                *v.value = previousIbar1;

            } else {
                auto v   = get_SetDataStorage_Ibar1();
                *v.value = Ibar1;
            }
        }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        CHIPFoamStrainEnergyMock elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        unsigned int elasticitySize = 9;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            elasticity = CHIPFoamStrainEnergyMock(this, elasticitySize, *getParameters());

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

    tardigradeHydra::floatType answer = -53.58;

    tardigradeHydra::floatType previousAnswer = -23.979999999999993;

    CHIPFoamStrainEnergyMock R(&hydra, 9, parameters);

    BOOST_TEST(answer == *R.get_WLB());

    BOOST_TEST(previousAnswer == *R.get_previousWLB());
}

BOOST_AUTO_TEST_CASE(test_CHIPFoamStrainEnergy_setWLB_derivatives, *boost::unit_test::tolerance(1e-6)) {
    class CHIPFoamStrainEnergyMock : public tardigradeHydra::CHIPFoamStrainEnergy {
       public:
        using tardigradeHydra::CHIPFoamStrainEnergy::CHIPFoamStrainEnergy;

        tardigradeHydra::floatType Je = 0.9;

        tardigradeHydra::floatType Ibar1 = 1.2;

        tardigradeHydra::floatType previousJe = 0.98;

        tardigradeHydra::floatType previousIbar1 = 2.2;

       protected:
        virtual void setJe(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousJe();
                *v.value = previousJe;

            } else {
                auto v   = get_SetDataStorage_Je();
                *v.value = Je;
            }
        }

        virtual void setIbar1(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousIbar1();
                *v.value = previousIbar1;

            } else {
                auto v   = get_SetDataStorage_Ibar1();
                *v.value = Ibar1;
            }
        }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        CHIPFoamStrainEnergyMock elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        unsigned int elasticitySize = 9;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            elasticity = CHIPFoamStrainEnergyMock(this, elasticitySize, *getParameters());

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

    CHIPFoamStrainEnergyMock R(&hydra, 9, parameters);

    {
        double                 eps      = 1e-6;
        constexpr unsigned int VAR_SIZE = 2;
        constexpr unsigned int OUT_SIZE = 1;
        std::vector<double>    x        = {R.Je, R.Ibar1};
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

            CHIPFoamStrainEnergyMock Rp(&hydrap, 9, parameters);
            CHIPFoamStrainEnergyMock Rm(&hydram, 9, parameters);

            Rp.Je    = xp[0];
            Rp.Ibar1 = xp[1];
            Rm.Je    = xm[0];
            Rm.Ibar1 = xm[1];

            auto rp = *Rp.get_WLB();
            auto rm = *Rm.get_WLB();

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                answer[VAR_SIZE * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_dWLBdD(), CHECK_PER_ELEMENT);
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int VAR_SIZE = 2;
        constexpr unsigned int OUT_SIZE = 2;
        std::vector<double>    x        = {R.Je, R.Ibar1};
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

            CHIPFoamStrainEnergyMock Rp(&hydrap, 9, parameters);
            CHIPFoamStrainEnergyMock Rm(&hydram, 9, parameters);

            Rp.Je    = xp[0];
            Rp.Ibar1 = xp[1];
            Rm.Je    = xm[0];
            Rm.Ibar1 = xm[1];

            auto rp = *Rp.get_dWLBdD();
            auto rm = *Rm.get_dWLBdD();

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                answer[VAR_SIZE * j + i] = (rp[j] - rm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_d2WLBdD2(), CHECK_PER_ELEMENT);
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int VAR_SIZE = 2;
        constexpr unsigned int OUT_SIZE = 1;
        std::vector<double>    x        = {R.previousJe, R.previousIbar1};
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

            CHIPFoamStrainEnergyMock Rp(&hydrap, 9, parameters);
            CHIPFoamStrainEnergyMock Rm(&hydram, 9, parameters);

            Rp.previousJe    = xp[0];
            Rp.previousIbar1 = xp[1];
            Rm.previousJe    = xm[0];
            Rm.previousIbar1 = xm[1];

            auto rp = *Rp.get_previousWLB();
            auto rm = *Rm.get_previousWLB();

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                answer[VAR_SIZE * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_dPreviousWLBdPreviousD(), CHECK_PER_ELEMENT);
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int VAR_SIZE = 2;
        constexpr unsigned int OUT_SIZE = 2;
        std::vector<double>    x        = {R.previousJe, R.previousIbar1};
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

            CHIPFoamStrainEnergyMock Rp(&hydrap, 9, parameters);
            CHIPFoamStrainEnergyMock Rm(&hydram, 9, parameters);

            Rp.previousJe    = xp[0];
            Rp.previousIbar1 = xp[1];
            Rm.previousJe    = xm[0];
            Rm.previousIbar1 = xm[1];

            auto rp = *Rp.get_dPreviousWLBdPreviousD();
            auto rm = *Rm.get_dPreviousWLBdPreviousD();

            for (unsigned int j = 0; j < OUT_SIZE; ++j) {
                answer[VAR_SIZE * j + i] = (rp[j] - rm[j]) / (2 * delta);
            }
        }

        BOOST_TEST(answer == *R.get_d2PreviousWLBdPreviousD2(), CHECK_PER_ELEMENT);
    }
}

BOOST_AUTO_TEST_CASE(test_CHIPFoamStrainEnergy_compute_f, *boost::unit_test::tolerance(1e-6)) {
    class CHIPFoamStrainEnergyMock : public tardigradeHydra::CHIPFoamStrainEnergy {
       public:
        using tardigradeHydra::CHIPFoamStrainEnergy::CHIPFoamStrainEnergy;

        tardigradeHydra::floatType Je = 0.9;

        tardigradeHydra::floatType Ibar1 = 1.2;

        tardigradeHydra::floatType previousJe = 0.98;

        tardigradeHydra::floatType previousIbar1 = 2.2;

       protected:
        virtual void setJe(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousJe();
                *v.value = previousJe;

            } else {
                auto v   = get_SetDataStorage_Je();
                *v.value = Je;
            }
        }

        virtual void setIbar1(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousIbar1();
                *v.value = previousIbar1;

            } else {
                auto v   = get_SetDataStorage_Ibar1();
                *v.value = Ibar1;
            }
        }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        CHIPFoamStrainEnergyMock elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        unsigned int elasticitySize = 9;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            elasticity = CHIPFoamStrainEnergyMock(this, elasticitySize, *getParameters());

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

    tardigradeHydra::floatType answer = 0.302232035116331;

    tardigradeHydra::floatType previousAnswer = 0.300078515500334;

    CHIPFoamStrainEnergyMock R(&hydra, 9, parameters);

    BOOST_TEST(answer == R.compute_f(R.Je));

    BOOST_TEST(previousAnswer == R.compute_f(R.previousJe));

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {R.Je};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_f(xp[i]);
            auto rm = R.compute_f(xm[i]);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_dfdJ(R.Je));
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {R.Je};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_dfdJ(xp[i]);
            auto rm = R.compute_dfdJ(xm[i]);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_d2fdJ2(R.Je));
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {R.previousJe};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_f(xp[i]);
            auto rm = R.compute_f(xm[i]);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_dfdJ(R.previousJe));
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {R.previousJe};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_dfdJ(xp[i]);
            auto rm = R.compute_dfdJ(xm[i]);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_d2fdJ2(R.previousJe));
    }
}

BOOST_AUTO_TEST_CASE(test_CHIPFoamStrainEnergy_compute_Jg, *boost::unit_test::tolerance(1e-6)) {
    class CHIPFoamStrainEnergyMock : public tardigradeHydra::CHIPFoamStrainEnergy {
       public:
        using tardigradeHydra::CHIPFoamStrainEnergy::CHIPFoamStrainEnergy;

        tardigradeHydra::floatType Je = 0.9;

        tardigradeHydra::floatType Ibar1 = 1.2;

        tardigradeHydra::floatType previousJe = 0.98;

        tardigradeHydra::floatType previousIbar1 = 2.2;

       protected:
        virtual void setJe(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousJe();
                *v.value = previousJe;

            } else {
                auto v   = get_SetDataStorage_Je();
                *v.value = Je;
            }
        }

        virtual void setIbar1(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousIbar1();
                *v.value = previousIbar1;

            } else {
                auto v   = get_SetDataStorage_Ibar1();
                *v.value = Ibar1;
            }
        }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        CHIPFoamStrainEnergyMock elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        unsigned int elasticitySize = 9;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            elasticity = CHIPFoamStrainEnergyMock(this, elasticitySize, *getParameters());

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

    tardigradeHydra::floatType answer = 0.875379939209726;

    tardigradeHydra::floatType Jbar = 0.94;

    CHIPFoamStrainEnergyMock R(&hydra, 9, parameters);

    BOOST_TEST(answer == R.compute_Jg(Jbar, R.Je));

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {Jbar};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_Jg(xp[i], R.Je);
            auto rm = R.compute_Jg(xm[i], R.Je);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_dJgdJbar(Jbar, R.Je));
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {R.Je};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_Jg(Jbar, xp[i]);
            auto rm = R.compute_Jg(Jbar, xm[i]);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_dJgdJe(Jbar, R.Je));
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {Jbar};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_dJgdJbar(xp[i], R.Je);
            auto rm = R.compute_dJgdJbar(xm[i], R.Je);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_d2JgdJbar2(Jbar, R.Je));
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {Jbar};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_dJgdJe(xp[i], R.Je);
            auto rm = R.compute_dJgdJe(xm[i], R.Je);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_d2JgdJedJbar(Jbar, R.Je));
    }
}

BOOST_AUTO_TEST_CASE(test_CHIPFoamStrainEnergy_compute_pg, *boost::unit_test::tolerance(1e-6)) {
    class CHIPFoamStrainEnergyMock : public tardigradeHydra::CHIPFoamStrainEnergy {
       public:
        using tardigradeHydra::CHIPFoamStrainEnergy::CHIPFoamStrainEnergy;

        tardigradeHydra::floatType Je = 0.9;

        tardigradeHydra::floatType Ibar1 = 1.2;

        tardigradeHydra::floatType previousJe = 0.98;

        tardigradeHydra::floatType previousIbar1 = 2.2;

       protected:
        virtual void setJe(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousJe();
                *v.value = previousJe;

            } else {
                auto v   = get_SetDataStorage_Je();
                *v.value = Je;
            }
        }

        virtual void setIbar1(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousIbar1();
                *v.value = previousIbar1;

            } else {
                auto v   = get_SetDataStorage_Ibar1();
                *v.value = Ibar1;
            }
        }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        CHIPFoamStrainEnergyMock elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        unsigned int elasticitySize = 9;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            elasticity = CHIPFoamStrainEnergyMock(this, elasticitySize, *getParameters());

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

    tardigradeHydra::floatType answer = 0.245792675727452;

    tardigradeHydra::floatType Jbar = 0.94;

    CHIPFoamStrainEnergyMock R(&hydra, 9, parameters);

    BOOST_TEST(answer == R.compute_pg(Jbar, R.Je));

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {Jbar};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_pg(xp[i], R.Je);
            auto rm = R.compute_pg(xm[i], R.Je);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_dpgdJbar(Jbar, R.Je));
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {R.Je};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_pg(Jbar, xp[i]);
            auto rm = R.compute_pg(Jbar, xm[i]);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_dpgdJe(Jbar, R.Je));
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {Jbar};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_dpgdJbar(xp[i], R.Je);
            auto rm = R.compute_dpgdJbar(xm[i], R.Je);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_d2pgdJbar2(Jbar, R.Je));
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {Jbar};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_dpgdJe(xp[i], R.Je);
            auto rm = R.compute_dpgdJe(xm[i], R.Je);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_d2pgdJedJbar(Jbar, R.Je));
    }
}

BOOST_AUTO_TEST_CASE(test_CHIPFoamStrainEnergy_compute_ptilde, *boost::unit_test::tolerance(1e-6)) {
    class CHIPFoamStrainEnergyMock : public tardigradeHydra::CHIPFoamStrainEnergy {
       public:
        using tardigradeHydra::CHIPFoamStrainEnergy::CHIPFoamStrainEnergy;

        tardigradeHydra::floatType Je = 0.9;

        tardigradeHydra::floatType Ibar1 = 1.2;

        tardigradeHydra::floatType previousJe = 0.98;

        tardigradeHydra::floatType previousIbar1 = 2.2;

       protected:
        virtual void setJe(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousJe();
                *v.value = previousJe;

            } else {
                auto v   = get_SetDataStorage_Je();
                *v.value = Je;
            }
        }

        virtual void setIbar1(const bool isPrevious) override {
            if (isPrevious) {
                auto v   = get_SetDataStorage_previousIbar1();
                *v.value = previousIbar1;

            } else {
                auto v   = get_SetDataStorage_Ibar1();
                *v.value = Ibar1;
            }
        }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        CHIPFoamStrainEnergyMock elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        unsigned int elasticitySize = 9;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            elasticity = CHIPFoamStrainEnergyMock(this, elasticitySize, *getParameters());

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

    tardigradeHydra::floatType answer = 24.8653575879148;

    tardigradeHydra::floatType Jbar = 0.94;

    CHIPFoamStrainEnergyMock R(&hydra, 9, parameters);

    BOOST_TEST(answer == R.compute_ptilde(Jbar, R.Je));

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {Jbar};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_ptilde(xp[i], R.Je);
            auto rm = R.compute_ptilde(xm[i], R.Je);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_dptildedJbar(Jbar, R.Je));
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {R.Je};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_ptilde(Jbar, xp[i]);
            auto rm = R.compute_ptilde(Jbar, xm[i]);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_dptildedJe(Jbar, R.Je));
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {Jbar};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_dptildedJbar(xp[i], R.Je);
            auto rm = R.compute_dptildedJbar(xm[i], R.Je);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_d2ptildedJbar2(Jbar, R.Je));
    }

    {
        double                 eps      = 1e-6;
        constexpr unsigned int NUM_VARS = 1;
        constexpr unsigned int NUM_OUTS = 1;
        std::vector<double>    x        = {Jbar};
        std::vector<double>    jacobian(NUM_VARS * NUM_OUTS, 0);

        for (unsigned int i = 0; i < NUM_VARS; ++i) {
            double delta = eps * std::fabs(x[i]) + eps;

            std::vector<double> xp = x;
            std::vector<double> xm = x;

            xp[i] += delta;
            xm[i] -= delta;

            auto rp = R.compute_dptildedJe(xp[i], R.Je);
            auto rm = R.compute_dptildedJe(xm[i], R.Je);

            for (unsigned int j = 0; j < NUM_OUTS; ++j) {
                jacobian[NUM_VARS * j + i] = (rp - rm) / (2 * delta);
            }
        }

        BOOST_TEST(jacobian[0] == R.compute_d2ptildedJedJbar(Jbar, R.Je));
    }
}
