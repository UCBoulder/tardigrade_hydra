/**
 * \file test_tardigrade_MooneyRivlinStrainEnergy.cpp
 *
 * Tests for tardigrade_MooneyRivlinStrainEnergy
 */

#include "tardigrade_MooneyRivlinStrainEnergy.h"
#include "tardigrade_SetDataStorage.h"

#define BOOST_TEST_MODULE test_tardigrade_MooneyRivlinStrainEnergy
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

BOOST_AUTO_TEST_CASE(test_MooneyRivlinStrainEnergy_setStrainEnergy, *boost::unit_test::tolerance(1e-5)) {
    class MooneyRivlinStrainEnergyMock : public tardigradeHydra::MooneyRivlinStrainEnergy {
       public:
        using tardigradeHydra::MooneyRivlinStrainEnergy::MooneyRivlinStrainEnergy;

        tardigradeHydra::floatVector Fe = {2, 1, 2, 3, 6, 3, 2, 1, 10};

        tardigradeHydra::floatVector previousFe = {1.1, 0.2, 0.3, 0.4, 1.5, 0.6, 0.7, 0.8, 1.9};

        tardigradeHydra::floatType public_compute_I1(const bool isPrevious) {
            return compute_I1<tardigradeHydra::floatType>(isPrevious);
        }

        void public_compute_dI1dFe(const bool isPrevious, std::vector<double> &dI1dFe) {
            compute_dI1dFe(isPrevious, std::begin(dI1dFe), std::end(dI1dFe));
        }

        void public_compute_d2I1dFe2(std::vector<double> &d2I1dFe2) {
            compute_d2I1dFe2(std::begin(d2I1dFe2), std::end(d2I1dFe2));
        }

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
        MooneyRivlinStrainEnergyMock elasticity;

        tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> remainder;

        unsigned int elasticitySize = 9;

        tardigradeHydra::floatVector elasticity_parameters = {1.23, 2.34, 3.45};

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            elasticity = MooneyRivlinStrainEnergyMock(this, elasticitySize, elasticity_parameters);

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

    tardigradeHydra::floatType answer = 17426.959549681244;

    tardigradeHydra::floatType previousAnswer = 11.594462065513707;

    MooneyRivlinStrainEnergyMock R(&hydra, 9, hydra.elasticity_parameters);

    BOOST_TEST(answer == *R.get_strainEnergy());

    BOOST_TEST(previousAnswer == *R.get_previousStrainEnergy());

    try {
        {
            double                 eps     = 2e-6;
            constexpr unsigned int NUM_VAR = 9;
            constexpr unsigned int NUM_OUT = 1;
            std::vector<double>    x       = R.Fe;
            std::vector<double>    jacobian(NUM_VAR * NUM_OUT, 0);

            for (unsigned int i = 0; i < NUM_VAR; ++i) {
                double delta = eps * std::fabs(x[i]) + eps;

                std::vector<double> xp = x;
                std::vector<double> xm = x;

                xp[i] += delta;
                xm[i] -= delta;

                hydraBaseMock hydrap(dof, model_configuration);
                hydraBaseMock hydram(dof, model_configuration);

                hydrap.initialize();
                hydram.initialize();

                MooneyRivlinStrainEnergyMock Rp(&hydrap, 9, hydra.elasticity_parameters);
                MooneyRivlinStrainEnergyMock Rm(&hydram, 9, hydra.elasticity_parameters);

                Rp.Fe = xp;
                Rm.Fe = xm;

                double vp = *Rp.get_strainEnergy();
                double vm = *Rm.get_strainEnergy();

                for (unsigned int j = 0; j < NUM_OUT; ++j) {
                    jacobian[NUM_VAR * j + i] = (vp - vm) / (2 * delta);
                }
            }

            BOOST_TEST(jacobian == *R.get_dStrainEnergydFe(), CHECK_PER_ELEMENT);
        }

        {
            double                 eps     = 1e-6;
            constexpr unsigned int NUM_VAR = 9;
            constexpr unsigned int NUM_OUT = 1;
            std::vector<double>    x       = R.previousFe;
            std::vector<double>    jacobian(NUM_VAR * NUM_OUT, 0);

            for (unsigned int i = 0; i < NUM_VAR; ++i) {
                double delta = eps * std::fabs(x[i]) + eps;

                std::vector<double> xp = x;
                std::vector<double> xm = x;

                xp[i] += delta;
                xm[i] -= delta;

                hydraBaseMock hydrap(dof, model_configuration);
                hydraBaseMock hydram(dof, model_configuration);

                hydrap.initialize();
                hydram.initialize();

                MooneyRivlinStrainEnergyMock Rp(&hydrap, 9, hydra.elasticity_parameters);
                MooneyRivlinStrainEnergyMock Rm(&hydram, 9, hydra.elasticity_parameters);

                Rp.previousFe = xp;
                Rm.previousFe = xm;

                double vp = *Rp.get_previousStrainEnergy();
                double vm = *Rm.get_previousStrainEnergy();

                for (unsigned int j = 0; j < NUM_OUT; ++j) {
                    jacobian[NUM_VAR * j + i] = (vp - vm) / (2 * delta);
                }
            }

            BOOST_TEST(jacobian == *R.get_dPreviousStrainEnergydPreviousFe(), CHECK_PER_ELEMENT);
        }

        {
            double                 eps     = 3e-6;
            constexpr unsigned int NUM_VAR = 9;
            constexpr unsigned int NUM_OUT = 9;
            std::vector<double>    x       = R.Fe;
            std::vector<double>    jacobian(NUM_VAR * NUM_OUT, 0);

            for (unsigned int i = 0; i < NUM_VAR; ++i) {
                double delta = eps * std::fabs(x[i]) + eps;

                std::vector<double> xp = x;
                std::vector<double> xm = x;

                xp[i] += delta;
                xm[i] -= delta;

                hydraBaseMock hydrap(dof, model_configuration);
                hydraBaseMock hydram(dof, model_configuration);

                hydrap.initialize();
                hydram.initialize();

                MooneyRivlinStrainEnergyMock Rp(&hydrap, 9, hydra.elasticity_parameters);
                MooneyRivlinStrainEnergyMock Rm(&hydram, 9, hydra.elasticity_parameters);

                Rp.Fe = xp;
                Rm.Fe = xm;

                std::vector<double> vp(NUM_OUT, 0);
                std::vector<double> vm(NUM_OUT, 0);

                vp = *Rp.get_dStrainEnergydFe();
                vm = *Rm.get_dStrainEnergydFe();

                for (unsigned int j = 0; j < NUM_OUT; ++j) {
                    jacobian[NUM_VAR * j + i] = (vp[j] - vm[j]) / (2 * delta);
                }
            }

            BOOST_TEST(jacobian == *R.get_d2StrainEnergydFe2(), CHECK_PER_ELEMENT);
        }

        {
            double                 eps     = 1e-6;
            constexpr unsigned int NUM_VAR = 9;
            constexpr unsigned int NUM_OUT = 9;
            std::vector<double>    x       = R.previousFe;
            std::vector<double>    jacobian(NUM_VAR * NUM_OUT, 0);

            for (unsigned int i = 0; i < NUM_VAR; ++i) {
                double delta = eps * std::fabs(x[i]) + eps;

                std::vector<double> xp = x;
                std::vector<double> xm = x;

                xp[i] += delta;
                xm[i] -= delta;

                hydraBaseMock hydrap(dof, model_configuration);
                hydraBaseMock hydram(dof, model_configuration);

                hydrap.initialize();
                hydram.initialize();

                MooneyRivlinStrainEnergyMock Rp(&hydrap, 9, hydra.elasticity_parameters);
                MooneyRivlinStrainEnergyMock Rm(&hydram, 9, hydra.elasticity_parameters);

                Rp.previousFe = xp;
                Rm.previousFe = xm;

                std::vector<double> vp(NUM_OUT, 0);
                std::vector<double> vm(NUM_OUT, 0);

                vp = *Rp.get_dPreviousStrainEnergydPreviousFe();
                vm = *Rm.get_dPreviousStrainEnergydPreviousFe();

                for (unsigned int j = 0; j < NUM_OUT; ++j) {
                    jacobian[NUM_VAR * j + i] = (vp[j] - vm[j]) / (2 * delta);
                }
            }

            BOOST_TEST(jacobian == *R.get_d2PreviousStrainEnergydPreviousFe2(), CHECK_PER_ELEMENT);
        }
    } catch (std::exception &e) {
        tardigradeErrorTools::printNestedExceptions(e);
        throw;
    }
}
