/**
 * \file test_tardigrade_ResidualBase.cpp
 *
 * Tests for tardigrade_ResidualBase
 */

#include "tardigrade_ResidualBase.h"
#include "tardigrade_SetDataStorage.h"

#define BOOST_TEST_MODULE test_tardigrade_ResidualBase
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

BOOST_AUTO_TEST_CASE(test_ResidualBase_ResidualBase, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class HydraMock {
       public:
        HydraMock() {

        };

        void addIterationData(tardigradeHydra::dataBase *) {

        };

        void addNLStepData(tardigradeHydra::dataBase *) {}
    };

    HydraMock hydra;

    unsigned int numEquations = 3;

    unsigned int numConstraints = 5;

    class ResidualBaseMock : public tardigradeHydra::ResidualBase<HydraMock> {
       public:
        using tardigradeHydra::ResidualBase<HydraMock>::ResidualBase;

        void public_setNumConstraints(const unsigned int &val) { setNumConstraints(val); }
    };

    ResidualBaseMock residual(&hydra, numEquations);

    residual.public_setNumConstraints(numConstraints);

    BOOST_CHECK(residual.hydra == &hydra);

    BOOST_CHECK(residual.getNumEquations() == numEquations);

    BOOST_CHECK(residual.getNumConstraints() == numConstraints);
}

BOOST_AUTO_TEST_CASE(test_ResidualBase_checkDefaults, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class HydraMock {
       public:
        HydraMock() {

        };

        void addIterationData(tardigradeHydra::dataBase *) {

        };

        void addNLStepData(tardigradeHydra::dataBase *) {}
    };

    HydraMock hydra;

    unsigned int numEquations = 3;

    tardigradeHydra::ResidualBase<HydraMock> residual(&hydra, numEquations);

#ifdef TARDIGRADE_ERROR_TOOLS_OPT
    BOOST_CHECK_THROW(residual.setResidual(), std::logic_error);

    BOOST_CHECK_THROW(residual.setJacobian(), std::logic_error);

    BOOST_CHECK_THROW(residual.setdRdF(), std::logic_error);

    BOOST_CHECK_THROW(residual.setdRdT(), std::logic_error);
#else
    BOOST_CHECK_THROW(residual.setResidual(), std::nested_exception);

    BOOST_CHECK_THROW(residual.setJacobian(), std::nested_exception);

    BOOST_CHECK_THROW(residual.setdRdF(), std::nested_exception);

    BOOST_CHECK_THROW(residual.setdRdT(), std::nested_exception);
#endif

    BOOST_CHECK_NO_THROW(residual.setAdditionalDerivatives());

    BOOST_CHECK(!residual.getUseProjection());
}

BOOST_AUTO_TEST_CASE(test_ResidualBase_setResidual, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class HydraMock {
       public:
        HydraMock() {

        };

        void addIterationData(tardigradeHydra::dataBase *) {

        };

        void addNLStepData(tardigradeHydra::dataBase *) {}
    };

    class ResidualBaseMock : public tardigradeHydra::ResidualBase<HydraMock> {
       public:
        using tardigradeHydra::ResidualBase<HydraMock>::setResidual;

        tardigradeHydra::floatVector residual = {1, 2, 3};

        ResidualBaseMock(HydraMock *hydra, unsigned int numEquations) : ResidualBase(hydra, numEquations) {}

        virtual void setResidual() { setResidual(residual); }
    };

    HydraMock hydra;

    unsigned int numEquations = 3;

    ResidualBaseMock residual(&hydra, numEquations);

    BOOST_TEST(*residual.getResidual() == residual.residual, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_ResidualBase_setJacobian, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class HydraMock {
       public:
        HydraMock() {

        };

        void addIterationData(tardigradeHydra::dataBase *) {

        };

        void addNLStepData(tardigradeHydra::dataBase *) {}
    };

    class ResidualBaseMock : public tardigradeHydra::ResidualBase<HydraMock> {
       public:
        using tardigradeHydra::ResidualBase<HydraMock>::setJacobian;

        tardigradeHydra::floatVector jacobian = {1, 2, 3, 4, 5, 6};

        ResidualBaseMock(HydraMock *hydra, unsigned int numEquations) : ResidualBase<HydraMock>(hydra, numEquations) {}

        virtual void setJacobian() { setJacobian(jacobian); }
    };

    HydraMock hydra;

    unsigned int numEquations = 3;

    ResidualBaseMock residual(&hydra, numEquations);

    BOOST_TEST(*residual.getJacobian() == residual.jacobian, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_ResidualBase_setdRdF, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class HydraMock {
       public:
        HydraMock() {

        };

        void addIterationData(tardigradeHydra::dataBase *) {

        };

        void addNLStepData(tardigradeHydra::dataBase *) {}
    };

    class ResidualBaseMock : public tardigradeHydra::ResidualBase<HydraMock> {
       public:
        using tardigradeHydra::ResidualBase<HydraMock>::setdRdF;

        tardigradeHydra::floatVector dRdF = {1, 2, 3, 4, 5, 6};

        ResidualBaseMock(HydraMock *hydra, unsigned int numEquations) : ResidualBase<HydraMock>(hydra, numEquations) {}

        virtual void setdRdF() { setdRdF(dRdF); }
    };

    HydraMock hydra;

    unsigned int numEquations = 3;

    ResidualBaseMock residual(&hydra, numEquations);

    BOOST_TEST(*residual.getdRdF() == residual.dRdF, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_ResidualBase_setdRdT, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class HydraMock {
       public:
        HydraMock() {

        };

        void addIterationData(tardigradeHydra::dataBase *) {

        };

        void addNLStepData(tardigradeHydra::dataBase *) {}
    };

    class ResidualBaseMock : public tardigradeHydra::ResidualBase<HydraMock> {
       public:
        using tardigradeHydra::ResidualBase<HydraMock>::setdRdT;

        tardigradeHydra::floatVector dRdT = {4, 5, 6};

        ResidualBaseMock(HydraMock *hydra, unsigned int numEquations) : ResidualBase<HydraMock>(hydra, numEquations) {}

        virtual void setdRdT() { setdRdT(dRdT); }
    };

    HydraMock hydra;

    unsigned int numEquations = 3;

    ResidualBaseMock residual(&hydra, numEquations);

    BOOST_TEST(*residual.getdRdT() == residual.dRdT, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_ResidualBase_setAdditionalDerivatives, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class HydraMock {
       public:
        HydraMock() {

        };

        void addIterationData(tardigradeHydra::dataBase *) {

        };

        void addNLStepData(tardigradeHydra::dataBase *) {}
    };

    class ResidualBaseMock : public tardigradeHydra::ResidualBase<HydraMock> {
       public:
        using tardigradeHydra::ResidualBase<HydraMock>::setAdditionalDerivatives;

        tardigradeHydra::floatVector additionalDerivatives = {4, 5, 6};

        ResidualBaseMock(HydraMock *hydra, unsigned int numEquations) : ResidualBase<HydraMock>(hydra, numEquations) {}

        virtual void setAdditionalDerivatives() { setAdditionalDerivatives(additionalDerivatives); }
    };

    HydraMock hydra;

    unsigned int numEquations = 3;

    ResidualBaseMock residual(&hydra, numEquations);

    BOOST_TEST(*residual.getAdditionalDerivatives() == residual.additionalDerivatives, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_ResidualBase_setStress, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class HydraMock {
       public:
        HydraMock() {

        };

        void addIterationData(tardigradeHydra::dataBase *) {

        };

        void addNLStepData(tardigradeHydra::dataBase *) {}
    };

    class ResidualBaseMock : public tardigradeHydra::ResidualBase<HydraMock> {
       public:
        using tardigradeHydra::ResidualBase<HydraMock>::setStress;

        tardigradeHydra::floatVector cauchyStress = {1, 2, 3, 4, 5, 6, 7, 8, 9};

        ResidualBaseMock(HydraMock *hydra, unsigned int numEquations) : ResidualBase<HydraMock>(hydra, numEquations) {}

        virtual void setStress() { setStress(cauchyStress); }
    };

    HydraMock hydra;

    unsigned int numEquations = 3;

    ResidualBaseMock residual(&hydra, numEquations);

    BOOST_TEST(*residual.getStress() == residual.cauchyStress, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_ResidualBase_setPreviousStress, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class HydraMock {
       public:
        HydraMock() {

        };

        void addIterationData(tardigradeHydra::dataBase *) {

        };

        void addNLStepData(tardigradeHydra::dataBase *) {}
    };

    class ResidualBaseMock : public tardigradeHydra::ResidualBase<HydraMock> {
       public:
        using tardigradeHydra::ResidualBase<HydraMock>::setPreviousStress;

        tardigradeHydra::floatVector previousCauchyStress = {1, 2, 3, 4, 5, 6, 7, 8, 9};

        ResidualBaseMock(HydraMock *hydra, unsigned int numEquations) : ResidualBase<HydraMock>(hydra, numEquations) {}

        virtual void setPreviousStress() { setPreviousStress(previousCauchyStress); }
    };

    HydraMock hydra;

    unsigned int numEquations = 3;

    ResidualBaseMock residual(&hydra, numEquations);

    BOOST_TEST(*residual.getPreviousStress() == residual.previousCauchyStress, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_ResidualBase_setCurrentAdditionalStateVariables,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class HydraMock {
       public:
        HydraMock() {}

        void addIterationData(tardigradeHydra::dataBase *) {}

        void addNLStepData(tardigradeHydra::dataBase *) {}
    };

    class ResidualBaseMock : public tardigradeHydra::ResidualBase<HydraMock> {
       public:
        using tardigradeHydra::ResidualBase<HydraMock>::setCurrentAdditionalStateVariables;

        tardigradeHydra::floatVector currentAdditionalStateVariables = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

        ResidualBaseMock(HydraMock *hydra, unsigned int numEquations) : ResidualBase<HydraMock>(hydra, numEquations) {}

        virtual void setCurrentAdditionalStateVariables() {
            setCurrentAdditionalStateVariables(currentAdditionalStateVariables);
        }
    };

    HydraMock hydra;

    unsigned int numEquations = 3;

    ResidualBaseMock residual(&hydra, numEquations);

    BOOST_TEST(*residual.getCurrentAdditionalStateVariables() == residual.currentAdditionalStateVariables,
               CHECK_PER_ELEMENT);
}
