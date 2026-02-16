/**
 * \file test_tardigrade_SubcyclerSolver.cpp
 *
 * Tests for tardigrade_SubcyclerSolver
 */

#include "tardigrade_SubcyclerSolver.h"
#include "tardigrade_hydra.h"

#define BOOST_TEST_MODULE test_tardigrade_SubcyclerSolver
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

BOOST_AUTO_TEST_CASE(test_SubcyclerSolver_callResidualPreSubcycler,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> {
       public:
        using tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>::ResidualBase;

        unsigned int numPreSubcyclerCalls = 0;

        virtual void preSubcycler() override { numPreSubcyclerCalls++; }

       protected:
    };

    class SubcyclerSolverMock : public tardigradeHydra::SubcyclerSolver {
       public:
        using tardigradeHydra::SubcyclerSolver::SubcyclerSolver;

        virtual void public_callResidualPreSubcycler() { callResidualPreSubcycler(); }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        residualMock r1;

        residualMock r2;

        residualMock r3;

        unsigned int s1 = 36;

        unsigned int s2 = 2;

        unsigned int s3 = 3;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            r1 = residualMock(this, s1);

            r2 = residualMock(this, s2);

            r3 = residualMock(this, s3);

            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            residuals[0] = &r1;

            residuals[1] = &r2;

            residuals[2] = &r3;

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

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079, 0.0,        0.1};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    hydraBaseMock hydra(dof,  previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables);

    SubcyclerSolverMock solver;
    solver.hydra = &hydra;
    hydra.solver = &solver;

    solver.public_callResidualPreSubcycler();

    BOOST_TEST(hydra.r1.numPreSubcyclerCalls == 1);

    BOOST_TEST(hydra.r2.numPreSubcyclerCalls == 1);

    BOOST_TEST(hydra.r3.numPreSubcyclerCalls == 1);
}

BOOST_AUTO_TEST_CASE(test_SubcyclerSolver_callResidualSubcyclerSuccess,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> {
       public:
        using tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>::ResidualBase;

        unsigned int numPostSubcyclerSuccessCalls = 0;

        virtual void postSubcyclerSuccess() override { numPostSubcyclerSuccessCalls++; }

       protected:
    };

    class SubcyclerSolverMock : public tardigradeHydra::SubcyclerSolver {
       public:
        using tardigradeHydra::SubcyclerSolver::SubcyclerSolver;

        virtual void public_callResidualPostSubcyclerSuccess() { callResidualPostSubcyclerSuccess(); }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        residualMock r1;

        residualMock r2;

        residualMock r3;

        unsigned int s1 = 36;

        unsigned int s2 = 2;

        unsigned int s3 = 3;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            r1 = residualMock(this, s1);

            r2 = residualMock(this, s2);

            r3 = residualMock(this, s3);

            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            residuals[0] = &r1;

            residuals[1] = &r2;

            residuals[2] = &r3;

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

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079, 0.0,        0.1};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    hydraBaseMock hydra(dof,  previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables);

    SubcyclerSolverMock solver;
    solver.hydra = &hydra;
    hydra.solver = &solver;

    solver.public_callResidualPostSubcyclerSuccess();

    BOOST_TEST(hydra.r1.numPostSubcyclerSuccessCalls == 1);

    BOOST_TEST(hydra.r2.numPostSubcyclerSuccessCalls == 1);

    BOOST_TEST(hydra.r3.numPostSubcyclerSuccessCalls == 1);
}

BOOST_AUTO_TEST_CASE(test_SubcyclerSolver_callResidualSubcyclerFailure,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> {
       public:
        using tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>::ResidualBase;

        unsigned int numPostSubcyclerFailureCalls = 0;

        virtual void postSubcyclerFailure() override { numPostSubcyclerFailureCalls++; }

       protected:
    };

    class SubcyclerSolverMock : public tardigradeHydra::SubcyclerSolver {
       public:
        using tardigradeHydra::SubcyclerSolver::SubcyclerSolver;

        virtual void public_callResidualPostSubcyclerFailure() { callResidualPostSubcyclerFailure(); }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        residualMock r1;

        residualMock r2;

        residualMock r3;

        unsigned int s1 = 36;

        unsigned int s2 = 2;

        unsigned int s3 = 3;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            r1 = residualMock(this, s1);

            r2 = residualMock(this, s2);

            r3 = residualMock(this, s3);

            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            residuals[0] = &r1;

            residuals[1] = &r2;

            residuals[2] = &r3;

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

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079, 0.0,        0.1};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    hydraBaseMock hydra(dof,  previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables);

    SubcyclerSolverMock solver;
    solver.hydra = &hydra;
    hydra.solver = &solver;

    solver.public_callResidualPostSubcyclerFailure();

    BOOST_TEST(hydra.r1.numPostSubcyclerFailureCalls == 1);

    BOOST_TEST(hydra.r2.numPostSubcyclerFailureCalls == 1);

    BOOST_TEST(hydra.r3.numPostSubcyclerFailureCalls == 1);
}
