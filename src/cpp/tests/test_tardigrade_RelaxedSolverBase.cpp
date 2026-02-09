/**
 * \file test_tardigrade_RelaxedSolverBase.cpp
 *
 * Tests for tardigrade_RelaxedSolverBase
 */

#include "tardigrade_RelaxedSolverBase.h"
#include "tardigrade_ResidualBase.h"
#include "tardigrade_hydra.h"

#define BOOST_TEST_MODULE test_tardigrade_RelaxedSolverBase
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

BOOST_AUTO_TEST_CASE(test_RelaxedSolverBase_callResidualRelaxedStepFailure,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> {
       public:
        residualMock() : tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>() {}

        residualMock(tardigradeHydra::hydraBase *h, unsigned int num_vals, unsigned int p)
            : tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>(h, num_vals), ncallsuccess(p) {}

        unsigned int numRelaxedStepFailureCalls = 0;

        virtual bool relaxedStepFailure() override {
            numRelaxedStepFailureCalls++;
            if (numRelaxedStepFailureCalls > ncallsuccess) {
                return true;
            }
            return false;
        }

       protected:
        unsigned int ncallsuccess;
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        residualMock r1;

        residualMock r2;

        residualMock r3;

        unsigned int s1 = 36;

        unsigned int s2 = 2;

        unsigned int s3 = 3;

        unsigned int p1 = 1;

        unsigned int p2 = 1;

        unsigned int p3 = 2;

        bool r1_return = false;

        bool r2_return = false;

        bool r3_return = false;

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() {
            r1 = residualMock(this, s1, p1);

            r2 = residualMock(this, s2, p2);

            r3 = residualMock(this, s3, p3);

            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            residuals[0] = &r1;

            residuals[1] = &r2;

            residuals[2] = &r3;

            setResidualClasses(residuals);
        }
    };

    class RelaxedSolverBaseMock : public tardigradeHydra::RelaxedSolverBase {
       public:
        using tardigradeHydra::RelaxedSolverBase::RelaxedSolverBase;

        virtual bool public_callResidualRelaxedStepFailure() { return callResidualRelaxedStepFailure(); }
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

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079, 0.0,        0.1};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra(time, deltaTime, temperature, previousTemperature, deformationGradient,
                        previousDeformationGradient, {}, {}, previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables, dimension);

    RelaxedSolverBaseMock solver;

    hydra.solver = &solver;
    solver.hydra = &hydra;

    BOOST_TEST(hydra.r1.numRelaxedStepFailureCalls == 0);

    BOOST_TEST(hydra.r2.numRelaxedStepFailureCalls == 0);

    BOOST_TEST(hydra.r3.numRelaxedStepFailureCalls == 0);

    BOOST_TEST(!solver.public_callResidualRelaxedStepFailure());

    BOOST_TEST(hydra.r1.numRelaxedStepFailureCalls == 1);

    BOOST_TEST(hydra.r2.numRelaxedStepFailureCalls == 1);

    BOOST_TEST(hydra.r3.numRelaxedStepFailureCalls == 1);

    BOOST_TEST(solver.public_callResidualRelaxedStepFailure());

    BOOST_TEST(hydra.r1.numRelaxedStepFailureCalls == 2);

    BOOST_TEST(hydra.r2.numRelaxedStepFailureCalls == 2);

    BOOST_TEST(hydra.r3.numRelaxedStepFailureCalls == 2);
}

BOOST_AUTO_TEST_CASE(test_RelaxedSolverBase_solve, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> {
       public:
        bool project_called = false;

        using tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>::ResidualBase;

        virtual void projectSuggestedX(std::vector<double> &trialX, const std::vector<double> &Xp) override {
            project_called = true;
        }
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        residualMock r1;

        residualMock r2;

        using tardigradeHydra::hydraBase::hydraBase;

        void setInitialX() { solver->initial_unknown = _mockInitialX; }

       protected:
        tardigradeHydra::floatVector _mockInitialX = {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2};

        virtual void updateUnknownVector(const tardigradeHydra::floatVector &newX) override {
            BOOST_TEST(solver->initial_unknown == newX, CHECK_PER_ELEMENT);
        }

        using tardigradeHydra::hydraBase::setResidualClasses;
    };

    class RelaxedSolverBaseMock : public tardigradeHydra::RelaxedSolverBase {
       public:
        using tardigradeHydra::RelaxedSolverBase::RelaxedSolverBase;

        unsigned int num_calls = 0;

        bool calledPerformRelaxedSolve = false;

        virtual void solve() override {
            num_calls++;

            performRelaxedSolve();
        }

       protected:
        virtual void performRelaxedSolve() override { calledPerformRelaxedSolve = true; }
    };

    class SolverBaseMock : public tardigradeHydra::SolverBase {
       public:
        using tardigradeHydra::SolverBase::SolverBase;

        unsigned int num_calls = 0;

        virtual void solve() override {
            num_calls++;

            throw tardigradeHydra::convergence_error("failure to converge");
        }
    };

    tardigradeHydra::floatType time = 1.1;

    tardigradeHydra::floatType deltaTime = 2.2;

    tardigradeHydra::floatType temperature = 5.3;

    tardigradeHydra::floatType previousTemperature = 23.4;

    tardigradeHydra::floatVector deformationGradient = {1.05, 0, 0, 0.00, 1, 0, 0.00, 1, 1};

    tardigradeHydra::floatVector previousDeformationGradient = {-0.21576496, -0.31364397, 0.45809941,
                                                                -0.12285551, -0.88064421, -0.20391149,
                                                                0.47599081,  -0.63501654, -0.64909649};

    tardigradeHydra::floatVector previousStateVariables = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    tardigradeHydra::floatVector parameters = {123.4, 56.7};

    unsigned int numConfigurations = 2;

    unsigned int numNonLinearSolveStateVariables = 0;

    unsigned int dimension = 3;

    hydraBaseMock hydra(time, deltaTime, temperature, previousTemperature, deformationGradient,
                        previousDeformationGradient, {}, {}, previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables, dimension);

    RelaxedSolverBaseMock solver;
    SolverBaseMock    internal_solver;
    solver.internal_solver = &internal_solver;
    internal_solver.hydra  = &hydra;

    hydra.solver = &solver;
    solver.hydra = &hydra;

    solver.solve();

    BOOST_TEST(!hydra.solver->step->getRankDeficientError());

    BOOST_TEST(solver.num_calls == 1);

    BOOST_TEST(solver.calledPerformRelaxedSolve);
}

BOOST_AUTO_TEST_CASE(test_RelaxedSolverBase_performRelaxedSolve, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> {
       public:
        using tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>::ResidualBase;

        unsigned int numSetupCalls = 0;

        unsigned int relaxedStepConverged = 0;

        unsigned int currentRelaxedStep = 0;

        virtual void setupRelaxedStep(const unsigned int &relaxedIteration) override {
            currentRelaxedStep = relaxedIteration;

            numSetupCalls++;
        }

        virtual bool checkRelaxedConvergence() override { return relaxedStepConverged <= currentRelaxedStep; }

        virtual void setConvergedRelaxedIncrement(const unsigned int value) { relaxedStepConverged = value; }

       protected:
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        residualMock r1;

        residualMock r2;

        residualMock r3;

        unsigned int s1 = 36;

        unsigned int s2 = 2;

        unsigned int s3 = 3;

        unsigned int i1 = 1;

        unsigned int i2 = 0;

        unsigned int i3 = 3;

        unsigned int numCallInitializeUnknownVector = 0;

        unsigned int numCallUpdateUnknownVector = 0;

        tardigradeHydra::floatVector baseX = {1, 2, 3};

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() override {
            r1 = residualMock(this, s1);

            r2 = residualMock(this, s2);

            r3 = residualMock(this, s3);

            r1.setConvergedRelaxedIncrement(i1);

            r2.setConvergedRelaxedIncrement(i2);

            r3.setConvergedRelaxedIncrement(i3);

            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            residuals[0] = &r1;

            residuals[1] = &r2;

            residuals[2] = &r3;

            setResidualClasses(residuals);
        }

        virtual std::vector<unsigned int> getNumSetupCalls() {
            std::vector<unsigned int> result = {r1.numSetupCalls, r2.numSetupCalls, r3.numSetupCalls};

            return result;
        }

       protected:
        virtual void updateUnknownVector(const tardigradeHydra::floatVector &X) override {
            numCallUpdateUnknownVector++;
        }

        virtual void initializeUnknownVector() override {
            setX(baseX);
            numCallInitializeUnknownVector++;
        }
    };

    class RelaxedSolverBaseMock : public tardigradeHydra::RelaxedSolverBase {
       public:
        using tardigradeHydra::RelaxedSolverBase::RelaxedSolverBase;
    };

    class SolverBaseMock : public tardigradeHydra::SolverBase {
       public:
        using tardigradeHydra::SolverBase::SolverBase;

        unsigned int numCallSolveNonLinearProblem = 0;
        unsigned int numCallReset                 = 0;

        virtual void solve() override { numCallSolveNonLinearProblem++; }
        virtual void reset() override { numCallReset++; };
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

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079, 0.0,        0.1};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra(time, deltaTime, temperature, previousTemperature, deformationGradient,
                        previousDeformationGradient, {}, {}, previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables, dimension);

    // Form the relaxed solver
    SolverBaseMock    internal_solver;
    RelaxedSolverBaseMock solver;
    solver.internal_solver = &internal_solver;
    internal_solver.hydra  = &hydra;

    hydra.solver = &solver;
    solver.hydra = &hydra;

    solver.performRelaxedSolve();

    std::vector<unsigned int> answer_1 = {4, 4, 4};

    BOOST_TEST(internal_solver.numCallSolveNonLinearProblem == 4);

    BOOST_TEST(hydra.numCallUpdateUnknownVector == 4);

    BOOST_TEST(hydra.numCallInitializeUnknownVector == 1);

    BOOST_TEST(answer_1 == hydra.getNumSetupCalls(), CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_RelaxedSolverBase_performRelaxedSolve2, *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class residualMock : public tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> {
       public:
        using tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase>::ResidualBase;

        unsigned int numSetupCalls = 0;

        unsigned int relaxedStepConverged = 0;

        unsigned int currentRelaxedStep = 0;

        unsigned int numRelaxedStepFailureCalls = 0;

        unsigned int numCheckRelaxedConvergenceCalls = 0;

        virtual void setupRelaxedStep(const unsigned int &relaxedIteration) override {
            currentRelaxedStep = relaxedIteration;

            numSetupCalls++;
        }

        virtual bool checkRelaxedConvergence() override {
            numCheckRelaxedConvergenceCalls++;
            return relaxedStepConverged <= currentRelaxedStep;
        }

        virtual void setConvergedRelaxedIncrement(const unsigned int value) { relaxedStepConverged = value; }

        virtual bool relaxedStepFailure() override {
            numRelaxedStepFailureCalls++;
            return true;
        }

       protected:
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        residualMock r1;

        residualMock r2;

        residualMock r3;

        unsigned int s1 = 36;

        unsigned int s2 = 2;

        unsigned int s3 = 3;

        unsigned int i1 = 1;

        unsigned int i2 = 0;

        unsigned int i3 = 3;

        unsigned int numCallInitializeUnknownVector = 0;

        unsigned int numCallUpdateUnknownVector = 0;

        tardigradeHydra::floatVector baseX = {1, 2, 3};

        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        virtual void setResidualClasses() override {
            r1 = residualMock(this, s1);

            r2 = residualMock(this, s2);

            r3 = residualMock(this, s3);

            r1.setConvergedRelaxedIncrement(i1);

            r2.setConvergedRelaxedIncrement(i2);

            r3.setConvergedRelaxedIncrement(i3);

            std::vector<tardigradeHydra::ResidualBase<tardigradeHydra::hydraBase> *> residuals(3);

            residuals[0] = &r1;

            residuals[1] = &r2;

            residuals[2] = &r3;

            setResidualClasses(residuals);
        }

        virtual std::vector<unsigned int> getNumSetupCalls() {
            std::vector<unsigned int> result = {r1.numSetupCalls, r2.numSetupCalls, r3.numSetupCalls};

            return result;
        }

        virtual std::vector<unsigned int> getNumRelaxedStepFailureCalls() {
            std::vector<unsigned int> result = {r1.numRelaxedStepFailureCalls, r2.numRelaxedStepFailureCalls,
                                                r3.numRelaxedStepFailureCalls};

            return result;
        }

        virtual std::vector<unsigned int> getNumCheckRelaxedConvergenceCalls() {
            std::vector<unsigned int> result = {r1.numCheckRelaxedConvergenceCalls, r2.numCheckRelaxedConvergenceCalls,
                                                r3.numCheckRelaxedConvergenceCalls};

            return result;
        }

       protected:
        virtual void updateUnknownVector(const tardigradeHydra::floatVector &X) override {
            numCallUpdateUnknownVector++;
        }

        virtual void initializeUnknownVector() override {
            setX(baseX);
            numCallInitializeUnknownVector++;
        }
    };

    class RelaxedSolverBaseMock : public tardigradeHydra::RelaxedSolverBase {
       public:
        using tardigradeHydra::RelaxedSolverBase::RelaxedSolverBase;
    };

    class SolverBaseMock : public tardigradeHydra::SolverBase {
       public:
        using tardigradeHydra::SolverBase::SolverBase;

        unsigned int numCallSolveNonLinearProblem = 0;

        virtual void solve() override {
            numCallSolveNonLinearProblem++;
            if (numCallSolveNonLinearProblem <= 1) {
                throw tardigradeHydra::convergence_error("failure to converge");
            }
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

    tardigradeHydra::floatVector previousStateVariables = {
        0.53155137, 0.53182759, 0.63440096, 0.84943179, 0.72445532, 0.61102351, 0.72244338, 0.32295891,
        0.36178866, 0.22826323, 0.29371405, 0.63097612, 0.09210494, 0.43370117, 0.43086276, 0.4936851,
        0.42583029, 0.31226122, 0.42635131, 0.89338916, 0.94416002, 0.50183668, 0.62395295, 0.1156184,
        0.31728548, 0.41482621, 0.86630916, 0.25045537, 0.48303426, 0.98555979, 0.51948512, 0.61289453,
        0.12062867, 0.8263408,  0.60306013, 0.54506801, 0.34276383, 0.30412079, 0.0,        0.1};

    tardigradeHydra::floatVector parameters = {1, 2, 3, 4, 5};

    unsigned int numConfigurations = 4;

    unsigned int numNonLinearSolveStateVariables = 5;

    unsigned int dimension = 3;

    hydraBaseMock hydra(time, deltaTime, temperature, previousTemperature, deformationGradient,
                        previousDeformationGradient, {}, {}, previousStateVariables, parameters, numConfigurations,
                        numNonLinearSolveStateVariables, dimension);

    // Form the relaxed solver
    SolverBaseMock    internal_solver;
    RelaxedSolverBaseMock solver;
    solver.internal_solver = &internal_solver;
    internal_solver.hydra  = &hydra;

    hydra.solver = &solver;
    solver.hydra = &hydra;

    solver.performRelaxedSolve();

    std::vector<unsigned int> answer_1 = {4, 4, 4};

    std::vector<unsigned int> answer_2 = {1, 1, 1};

    std::vector<unsigned int> answer_3 = {3, 3, 3};

    BOOST_TEST(internal_solver.numCallSolveNonLinearProblem == 4);

    BOOST_TEST(hydra.numCallUpdateUnknownVector == 4);

    BOOST_TEST(hydra.numCallInitializeUnknownVector == 1);

    BOOST_TEST(answer_1 == hydra.getNumSetupCalls(), CHECK_PER_ELEMENT);

    BOOST_TEST(answer_2 == hydra.getNumRelaxedStepFailureCalls(), CHECK_PER_ELEMENT);

    BOOST_TEST(answer_3 == hydra.getNumCheckRelaxedConvergenceCalls(), CHECK_PER_ELEMENT);
}
