/**
 * \file test_tardigrade_SequentialQuadraticProgrammingStep.cpp
 *
 * Tests for tardigrade_SequentialQuadraticProgrammingStep
 */

#include <tardigrade_SequentialQuadraticProgrammingStep.h>

#include "tardigrade_hydra.h"

#define BOOST_TEST_MODULE test_tardigrade_SequentialQuadraticProgrammingStep
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
            static void set_residual(hydraBase &hydra, const tardigradeHydra::floatVector &value) {
                hydra._residual.second = value;
                hydra._residual.first  = true;

                hydra.addIterationData(&hydra._residual);
            }

            static void set_unknownVector(hydraBase &hydra, const tardigradeHydra::floatVector &value) {
                hydra._X.second = value;
                hydra._X.first  = true;
            }

            static void set_flatJacobian(hydraBase &hydra, const tardigradeHydra::floatVector &value) {
                hydra._jacobian.second = value;
                hydra._jacobian.first  = true;

                hydra.addIterationData(&hydra._jacobian);
            }

            static void initializeUnknownVector(hydraBase &hydra) {
                BOOST_CHECK_NO_THROW(hydra.initializeUnknownVector());
            }
        };

        class StepDampingBaseTester {
           public:
            static void setMuk(StepDampingBase &damping, const tardigradeHydra::floatType &value) {
                damping.setMuk(value);
            }
        };

    }  // namespace unit_test

}  // namespace tardigradeHydra

BOOST_AUTO_TEST_CASE(test_SequentialQuadraticProgrammingStep_computeTrial,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class SequentialQuadraticProgrammingStepMock : public tardigradeHydra::SequentialQuadraticProgrammingStep {
       public:
        using tardigradeHydra::SequentialQuadraticProgrammingStep::SequentialQuadraticProgrammingStep;

        tardigradeHydra::floatVector initialUnknownVector = {2, 0};

       protected:
        virtual void assembleKKTRHSVector(const tardigradeHydra::floatVector &dx, tardigradeHydra::floatVector &RHS,
                                          const std::vector<bool> &active_constraints) override {
            RHS = tardigradeHydra::floatVector(7, 0);

            RHS[0] = 2 * (initialUnknownVector[0] + dx[0] - 1.0);
            RHS[1] = 2 * (initialUnknownVector[1] + dx[1] - 2.5);

            RHS[2 + 0] = (initialUnknownVector[0] + dx[0]) - 2 * (initialUnknownVector[1] + dx[1]) + 2;
            RHS[2 + 1] = -(initialUnknownVector[0] + dx[0]) - 2 * (initialUnknownVector[1] + dx[1]) + 6;
            RHS[2 + 2] = -(initialUnknownVector[0] + dx[0]) + 2 * (initialUnknownVector[1] + dx[1]) + 2;
            RHS[2 + 3] = initialUnknownVector[0] + dx[0];
            RHS[2 + 4] = initialUnknownVector[1] + dx[1];

            for (unsigned int i = 0; i < active_constraints.size(); i++) {
                if (!active_constraints[i]) {
                    RHS[2 + i] = 0;
                }
            }
        }

        virtual void assembleKKTMatrix(tardigradeHydra::floatVector &K,
                                       const std::vector<bool>      &active_constraints) override {
            const unsigned int numConstraints = getNumConstraints();

            K = tardigradeHydra::floatVector((2 + 5) * (2 + 5), 0);

            K[7 * 0 + 0] = 2;
            K[7 * 1 + 1] = 2;

            for (unsigned int i = 0; i < numConstraints; i++) {
                if (active_constraints[i]) {
                    K[7 * 0 + i + 2] = (*getConstraintJacobians())[2 * i + 0];
                    K[7 * 1 + i + 2] = (*getConstraintJacobians())[2 * i + 1];

                    K[7 * (i + 2) + 0] = (*getConstraintJacobians())[2 * i + 0];
                    K[7 * (i + 2) + 1] = (*getConstraintJacobians())[2 * i + 1];

                } else {
                    K[7 * (i + 2) + i + 2] = 1;
                }
            }
        }

        virtual void updateKKTMatrix(tardigradeHydra::floatVector &K,
                                     const std::vector<bool>      &active_constraints) override {
            const unsigned int numConstraints = getNumConstraints();

            for (unsigned int i = 0; i < numConstraints; i++) {
                if (active_constraints[i]) {
                    K[7 * 0 + i + 2] = (*getConstraintJacobians())[2 * i + 0];
                    K[7 * 1 + i + 2] = (*getConstraintJacobians())[2 * i + 1];

                    K[7 * (i + 2) + 0] = (*getConstraintJacobians())[2 * i + 0];
                    K[7 * (i + 2) + 1] = (*getConstraintJacobians())[2 * i + 1];

                    K[7 * (i + 2) + i + 2] = 0;

                } else {
                    K[7 * 0 + i + 2] = 0;
                    K[7 * 1 + i + 2] = 0;

                    K[7 * (i + 2) + 0] = 0;
                    K[7 * (i + 2) + 1] = 0;

                    K[7 * (i + 2) + i + 2] = 1;
                }
            }
        }

        virtual void initializeActiveConstraints(std::vector<bool> &active_constraints) override {
            active_constraints = {false, false, true, false, true};
        }
    };

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase {
       public:
        using tardigradeHydra::SolverStepBase::SolverStepBase;

        void setTrialStep(tardigradeHydra::TrialStepBase *ptr) { trial_step = ptr; }

       protected:
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        tardigradeHydra::floatVector initialUnknownVector = {2, 0};

       protected:
        virtual void setConstraints() override {
            auto constraints = get_SetDataStorage_constraints();

            *constraints.value = {2, 6, 2, 0, 0};

            for (unsigned int i = 0; i < 5; i++) {
                for (unsigned int j = 0; j < 2; j++) {
                    (*constraints.value)[i] += (*getConstraintJacobians())[2 * i + j] * initialUnknownVector[j];
                }
            }
        }

        virtual void setConstraintJacobians() override {
            auto constraintJacobians = get_SetDataStorage_constraintJacobians();

            *constraintJacobians.value = {1, -2, -1, -2, -1, 2, 1, 0, 0, 1};
        }

        virtual const unsigned int getNumUnknowns() override { return initialUnknownVector.size(); }

        virtual const unsigned int getNumConstraints() override { return 5; }
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

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::SolverBase            solver;
    SolverStepBaseMock                     step;
    SequentialQuadraticProgrammingStepMock trial_step;

    trial_step.step = &step;
    step.setTrialStep(&trial_step);

    solver.hydra = &hydra;
    solver.step  = &step;
    step.solver  = &solver;

    hydra.solver = &solver;

    tardigradeHydra::floatVector result = {0, 0};

    step.deltaX = result;

    tardigradeHydra::floatVector answer = {1.4, 1.7};

    trial_step.computeTrial();

    BOOST_TEST((step.deltaX + hydra.initialUnknownVector) == answer, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_SequentialQuadraticProgrammingStep_assembleKKTRHSVector,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class SequentialQuadraticProgrammingStepMock : public tardigradeHydra::SequentialQuadraticProgrammingStep {
       public:
        using tardigradeHydra::SequentialQuadraticProgrammingStep::SequentialQuadraticProgrammingStep;

        virtual void public_assembleKKTRHSVector(const tardigradeHydra::floatVector &dx,
                                                 tardigradeHydra::floatVector       &RHS,
                                                 const std::vector<bool>            &active_constraints) {
            assembleKKTRHSVector(dx, RHS, active_constraints);
        }
    };

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase {
       public:
        using tardigradeHydra::SolverStepBase::SolverStepBase;
    };

    class StepDampingBaseMock : public tardigradeHydra::StepDampingBase {
       public:
        using tardigradeHydra::StepDampingBase::StepDampingBase;
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        tardigradeHydra::floatVector initialUnknownVector = {2, 1};

        auto public_setMuk(const tardigradeHydra::floatType &value) {
            tardigradeHydra::unit_test::StepDampingBaseTester::setMuk(*(solver->step->damping), value);
        }

       protected:
        virtual void setConstraints() override {
            auto constraints = get_SetDataStorage_constraints();

            *constraints.value = {2, 6, 2, 0, 0};

            for (unsigned int i = 0; i < 5; i++) {
                for (unsigned int j = 0; j < 2; j++) {
                    (*constraints.value)[i] += (*getConstraintJacobians())[2 * i + j] * initialUnknownVector[j];
                }
            }
        }

        virtual void setConstraintJacobians() override {
            auto constraintJacobians = get_SetDataStorage_constraintJacobians();

            *constraintJacobians.value = {1, -2, -1, -2, -1, 2, 1, 0, 0, 1};
        }

        virtual void formNonLinearResidual() override {
            tardigradeHydra::floatVector residual = {1., 2.};

            tardigradeHydra::unit_test::hydraBaseTester::set_residual(*this, residual);
        }

        virtual void formNonLinearDerivatives() override {
            tardigradeHydra::floatVector jacobian = {std::pow(2, 0.5), 0.4, -0.1, std::pow(2, 0.5)};

            tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian(*this, jacobian);
        }

        virtual const unsigned int getNumUnknowns() override { return initialUnknownVector.size(); }

        virtual const unsigned int getNumConstraints() override { return 5; }
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

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::SolverBase            solver;
    SolverStepBaseMock                     step;
    StepDampingBaseMock                    damping;
    SequentialQuadraticProgrammingStepMock trial_step;

    step.trial_step = &trial_step;
    step.damping    = &damping;
    trial_step.step = &step;
    damping.step    = &step;

    hydra.solver = &solver;
    solver.hydra = &hydra;
    solver.step  = &step;
    step.solver  = &solver;

    tardigradeHydra::floatVector dx = {-0.2, 1.4};

    hydra.public_setMuk(0.1);

    tardigradeHydra::floatVector result_KKTRHSVector;
    std::vector<bool>            active_constraints(5, false);

    tardigradeHydra::floatVector answer1_KKTRHSVector = {1.38618326, 6.30757431, 0.0, 0.0, 0.0, 0.0, 0.0};

    tardigradeHydra::floatVector answer2_KKTRHSVector = {1.38618326, 6.30757431, 0.0, 0.0, 5., 0.0, 2.4};

    trial_step.public_assembleKKTRHSVector(dx, result_KKTRHSVector, active_constraints);

    BOOST_TEST(answer1_KKTRHSVector == result_KKTRHSVector, CHECK_PER_ELEMENT);

    active_constraints[2] = true;
    active_constraints[4] = true;

    result_KKTRHSVector.clear();

    trial_step.public_assembleKKTRHSVector(dx, result_KKTRHSVector, active_constraints);

    BOOST_TEST(answer2_KKTRHSVector == result_KKTRHSVector, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_SequentialQuadraticProgrammingStep_assembleKKTMatrix,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class SequentialQuadraticProgrammingStepMock : public tardigradeHydra::SequentialQuadraticProgrammingStep {
       public:
        using tardigradeHydra::SequentialQuadraticProgrammingStep::SequentialQuadraticProgrammingStep;

        virtual void public_assembleKKTMatrix(tardigradeHydra::floatVector &K,
                                              const std::vector<bool>      &active_constraints) {
            assembleKKTMatrix(K, active_constraints);
        }

        virtual void public_updateKKTMatrix(tardigradeHydra::floatVector &K,
                                            const std::vector<bool>      &active_constraints) {
            updateKKTMatrix(K, active_constraints);
        }
    };

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase {
       public:
        using tardigradeHydra::SolverStepBase::SolverStepBase;
    };

    class StepDampingBaseMock : public tardigradeHydra::StepDampingBase {
       public:
        using tardigradeHydra::StepDampingBase::StepDampingBase;
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

        tardigradeHydra::floatVector initialUnknownVector = {2, 1};

        auto public_setMuk(const tardigradeHydra::floatType &value) {
            tardigradeHydra::unit_test::StepDampingBaseTester::setMuk(*(solver->step->damping), value);
        }

       protected:
        virtual void setConstraints() override {
            auto constraints = get_SetDataStorage_constraints();

            *constraints.value = {2, 6, 2, 0, 0};

            for (unsigned int i = 0; i < 5; i++) {
                for (unsigned int j = 0; j < 2; j++) {
                    (*constraints.value)[i] += (*getConstraintJacobians())[2 * i + j] * initialUnknownVector[j];
                }
            }
        }

        virtual void setConstraintJacobians() override {
            auto constraintJacobians = get_SetDataStorage_constraintJacobians();

            *constraintJacobians.value = {1, -2, -1, -2, -1, 2, 1, 0, 0, 1};
        }

        virtual void formNonLinearDerivatives() override {
            tardigradeHydra::floatVector jacobian = {std::pow(2, 0.5), 0.4, -0.1, std::pow(2, 0.5)};

            tardigradeHydra::unit_test::hydraBaseTester::set_flatJacobian(*this, jacobian);
        }

        virtual const unsigned int getNumUnknowns() override { return initialUnknownVector.size(); }

        virtual const unsigned int getNumConstraints() override { return 5; }
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

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::SolverBase            solver;
    SolverStepBaseMock                     step;
    StepDampingBaseMock                    damping;
    SequentialQuadraticProgrammingStepMock trial_step;

    step.trial_step = &trial_step;
    step.damping    = &damping;
    trial_step.step = &step;
    damping.step    = &step;

    hydra.solver = &solver;

    solver.hydra = &hydra;
    solver.step  = &step;

    step.solver = &solver;

    hydra.public_setMuk(0.1);

    tardigradeHydra::floatVector result_KKT;
    std::vector<bool>            active_constraints(5, false);

    tardigradeHydra::floatVector answer1_KKTMatrix = {
        2.11, 0.42426407, 0., 0., 0., 0., 0., 0.42426407, 2.26, 0., 0., 0., 0., 0., 0., 0., 1.,
        0.,   0.,         0., 0., 0., 0., 0., 1.,         0.,   0., 0., 0., 0., 0., 0., 1., 0.,
        0.,   0.,         0., 0., 0., 0., 1., 0.,         0.,   0., 0., 0., 0., 0., 1.};

    tardigradeHydra::floatVector answer2_KKTMatrix = {
        2.11, 0.42426407, 0., 0., -1., 0., 0., 0.42426407, 2.26, 0., 0., 2.,  0., 1., 0., 0., 1.,
        0.,   0.,         0., 0., 0.,  0., 0., 1.,         0.,   0., 0., -1., 2., 0., 0., 0., 0.,
        0.,   0.,         0., 0., 0.,  0., 1., 0.,         0.,   1., 0., 0.,  0., 0., 0.};

    trial_step.public_assembleKKTMatrix(result_KKT, active_constraints);

    BOOST_TEST(answer1_KKTMatrix == result_KKT, CHECK_PER_ELEMENT);

    active_constraints[2] = true;
    active_constraints[4] = true;

    trial_step.public_updateKKTMatrix(result_KKT, active_constraints);

    BOOST_TEST(answer2_KKTMatrix == result_KKT, CHECK_PER_ELEMENT);

    result_KKT.clear();

    trial_step.public_assembleKKTMatrix(result_KKT, active_constraints);

    BOOST_TEST(answer2_KKTMatrix == result_KKT, CHECK_PER_ELEMENT);
}

BOOST_AUTO_TEST_CASE(test_SequentialQuadraticProgrammingStep_initializeActiveConstraints,
                     *boost::unit_test::tolerance(DEFAULT_TEST_TOLERANCE)) {
    class SequentialQuadraticProgrammingStepMock : public tardigradeHydra::SequentialQuadraticProgrammingStep {
       public:
        using tardigradeHydra::SequentialQuadraticProgrammingStep::SequentialQuadraticProgrammingStep;

        void public_initializeActiveConstraints(std::vector<bool> &active_constraints) {
            initializeActiveConstraints(active_constraints);
        }
    };

    class SolverStepBaseMock : public tardigradeHydra::SolverStepBase {
       public:
        using tardigradeHydra::SolverStepBase::SolverStepBase;
    };

    class hydraBaseMock : public tardigradeHydra::hydraBase {
       public:
        using tardigradeHydra::hydraBase::hydraBase;

        using tardigradeHydra::hydraBase::setResidualClasses;

       protected:
        virtual void setConstraints() override {
            auto constraints = get_SetDataStorage_constraints();

            *constraints.value = {2, 6, -2, 0.1, 2};
        }

        virtual const unsigned int getNumConstraints() override { return 5; }
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

    tardigradeHydra::ModelConfigurationBase model_configuration(previousStateVariables, parameters);

    hydraBaseMock hydra(dof, model_configuration, numConfigurations, numNonLinearSolveStateVariables);

    tardigradeHydra::SolverBase            solver;
    SolverStepBaseMock                     step;
    SequentialQuadraticProgrammingStepMock trial_step;

    hydra.solver = &solver;

    solver.hydra = &hydra;
    solver.step  = &step;

    step.solver = &solver;

    step.trial_step = &trial_step;
    trial_step.step = &step;

    std::vector<bool> answer = {false, false, true, false, false};

    std::vector<bool> result;

    trial_step.public_initializeActiveConstraints(result);

    BOOST_TEST(result == answer, CHECK_PER_ELEMENT);
}
