/**
 ******************************************************************************
 * \file tardigrade_hydra.h
 ******************************************************************************
 * A C++ library for constructing finite deformation constitutive models.
 ******************************************************************************
 */

#ifndef TARDIGRADE_HYDRA_H
#define TARDIGRADE_HYDRA_H

#include <functional>
#include <sstream>

#include "tardigrade_error_tools.h"
//! We will use the functions that depend on Eigen
#define USE_EIGEN
#include "tardigrade_vector_tools.h"

#ifdef TARDIGRADE_HYDRA_USE_LLXSMM
#include <libxsmm.h>
#endif

#include "tardigrade_CoreDefinitions.h"
#include "tardigrade_CustomErrors.h"
#include "tardigrade_MatrixMap.h"
#include "tardigrade_RelaxedSolver.h"
#include "tardigrade_ResidualBase.h"
#include "tardigrade_SetDataStorage.h"
#include "tardigrade_SolverBase.h"
#include "tardigrade_SubcyclerSolver.h"
// #include"tardigrade_PreconditionerBase.h"
// #include"tardigrade_SolverStepBase.h"

namespace tardigradeHydra {

    /*!
     * Specialization for the DataStorage zero function of the residuals
     */
    template <>
    inline void DataStorage<std::vector<ResidualBase<hydraBase> *> >::zero() {
        throw std::runtime_error("Zeroing the ResidualBase pointer vector is not allowed");
    }

    /*!
     * Specialization for the DataStorage function of the residuals
     *
     * \param size: The size of the DataStorage object
     */
    template <>
    inline void DataStorage<std::vector<ResidualBase<hydraBase> *> >::zero(const unsigned int size) {
        throw std::runtime_error("Zeroing the ResidualBase pointer vector is not allowed");
    }

    namespace unit_test {
        class hydraBaseTester;  //!< Friend class for hydraBase for unit testing
    }  // namespace unit_test

    /*!
     * hydraBase: A base class which can be used to construct finite deformation material models.
     *
     * The hydra class seeks to provide utilities for the construction of finite deformation constitutive models
     * more rapidly than would be possible previously. The user can define as many different configurations as desired
     * and provide a calculation of the Cauchy stress.
     *
     * A non-linear problem which is of the size ( dimension**2 * num_configurations + num_ISVs ) will be solved.
     */
    class hydraBase : public CachingDataBase {
       public:
        // Constructors
        //! Default constructor for hydraBase
        hydraBase() : _configuration_unknown_count(0) {}

        //! Main constructor for objects of type hydraBase. Sets all quantities required for most solves.
        hydraBase(const floatType &time, const floatType &deltaTime, const floatType &temperature,
                  const floatType &previousTemperature, const secondOrderTensor &deformationGradient,
                  const secondOrderTensor &previousDeformationGradient, const floatVector &additionalDOF,
                  const floatVector &previousAdditionalDOF, const floatVector &previousStateVariables,
                  const floatVector &parameters, const unsigned int numConfigurations,
                  const unsigned int numNonLinearSolveStateVariables, const unsigned int dimension = 3,
                  const unsigned int configuration_unknown_count = 9, const floatType tolr = 1e-9,
                  const floatType tola = 1e-9);

        virtual void initialize();

        // User defined functions

        // Setter functions

        const void setCurrentResidualIndexMeaningful(const bool &value);

        const void setCurrentResidualIndex(const unsigned int value);

        // Getter functions
        //! Get a reference to the number of unknowns in each configuration
        constexpr unsigned int getConfigurationUnknownCount() { return _configuration_unknown_count; }

        //! Get a reference to the number of components of the stress
        constexpr unsigned int getStressSize() { return _stress_size; }

        //! Get a reference to the current time
        const floatType getTime() { return getScaledTime(); }

        //! Get a reference to the change in time
        const floatType getDeltaTime() { return getScaledDeltaTime(); }

        //! Get a reference to the current temperature
        const floatType getTemperature() { return getScaledTemperature(); };

        //! Get a reference to the previous temperature
        const floatType getPreviousTemperature() { return _previousTemperature; };

        //! Get a reference to the deformation gradient
        const secondOrderTensor *getDeformationGradient() { return getScaledDeformationGradient(); }

        //! Get a reference to the previous deformation gradient
        const secondOrderTensor *getPreviousDeformationGradient() { return &_previousDeformationGradient; }

        //! Get a reference to the additional degrees of freedom
        const floatVector *getAdditionalDOF() { return getScaledAdditionalDOF(); }

        //! Get a reference to the previous additional degrees of freedom
        const floatVector *getPreviousAdditionalDOF() { return &_previousAdditionalDOF; }

        //! Get a reference to the previous values of the state variables
        const floatVector *getPreviousStateVariables() { return &_previousStateVariables; }

        //! Get a reference to the model parameters
        const floatVector *getParameters() { return &_parameters; }

        //! Get a reference to the number of configurations
        constexpr unsigned int getNumConfigurations() { return _numConfigurations; }

        //! Get a reference to the number of state variables involved in the non-linear solve
        constexpr unsigned int getNumNonLinearSolveStateVariables() { return _numNonLinearSolveStateVariables; }

        virtual const unsigned int getNumUnknowns();

        virtual const unsigned int getNumAdditionalDOF();

        //! Get the value of the number of constraint equations
        virtual const unsigned int getNumConstraints();

        //! Get the current residual index
        const unsigned int getCurrentResidualIndex();

        const unsigned int getCurrentResidualOffset();

        //! Get the dimension
        constexpr unsigned int getDimension() { return _dimension; }

        //! Get a second order tensor's dimension
        constexpr unsigned int getSOTDimension() { return _dimension * _dimension; }

        //! Get a third order tensor's dimension
        constexpr unsigned int getTOTDimension() { return _dimension * _dimension * _dimension; }

        //! Get a fourth order tensor's dimension
        constexpr unsigned int getFOTDimension() { return _dimension * _dimension * _dimension * _dimension; }

        //! Get the relative tolerance
        constexpr floatType getRelativeTolerance() { return _tolr; }

        //! Get the absolute tolerance
        constexpr floatType getAbsoluteTolerance() { return _tola; }

        secondOrderTensor getSubConfiguration(const floatVector &configurations, const unsigned int &lowerIndex,
                                              const unsigned int &upperIndex);

        secondOrderTensor getSubConfigurationJacobian(const floatVector &configurations, const unsigned int &lowerIndex,
                                                      const unsigned int &upperIndex);

        secondOrderTensor getSubConfiguration(const unsigned int &lowerIndex, const unsigned int &upperIndex);

        secondOrderTensor getPrecedingConfiguration(const unsigned int &index);

        secondOrderTensor getFollowingConfiguration(const unsigned int &index);

        secondOrderTensor getConfiguration(const unsigned int &index);

        secondOrderTensor getPreviousSubConfiguration(const unsigned int &lowerIndex, const unsigned int &upperIndex);

        secondOrderTensor getPreviousPrecedingConfiguration(const unsigned int &index);

        secondOrderTensor getPreviousFollowingConfiguration(const unsigned int &index);

        secondOrderTensor getPreviousConfiguration(const unsigned int &index);

        floatVector getSubConfigurationJacobian(const unsigned int &lowerIndex, const unsigned int &upperIndex);

        floatVector getPrecedingConfigurationJacobian(const unsigned int &index);

        floatVector getFollowingConfigurationJacobian(const unsigned int &index);

        floatVector getPreviousSubConfigurationJacobian(const unsigned int &lowerIndex, const unsigned int &upperIndex);

        floatVector getPreviousPrecedingConfigurationJacobian(const unsigned int &index);

        floatVector getPreviousFollowingConfigurationJacobian(const unsigned int &index);

        virtual void setResidualClasses();

        void setResidualClasses(std::vector<ResidualBase<hydraBase> *> &residualClasses);

        std::vector<ResidualBase<hydraBase> *> *getResidualClasses();

        const floatVector *getResidual();

        const floatVector *getFlatJacobian();

        floatMatrix getJacobian();

        const floatVector *getFlatdRdF();

        floatMatrix getdRdF();

        const floatVector *getdRdT();

        const floatVector *getFlatdRdAdditionalDOF();

        floatMatrix getdRdAdditionalDOF();

        const floatVector *getFlatAdditionalDerivatives();

        floatMatrix getAdditionalDerivatives();

        const floatVector *getUnknownVector();

        const floatVector *getStress();

        const floatVector *getPreviousStress();

        const floatVector *getPreviouslyConvergedStress();

        void setViscoplasticDamping(const floatType &factor);

        void clearViscoplasticDamping();

        //! Get the value of the viscoplastic damping
        const floatType getViscoplasticDamping() {return _viscoplastic_damping_factor;}

        const bool getViscoplasticDampingSet();

        virtual void evaluate();

        virtual void computeTangents();

        virtual void computedXdAdditionalDOF();

        const floatVector *getFlatdXdF();

        const floatVector *getFlatdXdT();

        const floatVector *getFlatdXdAdditionalDOF();

        //! Add data to the vector of values which will be cleared after each iteration
        virtual void addIterationData(dataBase *data) override { _iterationData.push_back(data); }

        //! Add data to the vector of values which will be cleared after each non-linear step
        virtual void addNLStepData(dataBase *data) override { _nlStepData.push_back(data); }

        void setFailureVerbosityLevel(const unsigned int &value);

        //! Set the failure output string to use scientific notation
        void setFailureOutputScientific() { _failure_output << std::scientific; }

        //! Get the verbosity level for failure outputs
        const unsigned int getFailureVerbosityLevel() { return _failure_verbosity_level; }

        void addToFailureOutput(const std::string &value, bool add_endline = false);

        void addToFailureOutput(const floatVector &value, bool add_endline = true);

        void addToFailureOutput(const std::vector<bool> &value, bool add_endline = true);

        void addToFailureOutput(const floatType &value, bool add_endline = true);

        /*!
         * Add a general non-iterable object to the output string
         *
         * \param &v: The value to add
         * \param add_endline: Whether to append a newline character after the value
         */
        template <class v_type>
        void addToFailureOutput(const v_type &v, bool add_endline = true){
            _failure_output << v;
            if (add_endline) {
                _failure_output << "\n";
            }
        }

        /*!
         * Add a general iterable object to the output string
         *
         * \param &v_begin: The starting iterator of the value vector
         * \param &v_end: The stopping iterator of the value vector
         * \param add_endline: Whether to append a newline character after the value
         */
        template <class v_iterator>
        void addToFailureOutput(const v_iterator &v_begin, const v_iterator &v_end, bool add_endline = true) {
            for (auto v = v_begin; v != v_end; ++v) {
                _failure_output << *v << ", ";
            }
            if (add_endline) {
                _failure_output << "\n";
            }
        }

        //! Get the failure output string
        const std::string getFailureOutput() { return _failure_output.str(); }

        //! Get a scale factor for the deformation
        const floatType getScaleFactor() { return _scale_factor; }

        virtual void setScaleFactor(const floatType &value);

        //! Get the value of the current scaled time
        const floatType getScaledTime() { return _scaled_time; }

        //! Get the vlaue of the scaled change in time
        const floatType getScaledDeltaTime() {
            return _scaled_deltaTime;
        }

        //! Get the value of the scaled current temperature
        const floatType getScaledTemperature() {return _scaled_temperature;}

        //! Get the value of the scaled current deformation gradient
        const floatVector *
        getScaledDeformationGradient() {
            return &_scaled_deformationGradient;
        }

        //! Get the value of the scaled current additional DOF
        const floatVector *getScaledAdditionalDOF() {
            return &_scaled_additionalDOF;
        }

        floatVector *getMutableResidual();

        floatVector *getMutableJacobian();

        floatVector *getMutabledRdT();

        floatVector *getMutabledRdF();

        floatVector *getMutabledRdAdditionalDOF();

        const bool currentResidualIndexMeaningful() {
            /*!
             * Return if the current residual index is meaningful or not
             */
            return _current_residual_index_set;
        }

        std::string getResidualParameterizationInfo();

        void updateAdditionalStateVariables() {
            /*!
             * Update the additional state variable vector
             */

            for (auto v = std::begin(*getResidualClasses()); v != std::end(*getResidualClasses()); ++v) {
                (*v)->updateAdditionalStateVariables(_additionalStateVariables.second);
            }
        }

        void setToleranceScaleFactor(floatType factor) {
            /*!
             * Loosen the convergence tolerance for the next iteration
             * Useful if a Residual's form is changing in a non-smooth way
             *
             * \param factor: The scale factor to be applied to the current tolerance
             */

            if (factor > _residual_scale_factor) {
                _residual_scale_factor = factor;
            }
        }

        const floatType getToleranceScaleFactor() {
            /*!
             * Get the scale factor for the tolerance
             */

            return _residual_scale_factor;
        }

        SolverBase *solver =
            &_solver;  //!< The class which performs the material point solve TODO: Make this an incoming pointer

       protected:
        SubcyclerSolver _solver;  //!< Default solver

        // Setters that the user may need to access but not override

        const void resetToleranceScaleFactor() {
            /*!
             * Reset the tolerance scale factor to 1
             */

            _residual_scale_factor = 1.0;
        }

        void setStress(const floatVector &stress);

        // Utility functions
        virtual void computeConfigurations(const floatVector *data_vector, const unsigned int start_index,
                                           const floatVector &total_transformation, floatVector &configurations,
                                           floatVector &inverseConfigurations, const bool add_eye = false);

        virtual void extractStress();

        virtual void updateConfigurationsFromUnknownVector();

        virtual void decomposeUnknownVector();

        virtual void decomposeStateVariableVector();

        virtual void formNonLinearResidual();

        virtual void formNonLinearDerivatives();

        virtual void initializeUnknownVector();

        virtual void updateUnknownVector(const floatVector &newUnknownVector);

        virtual void calculateFirstConfigurationJacobians(const floatVector &configurations, fourthOrderTensor &dC1dC,
                                                          floatVector &dC1dCn);

        //! Update the scaled quantities
        virtual void setScaledQuantities();

        void unexpectedError() {
            /*!
             * Function to throw for an unexpected error. A user should never get here!
             */

            TARDIGRADE_ERROR_TOOLS_CATCH(
                throw std::runtime_error("You shouldn't have gotten here. If you aren't developing the code then "
                                         "contact a developer with the stack trace."))
        }

        //! A pass through function that does nothing
        void passThrough() {}

        void setX(const floatVector &X) {
            /*!
             * Set the value of the unknown vector
             *
             * \param &X: The unknown vector
             */

            _X.second = X;

            _X.first = true;
        }

        std::string build_upper_index_out_of_range_error_string(const unsigned int upperIndex,
                                                                const unsigned int num_configurations);
        std::string build_lower_index_out_of_range_error_string(const unsigned int lowerIndex,
                                                                const unsigned int upperIndex);

        virtual tardigradeHydra::hydraBase::SetDataStorageIteration<secondOrderTensor> get_SetDataStorage_stress();

        virtual void setConstraints();

        virtual void setConstraintJacobians();

        virtual void resetProblem();

        void setAllowModifyGlobalResidual(const bool value) {
            /*!
             * Set a flag for if the global residual can be modified
             *
             * \param value: The updated value
             */
            _allow_modify_global_residual = value;
        }

        void setAllowModifyGlobalJacobian(const bool value) {
            /*!
             * Set a flag for if the global jacobian can be modified
             *
             * \param value: The updated value
             */
            _allow_modify_global_jacobian = value;
        }

        void setAllowModifyGlobaldRdT(const bool value) {
            /*!
             * Set a flag for if the global dRdT can be modified
             *
             * \param value: The updated value
             */
            _allow_modify_global_dRdT = value;
        }

        void setAllowModifyGlobaldRdF(const bool value) {
            /*!
             * Set a flag for if the global dRdF can be modified
             *
             * \param value: The updated value
             */
            _allow_modify_global_dRdF = value;
        }

        void setAllowModifyGlobaldRdAdditionalDOF(const bool value) {
            /*!
             * Set a flag for if the global dRdAdditionalDOF can be modified
             *
             * \param value: The updated value
             */
            _allow_modify_global_dRdAdditionalDOF = value;
        }

        void setPreviouslyConvergedStress(const floatVector &value) {
            /*!
             * Set the value of the previously converged stress.
             *
             * \param &value: The incoming value
             */

            _previouslyConvergedStress.second = value;
            _previouslyConvergedStress.first  = true;
        }

       private:
        // Friend classes
        friend class tardigradeHydra::SolverStepBase;  //!< The base class for the solver steps
        friend class tardigradeHydra::SolverBase;      //!< The base class for the solver
        friend class unit_test::hydraBaseTester;  //!< Friend class which allows modification of private variables. ONLY
                                                  //!< TO BE USED FOR TESTING!

        unsigned int _dimension;  //!< The spatial dimension of the problem

        unsigned int
            _configuration_unknown_count;  //!< The number of unknowns required for a configuration. Used to ensure that
                                           //!< the unknown and state variable vectors are the right size. Must be set
                                           //!< by all inheriting classes. For 3D classical continuum this will be 9,
                                           //!< for higher order theories this will change.

        unsigned int _stress_size;  //!< The number of terms in the stress measures. For 3D classical continuum this
                                    //!< will be 9, for higher order theories this will change.

        floatType _time;  //!< The current time

        floatType _deltaTime;  //!< The change in time

        floatType _temperature;  //!< The current temperature

        floatType _previousTemperature;  //!< The previous temperature

        secondOrderTensor _deformationGradient;  //!< The current deformation gradient

        secondOrderTensor _previousDeformationGradient;  //!< The previous deformation gradient

        floatVector _additionalDOF;  //!< The current additional degrees of freedom

        floatVector _previousAdditionalDOF;  //!< The previous additional degrees of freedom

        floatVector _previousStateVariables;  //!< The previous state variables

        floatVector _parameters;  //!< The model parameters

        floatType _scaled_time;  //!< The current time scaled by the scaling factor

        floatType _scaled_deltaTime;  //!< The change in time scaled by the scaling factor

        floatType _scaled_temperature;  //!< The current temperature scaled by the scaling factor

        secondOrderTensor
            _scaled_deformationGradient;  //!< The current deformation gradient scaled by the scaling factor

        floatVector _scaled_additionalDOF;  //!< The current additional degrees of freedom scaled by the scaling factor

        unsigned int _numConfigurations;  //!< The number of configurations

        unsigned int _numNonLinearSolveStateVariables;  //!< The number of state variables which will be solved in the
                                                        //!< Newton-Raphson loop

        floatType _residual_scale_factor =
            1;  //!< A scale factor for the residual which can be used to loosen the tolerance by a Residual

        floatType _tolr;  //!< The relative tolerance

        floatType _tola;  //!< The absolute tolerance

        bool _rank_deficient_error = true;  //!< Flag for whether a rank-deficient LHS will throw a convergence error

        std::vector<dataBase *>
            _iterationData;  //!< A vector of pointers to data which should be cleared at each iteration

        std::vector<dataBase *>
            _nlStepData;  //!< A vector of pointers to data which should be cleared after each nonlinear step

        DataStorage<std::vector<ResidualBase<hydraBase> *> >
            _residualClasses;  //!< A vector of classes which compute the terms in the residual equation

        DataStorage<floatVector> _residual;  //!< The residual vector for the global solve

        DataStorage<floatVector> _jacobian;  //!< The jacobian matrix in row-major form for the global solve

        DataStorage<floatVector> _previouslyConvergedStress;  //!< The previously converged stress

        DataStorage<floatVector> _dRdF;  //!< The gradient of the residual w.r.t. the deformation gradient in row-major
                                         //!< form for the global solve

        DataStorage<floatVector> _dRdT;  //!< The gradient of the residual w.r.t. the temperature for the global solve

        DataStorage<floatVector>
            _dRdAdditionalDOF;  //!< The derivatives of the residual w.r.t. the additional degrees of freedom

        DataStorage<floatVector> _additionalDerivatives;  //!< Additional derivatives of the residual

        DataStorage<floatVector> _X;  //!< The unknown vector { stress, F1, ..., Fn, xi1, ..., xim }

        DataStorage<floatVector>
            _stress;  //!< The stress in the current configuration as determined from the current state

        DataStorage<floatVector> _previousStress;  //!< The previous value of the stress in the current configuration as
                                                   //!< determined from the previous state

        DataStorage<floatVector>
            _flatdXdF;  //!< The total derivative of the unknown vector w.r.t. the deformation in row-major form

        DataStorage<floatVector> _flatdXdT;  //!< The total derivative of the unknown vector w.r.t. the temperature

        DataStorage<floatVector>
            _flatdXdAdditionalDOF;  //!< The total derivative of the unknown vector w.r.t. the additional DOF

        void setFirstConfigurationJacobians();

        void setPreviousFirstConfigurationJacobians();

        void resetIterationData();

        void resetNLStepData();

        unsigned int _failure_verbosity_level = 0;  //!< The verbosity level for failure.

        std::stringstream _failure_output;  //!< Additional failure output information

        floatType _scale_factor =
            1.0;  //!< A scale factor applied to the incoming loading (deformation, temperature, etc.)

        bool _allow_modify_global_residual = false;  //!< Flag for if the global residual can be modified

        bool _allow_modify_global_jacobian = false;  //!< Flag for if the global jacobian can be modified

        bool _allow_modify_global_dRdT = false;  //!< Flag for if the global dRdT can be modified

        bool _allow_modify_global_dRdF = false;  //!< Flag for if the global dRdF can be modified

        bool _allow_modify_global_dRdAdditionalDOF =
            false;  //!< Flag for if the global dRdAdditionalDOF can be modified

        bool _current_residual_index_set = false;  //!< Flag for whether the current residual index has been set

        int _current_residual_index = 0;  //!< The current residual index

        floatType _viscoplastic_damping_factor = 0;  //!< The fraction of the difference between the trial stress and
                                                     //!< the stress that will be suppressed to assist in convergence

        bool _viscoplastic_damping_set = false;  //!< Flag for whether the viscoplastic damping factor has been set

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, configurations, floatVector, passThrough)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousConfigurations, floatVector, passThrough)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, inverseConfigurations, floatVector, passThrough)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousInverseConfigurations, floatVector, passThrough)

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, nonLinearSolveStateVariables, floatVector, passThrough)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousNonLinearSolveStateVariables, floatVector,
                                                  passThrough)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, additionalStateVariables, floatVector, passThrough)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousAdditionalStateVariables, floatVector, passThrough)

        TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE(private, setConstraints, getConstraints, constraints,
                                                         floatVector, setConstraints)

        TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE(private, setConstraintJacobians, getConstraintJacobians,
                                                         constraintJacobians, floatVector, setConstraintJacobians)

        TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE(private, set_dF1dF, get_dF1dF, dF1dF, secondOrderTensor,
                                                         setFirstConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE(private, set_dF1dFn, get_dF1dFn, dF1dFn, floatVector,
                                                         setFirstConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_NAMED_PREVIOUS_STORAGE(private, set_previousdF1dF, get_previousdF1dF, previousdF1dF,
                                                        secondOrderTensor, setPreviousFirstConfigurationJacobians)

        TARDIGRADE_HYDRA_DECLARE_NAMED_PREVIOUS_STORAGE(private, set_previousdF1dFn, get_previousdF1dFn, previousdF1dFn,
                                                        floatVector, setPreviousFirstConfigurationJacobians)
    };

}  // namespace tardigradeHydra

#endif
