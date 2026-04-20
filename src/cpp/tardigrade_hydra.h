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
#include "tardigrade_DeformationBase.h"
#include "tardigrade_MatrixMap.h"
#include "tardigrade_ResidualBase.h"
#include "tardigrade_SetDataStorage.h"
#include "tardigrade_SolverBase.h"
#include "tardigrade_SubcyclerSolver.h"

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
        class HydraBaseTester;  //!< Friend class for HydraBase for unit testing
        class hydraBaseTester;  //!< Friend class for hydraBase for unit testing
    }  // namespace unit_test

    /*!
     * HydraConfigurationBase: A base class which defines the hydra configuration
     */
    class HydraConfigurationBase {
       public:
        //! The default constructor for HydraConfigurationBase
        HydraConfigurationBase() {}
        //! The number of unknowns in each configuration
        unsigned int configuration_unknown_count = 0;
        //! The relative tolerance
        floatType tolr                           = 1e-9;
        //! The absolute tolerance
        floatType tola                           = 1e-9;

       private:
        friend class hydraBase;
    };

    /*!
     * HydraClassicalConfiguration: A class which defines a classical deformation problem
     */
    class HydraClassicalConfiguration : public HydraConfigurationBase {
       public:
        HydraClassicalConfiguration() { configuration_unknown_count = 9; }
    };

    /*!
     * DOFStorageBase: A class which stores the degrees of freedom
     */
    class DOFStorageBase {
       public:
        /*!
         * Default constructor
         */
        DOFStorageBase()
            : _time(0),
              _deltaTime(0),
              _temperature(0),
              _previous_temperature(0),
              _deformation_gradient(floatVector(0, 0)),
              _previous_deformation_gradient(floatVector(0, 0)),
              _additional_dof(floatVector(0, 0)),
              _previous_additional_dof(floatVector(0, 0)) {}

        /*!
         * Constructor which sets the dof information
         *
         * TODO: We're eventually going to store the degrees of freedom into a single storage
         *       array but, for now, I'm leaving them as discrete just to help with the
         *       transition.
         *
         * \param &time: The current time
         * \param &deltaTime: The change in time from the previous time
         * \param &temperature: The current temperature
         * \param &previous_temperature: The previous temperature
         * \param &deformation_gradient: The deformation gradient
         * \param &previous_deformation_gradient: The previous deformation gradient
         * \param &additional_dof: The additional degrees of freedom
         * \param &previous_additional_dof: The previous additional degrees of freedom
         */
        DOFStorageBase(const floatType &time, const floatType &deltaTime, const floatType &temperature,
                       const floatType &previous_temperature, const floatVector &deformation_gradient,
                       const floatVector &previous_deformation_gradient, const floatVector &additional_dof,
                       const floatVector &previous_additional_dof)
            : _time(time),
              _deltaTime(deltaTime),
              _temperature(temperature),
              _previous_temperature(previous_temperature),
              _deformation_gradient(deformation_gradient),
              _previous_deformation_gradient(previous_deformation_gradient),
              _additional_dof(additional_dof),
              _previous_additional_dof(previous_additional_dof) {}

        //! The current time
        const floatType _time;

        //! The change in time from the previous timestep
        const floatType _deltaTime;

        //! The current temperature
        const floatType _temperature;

        //! The previous temperature
        const floatType _previous_temperature;

        //! The deformation gradient
        const floatVector _deformation_gradient;

        //! The previous deformation gradient
        const floatVector _previous_deformation_gradient;

        //! The additional degrees of freedom
        const floatVector _additional_dof;

        //! The previous additional degrees of freedom
        const floatVector _previous_additional_dof;

       protected:
    };

    /*!
     * A class which stores the model configuration
     */
    class ModelConfigurationBase {
       public:
        /*!
         * The default constructor
         */
        ModelConfigurationBase()
            : _previous_state_variables(floatVector(0, 0)),
              _parameters(floatVector(0, 0)),
              _num_configurations(0),
              _num_nonlinear_solve_state_variables(0) {}

        /*!
         * The constructor
         *
         * \param &previous_state_variables: The previous values of the internal state variables
         * \param &parameters: The model parameters
         * \param &num_configurations: The number of configurations
         * \param &num_nonlinear_solve_state_variables: The number of state variables required for the nonlinear solve
         */
        ModelConfigurationBase(const floatVector &previous_state_variables, const floatVector &parameters,
                               const unsigned int  num_configurations,
                               const unsigned int &num_nonlinear_solve_state_variables)
            : _previous_state_variables(previous_state_variables),
              _parameters(parameters),
              _num_configurations(num_configurations),
              _num_nonlinear_solve_state_variables(num_nonlinear_solve_state_variables) {}

        //! The previous values of the state variables
        const floatVector _previous_state_variables;

        //! The model parameters
        const floatVector _parameters;

        //! The number of configurations TODO: Maybe should be in the deformation?
        const unsigned int _num_configurations;

        //! Added the number of nonlinear solve state variables
        const unsigned int _num_nonlinear_solve_state_variables;
    };

    class HydraClassical3DConfiguration {
        public:
        protected:
        private:
    };

    /*!
     * HydraBase: A base class template which can be used to construct finite deformation material models.
     *
     * The hydra class seeks to provide utilities for the construction of finite deformation constitutive models
     * more rapidly than would be possible previously. The user can define as many different configurations as desired
     * and provide a calculation of the Cauchy stress.
     *
     * A non-linear problem which is of the size ( dimension**2 * num_configurations + num_ISVs ) will be solved.
     */
    template<class configuration>
    class HydraBase : public CachingDataBase { };

    /*!
     * A base class for 3D Classical continuum
     */
    class hydraBase : public HydraBase<HydraClassical3DConfiguration> {
       public:
        // Constructors
        //! Default constructor for hydraBase
        hydraBase() {}

        //! Main constructor for objects of type hydraBase. Sets all quantities required for most solves.
        hydraBase(const DOFStorageBase &DOFStorage, const ModelConfigurationBase &ModelConfiguration,
                  HydraConfigurationBase _hydra_configuration = HydraClassicalConfiguration());

        virtual void initialize();

        // User defined functions

        // Setter functions

        const void setCurrentResidualIndexMeaningful(const bool &value);

        const void setCurrentResidualIndex(const unsigned int value);

        // Getter functions
        //! Get a reference to the number of unknowns in each configuration
        constexpr unsigned int getConfigurationUnknownCount() {
            return hydra_configuration.configuration_unknown_count;
        }

        //! Get a reference to the number of components of the stress
        constexpr unsigned int getStressSize() { return _stress_size; }

        //! Get a reference to the current time
        const floatType getTime() { return getScaledTime(); }

        //! Get a reference to the change in time
        const floatType getDeltaTime() { return getScaledDeltaTime(); }

        //! Get a reference to the current temperature
        const floatType getTemperature() { return getScaledTemperature(); };

        //! Get a reference to the previous temperature
        const floatType getPreviousTemperature() { return dof->_previous_temperature; };

        //! Get a reference to the deformation gradient
        const secondOrderTensor *getDeformationGradient() { return getScaledDeformationGradient(); }

        //! Get a reference to the previous deformation gradient
        const secondOrderTensor *getPreviousDeformationGradient() { return &dof->_previous_deformation_gradient; }

        //! Get a reference to the additional degrees of freedom
        const floatVector *getAdditionalDOF() { return getScaledAdditionalDOF(); }

        //! Get a reference to the previous additional degrees of freedom
        const floatVector *getPreviousAdditionalDOF() { return &dof->_previous_additional_dof; }

        //! Get a reference to the previous values of the state variables
        const floatVector *getPreviousStateVariables() { return &model_configuration->_previous_state_variables; }

        //! Get a reference to the model parameters
        const floatVector *getParameters() { return &model_configuration->_parameters; }

        //! Get a reference to the number of configurations
        constexpr unsigned int getNumConfigurations() { return model_configuration->_num_configurations; }

        //! Get a reference to the number of state variables involved in the non-linear solve
        constexpr unsigned int getNumNonLinearSolveStateVariables() {
            return model_configuration->_num_nonlinear_solve_state_variables;
        }

        virtual const unsigned int getNumUnknowns();

        virtual const unsigned int getNumAdditionalDOF();

        //! Get the value of the number of constraint equations
        virtual const unsigned int getNumConstraints();

        //! Get the current residual index
        const unsigned int getCurrentResidualIndex();

        const unsigned int getCurrentResidualOffset();

        //! Get the dimension
        constexpr unsigned int getDimension() { return deformation->dimension; }

        //! Get a second order tensor's dimension
        constexpr unsigned int getSOTDimension() { return deformation->dimension * deformation->dimension; }

        //! Get a third order tensor's dimension
        constexpr unsigned int getTOTDimension() {
            return deformation->dimension * deformation->dimension * deformation->dimension;
        }

        //! Get a fourth order tensor's dimension
        constexpr unsigned int getFOTDimension() {
            return deformation->dimension * deformation->dimension * deformation->dimension * deformation->dimension;
        }

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
        const floatType getViscoplasticDamping() { return _viscoplastic_damping_factor; }

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

        template <class v_type>
        void addToFailureOutput(const v_type &v, bool add_endline = true);

        template <class v_iterator>
        void addToFailureOutput(const v_iterator &v_begin, const v_iterator &v_end, bool add_endline = true);

        //! Get the failure output string
        const std::string getFailureOutput() { return _failure_output.str(); }

        //! Get a scale factor for the deformation
        const floatType getScaleFactor() { return _scale_factor; }

        virtual void setScaleFactor(const floatType &value);

        //! Get the value of the current scaled time
        const floatType getScaledTime() { return _scaled_time; }

        //! Get the vlaue of the scaled change in time
        const floatType getScaledDeltaTime() { return _scaled_deltaTime; }

        //! Get the value of the scaled current temperature
        const floatType getScaledTemperature() { return _scaled_temperature; }

        //! Get the value of the scaled current deformation gradient
        const floatVector *getScaledDeformationGradient() { return &_scaled_deformationGradient; }

        //! Get the value of the scaled current additional DOF
        const floatVector *getScaledAdditionalDOF() { return &_scaled_additionalDOF; }

        floatVector *getMutableResidual();

        floatVector *getMutableJacobian();

        floatVector *getMutabledRdT();

        floatVector *getMutabledRdF();

        floatVector *getMutabledRdAdditionalDOF();

        const bool currentResidualIndexMeaningful();

        std::string getResidualParameterizationInfo();

        void updateAdditionalStateVariables();

        void setToleranceScaleFactor(floatType factor);

        //! Get the scale factor for the tolerance
        const floatType getToleranceScaleFactor() { return _residual_scale_factor; }

        //! The class which performs the material point solve TODO: Make this an incoming pointer
        SolverBase *solver = &_solver;

        //! The class which contains the deformation
        DeformationBase *deformation = &_deformation;

        //! The class which defines the hydra configuration
        HydraConfigurationBase hydra_configuration;

        //! The class which stores the degrees of freedom
        const DOFStorageBase *dof;

        //! The class which stores the model configuration
        const ModelConfigurationBase *model_configuration;

       protected:
        //! Default solver
        SubcyclerSolver _solver;

        //! Default deformation class
        DeformationBase _deformation = DeformationBase(this);

        // Setters that the user may need to access but not override

        //! Reset the tolerance scale factor to 1
        const void resetToleranceScaleFactor() { _residual_scale_factor = 1.0; }

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

        virtual void setScaledQuantities();

        void unexpectedError();

        //! A pass through function that does nothing
        void passThrough() {}

        void setX(const floatVector &X);

        std::string build_upper_index_out_of_range_error_string(const unsigned int upperIndex,
                                                                const unsigned int num_configurations);

        std::string build_lower_index_out_of_range_error_string(const unsigned int lowerIndex,
                                                                const unsigned int upperIndex);

        virtual tardigradeHydra::hydraBase::SetDataStorageIteration<secondOrderTensor> get_SetDataStorage_stress();

        virtual void setConstraints();

        virtual void setConstraintJacobians();

        virtual void resetProblem();

        void setAllowModifyGlobalResidual(const bool value);

        void setAllowModifyGlobalJacobian(const bool value);

        void setAllowModifyGlobaldRdT(const bool value);

        void setAllowModifyGlobaldRdF(const bool value);

        void setAllowModifyGlobaldRdAdditionalDOF(const bool value);

        void setPreviouslyConvergedStress(const floatVector &value);

       private:
        // Friend classes
        friend class tardigradeHydra::SolverBase;  //!< The base class for the solver
        friend class unit_test::hydraBaseTester;  //!< Friend class which allows modification of private variables. ONLY
                                                  //!< TO BE USED FOR TESTING!

        //! The number of terms in the stress measures
        unsigned int _stress_size;

        //! The current time scaled by the scaling factor
        floatType _scaled_time;

        //! The change in time scaled by the scaling factor
        floatType _scaled_deltaTime;

        //! The current temperature scaled by the scaling factor
        floatType _scaled_temperature;

        //! The current deformation gradient scaled by the scaling factor
        secondOrderTensor _scaled_deformationGradient;

        //! The current additional degrees of freedom scaled by the scaling factor
        floatVector _scaled_additionalDOF;

        //! A scale factor for the residual which can be used to loosen the tolerance
        floatType _residual_scale_factor = 1;

        //! A vector of pointers to data which should be cleared at each iteration
        std::vector<dataBase *> _iterationData;

        //! A vector of pointers to data which should be cleared after each nonlinear step
        std::vector<dataBase *> _nlStepData;

        //! A vector of classes which compute the terms in the residual equation
        DataStorage<std::vector<ResidualBase<hydraBase> *> > _residualClasses;

        //! The residual vector for the global solve
        DataStorage<floatVector> _residual;

        //! The jacobian matrix in row-major form for the global solve
        DataStorage<floatVector> _jacobian;

        //! The previously converged stress
        DataStorage<floatVector> _previouslyConvergedStress;

        //! The gradient of the residual w.r.t. the deformation gradient in row-major form
        DataStorage<floatVector> _dRdF;

        //! The gradient of the residual w.r.t. the temperature for the global solve
        DataStorage<floatVector> _dRdT;

        //! The derivatives of the residual w.r.t. the additional degrees of freedom
        DataStorage<floatVector> _dRdAdditionalDOF;

        //! Additional derivatives of the residual
        DataStorage<floatVector> _additionalDerivatives;

        //! The unknown vector { stress, F1, ..., Fn, xi1, ..., xim }
        DataStorage<floatVector> _X;

        //! The stress in the current configuration as determined from the current state
        DataStorage<floatVector> _stress;

        //! The previous value of the stress in the current configuration as determined from the previous state
        DataStorage<floatVector> _previousStress;

        //! The total derivative of the unknown vector w.r.t. the deformation in row-major form
        DataStorage<floatVector> _flatdXdF;

        //! The total derivative of the unknown vector w.r.t. the temperature
        DataStorage<floatVector> _flatdXdT;

        //! The total derivative of the unknown vector w.r.t. the additional dof
        DataStorage<floatVector> _flatdXdAdditionalDOF;

        void resetIterationData();

        void resetNLStepData();

        //! The verbosity level for failure
        unsigned int _failure_verbosity_level = 0;

        //! Additional failure output information
        std::stringstream _failure_output;

        //! A scale factor applied to the incoming loading (deformation, temperature, etc.)
        floatType _scale_factor = 1.0;

        //! Flag for if the global residual can be modified
        bool _allow_modify_global_residual = false;

        //! Flag for if the global jacobian can be modified
        bool _allow_modify_global_jacobian = false;

        //! Flag for if the global dRdT can be modified
        bool _allow_modify_global_dRdT = false;

        //! Flag for if the global dRdF can be modified
        bool _allow_modify_global_dRdF = false;

        //! Flag for if the global dRdAdditionalDOF can be modified
        bool _allow_modify_global_dRdAdditionalDOF = false;

        //! Flag for whether the current residual index has been set
        bool _current_residual_index_set = false;

        //! The current residual index
        int _current_residual_index = 0;

        //! The fraction of the difference between the trial and calculated stress that will be suppressed
        floatType _viscoplastic_damping_factor = 0;

        //! Flag for whether the viscoplastic damping factor has been set
        bool _viscoplastic_damping_set = false;

        TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE(private, nonLinearSolveStateVariables, floatVector, passThrough)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousNonLinearSolveStateVariables, floatVector,
                                                  passThrough)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, additionalStateVariables, floatVector, passThrough)

        TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(private, previousAdditionalStateVariables, floatVector, passThrough)

        TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE(private, setConstraints, getConstraints, constraints,
                                                         floatVector, setConstraints)

        TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE(private, setConstraintJacobians, getConstraintJacobians,
                                                         constraintJacobians, floatVector, setConstraintJacobians)
    };

}  // namespace tardigradeHydra

#include "tardigrade_hydra.tpp"

#endif
