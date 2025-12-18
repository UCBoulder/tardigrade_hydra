/**
  ******************************************************************************
  * \file tardigrade_hydra.h
  ******************************************************************************
  * A C++ library for constructing finite deformation constitutive models.
  ******************************************************************************
  */

#ifndef TARDIGRADE_HYDRA_H
#define TARDIGRADE_HYDRA_H

#include<sstream>
#include<functional>

#include"tardigrade_error_tools.h"
//!We will use the functions that depend on Eigen
#define USE_EIGEN
#include"tardigrade_vector_tools.h"

#ifdef TARDIGRADE_HYDRA_USE_LLXSMM
#include<libxsmm.h>
#endif

#include"tardigrade_CoreDefinitions.h"
#include"tardigrade_CustomErrors.h"
#include"tardigrade_SetDataStorage.h"
#include"tardigrade_MatrixMap.h"
#include"tardigrade_SolverBase.h"
#include"tardigrade_ResidualBase.h"
//#include"tardigrade_SolverStepBase.h"

namespace tardigradeHydra{

    template <>
    inline void DataStorage< std::vector< ResidualBase<hydraBase>* > >::zero( ){
                /*!
                 * The function to set the value to zero
                 */

        throw std::runtime_error( "Zeroing the ResidualBase pointer vector is not allowed" );

    }

    template <>
    inline void DataStorage< std::vector< ResidualBase<hydraBase>* > >::zero( const unsigned int size ){
                /*!
                 * The function to set the value to zero
                 */

        throw std::runtime_error( "Zeroing the ResidualBase pointer vector is not allowed" );

    }

    /*!
     * hydraBase: A base class which can be used to construct finite deformation material models.
     * 
     * The hydra class seeks to provide utilities for the construction of finite deformation constitutive models
     * more rapidly than would be possible previously. The user can define as many different configurations as desired
     * and provide a calculation of the Cauchy stress.
     * 
     * A non-linear problem which is of the size ( dimension**2 * num_configurations + num_ISVs ) will be solved.
     */
    class hydraBase{

        public:

            // Constructors
            //! Default constructor for hydraBase
            hydraBase( ) : _configuration_unknown_count( 0 ){ }

            //! Main constructor for objects of type hydraBase. Sets all quantities required for most solves.
            hydraBase( const floatType &time, const floatType &deltaTime,
                       const floatType &temperature, const floatType &previousTemperature,
                       const secondOrderTensor &deformationGradient, const secondOrderTensor &previousDeformationGradient,
                       const floatVector &additionalDOF, const floatVector &previousAdditionalDOF,
                       const floatVector &previousStateVariables, const floatVector &parameters,
                       const unsigned int numConfigurations, const unsigned int numNonLinearSolveStateVariables,
                       const unsigned int dimension=3, const unsigned int configuration_unknown_count=9,
                       const floatType tolr=1e-9, const floatType tola=1e-9,
                       const unsigned int maxIterations=20, const unsigned int maxLSIterations=5, const floatType lsAlpha=1e-4,
                       const bool use_preconditioner=false, const unsigned int preconditioner_type=0 );

            // User defined functions

            // Setter functions

            // Getter functions
            //! Get a reference to the number of unknowns in each configuration
            constexpr unsigned int getConfigurationUnknownCount( ){ return _configuration_unknown_count; }

            //! Get a reference to the number of components of the stress
            constexpr unsigned int getStressSize( ){ return _stress_size; }

            //! Get a reference to the current time
            const floatType getTime( ){ return getScaledTime( ); }

            //! Get a reference to the change in time
            const floatType getDeltaTime( ){ return getScaledDeltaTime( ); }

            //! Get a reference to the current temperature
            const floatType getTemperature( ){ return getScaledTemperature( ); };

            //! Get a reference to the previous temperature
            const floatType getPreviousTemperature( ){ return _previousTemperature; };

            //! Get a reference to the deformation gradient
            const secondOrderTensor* getDeformationGradient( ){ return getScaledDeformationGradient( ); }

            //! Get a reference to the previous deformation gradient
            const secondOrderTensor* getPreviousDeformationGradient( ){ return &_previousDeformationGradient; }

            //! Get a reference to the additional degrees of freedom
            const floatVector* getAdditionalDOF( ){ return getScaledAdditionalDOF( ); }

            //! Get a reference to the previous additional degrees of freedom
            const floatVector* getPreviousAdditionalDOF( ){ return &_previousAdditionalDOF; }

            //! Get a reference to the previous values of the state variables
            const floatVector* getPreviousStateVariables( ){ return &_previousStateVariables; }

            //! Get a reference to the model parameters
            const floatVector* getParameters( ){ return &_parameters; }

            //! Get a reference to the number of configurations
            constexpr unsigned int getNumConfigurations( ){ return _numConfigurations; }

            //! Get a reference to the number of state variables involved in the non-linear solve
            constexpr unsigned int getNumNonLinearSolveStateVariables( ){ return _numNonLinearSolveStateVariables; }

            //! Get the number of terms in the unknown vector
            virtual const unsigned int getNumUnknowns( ){ return getNumConfigurations( ) * getConfigurationUnknownCount( ) + getNumNonLinearSolveStateVariables( ); }

            //! Get the number of additional degrees of freedom
            virtual const unsigned int getNumAdditionalDOF( ){ return getAdditionalDOF( )->size( ); }

            //! Get the value of the number of constraint equations
            virtual const unsigned int getNumConstraints( );

            //! Get the current residual index
            const unsigned int getCurrentResidualIndex( ){

                TARDIGRADE_ERROR_TOOLS_CHECK( currentResidualIndexMeaningful( ), "The current residual index isn't meaningful" );
                return _current_residual_index;

            }

            const unsigned int getCurrentResidualOffset( );

            //! Get the dimension
            constexpr unsigned int getDimension( ){ return _dimension; }

            //! Get a second order tensor's dimension
            constexpr unsigned int getSOTDimension( ){ return _dimension * _dimension; }

            //! Get a third order tensor's dimension
            constexpr unsigned int getTOTDimension( ){ return _dimension * _dimension * _dimension; }

            //! Get a fourth order tensor's dimension
            constexpr unsigned int getFOTDimension( ){ return _dimension * _dimension * _dimension * _dimension; }

            //! Get the relative tolerance
            constexpr floatType getRelativeTolerance( ){ return _tolr; }

            //! Get the absolute tolerance
            constexpr floatType getAbsoluteTolerance( ){ return _tola; }

            //! Get the line-search alpha
            constexpr floatType getLSAlpha( ){ return _lsAlpha; }

            //! Get whether to use a preconditioner
            constexpr bool getUsePreconditioner( ){ return _use_preconditioner; }

            //! Get the preconditioner type
            constexpr unsigned int getPreconditionerType( ){ return _preconditioner_type; }

            //! Get whether the preconditioner is diagonal or not
            constexpr bool getPreconditionerIsDiagonal( ){ return _preconditioner_is_diagonal; }

            //!< Get the gradient descent sigma parameter
            const floatType getGradientSigma( ){ return _gradientSigma; }

            //!< Get the gradient descent beta parameter
            const floatType getGradientBeta( ){ return _gradientBeta; }

            //!< Get the max allowable number of gradient iterations
            const unsigned int getMaxGradientIterations( ){ return _maxGradientIterations; }

            //!< Get the gradient descent rho parameter
            const floatType getGradientRho( ){ return _gradientRho; }

            //!< Get the gradient descent p parameter
            const floatType getGradientP( ){ return _gradientP; }

            //!< Get the Newton step should be a LevenbergMarquardt step
            const bool getUseLevenbergMarquardt( ){ return _use_LM_step; }

            //!< Get a reference to whether the Newton step should be a relaxed solve
            const bool getUseRelaxedSolve( ){ return _use_relaxed_solve; }

            //!< Get a reference to whether Gradient descent is allowed
            const bool getUseGradientDescent( ){ return _use_gradient_descent; }

            //!< Get a reference to the flag for whether to throw an error if the LHS matrix is rank-deficient
            const bool getRankDeficientError( ){ return _rank_deficient_error; }

            //!< Set the gradient descent sigma parameter
            void setGradientSigma( const floatType &value ){
               /*!
                * Set the value of the sigma parameter for gradient descent steps
                *
                * \param &value: The value of the parameter
                */

                _gradientSigma = value;

            }

            //!< Set the gradient descent beta parameter
            void setGradientBeta( const floatType &value ){
               /*!
                * Set the value of the beta parameter for gradient descent steps
                *
                * \param &value: The value of the parameter
                */
 
                _gradientBeta = value;

            }

            //!< Set the max allowable number of gradient iterations
            void setMaxGradientIterations( const unsigned int &value ){
               /*!
                * Set the value of the maximum number of iterations for gradient descent steps
                *
                * \param &value: The value of the parameter
                */

                _maxGradientIterations = value;

            }

            //!< Set the gradient descent rho parameter
            void setGradientRho( const floatType &value ){
               /*!
                * Set the value of the rho parameter for gradient descent steps
                *
                * \param &value: The value of the parameter
                */
 
                _gradientRho = value;

            }

            //!< Set the gradient descent p parameter
            void setGradientP( const floatType &value ){
               /*!
                * Set the value of the p parameter for gradient descent steps
                *
                * \param &value: The value of the parameter
                */
 
                _gradientP = value;

            }

            //!< Set whether to attempt a Levenberg-Marquardt step
            void setUseLevenbergMarquardt( const bool &value ){
                /*!
                 * Set whether to attempt a Levenberg-Marquardt step
                 * 
                 * \param &value: The value of the parameter
                 */

                setUseGradientDescent( value );

                _use_LM_step = value;

            }

            //!< Set whether to attempt a Relaxed-solve
            void setUseRelaxedSolve( const bool &value ){
                /*!
                 * Set whether to attempt a relaxed solve
                 * 
                 * \param &value: The value of the parameter
                 */

                _use_relaxed_solve = value;

            }

            //!< Set whether to use gradient descent
            void setUseGradientDescent( const bool &value ){
                /*!
                 * Set whether to attempt a gradient descent step
                 * 
                 * \param &value: The value of the parameter
                 */

                _use_gradient_descent = value;

            }

            //!< Set whether rank deficiency is a reason to throw an error
            void setRankDeficientError( const bool &value ){
                /*!
                 * Set whether a rank-deficient LHS will cause an error
                 * 
                 * \param &value: The value of the parameter
                 */

                _rank_deficient_error = value;

            }

            secondOrderTensor getSubConfiguration( const floatVector &configurations, const unsigned int &lowerIndex, const unsigned int &upperIndex );

            secondOrderTensor getSubConfigurationJacobian( const floatVector &configurations, const unsigned int &lowerIndex, const unsigned int &upperIndex );

            secondOrderTensor getSubConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            secondOrderTensor getPrecedingConfiguration( const unsigned int &index );

            secondOrderTensor getFollowingConfiguration( const unsigned int &index );

            secondOrderTensor getConfiguration( const unsigned int &index );

            secondOrderTensor getPreviousSubConfiguration( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            secondOrderTensor getPreviousPrecedingConfiguration( const unsigned int &index );

            secondOrderTensor getPreviousFollowingConfiguration( const unsigned int &index );

            secondOrderTensor getPreviousConfiguration( const unsigned int &index );

            floatVector getSubConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getPrecedingConfigurationJacobian( const unsigned int &index );

            floatVector getFollowingConfigurationJacobian( const unsigned int &index );

            floatVector getPreviousSubConfigurationJacobian( const unsigned int &lowerIndex, const unsigned int &upperIndex );

            floatVector getPreviousPrecedingConfigurationJacobian( const unsigned int &index );

            floatVector getPreviousFollowingConfigurationJacobian( const unsigned int &index );

            const floatType* getLSResidualNorm( );

            virtual void setResidualClasses( );

            void setResidualClasses( std::vector< ResidualBase<hydraBase>* > &residualClasses );

            std::vector< ResidualBase<hydraBase>* >* getResidualClasses( );

            const floatVector* getResidual( );

            const floatVector* getFlatJacobian( );

            virtual const floatVector* getFlatNonlinearLHS( );

            virtual const floatVector* getNonlinearRHS( );

            const floatVector* getFlatPreconditioner( );

            floatMatrix getJacobian( );

            const floatVector* getFlatdRdF( );

            floatMatrix getdRdF( );

            const floatVector* getdRdT( );

            const floatVector* getFlatdRdAdditionalDOF( );

            floatMatrix getdRdAdditionalDOF( );

            const floatVector* getFlatAdditionalDerivatives( );

            floatMatrix getAdditionalDerivatives( );

            const floatVector* getUnknownVector( );

            const floatVector* getTolerance( );

            virtual bool checkConvergence( );

            virtual bool checkLSConvergence( );

            virtual bool checkGradientConvergence( const floatVector &X0 );

            virtual bool checkDescentDirection( const floatVector &dx );

            const floatVector* getStress( );

            const floatVector* getPreviousStress( );

            const floatVector* getPreviouslyConvergedStress( );

            void setViscoplasticDamping( const floatType &factor );

            void clearViscoplasticDamping( );

            const floatType getViscoplasticDamping( ){ /*! Get the value of the viscoplastic damping */ return _viscoplastic_damping_factor; }

            const bool getViscoplasticDampingSet( );

            virtual void evaluate( const bool &use_subcycler = false );

            virtual void computeTangents( );

            virtual void computedXdAdditionalDOF( );

            const floatVector *getFlatdXdF( );

            const floatVector *getFlatdXdT( );

            const floatVector *getFlatdXdAdditionalDOF( );

            //! Return the flag which indicates whether hydra should initialize the unknown vector
            const bool getInitializeUnknownVector( ){ return _initializeUnknownVector; }

            //! Add data to the vector of values which will be cleared after each iteration
            void addIterationData( dataBase *data ){ _iterationData.push_back( data ); }

            //! Add data to the vector of values which will be cleared after each non-linear step
            void addNLStepData( dataBase *data ){ _nlStepData.push_back( data ); }

            unsigned int getNumNewton( ){ /*! Get the number of newton steps performed */  return _NUM_NEWTON; }

            unsigned int getNumLS( ){ /*! Get the number of line search steps performed */ return _NUM_LS; }

            unsigned int getNumGrad( ){ /*! Get the number of gradient descent steps performed */  return _NUM_GRAD; }

            const bool getUseSQPSolver( ){ /*! Return a flag for whether to use the SQP solver */ return _useSQPSolver; }

            const void setMaxRelaxedIterations( const unsigned int &value ){
                /*!
                 * Set the maximum allowable number of relaxed iterations
                 * 
                 * \param &value: The number of relaxed iterations
                 */

                _maxRelaxedIterations = value;

            }

            const unsigned int getMaxRelaxedIterations( ){
                /*!
                 * Get the maximum number of relaxed iterations
                 */

                return _maxRelaxedIterations;

            }

            void setFailureVerbosityLevel( const unsigned int &value ){
                /*!
                 * Set the verbosity level for failures
                 * 
                 * \param &value: The verbosity level of the failure (defaults to zero)
                 */

                _failure_verbosity_level = value;

            }

            //! Set the failure output string to use scientific notation
            void setFailureOutputScientific( ){ _failure_output << std::scientific; }

            //! Get the verbosity level for failure outputs
            const unsigned int getFailureVerbosityLevel( ){ return _failure_verbosity_level; }

            //! Add a string to the failure output string
            void addToFailureOutput( const std::string &additional ){ _failure_output << additional; }

            //! Add a floating point vector to the output string
            void addToFailureOutput( const floatVector &value, bool add_endline = true ){ for ( auto v = value.begin( ); v != value.end( ); v++ ){ _failure_output << *v << ", "; } if ( add_endline ){_failure_output << "\n"; } }

            //! Add a boolean vector to the output string
            void addToFailureOutput( const std::vector<bool> &value, bool add_endline = true ){ for ( auto v = value.begin( ); v != value.end( ); v++ ){ _failure_output << *v << ", "; } if ( add_endline ){_failure_output << "\n"; } }

            //! Add a floating point value to the output string
            void addToFailureOutput( const floatType &value, bool add_endline = true ){ _failure_output << value; if ( add_endline ){_failure_output << "\n"; } }

            //! Add a general iterable object to the output string
            template< class v_iterator >
            void addToFailureOutput( const v_iterator &v_begin, const v_iterator &v_end, bool add_endline = true ){ for ( auto v = v_begin; v != v_end; ++v ){ _failure_output << *v << ", "; } if ( add_endline ){_failure_output << "\n"; } }

            //! Get the failure output string
            const std::string getFailureOutput( ){ return _failure_output.str( ); }

            //! Get a scale factor for the deformation
            const floatType getScaleFactor( ){ return _scale_factor; }

            //! Set the value of the scale factor. Will automatically re-calculate the deformation and trial stresses
            virtual void setScaleFactor( const floatType &value );

            const floatType   getScaledTime( ){ /*! Get the value of the scaled current time */ return _scaled_time; }

            const floatType   getScaledDeltaTime( ){ /*! Get the value of the scaled changed in time */ return _scaled_deltaTime; }

            const floatType   getScaledTemperature( ){ /*! Get the value of the scaled current temperature */ return _scaled_temperature; }

            const floatVector *getScaledDeformationGradient( ){ /*! Get the value of the scaled current deformation gradient */ return &_scaled_deformationGradient; }

            const floatVector *getScaledAdditionalDOF( ){ /*! Get the value of the scaled current additional DOF */ return &_scaled_additionalDOF; }

            const floatType getCutbackFactor( ){ /*! Get the value of the cutback factor */ return _cutback_factor; }

            const unsigned int getNumGoodControl( ){ /*! Get the number of good iterations we need to have before increasing the timestep */ return _num_good_control; }

            const floatType getGrowthFactor( ){ /*! Get the growth factor for the timestep increase */ return _growth_factor; }

            const floatType getMinDS( ){ /*! Get the minimum allowable ratio of the total timestep to the cutback timestep */ return _minDS; }

            void setCutbackFactor( const floatType &value ){ /*! Get the current value of the cutback factor. \param &value: The value of the cutback */  _cutback_factor = value; }

            void setNumGoodControl( const unsigned int &value ){ /*! Set the number of good iterations that need to happen before the timestep increases. \param &value: The value of the number of good iterations prior to increasing the relative timestep */  _num_good_control = value; }

            void setGrowthFactor( const floatType &value ){ /*! Set the relative growth factor for the local timestep increase \param &value: The new value */ _growth_factor = value; }

            void setMinDS( const floatType &value ){ /*! Set the minimum value of the relative cutback timestep \param &value: The new value */  _minDS = value; }

            const bool allowStepGrowth( const unsigned int &num_good );

            floatVector *getMutableResidual( ){
                /*! Get a reference to the full residual that is mutable. Returns NULL if it's not allowed.
                 *
                 * This should only be called in residual classes that need to modify the full residual in their modifyGlobalResidual methods.
                 *
                 * Be careful!
                 */
                
                if ( _allow_modify_global_residual ){
                    return &_residual.second;
                }

                return NULL;

            }

            floatVector *getMutableJacobian( ){
                /*! Get a reference to the full jacobian that is mutable. Returns NULL if it's not allowed.
                 *
                 * This should only be called in residual classes that need to modify the full residual in their modifyGlobalJacobian methods.
                 *
                 * Be careful!
                 */
                
                if ( _allow_modify_global_jacobian ){
                    return &_jacobian.second;
                }

                return NULL;

            }

            floatVector *getMutabledRdT( ){
                /*! Get a reference to the full dRdT that is mutable. Returns NULL if it's not allowed.
                 *
                 * This should only be called in residual classes that need to modify the full residual in their modifyGlobaldRdT methods.
                 *
                 * Be careful!
                 */
                
                if ( _allow_modify_global_dRdT ){
                    return &_dRdT.second;
                }

                return NULL;

            }

            floatVector *getMutabledRdF( ){
                /*! Get a reference to the full dRdF that is mutable. Returns NULL if it's not allowed.
                 *
                 * This should only be called in residual classes that need to modify the full residual in their modifyGlobaldRdF methods.
                 *
                 * Be careful!
                 */
                
                if ( _allow_modify_global_dRdF ){
                    return &_dRdF.second;
                }

                return NULL;

            }

            floatVector *getMutabledRdAdditionalDOF( ){
                /*! Get a reference to the full dRdAdditionalDOF that is mutable. Returns NULL if it's not allowed.
                 *
                 * This should only be called in residual classes that need to modify the full residual in their modifyGlobaldRdAdditionalDOF methods.
                 *
                 * Be careful!
                 */
                
                if ( _allow_modify_global_dRdAdditionalDOF ){
                    return &_dRdAdditionalDOF.second;
                }

                return NULL;

            }

            const bool currentResidualIndexMeaningful( ){
                /*!
                 * Return if the current residual index is meaningful or not
                 */
                return _current_residual_index_set;
            }

            std::string getResidualParameterizationInfo( );

            void updateAdditionalStateVariables( ){
                /*!
                 * Update the additional state variable vector
                 */

                for ( auto v = std::begin( *getResidualClasses( ) ); v != std::end( *getResidualClasses( ) ); ++v ){

                    ( *v )->updateAdditionalStateVariables( _additionalStateVariables.second );

                }

            }

            void setToleranceScaleFactor( floatType factor ){
                /*!
                 * Loosen the convergence tolerance for the next iteration
                 * Useful if a Residual's form is changing in a non-smooth way
                 *
                 * \param factor: The scale factor to be applied to the current tolerance
                 */

                if ( factor > _residual_scale_factor ){

                    _residual_scale_factor = factor;

                }

            }

            const floatType getToleranceScaleFactor( ){
                /*!
                 * Get the scale factor for the tolerance
                 */

                return _residual_scale_factor;

            }

            const unsigned int getRelaxedIteration( ){
                /*!
                 * Get the current relaxed iteration
                 */

                return _relaxedIteration;
            }

            // TEMP REMOVE THESE

            //!< Get the current value of mu_k
            const floatType getMuk( ){ return _solver._step._mu_k; }

            //!< Set the Levenberg-Marquardt mu_k
            void setMuk( const floatType &value ){
               /*!
                * Set the value of the mu_k parameter for Levenberg-Marquardt steps
                *
                * \param &value: The value of the parameter
                */
 
                _solver._step.setMuk( value );

            }

            //!< Get the Levenberg-Marquardt mu parameter
            const floatType getLMMu( ){ return _solver._step._lm_mu; }

            //!< Set the Levenberg-Marquardt mu parameter
            void setLMMu( const floatType &value ){
               /*!
                * Set the value of the mu parameter for Levenberg-Marquardt steps
                *
                * \param &value: The value of the parameter
                */
 
                _solver._step.setLMMu( value );

            }

            virtual void setBaseQuantities( ){ _solver._step.setBaseQuantities( ); }

            // END TEMP

        protected:

            SolverBase _solver; //!< Temporary
            SolverBase *solver = &_solver; //!< The class which performs the material point solve TODO: Make this an incoming pointer

            // Setters that the user may need to access but not override

            const void setCurrentResidualIndexMeaningful( const bool &value ){
                /*!
                 * Set if the current residual index is meaningful
                 * 
                 * \param &value: Set if the current residual index is meaningful or not
                 */

                _current_residual_index_set = value;

            }

            const void setCurrentResidualIndex( const unsigned int value ){
                /*!
                 * Set if the current residual index is meaningful
                 * 
                 * \param value: Set the value of the current residual index
                 */

                _current_residual_index = value;

            }

            const void resetToleranceScaleFactor( ){
                /*!
                 * Reset the tolerance scale factor to 1
                 */

                _residual_scale_factor = 1.0;

            }

            void setStress( const floatVector &stress );

            // Utility functions
            virtual void computeConfigurations( const floatVector *data_vector, const unsigned int start_index,
                                                const floatVector &total_transformation,
                                                floatVector &configurations, floatVector &inverseConfigurations,
                                                const bool add_eye=false );

            virtual void extractStress( );

            virtual void updateConfigurationsFromUnknownVector( );

            virtual void decomposeUnknownVector( );

            virtual void decomposeStateVariableVector( );

            virtual void formNonLinearResidual( );

            virtual void formNonLinearDerivatives( );

            virtual void formPreconditioner( );

            virtual void solveNonLinearProblem( );

            virtual void formMaxRowPreconditioner( );

            virtual void initializeUnknownVector( );

            virtual void setTolerance( );

            virtual void initializePreconditioner( );

            virtual void evaluateInternal( );

            //! Update the line-search lambda parameter
            virtual void updateLambda( ){ _lambda *= 0.5; }

            virtual void updateUnknownVector( const floatVector &newUnknownVector );

            virtual void calculateFirstConfigurationJacobians( const floatVector &configurations, fourthOrderTensor &dC1dC, floatVector &dC1dCn );

            virtual void performArmijoTypeLineSearch( const floatVector &X0, const floatVector &deltaX );

            virtual void performGradientStep( const floatVector &X0 );

            //! Update the scaled quantities
            virtual void setScaledQuantities( );

            const floatType *get_baseResidualNorm( );

            const floatVector *get_basedResidualNormdX( );

            DataStorage< floatType > _baseResidualNorm; //!< The base value of the norm of the residual

            DataStorage< floatVector > _basedResidualNormdX; //!< The base value of the derivative of the norm of the residual w.r.t. the unknown vector

            void set_baseResidualNorm( const floatType &value ){ /*! Set the base value of the residual norm \param &value: The new value */  solver->step->set_baseResidualNorm( value ); }

            void set_basedResidualNormdX( const floatVector &value ){ /*! Set the base derivative of the residual norm w.r.t. the unknown vector \param &value: The new value */  solver->step->set_basedResidualNormdX( value ); }

            void resetNumNewton( ){ /*! Reset the number of newton steps */ _NUM_NEWTON = 0; }

            void resetNumLS( ){ /*! Reset the number of line search steps */ _NUM_LS = 0; }

            void resetNumGrad( ){ /*! Reset the number of gradient descent steps */ _NUM_GRAD = 0; }

            void incrementNumNewton( ){ /*! Reset the number of newton steps */ _NUM_NEWTON++; }

            void incrementNumLS( ){ /*! Reset the number of line search steps */ _NUM_LS++; }

            void incrementNumGrad( ){ /*! Reset the number of gradient descent steps */ _NUM_GRAD++; }

            void resetIterations( ){ /*! Reset the number of iterations */ _iteration = 0; }

            virtual void callResidualSuccessfulNLStep( );

            virtual void callResidualPreNLSolve( );

            virtual void callResidualPostNLSolve( );

            virtual void callResidualPreSubcycler( );

            virtual void callResidualPostSubcyclerSuccess( );

            virtual void callResidualPostSubcyclerFailure( );

            virtual bool callResidualRelaxedStepFailure( );

            template<class T>
            void setIterationData( const T &data, DataStorage<T> &storage ){
                /*!
                 * Template function for adding iteration data. These values are cleared
                 * every time the unknown vector is updated.
                 *
                 * \param &data: The data to be added
                 * \param &storage: The storage to add the data to
                 */

                storage.second = data;

                storage.first = true;

                addIterationData( &storage );

            }

            template<class T>
            void setNLStepData( const T &data, DataStorage<T> &storage ){
                /*!
                 * Template function for adding nonlinear step data. These values are cleared
                 * every time the nonlinear step is advanced.
                 *
                 * \param &data: The data to be added
                 * \param &storage: The storage to add the data to
                 */

                storage.second = data;

                storage.first = true;

                addNLStepData( &storage );

            }

            template<class T>
            void setPreviousData( const T &data, DataStorage<T> &storage ){
                /*!
                 * Template function for adding previous data
                 * 
                 * \param &data: The data to be added
                 * \param &storage: The storage to add the data to
                 */

                storage.second = data;

                storage.first = true;

            }

            template<class T>
            void setConstantData( const T &data, DataStorage<T> &storage ){
                /*!
                 * Template function for adding constant data
                 * 
                 * \param &data: The data to be added
                 * \param &storage: The storage to add the data to
                 */

                storage.second = data;

                storage.first = true;

            }

            void unexpectedError( ){
                /*!
                 * Function to throw for an unexpected error. A user should never get here!
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "You shouldn't have gotten here. If you aren't developing the code then contact a developer with the stack trace." ) )

            }

            //!A pass through function that does nothing
            void passThrough( ){ }

            void setX( const floatVector &X ){
                /*!
                 * Set the value of the unknown vector
                 * 
                 * \param &X: The unknown vector
                 */

                _X.second = X;

                _X.first = true;

            }

            virtual void setResidualNorm( );

            virtual void setdResidualNormdX( );

            std::string build_upper_index_out_of_range_error_string( const unsigned int upperIndex, const unsigned int num_configurations );
            std::string build_lower_index_out_of_range_error_string( const unsigned int lowerIndex, const unsigned int upperIndex );

            floatVector _prerelaxed_initialX; //!< The initial value of the unknown vector prior to a relaxed step

            floatVector _initialX; //!< The initial value of the unknown vector

            //! A data storage class that resets at every iteration
            template< typename T >
            class SetDataStorageIteration : public SetDataStorageIterationBase< hydraBase, T > {

                using tardigradeHydra::SetDataStorageIterationBase<hydraBase,T>::SetDataStorageIterationBase;

            };

            /*!
             * Class which defines setting values defined at the previous timestep
             */
            template< typename T >
            class SetDataStoragePrevious : public SetDataStorageBase< T > {

                public:

                    //! Create a data storage object that will be reset whenever the previous value gets reset
                    SetDataStoragePrevious( DataStorage< T > *ds ) : SetDataStorageBase< T >( ds ){
                        /*!
                         * Constructor for data storage objects for temporally previous objects
                         * 
                         * \param *ds: The data storage object to modify
                         */
                    }

            };

            /*!
             * Class which defines setting constant values regardless of the timestep
             */
            template< typename T >
            class SetDataStorageConstant : public SetDataStorageBase< T > {

                public:

                    SetDataStorageConstant( DataStorage< T > *ds ) : SetDataStorageBase< T >( ds ){
                        /*!
                         * Constructor for constant data storage objects
                         * 
                         * \param *ds: The data storage object
                         */
                    }

            };

            virtual tardigradeHydra::hydraBase::SetDataStorageConstant<floatVector> get_SetDataStorage_tolerance( );

            virtual tardigradeHydra::hydraBase::SetDataStorageIteration<secondOrderTensor> get_SetDataStorage_stress( );

            virtual void assembleKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints );

            virtual void updateKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints );

            virtual void assembleKKTRHSVector( const floatVector &dx, floatVector &KKTRHSVector, const std::vector< bool > &active_constraints );

            virtual void solveConstrainedQP( floatVector &dx, const unsigned int kmax=100 );

            virtual void setConstraints( );

            virtual void setConstraintJacobians( );

            virtual void initializeActiveConstraints( std::vector< bool > &active_constraints );

            void setUseSQPSolver( const unsigned int &value ){ /*! Set whether to use the SQP solver \param &value: The updated value */ _useSQPSolver = value; }

            void setInitializeUnknownVector( const bool &value ){
                /*!
                 * Set the initialize unknown vector flag
                 * 
                 * \param &value: The value of the flag
                 */

                _initializeUnknownVector = value;

            }

            virtual void resetProblem( );

            virtual void performRelaxedSolve( );

            void setAllowModifyGlobalResidual( const bool value){
                /*!
                 * Set a flag for if the global residual can be modified
                 *
                 * \param value: The updated value
                 */
                _allow_modify_global_residual = value;
            }

            void setAllowModifyGlobalJacobian( const bool value){
                /*!
                 * Set a flag for if the global jacobian can be modified
                 *
                 * \param value: The updated value
                 */
                _allow_modify_global_jacobian = value;
            }

            void setAllowModifyGlobaldRdT( const bool value){
                /*!
                 * Set a flag for if the global dRdT can be modified
                 *
                 * \param value: The updated value
                 */
                _allow_modify_global_dRdT = value;
            }

            void setAllowModifyGlobaldRdF( const bool value){
                /*!
                 * Set a flag for if the global dRdF can be modified
                 *
                 * \param value: The updated value
                 */
                _allow_modify_global_dRdF = value;

            }

            void setAllowModifyGlobaldRdAdditionalDOF( const bool value){
                /*!
                 * Set a flag for if the global dRdAdditionalDOF can be modified
                 *
                 * \param value: The updated value
                 */
                _allow_modify_global_dRdAdditionalDOF = value;

            }

            void setPreviouslyConvergedStress( const floatVector &value ){
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
            friend class tardigradeHydra::SolverStepBase; //!< The base class for the solver steps
            friend class tardigradeHydra::SolverBase; //!< The base class for the solver
            friend class unit_test::hydraBaseTester; //!< Friend class which allows modification of private variables. ONLY TO BE USED FOR TESTING!

            unsigned int _dimension; //!< The spatial dimension of the problem

            unsigned int _configuration_unknown_count; //!< The number of unknowns required for a configuration. Used to ensure that the unknown and state variable vectors are the right size. Must be set by all inheriting classes. For 3D classical continuum this will be 9, for higher order theories this will change.

            unsigned int _stress_size; //!< The number of terms in the stress measures. For 3D classical continuum this will be 9, for higher order theories this will change.

            floatType _time; //!< The current time

            floatType _deltaTime; //!< The change in time

            floatType _temperature; //!< The current temperature

            floatType _previousTemperature; //!< The previous temperature

            secondOrderTensor _deformationGradient; //!< The current deformation gradient

            secondOrderTensor _previousDeformationGradient; //!< The previous deformation gradient

            floatVector _additionalDOF; //!< The current additional degrees of freedom

            floatVector _previousAdditionalDOF; //!< The previous additional degrees of freedom

            floatVector _previousStateVariables; //!< The previous state variables

            floatVector _parameters; //!< The model parameters

            floatType _scaled_time; //!< The current time scaled by the scaling factor

            floatType _scaled_deltaTime; //!< The change in time scaled by the scaling factor

            floatType _scaled_temperature; //!< The current temperature scaled by the scaling factor

            secondOrderTensor _scaled_deformationGradient; //!< The current deformation gradient scaled by the scaling factor

            floatVector _scaled_additionalDOF; //!< The current additional degrees of freedom scaled by the scaling factor

            unsigned int _numConfigurations; //!< The number of configurations

            unsigned int _numNonLinearSolveStateVariables; //!< The number of state variables which will be solved in the Newton-Raphson loop

            floatType _residual_scale_factor = 1; //!< A scale factor for the residual which can be used to loosen the tolerance by a Residual

            floatType _tolr; //!< The relative tolerance

            floatType _tola; //!< The absolute tolerance

            unsigned int _maxIterations; //!< The maximum number of allowable iterations

            unsigned int _maxLSIterations; //!< The maximum number of line-search iterations

            unsigned int _maxGradientIterations = 10; //!< The maximum number of gradient iterations

            unsigned int _NUM_NEWTON = 0; //!< The number of Newton steps performed

            unsigned int _NUM_LS = 0; //!< The number of line search steps performed

            unsigned int _NUM_GRAD = 0; //!< The number of gradient descent steps performed

            bool _use_LM_step = false; //!< Flag for whether to attempt a Levenberg-Marquardt step

            bool _use_relaxed_solve = true; //!< Flag for whether to attempt a relaxed solve in case of failure

            bool _use_gradient_descent = false; //!< Flag for whether to attempt a gradient descent step

            bool _rank_deficient_error = true; //!< Flag for whether a rank-deficient LHS will throw a convergence error

            bool _initializeUnknownVector = true; //!< Flag for whether to initialize the unknown vector in the non-linear solve

            unsigned int _maxRelaxedIterations = 5; //!< The number of allowed relaxed iterations

            floatType _lsAlpha; //!< The line-search alpha value i.e., the term by which it is judged that the line-search is converging

            std::vector< dataBase* > _iterationData; //!< A vector of pointers to data which should be cleared at each iteration

            std::vector< dataBase* > _nlStepData; //!< A vector of pointers to data which should be cleared after each nonlinear step

            DataStorage< std::vector< ResidualBase<hydraBase>* > > _residualClasses; //!< A vector of classes which compute the terms in the residual equation

            DataStorage< floatVector > _residual; //!< The residual vector for the global solve

            DataStorage< floatVector > _jacobian; //!< The jacobian matrix in row-major form for the global solve

            DataStorage< floatVector > _nonlinearRHS; //!< The right hand side vector for the Newton solve

            DataStorage< floatVector > _flatNonlinearLHS; //!< The left hand side vector for the Newton solve

            DataStorage< floatVector > _preconditioner; //!< The pre-conditioner matrix in row-major form for the global solve

            DataStorage< floatVector > _previouslyConvergedStress; //!< The previously converged stress

            bool _use_preconditioner; //!< Flag for whether to pre-condition the Jacobian or not

            unsigned int _preconditioner_type; //<! The type of preconditioner to use

            bool _preconditioner_is_diagonal; //!< Flag for if the pre-conditioner only stores the diagonal elements

            floatType _gradientSigma = 1e-4; //!< The sigma parameter for the gradient descent step

            floatType _gradientBeta  = 0.9; //!< The beta parameter for the gradient descent step

            floatType _gradientRho   = 1e-8; //!< The rho parameter for the gradient descent step

            floatType _gradientP     = 2.1; //!< The p parameter for the gradient descent step

            DataStorage< floatVector > _dRdF; //!< The gradient of the residual w.r.t. the deformation gradient in row-major form for the global solve

            DataStorage< floatVector > _dRdT; //!< The gradient of the residual w.r.t. the temperature for the global solve

            DataStorage< floatVector > _dRdAdditionalDOF; //!< The derivatives of the residual w.r.t. the additional degrees of freedom

            DataStorage< floatVector > _additionalDerivatives; //!< Additional derivatives of the residual

            DataStorage< floatVector > _X; //!< The unknown vector { stress, F1, ..., Fn, xi1, ..., xim }

            DataStorage< floatVector > _tolerance; //!< The tolerance vector for the non-linear solve

            DataStorage< floatType > _lsResidualNorm; //!< The reference residual norm for the line-search convergence criteria

            DataStorage< floatVector > _stress; //!< The stress in the current configuration as determined from the current state

            DataStorage< floatVector > _previousStress; //!< The previous value of the stress in the current configuration as determined from the previous state

            DataStorage< floatVector > _flatdXdF; //!< The total derivative of the unknown vector w.r.t. the deformation in row-major form

            DataStorage< floatVector > _flatdXdT; //!< The total derivative of the unknown vector w.r.t. the temperature

            DataStorage< floatVector > _flatdXdAdditionalDOF; //!< The total derivative of the unknown vector w.r.t. the additional DOF

            unsigned int _iteration = 0; //!< The current iteration of the non-linear problem

            unsigned int _LSIteration = 0; //!< The current line search iteration of the non-linear problem

            unsigned int _gradientIteration = 0; //!< The current gradient iteration of the non-linear problem

            unsigned int _relaxedIteration = 0; //!< The current relaxed iteration of the non-linear problem

            floatType _lambda = 1;

            bool _useSQPSolver = false;

            void setFirstConfigurationJacobians( );

            void setPreviousFirstConfigurationJacobians( );

            void performPreconditionedSolve( floatVector &deltaX_tr );

            void solveNewtonUpdate( floatVector &deltaX_tr );

            void setTolerance( const floatVector &tolerance );

            void incrementIteration( ){ _iteration++; resetLSIteration( ); }

            void incrementLSIteration( ){ _LSIteration++; }

            void incrementGradientIteration( ){ _gradientIteration++; }

            void resetLSIteration( ){ _LSIteration = 0; _lambda = 1.0;
                                      _lsResidualNorm.second = tardigradeVectorTools::l2norm( *getResidual( ) );
                                      _lsResidualNorm.first = true; }

            void resetGradientIteration( ){ _gradientIteration = 0; }

            const floatType getLambda( ){ return _lambda; }

            bool checkIteration( ){ return _iteration < _maxIterations; }

            bool checkLSIteration( ){ return _LSIteration < _maxLSIterations; }

            bool checkGradientIteration( ){ return _gradientIteration < _maxGradientIterations; }

            void resetIterationData( );

            void resetNLStepData( );

            void setRelaxedIteration( const unsigned int &value ){
                /*! Set the relaxed iteration number
                 * \param &value: The incoming value
                 */

                _relaxedIteration = value;

            }

            unsigned int _failure_verbosity_level = 0; //!< The verbosity level for failure.

            std::stringstream _failure_output; //!< Additional failure output information

            floatType _scale_factor = 1.0; //!< A scale factor applied to the incoming loading (deformation, temperature, etc.)

            floatType _cutback_factor = 0.5; //!< The factor by which the pseudo-time will be scaled if a solve fails

            floatType _growth_factor  = 1.2; //!< The factor by which the pseudo-time will be scaled if we can grow the pseudo-timestep

            unsigned int _num_good_control = 2; //!< The number of good iterations we need to have before we try and increase the timestep

            floatType _minDS = 1e-2; //!< The minimum allowable pseudo-timestep

            bool _allow_modify_global_residual = false; //!< Flag for if the global residual can be modified

            bool _allow_modify_global_jacobian = false; //!< Flag for if the global jacobian can be modified

            bool _allow_modify_global_dRdT = false; //!< Flag for if the global dRdT can be modified

            bool _allow_modify_global_dRdF = false; //!< Flag for if the global dRdF can be modified

            bool _allow_modify_global_dRdAdditionalDOF = false; //!< Flag for if the global dRdAdditionalDOF can be modified

            bool _current_residual_index_set = false; //!< Flag for whether the current residual index has been set

            int _current_residual_index = 0; //!< The current residual index

            floatType _viscoplastic_damping_factor = 0; //!< The fraction of the difference between the trial stress and the stress that will be suppressed to assist in convergence

            bool _viscoplastic_damping_set = false; //!< Flag for whether the viscoplastic damping factor has been set

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, configurations,                       floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousConfigurations,               floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, inverseConfigurations,                floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousInverseConfigurations,        floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, nonLinearSolveStateVariables,         floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousNonLinearSolveStateVariables, floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE( private, additionalStateVariables,             floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_PREVIOUS_STORAGE(  private, previousAdditionalStateVariables,     floatVector, passThrough )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, residualNorm,       floatType,          setResidualNorm )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dResidualNormdX,    floatVector,        setdResidualNormdX )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setConstraints,         getConstraints,         constraints,         floatVector,       setConstraints )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setConstraintJacobians, getConstraintJacobians, constraintJacobians, floatVector,       setConstraintJacobians )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, set_dF1dF,              get_dF1dF,              dF1dF,               secondOrderTensor, setFirstConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, set_dF1dFn,             get_dF1dFn,             dF1dFn,              floatVector,       setFirstConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_NAMED_PREVIOUS_STORAGE(  private, set_previousdF1dF,      get_previousdF1dF,      previousdF1dF,       secondOrderTensor, setPreviousFirstConfigurationJacobians )

            TARDIGRADE_HYDRA_DECLARE_NAMED_PREVIOUS_STORAGE(  private, set_previousdF1dFn,     get_previousdF1dFn,     previousdF1dFn,      floatVector,       setPreviousFirstConfigurationJacobians )

    };

    //! A data storage class that updates at every iteration
    template< typename T >
    class SetDataStorageIteration : public SetDataStorageIterationBase< ResidualBase<hydraBase>, T > {

        using tardigradeHydra::SetDataStorageIterationBase<ResidualBase<hydraBase>,T>::SetDataStorageIterationBase;

    };

    //! A data storage class that updates whenever the previous values change
    template< typename T >
    class SetDataStoragePrevious : public SetDataStorageBase< T > {

        public:

            SetDataStoragePrevious( DataStorage< T > *ds ) : SetDataStorageBase< T >( ds ){
                /*!
                 * Constructor for data storage objects for temporally previous objects
                 * 
                 * \param *ds: The data storage object to modify
                 */
            }

    };

    //! A data storage class that is constant
    template< typename T >
    class SetDataStorageConstant : public SetDataStorageBase< T > {

        public:

            SetDataStorageConstant( DataStorage< T > *ds ) : SetDataStorageBase< T >( ds ){
                /*!
                 * Constructor for constant data storage objects
                 * 
                 * \param *ds: The data storage object
                 */
            }

    };

}

#endif
