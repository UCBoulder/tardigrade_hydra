/**
 ******************************************************************************
 * \file tardigrade_ResidualBase.h
 ******************************************************************************
 * The base class for residuals
 ******************************************************************************
 */

#ifndef TARDIGRADE_RESIDUALBASE
#define TARDIGRADE_RESIDUALBASE

#include"tardigrade_CoreDefinitions.h"
#include"tardigrade_SetDataStorage.h"

namespace tardigradeHydra{

    /*!
     * A class to contain the residual computations associated with some part of a non-linear solve
     */
    class ResidualBase{

        public:

            /*!
             * Default residual
             */
            ResidualBase( ) : hydra( NULL ), _numEquations( 0 ){ };

            /*!
             * Main utilization constructor
             * 
             * \param *_hydra: A pointer to a hydraBase object
             * \param &_numEquations: The number of equations defined by the residual
             */
            ResidualBase( hydraBase *_hydra, const unsigned int &_numEquations ) : hydra( _hydra ), _numEquations( _numEquations ){ }

            /*!
             * Copy constructor
             * 
             * \param &r: The residual to be copied
             */
            ResidualBase( ResidualBase &r ) : hydra( r.hydra ), _numEquations( r.getNumEquations( ) ), _numConstraints( r.getNumConstraints( ) ){ }

            hydraBase* hydra; //!< The hydra class which owns the ResidualBase object

            // User defined setter functions

            virtual void setResidual( ){
                /*!
                 * The user-defined residual equation. Must have a size of numEquations
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The residual is not implemented" ) );

            }

            virtual void setJacobian( ){
                /*!
                 * The user-defined jacobian equation. Must have a size of numEquations x numUnknowns
                 * 
                 * The order of the unknowns are the cauchy stress, the configurations in order (minus the first one),
                 * and the state variables solved for in the non-linear solve.
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The jacobian is not implemented" ) );

            }

            virtual void setdRdF( ){
                /*!
                 * The user-defined derivative of the residual w.r.t. the deformation gradient.
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The derivative of the residual w.r.t. the deformation gradient is not implemented" ) );
 
            }

            virtual void setdRdT( ){
                /*!
                 * The user-defined derivative of the residual w.r.t. the temperature
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The derivative of the residual w.r.t. the temperature is not implemented" ) );
 
            }

            virtual void setdRdAdditionalDOF( ){
                /*!
                 * The user-defined derivative of the residual w.r.t. the additional DOF
                 */

            }

            virtual void setAdditionalDerivatives( ){
                /*!
                 * The user-defined derivative of the residual w.r.t. additional values
                 */

            }

            virtual void setStress( ){
                /*!
                 * Compute the current stress
                 * 
                 * Only needs to be defined for the first residual
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The calculation of the stress is not implemented" ) );

            }

            virtual void setPreviousStress( ){
                /*!
                 * Compute the previous stress
                 * 
                 * Only needs to be defined for the first residual
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The calculation of the previous stress is not implemented" ) );

            }

            virtual void setConstraints( ){
                /*!
                 * Compute the contraints
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The calculation of the constraints is not implemented" ) );

            }

            virtual void setConstraintJacobians( ){
                /*!
                 * Compute the contraint Jacobians
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::logic_error( "The calculation of the constraint jacobians is not implemented" ) );

            }

            virtual void setCurrentAdditionalStateVariables( ){
                /*!
                 * Set the current additional state variables
                 * 
                 * Doesn't need to be defined
                 */

                setCurrentAdditionalStateVariables( floatVector( 0, 0 ) );

                return;

            }

            virtual void suggestInitialIterateValues( std::vector< unsigned int >   &indices,
                                                      std::vector< floatType > &values ){

                /*!
                 * Function which is called which allows the residual to suggest initial values for given
                 * configurations. This is called when the unknown vector is being initialized. If more than
                 * one residual attempts to set the initial vector the last residual will override all of the others.
                 *
                 * After the initial iterate has been suggested, the iteration data is cleared so that the residual
                 * starts the iteration in a clean state.
                 * 
                 * \param &indices: The indices of the unknown vector to set
                 * \param &values:  The values to be set in the unknown vector
                 */

                indices.clear( );
                values.clear( );

            }

            virtual void projectSuggestedX( std::vector< floatType > &trialX,
                                            const std::vector< floatType > &Xp ){
                /*!
                 * Project the suggested unknown vector to the allowable space
                 * 
                 * Called whenever hydra calls updateUnknownVector. It is assumed that the
                 * initial value as suggested by `residual::suggestInitialIterationValues` is
                 * in the allowable space.
                 * 
                 * \param &trialX: The trial value of X
                 * \param &Xp: The previously accepted value of X
                 */

            }

            void setUseProjection( const bool &value ){
                /*!
                 * Set whether to use the projection or not
                 * 
                 * \param &value: The value of the parameter
                 */

                _useProjection = value;

            }

            virtual bool checkRelaxedConvergence( ){
                /*!
                 * When performing a relaxed solve the residuals must return if they are converged or not.
                 * This function returns true or false (default true)
                 */

                return true;

            }

            virtual void modifyGlobalResidual( ){
                /*!
                 * Function that is called to modify the global residual.
                 *
                 * Called after all of the residuals are agglomerated into the whole.
                 *
                 * A mutable version of the global residual is accessable with hydra->getMutableResidual( )
                 */
            }

            virtual void modifyGlobalJacobian( ){
                /*!
                 * Function that is called to modify the global jacobian.
                 *
                 * Called after all of the residuals are agglomerated into the whole
                 *
                 * A mutable version of the global jacobian is accessable with hydra->getMutableJacobian( )
                 */
            }

            virtual void modifyGlobaldRdT( ){
                /*!
                 * Function that is called to modify the global derivative of the residual w.r.t. the temperature
                 *
                 * Called after all of the residuals are agglomerated into the whole
                 *
                 * A mutable version of the global dRdT is accessable with hydra->getMutabledRdT( )
                 */
            }

            virtual void modifyGlobaldRdF( ){
                /*!
                 * Function that is called to modify the global derivative of the residual w.r.t. the deformation
                 *
                 * Called after all of the residuals are agglomerated into the whole
                 *
                 * A mutable version of the global dRdF is accessable with hydra->getMutabledRdF( )
                 */
            }

            virtual void modifyGlobaldRdAdditionalDOF( ){
                /*!
                 * Function that is called to modify the global derivative of the residual w.r.t. the additional DOF
                 *
                 * Called after all of the residuals are agglomerated into the whole
                 *
                 * A mutable version of the global dRdAdditionalDOF is accessable with hydra->getMutabledRdAdditionalDOF( )
                 */
            }

            virtual void preNLSolve( ){
                /*!
                 * Function that is called prior to a nonlinear solve
                 */
            };

            virtual void postNLSolve( ){
                /*!
                 * Function that is called after a nonlinear solve
                 */
            };

            virtual void successfulNLStep( ){
                /*!
                 * Function that is called whenever a successful nonlinear step is taken
                 */
            };

            virtual void preSubcycler( ){
                /*!
                 * Function that is called prior to entering the subcycler
                 */
            }

            virtual void postSubcyclerSuccess( ){
                /*!
                 * Function that is called whenever the subcycler succeeds
                 */
            }

            virtual void postSubcyclerFailure( ){
                /*!
                 * Function that is called whenever the subcycler fails
                 */
            }

            virtual bool relaxedStepFailure( ){
                /*!
                 * Function that is called whenever the relaxed solve fails in a step
                 *
                 * Should return true if it is desirable to try another relaxed step
                 */

                return false;

            }

            virtual void setupRelaxedStep( const unsigned int &relaxedStep );

            //! Get the flag for whether to use the projection or not
            const bool getUseProjection( ){ return _useProjection; }

            // Getter functions

            //! Get the number of equations the residual defined
            const unsigned int getNumEquations( ){ return _numEquations; }

            //! Get the number of constraints the residual defined
            const unsigned int getNumConstraints( ){ return _numConstraints; }

            void addIterationData( dataBase *data );

            void addNLStepData( dataBase *data );

            template<class T>
            void setIterationData( const T &data, DataStorage<T> &storage ){
                /*!
                 * Template function for adding iteration data
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
                 * Template function for adding nonlinear step data
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

            virtual void addParameterizationInfo( std::string &parameterization_info ){
                /*!
                 * Add parameterization information to the provided docstring
                 * 
                 * Information on the current parameters and their values should be added to the
                 * incoming string.
                 * 
                 * \param &parameterization_info: The parameterization information string. Append
                 *     this class' information
                 */

                parameterization_info += "NO INFORMATION DEFINED FOR CLASS\n";

            }

            virtual void updateAdditionalStateVariables( floatVector &additionalStateVariables ){
                /*!
                 * Update the additional state variable vector
                 *
                 * \param &additionalStateVariables: The additional state variable vector
                 */

            }

        protected:

            void setNumConstraints( const unsigned int numConstraints ){
                /*!
                 * Set the number of constraints for the solve
                 * 
                 * \param numConstraints: The number of constraints
                 */

                _numConstraints = numConstraints;

            }

            void setPenaltyIndices( const std::vector< unsigned int > &indices ){
                /*!
                 * Set the indices where the penalties are defined
                 * 
                 * \param &indices: The indices where the penalty is defined
                 */
                _penalty_indices = indices;
            }

            void unexpectedError( ){
                /*!
                 * Function to throw for an unexpected error. A user should never get here!
                 */

                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "You shouldn't have gotten here. If you aren't developing the code then contact a developer with the stack trace." ) )

            }

            //! Class which defines data storage objects which are reset at each iteration
            template< typename T >
            class SetDataStorageIteration : public SetDataStorageIterationBase< ResidualBase, T > {

              public:

                  using tardigradeHydra::SetDataStorageIterationBase<ResidualBase,T>::SetDataStorageIterationBase;

            };

            //! Class which defines data storage objects which are reset at each nonlinear step
            template< typename T >
            class SetDataStorageNLStep : public SetDataStorageNLStepBase< ResidualBase, T > {

              public:

                  using tardigradeHydra::SetDataStorageNLStepBase<ResidualBase,T>::SetDataStorageNLStepBase;

            };

            //! Class which defines data storage objects for values defined at the previous timestep
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

            /*!
             * Class that is a constant data storage object
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

        private:

            unsigned int _numEquations; //!< The number of residual equations

            unsigned int _numConstraints = 0; //!< The number of constraint equations

            bool _useProjection = false; //!< Flag for whether to use the projection or not

            std::vector< unsigned int > _penalty_indices; //!< The indices of the variables which should be penalized for negative values

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setResidual,                         getResidual,                        residual,                        floatVector, setResidual )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setJacobian,                         getJacobian,                        jacobian,                        floatVector, setJacobian )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setdRdF,                             getdRdF,                            dRdF,                            floatVector, setdRdF )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setdRdT,                             getdRdT,                            dRdT,                            floatVector, setdRdT )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setdRdAdditionalDOF,                 getdRdAdditionalDOF,                dRdAdditionalDOF,                floatVector, setdRdAdditionalDOF )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setAdditionalDerivatives,            getAdditionalDerivatives,           additionalDerivatives,           floatVector, setAdditionalDerivatives )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setStress,                           getStress,                          stress,                          floatVector, setStress )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setPreviousStress,                   getPreviousStress,                  previousStress,                  floatVector, setPreviousStress )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setConstraints,                      getConstraints,                     constraints,                     floatVector, setConstraints )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setConstraintJacobians,              getConstraintJacobians,             constraintJacobians,             floatVector, setConstraintJacobians )

            TARDIGRADE_HYDRA_DECLARE_NAMED_ITERATION_STORAGE( private, setCurrentAdditionalStateVariables,  getCurrentAdditionalStateVariables, currentAdditionalStateVariables, floatVector, setCurrentAdditionalStateVariables )

    };

}

#endif
