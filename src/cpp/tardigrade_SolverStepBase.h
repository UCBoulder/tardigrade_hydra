/**
 ******************************************************************************
 * \file tardigrade_SolverStepBase.h
 ******************************************************************************
 * The base class for solver steps
 ******************************************************************************
 */

#ifndef TARDIGRADE_SOLVERSTEPBASE
#define TARDIGRADE_SOLVERSTEPBASE

#include"tardigrade_CoreDefinitions.h"
#include"tardigrade_SetDataStorage.h"

namespace tardigradeHydra{

    /*!
     * Base class for Solver Steps
     */
    class SolverStepBase{

        public:

            SolverStepBase( ) : solver(NULL){
                /*!
                 * Constructor for NonlinearStepBase
                 */
            }

            SolverStepBase( SolverBase *_solver ) : solver(_solver){
                /*!
                 * Constructor for NonlinearStepBase
                 *
                 * \param *_solver: The containing solver object
                 */
            }

            void incrementSolution( );

            floatVector X0; //!< The initial value of the unknown vector

            void setSolver( SolverBase *_solver ){
                /*! Set the containing solver object
                 * \param *_solver: The containing solver object
                 */
                solver = _solver;
            }

            // CACHED DATA STORAGE OPERATIONS
            void addIterationData( dataBase *data );

            void addNLStepData( dataBase *data );
            // END CACHED DATA STORAGE OPERATIONS

            // LEVENBERG-MARQUARDT FUNCTIONS (MOVE TO OWN CLASS)

            const floatType *get_baseResidualNorm( );

            const floatVector *get_basedResidualNormdX( );

            //! Get the current value of mu_k
            const floatType getMuk( ){ return _mu_k; }

            //! Get the Levenberg-Marquardt mu parameter
            const floatType getLMMu( ){ return _lm_mu; }

            // END LEVENBERG-MARQUARDT FUNCTIONS

            // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            const bool getUseSQPSolver( ){ /*! Return a flag for whether to use the SQP solver */ return _useSQPSolver; }

            // END SQP SOLVER FUNCTIONS
        protected:

            SolverBase *solver; //!< Pointer to the containing SolverBase object

            // CACHED DATA STORAGE OPERATIONS

            //! Class which defines data storage objects which are reset at each iteration
            template< typename T >
            class SetDataStorageIteration : public SetDataStorageIterationBase< SolverStepBase, T > {

              public:

                  using tardigradeHydra::SetDataStorageIterationBase<SolverStepBase,T>::SetDataStorageIterationBase;

            };

            //! Class which defines data storage objects which are reset at each nonlinear step
            template< typename T >
            class SetDataStorageNLStep : public SetDataStorageNLStepBase< SolverStepBase, T > {

              public:

                  using tardigradeHydra::SetDataStorageNLStepBase<SolverStepBase,T>::SetDataStorageNLStepBase;

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

            // END CACHED DATA STORAGE OPERATIONS

            // LEVENBERG-MARQUARDT FUNCTIONS (MOVE TO OWN CLASS)

            virtual void setResidualNorm( );

            virtual void setdResidualNormdX( );

            void set_baseResidualNorm( const floatType &value );

            void set_basedResidualNormdX( const floatVector &value );

            virtual void setBaseQuantities( );

            //!< Set the Levenberg-Marquardt mu_k
            void setMuk( const floatType &value ){
               /*!
                * Set the value of the mu_k parameter for Levenberg-Marquardt steps
                *
                * \param &value: The value of the parameter
                */
 
                _mu_k = value;

            }

            //!< Set the Levenberg-Marquardt mu parameter
            void setLMMu( const floatType &value ){
               /*!
                * Set the value of the mu parameter for Levenberg-Marquardt steps
                *
                * \param &value: The value of the parameter
                */
 
                _lm_mu = value;

            }

            floatType _mu_k = -1; //!< The Levenberg-Marquardt scaling parameter

            floatType _lm_mu = 1e-8; //!< The mu parameter for Levenberg-Marquardt iterations

            // END LEVENBERG-MARQUARDT FUNCTIONS

            // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            void setUseSQPSolver( const unsigned int &value ){ /*! Set whether to use the SQP solver \param &value: The updated value */ _useSQPSolver = value; }

            // END SQP SOLVER FUNCTIONS

        private:

            friend class tardigradeHydra::hydraBase; //!< TEMP REMOVE THIS
            friend class tardigradeHydra::unit_test::SolverStepBaseTester; //!< The unit tester for the class
            DataStorage< floatType > _baseResidualNorm; //!< The base value of the norm of the residual

            DataStorage< floatVector > _basedResidualNormdX; //!< The base value of the derivative of the norm of the residual w.r.t. the unknown vector

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, residualNorm,       floatType,          setResidualNorm )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dResidualNormdX,    floatVector,        setdResidualNormdX )

            // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            bool _useSQPSolver = false; //!< The flag for whether to use the SQP solver

            // END SQP SOLVER FUNCTIONS

    };

}

#endif
