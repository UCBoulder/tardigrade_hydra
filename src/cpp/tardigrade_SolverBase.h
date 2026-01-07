/**
  ******************************************************************************
  * \file tardigrade_SolverBase.h
  ******************************************************************************
  * A C++ library for the base classes for solvers
  ******************************************************************************
  */

#ifndef TARDIGRADE_SOLVERBASE_H
#define TARDIGRADE_SOLVERBASE_H

#include"tardigrade_CoreDefinitions.h"
#include"tardigrade_SetDataStorage.h"
#include"tardigrade_SolverStepBase.h"
#include"tardigrade_PreconditionerBase.h"

namespace tardigradeHydra{

    /*!
     * Base Solver class
     */
    class SolverBase : public CachingDataBase {

        public:

            SolverBase( ) : hydra(NULL){//, step(NULL){ ///!< TODO: Re-enable this

            }

            /*!
             * Constructor for NonlinearStepBase
             *
             * \param *_hydra: The containing hydraBase object
             */
            SolverBase( hydraBase * _hydra ) : hydra( _hydra ){ }

            hydraBase *hydra; //!< Pointer to the containing hydra object

            SolverStepBase _step; //!< Temporary object
            SolverStepBase *step = &_step; //!< The object that defines the step to be taken by the solver TODO: Make this an incoming pointer

            PreconditionerBase _preconditioner; //!< Temporary object
            PreconditionerBase *preconditioner = &_preconditioner; //!< The object that defines the preconditioner TODO: Make this an incoming pointer

            virtual void solve( );

            const bool getRankDeficientError( );

            void setRankDeficientError( const bool &value );

            // CACHED DATA STORAGE OPERATIONS
            virtual void addIterationData( dataBase *data ) override;

            virtual void addNLStepData( dataBase *data ) override;
            // END CACHED DATA STORAGE OPERATIONS

            // NONLINEAR FUNCTIONS (MOVE TO OWN CLASS)

            //! Return the flag which indicates whether hydra should initialize the unknown vector
            const bool getInitializeUnknownVector( ){ return _initializeUnknownVector; }

            void resetNLStepData( );

            //! Return the maximum number of allowable iterations
            const unsigned int getMaxIterations( ){ return _maxIterations; }

            void setMaxIterations( const unsigned int &value );

            //! Get the current nonlinear iteration number
            const unsigned int getIteration( ){ return _iteration; }

            //! Check if the number of nonlinear iterations has exceeded the allowable count
            bool checkIteration( ){ return getIteration( ) < getMaxIterations( ); }

            virtual bool checkConvergence( );

            const floatVector* getTolerance( );

            // END NONLINEAR FUNCTIONS

            // Pass-through functions
            const floatType getRelativeTolerance( ); //TODO: Want to allow this to be constexpr

            const floatType getAbsoluteTolerance( ); //TODO: Want to allow this to be constexpr

            const unsigned int getNumUnknowns( ); //TODO: Want to allow this to be constexpr

            const floatVector *getUnknownVector( ); //TODO: Want to generalize this

            void updateUnknownVector( const floatVector &value ); //TODO: Want to generalize this

            const floatVector *getResidual( ); //TODO: Want to generalize this

            const floatVector *getFlatJacobian( ); //TODO: Want to generalize this

            const unsigned int getNumConstraints( ); //TODO: Want to allow this to be constexpr

            const floatVector *getConstraints( ); //TODO: Want to generalize this

            const floatVector *getConstraintJacobians( ); //TODO: Want to generalize this

            const unsigned int getFailureVerbosityLevel( );

            void addToFailureOutput( const std::string &string );

            void addToFailureOutput( const floatVector &value, bool add_endline = true );

            void addToFailureOutput( const std::vector<bool> &value, bool add_endline = true );

            void addToFailureOutput( const floatType &value, bool add_endline = true );

            const floatType getToleranceScaleFactor( );

            void resetToleranceScaleFactor( );

            const floatVector *getFlatNonlinearLHS( ); //TODO: Want to generalize this

            /*!
             * Add a general iterable object to the output string
             * 
             * \param &v_begin: The starting iterator
             * \param &v_end: The stopping iterator
             * \param add_endline: Whether to add an endline to the string or not
             */
            template< class v_iterator >
            void addToFailureOutput( const v_iterator &v_begin, const v_iterator &v_end, bool add_endline = true ){

                std::stringstream failure_output;

                for ( auto v = v_begin; v != v_end; ++v ){ failure_output << *v << ", "; }

                if ( add_endline ){ failure_output << "\n"; }

                addToFailureOutput( failure_output.str( ) );

            }

        protected:

            // NONLINEAR FUNCTIONS (MOVE TO OWN CLASS)

            void setInitializeUnknownVector( const bool &value );

            virtual void callResidualSuccessfulNLStep( );

            virtual void callResidualPreNLSolve( );

            virtual void callResidualPostNLSolve( );

            //! Reset the number of iterations
            void resetIterations( ){ _iteration = 0; }

            void incrementIteration( );

            virtual void setTolerance( );

            void setTolerance( const floatVector &tolerance );

            virtual tardigradeHydra::SolverBase::SetDataStorageConstant<floatVector> get_SetDataStorage_tolerance( );

            // END NONLINEAR FUNCTIONS

        private:

            bool _rank_deficient_error = false; //!< Flag for whether a rank-deficient Jacobian should cause an error

            friend class tardigradeHydra::hydraBase; //!< TEMP REMOVE THIS
            friend class tardigradeHydra::unit_test::SolverBaseTester; //!< The unit tester for the class

            // NONLINEAR FUNCTIONS (MOVE TO OWN CLASS)

            bool _initializeUnknownVector = true; //!< Flag for whether to initialize the unknown vector in the non-linear solve

            unsigned int _maxIterations = 20; //!< The maximum number of allowable iterations

            unsigned int _iteration = 0; //!< The current iteration of the non-linear problem

            DataStorage< floatVector > _tolerance; //!< The tolerance vector for the non-linear solve

            // END NONLINEAR FUNCTIONS

    };

}

#endif
