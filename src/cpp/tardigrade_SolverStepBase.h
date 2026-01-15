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
#include"tardigrade_TrialStepBase.h"
#include"tardigrade_StepDampingBase.h"

namespace tardigradeHydra{

    /*!
     * Base class for Solver Steps
     */
    class SolverStepBase : public CachingDataBase {

        public:

            /*!
             * Constructor for NonlinearStepBase
             */
            SolverStepBase( ) : solver(NULL){

                initializeDefaults( );

            }

            /*!
             * Constructor for NonlinearStepBase
             *
             * \param *_solver: The containing solver object
             */
            SolverStepBase( SolverBase *_solver ) : solver(_solver){

                initializeDefaults( );

            }

            virtual void reset( );

            void incrementSolution( );

            floatVector X0; //!< The initial value of the unknown vector

            floatVector deltaX; //!< The change in the unknown vector

            /*! Set the containing solver object
             * \param *_solver: The containing solver object
             */
            void setSolver( SolverBase *_solver ){
                solver = _solver;
            }

            // CACHED DATA STORAGE OPERATIONS
            virtual void addIterationData( dataBase *data ) override;

            virtual void addNLStepData( dataBase *data ) override;
            // END CACHED DATA STORAGE OPERATIONS

            // PASS-THROUGH functions

            const floatType getRelativeTolerance( );

            const floatType getAbsoluteTolerance( );

            const unsigned int getIteration( );

            const floatVector *getResidual( );

            const unsigned int getNumUnknowns( );

            const floatVector *getUnknownVector( );

            void updateUnknownVector( const floatVector &value );

            const floatVector *getFlatJacobian( );

            const unsigned int getNumConstraints( );

            const floatVector *getConstraints( );

            const floatVector *getConstraintJacobians( );

            const floatType getToleranceScaleFactor( );

            void resetToleranceScaleFactor( );

            bool getRankDeficientError( );

            const unsigned int getFailureVerbosityLevel( );

            void addToFailureOutput( const std::string &string );

            void addToFailureOutput( const floatVector &value, bool add_endline = true );

            void addToFailureOutput( const std::vector<bool> &value, bool add_endline = true );

            void addToFailureOutput( const floatType &value, bool add_endline = true );

            // END PASS-THROUGH FUNCTIONS

            // NONLINEAR FUNCTIONS (MOVE TO OWN CLASS)
            virtual const floatVector* getNonlinearRHS( );

            virtual const floatVector* getFlatNonlinearLHS( );

            //! Get the number of Newton steps performed
            unsigned int getNumNewton( ){ return _NUM_NEWTON; }

            // END NONLINEAR FUNCTIONS

            // LEVENBERG-MARQUARDT FUNCTIONS (MOVE TO OWN CLASS)

            //! Get if the Newton step should be a LevenbergMarquardt step
            const bool getUseLevenbergMarquardt( ){ return _use_LM_step; }

            void setUseLevenbergMarquardt( const bool &value );

            const floatType *get_baseResidualNorm( );

            const floatVector *get_basedResidualNormdX( );

            //! Get the current value of mu_k
            const floatType getMuk( ){ return _mu_k; }

            //! Get the Levenberg-Marquardt mu parameter
            const floatType getLMMu( ){ return _lm_mu; }

            // END LEVENBERG-MARQUARDT FUNCTIONS

            // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            //! Return a flag for whether to use the SQP solver
            const bool getUseSQPSolver( ){ return _useSQPSolver; }

            // END SQP SOLVER FUNCTIONS

            void performPreconditionedSolve( floatVector &deltaX_tr ); // TEMP REMOVE THIS

            // GRADIENT DESCENT FUNCTIONS (MOVE TO OWN CLASS)

            virtual bool checkDescentDirection( const floatVector &dx );

            virtual void performGradientStep( const floatVector &X0 );

            virtual bool checkGradientConvergence( const floatVector &X0 );

            // END GRADIENT DESCENT FUNCTIONS

            TrialStepBase *trial_step; //!< The trial step class which proposes a step to reduce the residual
            StepDampingBase *damping; //!< The damping class which reduces the proposed step to improve stability

        protected:

            SolverBase *solver; //!< Pointer to the containing SolverBase object

            StepDampingBase _damping; //!< The default step damping
            TrialStepBase _trial_step; //!< The default trial step

            void initializeDefaults( );

            // NONLINEAR FUNCTIONS (MOVE TO OWN CLASS)

            //! Reset the number of Newton steps
            void resetNumNewton( ){ _NUM_NEWTON = 0; }

            //! Increment the number of Newton steps
            void incrementNumNewton( ){ _NUM_NEWTON++; }

            // END NONLINEAR FUNCTIONS

            // LEVENBERG-MARQUARDT FUNCTIONS (MOVE TO OWN CLASS)

            virtual void setResidualNorm( );

            virtual void setdResidualNormdX( );

            void set_baseResidualNorm( const floatType &value );

            void set_basedResidualNormdX( const floatVector &value );

            virtual void setBaseQuantities( );

            void setMuk( const floatType &value );

            void setLMMu( const floatType &value );

            floatType _mu_k = -1; //!< The Levenberg-Marquardt scaling parameter

            floatType _lm_mu = 1e-8; //!< The mu parameter for Levenberg-Marquardt iterations

            // END LEVENBERG-MARQUARDT FUNCTIONS

            // NEWTON SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            void solveNewtonUpdate( floatVector &deltaX_tr );

            // END NEWTON SOLVER FUNCTIONS

            // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            /*!
             * Set whether to use the SQP solver
             *
             * \param &value: The updated value
             */
            void setUseSQPSolver( const unsigned int &value ){ _useSQPSolver = value; }

            virtual void solveConstrainedQP( floatVector &dx, const unsigned int kmax=100 );

            virtual void initializeActiveConstraints( std::vector< bool > &active_constraints );

            virtual void assembleKKTRHSVector( const floatVector &dx, floatVector &KKTRHSVector, const std::vector< bool > &active_constraints );

            virtual void assembleKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints );

            virtual void updateKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints );

            // END SQP SOLVER FUNCTIONS

        private:

            friend class tardigradeHydra::hydraBase; //!< TEMP REMOVE THIS
            friend class tardigradeHydra::unit_test::SolverStepBaseTester; //!< The unit tester for the class
            DataStorage< floatType > _baseResidualNorm; //!< The base value of the norm of the residual

            DataStorage< floatVector > _basedResidualNormdX; //!< The base value of the derivative of the norm of the residual w.r.t. the unknown vector

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, residualNorm,       floatType,          setResidualNorm )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dResidualNormdX,    floatVector,        setdResidualNormdX )

            // NONLINEAR DATA STORAGE

            unsigned int _NUM_NEWTON = 0; //!< The number of Newton steps performed

            // END NONLINEAR DATA STORAGE

            // LM Functions (MOVE TO OWN CLASS)

            bool _use_LM_step = false; //!< Flag for whether to attempt a Levenberg-Marquardt step

            // END LM Functions

            // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            bool _useSQPSolver = false; //!< The flag for whether to use the SQP solver

            // END SQP SOLVER FUNCTIONS

    };

}

#endif
