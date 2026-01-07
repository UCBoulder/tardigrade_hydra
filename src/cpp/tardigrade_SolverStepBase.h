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
    class SolverStepBase : public CachingDataBase {

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

            floatVector deltaX; //!< The change in the unknown vector

            void setSolver( SolverBase *_solver ){
                /*! Set the containing solver object
                 * \param *_solver: The containing solver object
                 */
                solver = _solver;
            }

            // CACHED DATA STORAGE OPERATIONS
            virtual void addIterationData( dataBase *data ) override;

            virtual void addNLStepData( dataBase *data ) override;
            // END CACHED DATA STORAGE OPERATIONS

            // NONLINEAR FUNCTIONS (MOVE TO OWN CLASS)
            virtual const floatVector* getNonlinearRHS( );

            virtual const floatVector* getFlatNonlinearLHS( );
            // END NONLINEAR FUNCTIONS

            // GRADIENT FUNCTIONS (MOVE TO OWN CLASS)

            //!< Get whether Gradient descent is allowed
            const bool getUseGradientDescent( ){ return _use_gradient_descent; }

            void setUseGradientDescent( const bool &value );

            // END GRADIENT FUNCTIONS

            // LEVENBERG-MARQUARDT FUNCTIONS (MOVE TO OWN CLASS)

            //!< Get the Newton step should be a LevenbergMarquardt step
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

            const bool getUseSQPSolver( ){ /*! Return a flag for whether to use the SQP solver */ return _useSQPSolver; }

            // END SQP SOLVER FUNCTIONS

            void performPreconditionedSolve( floatVector &deltaX_tr ); // TEMP REMOVE THIS

            // LINESEARCH FUNCTIONS (MOVE TO OWN CLASS)

            virtual bool checkLSConvergence( );

            void setLSAlpha( const floatType &value );

            //! Get the maximum number of line-search iterations
            constexpr unsigned int getMaxLSIterations( ){ return _maxLSIterations; }

            void setMaxLSIterations( const unsigned int &value );

            //! Get the line-search alpha parameter
            constexpr floatType getLSAlpha( ){ return _lsAlpha; }

            virtual void performArmijoTypeLineSearch( const floatVector &X0, const floatVector &deltaX );

            //! Get the current value of the line-search iteration
            const unsigned int getLSIteration( ){ return _LSIteration; }

            void resetLSIteration( );

            bool checkLSIteration( );

            //! Get the linesearch lambda parameter
            const floatType getLambda( ){ return _lambda; }

            unsigned int getNumLS( ){ /*! Get the number of line search steps performed */ return _NUM_LS; }

            const floatType* getLSResidualNorm( );
            // END LINESEARCH FUNCTIONS

            // GRADIENT DESCENT FUNCTIONS (MOVE TO OWN CLASS)

            //! Get the gradient descent rho parameter
            const floatType getGradientRho( ){ return _gradientRho; }

            //! Get the gradient descent p parameter
            const floatType getGradientP( ){ return _gradientP; }

            //! Get the gradient descent beta parameter
            const floatType getGradientBeta( ){ return _gradientBeta; }

            //! Get the gradient descent sigma parameter
            const floatType getGradientSigma( ){ return _gradientSigma; }

            void setGradientRho( const floatType &value );

            void setGradientP( const floatType &value );

            void setGradientBeta( const floatType &value );

            void setGradientSigma( const floatType &value );

            virtual bool checkDescentDirection( const floatVector &dx );

            //! Get the max allowable number of gradient iterations
            const unsigned int getMaxGradientIterations( ){ return _maxGradientIterations; }

            //! Get the current gradient iteration
            const unsigned int getGradientIteration( ){ return _gradientIteration; }

            void setMaxGradientIterations( const unsigned int &value );

            bool checkGradientIteration( );

            virtual bool checkGradientConvergence( const floatVector &X0 );

            // END GRADIENT DESCENT FUNCTIONS

        protected:

            SolverBase *solver; //!< Pointer to the containing SolverBase object

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

            // NEWTON SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            void solveNewtonUpdate( floatVector &deltaX_tr );

            // END NEWTON SOLVER FUNCTIONS

            // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            void setUseSQPSolver( const unsigned int &value ){ /*! Set whether to use the SQP solver \param &value: The updated value */ _useSQPSolver = value; }

            virtual void solveConstrainedQP( floatVector &dx, const unsigned int kmax=100 );

            virtual void initializeActiveConstraints( std::vector< bool > &active_constraints );

            virtual void assembleKKTRHSVector( const floatVector &dx, floatVector &KKTRHSVector, const std::vector< bool > &active_constraints );

            virtual void assembleKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints );

            virtual void updateKKTMatrix( floatVector &KKTMatrix, const std::vector< bool > &active_constraints );

            // END SQP SOLVER FUNCTIONS

            // LINESEARCH PARAMETERS (MOVE TO OWN CLASS)

            //! Update the line-search lambda parameter
            virtual void updateLambda( ){ _lambda *= 0.5; }

            //! Reset the number of line search steps
            void resetNumLS( ){ _NUM_LS = 0; }

            //! Increment the number of line search steps
            void incrementNumLS( ){ _NUM_LS++; }

            // END LINESEARCH PARAMETERS

            // GRADIENT DESCENT FUNCTIONS (MOVE TO OWN CLASS)

            //! Increment the number of gradient descent steps
            void incrementGradientIteration( ){ _gradientIteration++; }

            //! Reset the number of gradient descent steps
            void resetGradientIteration( ){ _gradientIteration = 0; }

            // END GRADIENT DESCENT FUNCTIONS

        private:

            friend class tardigradeHydra::hydraBase; //!< TEMP REMOVE THIS
            friend class tardigradeHydra::unit_test::SolverStepBaseTester; //!< The unit tester for the class
            DataStorage< floatType > _baseResidualNorm; //!< The base value of the norm of the residual

            DataStorage< floatVector > _basedResidualNormdX; //!< The base value of the derivative of the norm of the residual w.r.t. the unknown vector

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, residualNorm,       floatType,          setResidualNorm )

            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dResidualNormdX,    floatVector,        setdResidualNormdX )

            // NONLINEAR DATA STORAGE

            DataStorage< floatVector > _nonlinearRHS; //!< The right hand side vector for the Newton solve

            DataStorage< floatVector > _flatNonlinearLHS; //!< The left hand side vector for the Newton solve

            // END NONLINEAR DATA STORAGE

            // GRADIENT FUNCTIONS (MOVE TO OWN CLASS)

            bool _use_gradient_descent = false; //!< Flag for whether to attempt a gradient descent step

            // END GRADIENT FUNCTIONS

            // LM Functions (MOVE TO OWN CLASS)

            bool _use_LM_step = false; //!< Flag for whether to attempt a Levenberg-Marquardt step

            // END LM Functions

            // SQP SOLVER FUNCTIONS (MOVE TO OWN CLASS)

            bool _useSQPSolver = false; //!< The flag for whether to use the SQP solver

            // END SQP SOLVER FUNCTIONS

            // LS Functions (MOVE TO OWN CLASS)

            floatType _lsAlpha; //!< The line-search alpha value i.e., the term by which it is judged that the line-search is converging

            unsigned int _maxLSIterations; //!< The maximum number of line-search iterations

            unsigned int _LSIteration = 0; //!< The current line search iteration of the non-linear problem

            void incrementLSIteration( ){ _LSIteration++; }

            floatType _lambda = 1;

            DataStorage< floatType > _lsResidualNorm; //!< The reference residual norm for the line-search convergence criteria

            unsigned int _NUM_LS = 0; //!< The number of line search steps performed

            // END LS FUNCTIONS

            // GRADIENT DESCENT FUNCTIONS (MOVE TO OWN CLASS)

            floatType _gradientRho   = 1e-8; //!< The rho parameter for the gradient descent step

            floatType _gradientP     = 2.1; //!< The p parameter for the gradient descent step

            floatType _gradientBeta  = 0.9; //!< The beta parameter for the gradient descent step

            floatType _gradientSigma = 1e-4; //!< The sigma parameter for the gradient descent step

            unsigned int _gradientIteration = 0; //!< The current gradient iteration of the non-linear problem

            unsigned int _maxGradientIterations = 10; //!< The maximum number of gradient iterations

            // END GRADIENT DESCENT FUNCTIONS
    };

}

#endif
