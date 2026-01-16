/**
 ******************************************************************************
 * \file tardigrade_ArmijoGradientDamping.h
 ******************************************************************************
 * A class for combined line-search gradient damping
 ******************************************************************************
 */

#ifndef TARDIGRADE_ARMIJOGRADIENTDAMPING
#define TARDIGRADE_ARMIJOGRADIENTDAMPING

#include"tardigrade_CoreDefinitions.h"
#include"tardigrade_StepDampingBase.h"

namespace tardigradeHydra{

    namespace unit_test{

        class ArmijoGradientDampingTester; //!< Forward declaration of the unit tester for the class

    };

    /*!
     * The base class for step damping operations to improve
     * stability
     */
    class ArmijoGradientDamping : public StepDampingBase{

        public:

            using tardigradeHydra::StepDampingBase::StepDampingBase;

            virtual void setBaseQuantities( ) override;

            virtual void resetCounts( ) override;

            virtual void reset( ) override;

            const bool applyDamping( ) override;

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

            //! Get the number of line search steps performed
            unsigned int getNumLS( ){ return _NUM_LS; }

            const floatType* getLSResidualNorm( );
            // END LINESEARCH FUNCTIONS
//
//            // GRADIENT DESCENT FUNCTIONS (MOVE TO OWN CLASS)
//
//            //!< Get whether Gradient descent is allowed
//            const bool getUseGradientDescent( ){ return _use_gradient_descent; }
//
//            void setUseGradientDescent( const bool &value );
//
//            //! Get the gradient descent rho parameter
//            const floatType getGradientRho( ){ return _gradientRho; }
//
//            //! Get the gradient descent p parameter
//            const floatType getGradientP( ){ return _gradientP; }
//
//            //! Get the gradient descent beta parameter
//            const floatType getGradientBeta( ){ return _gradientBeta; }
//
//            //! Get the gradient descent sigma parameter
//            const floatType getGradientSigma( ){ return _gradientSigma; }
//
//            void setGradientRho( const floatType &value );
//
//            void setGradientP( const floatType &value );
//
//            void setGradientBeta( const floatType &value );
//
//            void setGradientSigma( const floatType &value );
//
//            virtual bool checkDescentDirection( const floatVector &dx );
//
//            //! Get the max allowable number of gradient iterations
//            const unsigned int getMaxGradientIterations( ){ return _maxGradientIterations; }
//
//            //! Get the current gradient iteration
//            const unsigned int getGradientIteration( ){ return _gradientIteration; }
//
//            //! Reset the number of gradient descent steps
//            void resetGradientIteration( ){ _gradientIteration = 0; }
//
//            void setMaxGradientIterations( const unsigned int &value );
//
//            bool checkGradientIteration( );
//
//            virtual bool checkGradientConvergence( const floatVector &X0 );
//
//            //! Get the number of gradient descent steps performed
//            unsigned int getNumGrad( ){ return _NUM_GRAD; }
//
//            virtual void performGradientStep( const floatVector &X0 );
//
//            //! Get the current value of mu_k
//            const floatType getMuk( ){ return _mu_k; }
//
//            //! Get the Levenberg-Marquardt mu parameter
//            const floatType getLMMu( ){ return _lm_mu; }
//
//            const floatType *get_baseResidualNorm( );
//
//            const floatVector *get_basedResidualNormdX( );
//
//            // END GRADIENT DESCENT FUNCTIONS
        protected:

            // LINESEARCH PARAMETERS (MOVE TO OWN CLASS)

            //! Update the line-search lambda parameter
            virtual void updateLambda( ){ _lambda *= 0.5; }

            void setLSResidualNorm( );

            //! Reset the number of line search steps
            void resetNumLS( ){ _NUM_LS = 0; }

            //! Increment the number of line search steps
            void incrementNumLS( ){ _NUM_LS++; }

            // END LINESEARCH PARAMETERS

//            // GRADIENT DESCENT FUNCTIONS (MOVE TO OWN CLASS)
//
//            //! Increment the number of gradient descent steps
//            void incrementGradientIteration( ){ _gradientIteration++; }
//
//            //! Reset the number of gradient descent steps
//            void resetNumGrad( ){ _NUM_GRAD = 0; }
//
//            //! Increment the number of gradient descent steps
//            void incrementNumGrad( ){ _NUM_GRAD++; }
//
//            virtual void setResidualNorm( );
//
//            virtual void setdResidualNormdX( );
//
//            void setMuk( const floatType &value );
//
//            void setLMMu( const floatType &value );
//
//            void set_baseResidualNorm( const floatType &value );
//
//            void set_basedResidualNormdX( const floatVector &value );
//
//            // END GRADIENT DESCENT FUNCTIONS

        private:

            friend class tardigradeHydra::unit_test::ArmijoGradientDampingTester; //!< The unit tester for the class
            friend class tardigradeHydra::SolverStepBase; //!< TODO: REMOVE THIS
            // LS Functions (MOVE TO OWN CLASS)

            floatType _lsAlpha; //!< The line-search alpha value i.e., the term by which it is judged that the line-search is converging

            unsigned int _maxLSIterations; //!< The maximum number of line-search iterations

            unsigned int _LSIteration = 0; //!< The current line search iteration of the non-linear problem

            void incrementLSIteration( ){ _LSIteration++; }

            floatType _lambda = 1;

            DataStorage< floatType > _lsResidualNorm; //!< The reference residual norm for the line-search convergence criteria

            unsigned int _NUM_LS = 0; //!< The number of line search steps performed

            // END LS FUNCTIONS
//
//            // GRADIENT DESCENT FUNCTIONS (MOVE TO OWN CLASS)
//
//            bool _use_gradient_descent = false; //!< Flag for whether to attempt a gradient descent step
//
//            floatType _gradientRho   = 1e-8; //!< The rho parameter for the gradient descent step
//
//            floatType _gradientP     = 2.1; //!< The p parameter for the gradient descent step
//
//            floatType _gradientBeta  = 0.9; //!< The beta parameter for the gradient descent step
//
//            floatType _gradientSigma = 1e-4; //!< The sigma parameter for the gradient descent step
//
//            unsigned int _gradientIteration = 0; //!< The current gradient iteration of the non-linear problem
//
//            unsigned int _maxGradientIterations = 10; //!< The maximum number of gradient iterations
//
//            unsigned int _NUM_GRAD = 0; //!< The number of gradient descent steps performed
//
//            floatType _mu_k = -1; //!< The Gradient-descent scaling parameter
//
//            floatType _lm_mu = 1e-8; //!< The mu parameter for Levenberg-Marquardt iterations
//
//            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, residualNorm,       floatType,          setResidualNorm )
//
//            TARDIGRADE_HYDRA_DECLARE_ITERATION_STORAGE( private, dResidualNormdX,    floatVector,        setdResidualNormdX )
//
//            DataStorage< floatType > _baseResidualNorm; //!< The base value of the norm of the residual
//
//            DataStorage< floatVector > _basedResidualNormdX; //!< The base value of the derivative of the norm of the residual w.r.t. the unknown vector
//            // END GRADIENT DESCENT FUNCTIONS
    };

}

#endif
