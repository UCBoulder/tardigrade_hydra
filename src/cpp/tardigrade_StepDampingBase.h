/**
 ******************************************************************************
 * \file tardigrade_StepDampingBase.h
 ******************************************************************************
 * The base class for step damping operations
 ******************************************************************************
 */

#ifndef TARDIGRADE_STEPDAMPINGBASE
#define TARDIGRADE_STEPDAMPINGBASE

#include"tardigrade_CoreDefinitions.h"
#include"tardigrade_SetDataStorage.h"

namespace tardigradeHydra{

    namespace unit_test{

        class StepDampingBaseTester; //!< Forward declaration of the unit tester for the class

    }
    /*!
     * The base class for step damping operations to improve
     * stability
     */
    class StepDampingBase : public CachingDataBase{

        public:

            /*!
             * The constructor for StepDampingBase
             */
            StepDampingBase( ) : step(NULL){

            }

            /*!
             * The constructor for StepDampingBase
             *
             * \param *_step: The containing step object
             */
            StepDampingBase( SolverStepBase *_step ) : step( _step ){

            }

            SolverStepBase *step; //!< The containing step class

            /*!
             * Apply the damping to the proposed step
             */
            void applyDamping( );

            // PASS-THROUGH functions

            const floatVector *getResidual( );

            // END PASS-THROUGH FUNCTIONS

            // LINESEARCH FUNCTIONS (MOVE TO OWN CLASS)

            virtual bool checkLSConvergence( );

            void setLSAlpha( const floatType &value );

            //! Get the maximum number of line-search iterations
            constexpr unsigned int getMaxLSIterations( ){ return _maxLSIterations; }

            void setMaxLSIterations( const unsigned int &value );

            //! Get the line-search alpha parameter
            constexpr floatType getLSAlpha( ){ return _lsAlpha; }

//            virtual void performArmijoTypeLineSearch( const floatVector &X0, const floatVector &deltaX );

            //! Get the current value of the line-search iteration
            const unsigned int getLSIteration( ){ return _LSIteration; }

            void resetLSIteration( );

//            bool checkLSIteration( );

            //! Get the linesearch lambda parameter
            const floatType getLambda( ){ return _lambda; }

//            //! Get the number of line search steps performed
//            unsigned int getNumLS( ){ return _NUM_LS; }

            const floatType* getLSResidualNorm( );

            void setLSResidualNorm( ); //TODO: Move this to protected
            // END LINESEARCH FUNCTIONS

        protected:

            // LINESEARCH PARAMETERS (MOVE TO OWN CLASS)

            //! Update the line-search lambda parameter
            virtual void updateLambda( ){ _lambda *= 0.5; }

//            //! Reset the number of line search steps
//            void resetNumLS( ){ _NUM_LS = 0; }
//
//            //! Increment the number of line search steps
//            void incrementNumLS( ){ _NUM_LS++; }
//
//            // END LINESEARCH PARAMETERS

        private:

            friend class tardigradeHydra::unit_test::StepDampingBaseTester; //!< The unit tester for the class
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

    };

}

#endif
