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

            const floatType *get_baseResidualNorm( );

            const floatVector *get_basedResidualNormdX( );

            //!< Get the current value of mu_k
            const floatType getMuk( ){ return _mu_k; }

            //!< Set the Levenberg-Marquardt mu_k
            void setMuk( const floatType &value ){
               /*!
                * Set the value of the mu_k parameter for Levenberg-Marquardt steps
                *
                * \param &value: The value of the parameter
                */
 
                _mu_k = value;

            }

            //!< Get the Levenberg-Marquardt mu parameter
            const floatType getLMMu( ){ return _lm_mu; }

            //!< Set the Levenberg-Marquardt mu parameter
            void setLMMu( const floatType &value ){
               /*!
                * Set the value of the mu parameter for Levenberg-Marquardt steps
                *
                * \param &value: The value of the parameter
                */
 
                _lm_mu = value;

            }
        protected:

            SolverBase *solver; //!< Pointer to the containing SolverBase object

            void set_baseResidualNorm( const floatType &value );

            void set_basedResidualNormdX( const floatVector &value );

            floatType _mu_k = -1; //!< The Levenberg-Marquardt scaling parameter

            floatType _lm_mu = 1e-8; //!< The mu parameter for Levenberg-Marquardt iterations

        private:

            friend class tardigradeHydra::hydraBase; //!< TEMP REMOVE THIS
            friend class tardigradeHydra::unit_test::SolverStepBaseTester; //!< The unit tester for the class
            DataStorage< floatType > _baseResidualNorm; //!< The base value of the norm of the residual

            DataStorage< floatVector > _basedResidualNormdX; //!< The base value of the derivative of the norm of the residual w.r.t. the unknown vector

    };

}

#include"tardigrade_SolverStepBase.cpp"

#endif
