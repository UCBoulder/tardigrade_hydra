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

        protected:

        private:

    };

}

#endif
