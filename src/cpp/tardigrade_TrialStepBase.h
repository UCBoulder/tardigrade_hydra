/**
 ******************************************************************************
 * \file tardigrade_TrialStepBase.h
 ******************************************************************************
 * The base class for trial step proposal operations
 ******************************************************************************
 */

#ifndef TARDIGRADE_TRIALSTEPBASE
#define TARDIGRADE_TRIALSTEPBASE

#include"tardigrade_CoreDefinitions.h"
#include"tardigrade_SetDataStorage.h"

namespace tardigradeHydra{

    /*!
     * The base class for step damping operations to improve
     * stability
     */
    class TrialStepBase : public CachingDataBase{

        public:

            /*!
             * The constructor for TrialStepBase
             */
            TrialStepBase( ) : step(NULL){

            }

            /*!
             * The constructor for TrialStepBase
             *
             * \param *_step: The containing step object
             */
            TrialStepBase( SolverStepBase *_step ) : step( _step ){

            }

            SolverStepBase *step; //!< The containing step class

            virtual void resetCounts( );

            virtual void reset( );

            virtual void computeTrial( );

        protected:

        private:

    };

}

#endif
