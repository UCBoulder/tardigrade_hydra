/**
 ******************************************************************************
 * \file tardigrade_ResidualBase.cpp
 ******************************************************************************
 * The base class for residuals
 ******************************************************************************
 */

#include"tardigrade_ResidualBase.h"
#include"tardigrade_hydra.h"

namespace tardigradeHydra{

    void residualBase::addIterationData( dataBase *data ){
        /*!
         * Add data to the vector of values which will be cleared after each iteration
         * 
         * \param *data: The dataBase object to be cleared
         */

        hydra->addIterationData( data );

    }

    void residualBase::addNLStepData( dataBase *data ){
        /*!
         * Add data to the vector of values which will be cleared after each nonlinear step
         * 
         * \param *data: The dataBase object to be cleared
         */

        hydra->addNLStepData( data );

    }

    void residualBase::setupRelaxedStep( const unsigned int &relaxedStep ){
        /*!
         * When performing a relaxed iteration this function is called prior to the solution of the non-linear
         * problem. Users can use this function to dynamically adjust parameters or perform other tuning tasks.
         *
         * \param &relaxedStep: The current relaxed step.
         */
    }

}
