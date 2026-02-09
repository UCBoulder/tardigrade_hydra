/**
 ******************************************************************************
 * \file tardigrade_DeformationBase.cpp
 ******************************************************************************
 * A C++ library for defining deformations
 ******************************************************************************
 */

#include "tardigrade_DeformationBase.h"
#include "tardigrade_hydra.h"

namespace tardigradeHydra {

    /*!
     * Add data to the vector of values which will be cleared after each iteration
     *
     * \param *data: The dataBase object to be cleared
     */
    void DeformationBase::addIterationData(dataBase *data) {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        hydra->addIterationData(data);
    }

    /*!
     * Add data to the vector of values which will be cleared after each nonlinear step
     *
     * \param *data: The dataBase object to be cleared
     */
    void DeformationBase::addNLStepData(dataBase *data) {
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "Hydra has not been defined");
        hydra->addNLStepData(data);
    }

}
