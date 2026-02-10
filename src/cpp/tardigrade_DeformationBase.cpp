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

    /*!
     * Get a sub-configuration \f$\bf{F}^{sc}\f$ defined as
     *
     * \f$ F^{sc}_{iI} = F^{\text{lowerIndex}}_{i\hat{I}} F^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots
     * F^{\text{upperIndex-1}}_{\bar{I}I} \f$ \param &lowerIndex: The index of the lower configuration (starts at 0
     * and goes to numConfigurations - 1) \param &upperIndex: The index of the upper configuration (starts at 0 and
     * goes to numConfigurations) Note, the configuration indicated by the index is NOT included in the
     * sub-configuration
     */
    secondOrderTensor DeformationBase::getSubConfiguration(const unsigned int &lowerIndex, const unsigned int &upperIndex) {
        return getSubConfiguration<3,3,3>(*get_configurations(), lowerIndex, upperIndex);
    }

}
