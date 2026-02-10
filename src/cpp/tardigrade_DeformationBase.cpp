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

    /*!
     * Get the sub-configuration preceding but not including the index
     *
     * \param &index: The index of the configuration immediately following the sub-configuration
     */
    secondOrderTensor DeformationBase::getPrecedingConfiguration(const unsigned int &index) {
        return getSubConfiguration(0, index);
    }

    /*!
     * Get the sub-configuration following but not including the index
     *
     * \param &index: The index of the current configuration immediately before the sub-configuration
     */
    secondOrderTensor DeformationBase::getFollowingConfiguration(const unsigned int &index) {
        return getSubConfiguration(index + 1, getNumConfigurations());
    }

    /*!
     * Get the configuration indicated by the provided index
     *
     * \param &index: The index of the current configuration to be extracted
     */
    secondOrderTensor DeformationBase::getConfiguration(const unsigned int &index) {
        return getSubConfiguration(index, index + 1);
    }

    /*!
     * Get a previous sub-configuration \f$\bf{F}^{sc}\f$ defined as
     *
     * \f$ F^{sc}_{iI} = F^{\text{lowerIndex}}_{i\hat{I}} F^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots
     * F^{\text{upperIndex-1}}_{\bar{I}I} \f$ \param &lowerIndex: The index of the lower configuration (starts at 0
     * and goes to numConfigurations - 1) \param &upperIndex: The index of the upper configuration (starts at 0 and
     * goes to numConfigurations) Note, the configuration indicated by the index is NOT included in the
     * sub-configuration
     */
    secondOrderTensor DeformationBase::getPreviousSubConfiguration(const unsigned int &lowerIndex,
                                                                   const unsigned int &upperIndex) {
        return getSubConfiguration<3,3,3>(*get_previousConfigurations(), lowerIndex, upperIndex);
    }

    /*!
     * Get the previous sub-configuration preceding but not including the index
     *
     * \param &index: The index of the configuration immediately following the sub-configuration
     */
    secondOrderTensor DeformationBase::getPreviousPrecedingConfiguration(const unsigned int &index) {
        return getPreviousSubConfiguration(0, index);
    }

    /*!
     * Get the previous sub-configuration following but not including the index
     *
     * \param &index: The index of the current configuration immediately before the sub-configuration
     */
    secondOrderTensor DeformationBase::getPreviousFollowingConfiguration(const unsigned int &index) {
        return getPreviousSubConfiguration(index + 1, getNumConfigurations());
    }

    /*!
     * Get the number of configurations
     */
    unsigned int DeformationBase::getNumConfigurations(){
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "The containing hydraBase class has not been set")
        return hydra->getNumConfigurations();
    }
}
