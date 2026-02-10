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
     * Get the previous configuration indicated by the provided index
     *
     * \param &index: The index of the current configuration to be extracted
     */
    secondOrderTensor DeformationBase::getPreviousConfiguration(const unsigned int &index) {
        return getPreviousSubConfiguration(index, index + 1);
    }

    /*!
     * Get the jacobian of a sub-configuration \f$\bf{F}^{sc}\f$ defined as
     *
     * \f$ F^{sc}_{iI} = F^{\text{lowerIndex}}_{i\hat{I}} F^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots
     * F^{\text{upperIndex-1}}_{\bar{I}I} \f$
     *
     * with respect to the current configurations.
     *
     * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
     * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
     *   Note, the configuration indicated by the index is NOT included in the sub-configuration
     */
    floatVector DeformationBase::getSubConfigurationJacobian(const unsigned int &lowerIndex, const unsigned int &upperIndex) {
        return getSubConfigurationJacobian<3,3,3>(*get_configurations(), lowerIndex, upperIndex);
    }

    /*!
     * Get the jacobian of the sub-configuration preceding but not including the index with respect to the current
     * configurations.
     *
     * \param &index: The index of the configuration immediately following the sub-configuration
     */
    floatVector DeformationBase::getPrecedingConfigurationJacobian(const unsigned int &index) {
        return getSubConfigurationJacobian(0, index);
    }

    /*!
     * Get the jacobian of the sub-configuration following but not including the index with respect to the current
     * configurations.
     *
     * \param &index: The index of the current configuration immediately before the sub-configuration
     */
    floatVector DeformationBase::getFollowingConfigurationJacobian(const unsigned int &index) {
        return getSubConfigurationJacobian(index + 1, getNumConfigurations());
    }

    /*!
     * Get the number of configurations
     */
    unsigned int DeformationBase::getNumConfigurations(){
        TARDIGRADE_ERROR_TOOLS_CHECK(hydra != nullptr, "The containing hydraBase class has not been set")
        return hydra->getNumConfigurations();
    }

    /*!
     * Get the jacobian of a previous sub-configuration \f$\bf{F}^{sc}\f$ defined as
     *
     * \f$ F^{sc}_{iI} = F^{\text{lowerIndex}}_{i\hat{I}} F^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots
     * F^{\text{upperIndex-1}}_{\bar{I}I} \f$
     *
     * with respect to the previous configurations.
     *
     * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
     * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
     *   Note, the configuration indicated by the index is NOT included in the sub-configuration
     */
    floatVector DeformationBase::getPreviousSubConfigurationJacobian(const unsigned int &lowerIndex,
                                                               const unsigned int &upperIndex) {
        return getSubConfigurationJacobian<3,3,3>(*get_previousConfigurations(), lowerIndex, upperIndex);
    }

    /*!
     * Get the jacobian of the previous sub-configuration following but not including the index with
     * respect to the previous configurations
     *
     * \param &index: The index of the current configuration immediately before the sub-configuration
     */
    floatVector DeformationBase::getPreviousFollowingConfigurationJacobian(const unsigned int &index) {
        return getPreviousSubConfigurationJacobian(index + 1, getNumConfigurations());
    }

    /*!
     * Get the Jacobian of the first configuration w.r.t. the total mapping and the remaining configurations.
     *
     * \param &configurations: The configurations which describe the mapping from the current to the reference
     * configuration \param &dC1dC: The Jacobian of the first entry w.r.t. the total \param &dC1dCn: The Jacobian of
     * the first entry w.r.t. the remaining terms
     *
     * whre \f$C^n = C^2, C^3, \cdots \f$
     */
    void DeformationBase::calculateFirstConfigurationJacobians(const floatVector &configurations, fourthOrderTensor &dC1dC,
                                                         floatVector &dC1dCn) {
        constexpr unsigned int dim         = 3;
        constexpr unsigned int sot_dim     = dim * dim;
        auto                   num_configs = getNumConfigurations();

        secondOrderTensor fullConfiguration = getSubConfiguration<3,3,3>(configurations, 0, num_configs);

        dC1dC  = secondOrderTensor(sot_dim * sot_dim, 0);
        dC1dCn = floatVector(sot_dim * (num_configs - 1) * sot_dim, 0);

        DeformationDecompositionBase<dim, dim, dim> decomposition;
        decomposition.solveForLeadingConfigurationTotalConfigurationJacobian(std::begin(fullConfiguration),
                                                                             std::end(fullConfiguration),
                                                                             std::begin(configurations) + sot_dim,
                                                                             std::end(configurations),
                                                                             std::begin(dC1dC), std::end(dC1dC));
        decomposition.solveForLeadingConfigurationConfigurationJacobian(std::begin(fullConfiguration),
                                                                        std::end(fullConfiguration),
                                                                        std::begin(configurations) + sot_dim,
                                                                        std::end(configurations), std::begin(dC1dCn),
                                                                        std::end(dC1dCn));
    }

    /*!
     * Set the Jacobians of the first configuration w.r.t. the total configuration and the remaining
     * sub-configurations
     */
    void DeformationBase::setFirstConfigurationJacobians() {
        auto dF1dF = get_SetDataStorage_dF1dF();

        auto dF1dFn = get_SetDataStorage_dF1dFn();

        calculateFirstConfigurationJacobians(*get_configurations(), *dF1dF.value, *dF1dFn.value);
    }

    /*!
     * Set the Jacobians of the previous first configuration w.r.t. the total configuration and the remaining
     * sub-configurations
     */
    void DeformationBase::setPreviousFirstConfigurationJacobians() {
        auto dF1dF = get_SetDataStorage_previousdF1dF();

        auto dF1dFn = get_SetDataStorage_previousdF1dFn();

        calculateFirstConfigurationJacobians(*get_previousConfigurations(), *dF1dF.value, *dF1dFn.value);
    }

    /*!
     * Get the jacobian of the previous sub-configuration preceding but not including the index with
     * respect to the previous configurations.
     *
     * \param &index: The index of the configuration immediately following the sub-configuration
     */
    floatVector DeformationBase::getPreviousPrecedingConfigurationJacobian(const unsigned int &index) {
        return getPreviousSubConfigurationJacobian(0, index);
    }

}
