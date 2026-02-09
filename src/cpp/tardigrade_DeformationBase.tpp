/**
 ******************************************************************************
 * \file tardigrade_DeformationBase.tpp
 ******************************************************************************
 * A C++ library for defining deformations (template definitions)
 ******************************************************************************
 */

#include "tardigrade_DeformationBase.h"
#include "tardigrade_DeformationDecompositionBase.h"

namespace tardigradeHydra {
    /*!
     * Get a sub-configuration \f$\bf{F}^{sc}\f$ defined as
     *
     * TODO: Determine if returning identity for an equal lower and upper configuration is desirable
     * TODO: Generalize for a differently shaped leading configuration
     *
     * \f$ F^{sc}_{iI} = F^{\text{lowerIndex}}_{i\hat{I}} F^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots
     * F^{\text{upperIndex-1}}_{\bar{I}I} \f$ \param &configurations: The configurations to operate on \param
     * &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1) \param
     * &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations) Note, the
     * configuration indicated by the index is NOT included in the sub-configuration
     */
    template<unsigned int leading_rows, unsigned int size, unsigned int dim>
    floatVector DeformationBase::getSubConfiguration(const floatVector &configurations, const unsigned int &lowerIndex,
                                                     const unsigned int &upperIndex) {
        constexpr unsigned int sot_dim = size * size;
        TARDIGRADE_ERROR_TOOLS_CHECK(leading_rows == size, "leading_rows and size must be the same")
        TARDIGRADE_ERROR_TOOLS_CHECK(leading_rows == dim, "leading_rows and dim must be the same")

        secondOrderTensor Fsc(dim * dim, 0);

        if (lowerIndex == upperIndex) {
            for (unsigned int i = 0; i < size; ++i) {
                Fsc[dim * i + i] = 1;
            }

        } else {
            DeformationDecompositionBase<leading_rows, size, dim> decomposition;
            decomposition.getNetConfiguration(std::begin(configurations) + sot_dim * lowerIndex,
                                              std::begin(configurations) + sot_dim * upperIndex, std::begin(Fsc),
                                              std::end(Fsc));
        }

        return Fsc;
    }

    /*!
     * Get the jacobian of a sub-configuration \f$\bf{F}^{sc}\f$ defined as
     *
     * \f$ F^{sc}_{iI} = F^{\text{lowerIndex}}_{i\hat{I}} F^{\text{lowerIndex} + 1}_{\hat{I}\breve{I}} \cdots
     * F^{\text{upperIndex-1}}_{\bar{I}I} \f$
     *
     * with respect to all of the configurations. The returned matrix will be of size ( dimensions**2,
     * configurations.size( ) * dimensions**2 )
     *
     * \param &configurations: The configurations to operate on
     * \param &lowerIndex: The index of the lower configuration (starts at 0 and goes to numConfigurations - 1)
     * \param &upperIndex: The index of the upper configuration (starts at 0 and goes to numConfigurations)
     *   Note, the configuration indicated by the index is NOT included in the sub-configuration
     */
    template<unsigned int leading_rows, unsigned int size, unsigned int dim>
    floatVector DeformationBase::getSubConfigurationJacobian(const floatVector  &configurations,
                                                       const unsigned int &lowerIndex, const unsigned int &upperIndex) {
        constexpr unsigned int sot_dim = size * size;
        TARDIGRADE_ERROR_TOOLS_CHECK(leading_rows == size, "leading_rows and size must be the same")
        TARDIGRADE_ERROR_TOOLS_CHECK(leading_rows == dim, "leading_rows and dim must be the same")
        const unsigned int     num_incoming_configs = configurations.size() / sot_dim;

        floatVector gradient(sot_dim * sot_dim * num_incoming_configs, 0);

        if (lowerIndex != upperIndex) {
            DeformationDecompositionBase<leading_rows,size,dim> decomposition;
            decomposition.getNetConfigurationJacobian(std::begin(configurations) + sot_dim * lowerIndex,
                                                      std::begin(configurations) + sot_dim * upperIndex,
                                                      std::begin(gradient), std::end(gradient), lowerIndex * sot_dim,
                                                      num_incoming_configs * sot_dim);
        }

        return gradient;
    }

}
