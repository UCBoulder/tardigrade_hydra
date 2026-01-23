/**
 ******************************************************************************
 * \file tardigrade_GradientDamping.cpp
 ******************************************************************************
 * A C++ library for the Gradient damping
 ******************************************************************************
 */

#include "tardigrade_GradientDamping.h"

namespace tardigradeHydra {

    /*!
     * Set the base quantities required for gradient steps
     */
    void GradientDamping::setBaseQuantities() {
        set_baseResidualNorm(*get_residualNorm());

        set_basedResidualNormdX(*get_dResidualNormdX());

        if (getMuk() < 0) {
            setMuk(0.5 * getLMMu() * (*get_baseResidualNorm()));

        } else {
            setMuk(std::fmin(getMuk(), (*get_baseResidualNorm())));
        }
    }

}  // namespace tardigradeHydra
