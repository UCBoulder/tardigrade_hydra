/**
 ******************************************************************************
 * \file tardigrade_GradientDeformationEvolutionBase.h
 ******************************************************************************
 * The base class for defining an evolution equation for the gradient of the
 * deformation
 ******************************************************************************
 */

#ifndef TARDIGRADE_GRADIENTDEFORMATIONEVOLUTIONBASE_H
#define TARDIGRADE_GRADIENTDEFORMATIONEVOLUTIONBASE_H

#include "tardigrade_DeformationEvolutionBase.h"

namespace tardigradeHydra {

    /*!
     * Base class for defining the evolution of a gradient of the deformation
     *
     */
    template <class container, int size>
    class DeformationEvolutionBase : public ResidualBase<container> {
       public:
        using tardigradeHydra::DeformationEvolutionBase<container, size>::DeformationEvolutionBase;

        //! The spatial dimension
        using tardigradeHydra::ResidualBase<container>::dimension;

       protected:
       private:
    };

}  // namespace tardigradeHydra

#include "tardigrade_GradientDeformationEvolutionBase.tpp"

#endif
