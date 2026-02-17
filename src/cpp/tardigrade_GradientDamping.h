/**
 ******************************************************************************
 * \file tardigrade_GradientDamping.h
 ******************************************************************************
 * A C++ library for gradient damping
 ******************************************************************************
 */

#ifndef TARDIGRADE_GRADIENTDAMPING_H
#define TARDIGRADE_GRADIENTDAMPING_H

#include "tardigrade_CoreDefinitions.h"
#include "tardigrade_StepDampingBase.h"

namespace tardigradeHydra {

    namespace unit_test {

        class GradientDampingTester;

    }

    /*!
     * The Gradient step class
     */
    class GradientDamping : virtual public StepDampingBase {
       public:
        using tardigradeHydra::StepDampingBase::StepDampingBase;

       protected:
        virtual void setBaseQuantities() override;

       private:
        friend tardigradeHydra::unit_test::GradientDampingTester;
    };

}  // namespace tardigradeHydra

#include "tardigrade_GradientDamping.tpp"

#endif
