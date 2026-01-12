/**
  ******************************************************************************
  * \file tardigrade_GradientStep.h
  ******************************************************************************
  * A C++ library for the gradient step
  ******************************************************************************
  */

#ifndef TARDIGRADE_GRADIENTSTEP_H
#define TARDIGRADE_GRADIENTSTEP_H

#include"tardigrade_CoreDefinitions.h"
#include"tardigrade_SolverBase.h"

namespace tardigradeHydra{

    /*!
     * The Gradient step class
     */
    class GradientStep : public SolverStepBase {

        public:

            using tardigradeHydra::SolverStepBase::SolverStepBase;

        protected:

            virtual void setBaseQuantities( ) override;

        private:

    };

}

#endif
