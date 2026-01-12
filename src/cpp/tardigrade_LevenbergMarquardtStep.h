/**
  ******************************************************************************
  * \file tardigrade_LevenbergMarquardtStep.h
  ******************************************************************************
  * A C++ library for the Levenberg Marquardt step
  ******************************************************************************
  */

#ifndef TARDIGRADE_LEVENBERGMARQUARDTSTEP_H
#define TARDIGRADE_LEVENBERGMARQUARDTSTEP_H

#include"tardigrade_CoreDefinitions.h"
#include"tardigrade_GradientStep.h"

namespace tardigradeHydra{

    /*!
     * The Levenberg Marquardt step class
     */
    class LevenbergMarquardtStep : public GradientStep {

        public:

            using tardigradeHydra::GradientStep::GradientStep;

            virtual const floatVector* getNonlinearRHS( ) override;

            virtual const floatVector* getFlatNonlinearLHS( ) override;

        protected:

        private:

            DataStorage< floatVector > _nonlinearRHS; //!< The right hand side vector for the Newton solve

            DataStorage< floatVector > _flatNonlinearLHS; //!< The left hand side vector for the Newton solve


    };

}

#endif
