/**
  ******************************************************************************
  * \file tardigrade_LevenbergMarquardtStep.cpp
  ******************************************************************************
  * A C++ library for the Levenberg Marquardt step
  ******************************************************************************
  */

#include"tardigrade_LevenbergMarquardtStep.h"
#include"tardigrade_GradientDamping.h"
#include"tardigrade_ResidualBase.h"

namespace tardigradeHydra{

    /*!
     * The constructor for LevenbergMarquardtStep
     */
    LevenbergMarquardtStep::LevenbergMarquardtStep() : NewtonStep() { }

    /*!
     * The constructor for LevenbergMarquardtStep
     *
     * \param *_step: The containing step object
     */
    LevenbergMarquardtStep::LevenbergMarquardtStep(SolverStepBase *_step) : NewtonStep(_step) {
         step->setRankDeficientError(false);
    }

    /*!
     * The constructor for LevenbergMarquardtStep
     *
     * \param *_step: The containing step object
     * \param *_preconditioner_ptr: The preconditioner object used by the trial step
     */
    LevenbergMarquardtStep::LevenbergMarquardtStep(SolverStepBase *_step, PreconditionerBase *_preconditioner_ptr)
        : NewtonStep(_step,_preconditioner_ptr) {
         step->setRankDeficientError(false);
    }

    /*!
     * Get the RHS vector for the non-linear problem
     */
    const floatVector *LevenbergMarquardtStep::getNonlinearRHS( ){

        if ( !_nonlinearRHS.first ){

            const unsigned int xsize = getNumUnknowns( );

            _nonlinearRHS.first = true;

            _nonlinearRHS.second = floatVector( xsize, 0 );

            const floatVector *residual = getResidual( );

            const floatVector *jacobian = getFlatJacobian( );

            auto Jmap = tardigradeHydra::getDynamicSizeMatrixMap( jacobian->data( ), xsize, xsize );

            auto Rmap = tardigradeHydra::getDynamicSizeVectorMap( residual->data( ), xsize );

            auto RHSmap = tardigradeHydra::getDynamicSizeVectorMap( _nonlinearRHS.second.data( ), xsize );

            RHSmap = ( Jmap.transpose( ) * Rmap ).eval( );

            addIterationData( &_nonlinearRHS );

        }

        return &_nonlinearRHS.second;

    }

    /*!
     * Get the LHS matrix for the non-linear problem
     */
    const floatVector* LevenbergMarquardtStep::getFlatNonlinearLHS( ){

        TARDIGRADE_ERROR_TOOLS_CHECK(getDamping() != nullptr, "The damping has not been defined" );

        auto local_damping = dynamic_cast<tardigradeHydra::GradientDamping*>(getDamping());

        TARDIGRADE_ERROR_TOOLS_CHECK(local_damping != nullptr, "The damping method must be GradientDamping to work with Levenberg-Marquardt")

        if ( !_flatNonlinearLHS.first ){

            const unsigned int xsize = getNumUnknowns( );

            _flatNonlinearLHS.first = true;

            _flatNonlinearLHS.second = floatVector( xsize * xsize, 0 );

            const floatVector * jacobian = getFlatJacobian( );

            auto Jmap = tardigradeHydra::getDynamicSizeMatrixMap( jacobian->data( ), xsize, xsize );

            auto LHSmap = tardigradeHydra::getDynamicSizeMatrixMap( _flatNonlinearLHS.second.data( ), xsize, xsize );

            LHSmap = ( Jmap.transpose( ) * Jmap ).eval( );

            for ( unsigned int i = 0; i < xsize; i++ ){
                _flatNonlinearLHS.second[ xsize * i + i ] += local_damping->getMuk( );
            }

            addIterationData( &_flatNonlinearLHS );

        }

        return &_flatNonlinearLHS.second;

    }

    /*!
     * Enable the residual classes to project the unknown vector into the allowable space
     */
    void LevenbergMarquardtStep::enableProjection() {
        setCurrentResidualIndexMeaningful(true);

        for (auto residual_ptr = getResidualClasses()->begin(); residual_ptr != getResidualClasses()->end();
             ++residual_ptr) {
            setCurrentResidualIndex(residual_ptr - getResidualClasses()->begin());

            (*residual_ptr)->setUseProjection(true);
        }

        setCurrentResidualIndexMeaningful(false);
    }

}
