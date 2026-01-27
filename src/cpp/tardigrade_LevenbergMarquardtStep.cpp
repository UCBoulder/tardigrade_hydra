/**
  ******************************************************************************
  * \file tardigrade_LevenbergMarquardtStep.cpp
  ******************************************************************************
  * A C++ library for the Levenberg Marquardt step
  ******************************************************************************
  */

#include"tardigrade_LevenbergMarquardtStep.h"
#include"tardigrade_GradientDamping.h"

namespace tardigradeHydra{

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


}
