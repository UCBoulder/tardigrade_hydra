/**
  ******************************************************************************
  * \file tardigrade_GradientStep.cpp
  ******************************************************************************
  * A C++ library for the Gradient step
  ******************************************************************************
  */

#include"tardigrade_GradientStep.h"

namespace tardigradeHydra{

    /*!
     * Set the base quantities required for gradient steps
     */
    void GradientStep::setBaseQuantities( ){

        set_baseResidualNorm( *get_residualNorm( ) );

        set_basedResidualNormdX( *get_dResidualNormdX( ) );

        if ( _mu_k < 0 ){

            setMuk( 0.5 * getLMMu( ) * ( *get_baseResidualNorm( ) ) );

        }
        else{

            setMuk( std::fmin( _mu_k, ( *get_baseResidualNorm( ) ) ) );

        }

    }

}
