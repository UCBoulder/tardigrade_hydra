/**
  ******************************************************************************
  * \file tardigrade_hydraMicromorphicRadialReturnLinearElasticity.cpp
  ******************************************************************************
  * An implementation of micromorphic linear elasticity using the hydra
  * framework which is compatible with radial return algorithms.
  ******************************************************************************
  */

#include<tardigrade_hydraMicromorphicRadialReturnLinearElasticity.h>

namespace tardigradeHydra{

    namespace micromorphicRadialReturnLinearElasticity{

        void residual::setTrialStress( ){
            /*!
             * Set the value of the trial stress
             */

            auto trialStress = get_setDataStorage_trialStress( );
            trialStress.zero( ( unsigned int )( std::end( *getStress( ) ) - std::begin( *getStress( ) ) ) ); 

            std::copy(
                    std::cbegin( *getStress( ) ),
                    std::cend( *getStress( ) ),
                    std::begin( *trialStress.value )
            );

        }

        void residual::setdTrialStressdD( ){
            /*!
             * Set the value of the derivative of the trial stress w.r.t. the deformation
             */

            auto dTrialStressdD = get_setDataStorage_dTrialStressdD( );

        }

        void residual::setResidual( ){
            /*!
             * Set the value of the residual
             */

            auto residual = get_setDataStorage_residual( );

            residual.zero( *getNumEquations( ) );

            std::transform(
                std::begin( *getStress( ) ),
                std::end(   *getStress( ) ),
                std::begin( *get_trialStress( ) ),
                residual.begin( ),
                std::minus<>( )
            );

        }

        void residual::setJacobian( ){
            /*!
             * Set the value of the Jacobian
             */
        }

        void residual::setdRdD( ){
            /*!
             * Set the value of the derivative of the residual w.r.t. the deformation
             */
        }

    }

}
