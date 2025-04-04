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

            getStress( );

            trialStress.zero( ( unsigned int )( std::end( *getStress( ) ) - std::begin( *getStress( ) ) ) ); 

            std::copy(
                    std::cbegin( *getStress( ) ),
                    std::cend( *getStress( ) ),
                    trialStress.begin( )
            );

            setdTrialStressdD( ); // Do this so the trial stress and its Jacobian are consistent

        }

        void residual::setdTrialStressdD( ){
            /*!
             * Set the value of the derivative of the trial stress w.r.t. the deformation
             */

            auto dTrialStressdD = get_setDataStorage_dTrialStressdD( );

            auto sot_dim = hydra->getSOTDimension( );

            auto tot_dim = hydra->getTOTDimension( );

            auto num_equations = *getNumEquations( );

            auto num_configurationUnknowns = *hydra->getConfigurationUnknownCount( );

            dTrialStressdD.zero( num_equations * num_configurationUnknowns );

            //Get references to the stress Jacobians. Doing it this way to allow changing the residual to the current configuration in the future.
            auto dS1dF       = get_dPK2dF( );

            auto dS1dChi     = get_dPK2dChi( );

            auto dS1dGradChi = get_dPK2dGradChi( );

            auto dS2dF       = get_dSIGMAdF( );

            auto dS2dChi     = get_dSIGMAdChi( );

            auto dS2dGradChi = get_dSIGMAdGradChi( );

            auto dS3dF       = get_dMdF( );

            auto dS3dChi     = get_dMdChi( );

            auto dS3dGradChi = get_dMdGradChi( );

            constexpr unsigned int numStresses = 3;
            constexpr unsigned int numDeformationMeasures = 3;

            std::array< const floatVector *, numStresses * numDeformationMeasures > stressReferences = { dS1dF, dS1dChi, dS1dGradChi,
                                                                                                         dS2dF, dS2dChi, dS2dGradChi,
                                                                                                         dS3dF, dS3dChi, dS3dGradChi };

            const std::array< unsigned int, numDeformationMeasures > dims = { sot_dim, sot_dim, tot_dim };

            unsigned int row = 0;
            for ( unsigned int S = 0; S < numStresses; ++S ){

                // Loop through the stress Jacobians

                for ( unsigned int i = 0; i < dims[ S ]; ++i ){

                    unsigned int col = 0;

                    // Jacobians w.r.t. the deformation

                    for ( auto Sn = std::begin( stressReferences ) + numDeformationMeasures * S; Sn != std::begin( stressReferences ) + numDeformationMeasures * ( S + 1 ); ++Sn ){

                        unsigned int dof_index = ( unsigned int )( Sn - ( std::begin( stressReferences ) + numDeformationMeasures * S ) );

                        for ( unsigned int j = 0; j < dims[ dof_index ]; ++j ){

                            ( *dTrialStressdD.value )[ num_configurationUnknowns * row + col ] += ( **Sn )[ dims[ dof_index ] * i + j ];

                            col++;

                        }

                    }

                    row++;

                }

            }

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

            tardigradeHydra::micromorphicLinearElasticity::residual::setJacobian( );

            auto jacobian = get_setDataStorage_jacobian( );

            std::transform(
                jacobian.begin( ),
                jacobian.end( ),
                jacobian.begin( ),
                std::negate<>( )
            );

        }

        void residual::setdRdD( ){
            /*!
             * Set the derivative of the residual w.r.t. the deformation
             */

            auto dRdD = get_setDataStorage_dRdD( );

            auto sot_dim = hydra->getSOTDimension( );

            auto tot_dim = hydra->getTOTDimension( );

            auto num_equations = *getNumEquations( );

            auto num_configurationUnknowns = *hydra->getConfigurationUnknownCount( );

            dRdD.zero( num_equations * num_configurationUnknowns );

            //Get references to the stress Jacobians. Doing it this way to allow changing the residual to the current configuration in the future.
            auto dS1dF       = get_dPK2dF( );

            auto dS1dChi     = get_dPK2dChi( );

            auto dS1dGradChi = get_dPK2dGradChi( );

            auto dS2dF       = get_dSIGMAdF( );

            auto dS2dChi     = get_dSIGMAdChi( );

            auto dS2dGradChi = get_dSIGMAdGradChi( );

            auto dS3dF       = get_dMdF( );

            auto dS3dChi     = get_dMdChi( );

            auto dS3dGradChi = get_dMdGradChi( );

            constexpr unsigned int numStresses = 3;
            constexpr unsigned int numDeformationMeasures = 3;

            std::array< const floatVector *, numStresses * numDeformationMeasures > stressReferences = { dS1dF, dS1dChi, dS1dGradChi,
                                                                                                         dS2dF, dS2dChi, dS2dGradChi,
                                                                                                         dS3dF, dS3dChi, dS3dGradChi };

            const std::array< unsigned int, numDeformationMeasures > dims = { sot_dim, sot_dim, tot_dim };

            unsigned int row = 0;
            for ( unsigned int S = 0; S < numStresses; ++S ){

                // Loop through the stress Jacobians

                for ( unsigned int i = 0; i < dims[ S ]; ++i ){

                    unsigned int col = 0;

                    // Jacobians w.r.t. the deformation

                    for ( auto Sn = std::begin( stressReferences ) + numDeformationMeasures * S; Sn != std::begin( stressReferences ) + numDeformationMeasures * ( S + 1 ); ++Sn ){

                        for ( unsigned int j = 0; j < dims[ ( unsigned int )( Sn - std::begin( stressReferences ) + numDeformationMeasures * S ) ]; ++j ){

                            ( *dRdD.value )[ num_configurationUnknowns * row + col ] += ( **Sn )[ dims[ S ] * i + j ];

                            col++;

                        }

                    }

                    row++;

                }

            }
            
            std::transform(
                dRdD.begin( ),
                dRdD.end( ),
                std::begin( *get_dTrialStressdD( ) ),
                dRdD.begin( ),
                std::minus<>( )
            );

        }

    }

}
