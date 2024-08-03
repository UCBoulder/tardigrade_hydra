/**
  ******************************************************************************
  * \file tardigrade_hydraMicromorphicLinearElasticity.cpp
  ******************************************************************************
  * An implementation of micromorphic linear elasticity using the hydra
  * framework. Used as an example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraMicromorphicLinearElasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_micromorphic_tools.h>

namespace tardigradeHydra{

    namespace micromorphicLinearElasticity{

        errorOut linearElasticity( const variableVector &deformationGradient, const variableVector &microDeformation,
                                   const variableVector &gradientMicroDeformation,
                                   const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                   const parameterVector &D,
                                   variableVector &cauchyStress, variableVector &microStress,
                                   variableVector &higherOrderStress ){
            /*!
             * Compute the stress measures in the current configuration for micromorphic linear elasticity based off 
             * of a quadratic decomposition of the energy.
             *
             * :param const variableVector &deformationGradient: The deformation gradient
             * :param const variableVector &microDeformation: The micro-deformation
             * :param const variableVector &gradientMicroDeformation: The spatial gradient of the micro-deformation
             * :param const parameterVector &A: The A stiffness matrix.
             * :param const parameterVector &B: The B stiffness matrix.
             * :param const parameterVector &C: The C stiffness matrix.
             * :param const parameterVector &D: The D stiffness matrix.
             * :param variableVector &cauchyStress: The Cauchy stress.
             * :param variableVector &microStress: The symmetric micro-stress.
             * :param variableVector &higherOrderStress: The higher-order stress.
             */
    
            variableVector PK2Stress, referenceMicroStress, referenceHigherOrderStress;
            TARDIGRADE_ERROR_TOOLS_CATCH( linearElasticityReference( deformationGradient, microDeformation, gradientMicroDeformation,
                                                                     A, B, C, D,
                                                                     PK2Stress, referenceMicroStress, referenceHigherOrderStress ) );
    
            TARDIGRADE_ERROR_TOOLS_CATCH( mapStressMeasuresToCurrent( deformationGradient, microDeformation, PK2Stress,
                                                                      referenceMicroStress, referenceHigherOrderStress,
                                                                      cauchyStress, microStress, higherOrderStress ) );
    
            return NULL;
        }
    
        errorOut linearElasticity( const variableVector &deformationGradient, const variableVector &microDeformation,
                                   const variableVector &gradientMicroDeformation,
                                   const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                   const parameterVector &D,
                                   variableVector &cauchyStress, variableVector &microStress,
                                   variableVector &higherOrderStress,
                                   variableVector &dCauchyStressdF, variableVector &dCauchyStressdChi, variableVector &dCauchyStressdGradChi,
                                   variableVector &dMicroStressdF, variableVector &dMicroStressdChi, variableVector &dMicroStressdGradChi,
                                   variableVector &dHigherOrderStressdF, variableVector &dHigherOrderStressdChi,
                                   variableVector &dHigherOrderStressdGradChi ){
            /*!
             * Compute the stress measures in the current configuration for micromorphic linear elasticity based off 
             * of a quadratic decomposition of the energy.
             *
             * Also compute the Jacobians
             *
             * \param &deformationGradient: The deformation gradient
             * \param &microDeformation: The micro-deformation
             * \param &gradientMicroDeformation: The spatial gradient of the micro-deformation
             * \param &A: The A stiffness matrix.
             * \param &B: The B stiffness matrix.
             * \param &C: The C stiffness matrix.
             * \param &D: The D stiffness matrix.
             * \param &cauchyStress: The Cauchy stress.
             * \param &microStress: The symmetric micro-stress.
             * \param &higherOrderStress: The higher-order stress.
             * \param &dCauchyStressdF: The Jacobian of the Cauchy stress w.r.t. the deformation gradient
             * \param &dCauchyStressdChi: The Jacobian of the Cauchy stress w.r.t. the micro deformation.
             * \param &dCauchyStressdGradChi: The Jacobian of the Cauchy stress w.r.t. the gradient of the 
             * \   micro-deformation.
             * \param &dMicroStressdF: The Jacobian of the Micro stress w.r.t. the deformation gradient
             * \param &dMicroStressdChi: The Jacobian of the Micro stress w.r.t. the micro deformation.
             * \param &dMicroStressdGradChi: The Jacobian of the Micro stress w.r.t. the gradient of the 
             * \   micro-deformation.
             * \param &dHigherOrderStressdF: The Jacobian of the Higher Order stress w.r.t. the deformation gradient
             * \param &dHigherOrderStressdChi: The Jacobian of the Higher Order stress w.r.t. the micro deformation.
             * \param &dHigherOrderStressdGradChi: The Jacobian of the Higher Order stress w.r.t. the gradient of the 
             *     micro-deformation.
             */

            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;
 
            variableVector PK2Stress, referenceMicroStress, referenceHigherOrderStress;
    
            variableVector dPK2StressdF, dPK2StressdChi, dPK2StressdGradChi;
            variableVector dReferenceMicroStressdF, dReferenceMicroStressdChi, dReferenceMicroStressdGradChi;
            variableVector dReferenceHigherOrderStressdF, dReferenceHigherOrderStressdGradChi;
    
            TARDIGRADE_ERROR_TOOLS_CATCH( linearElasticityReference( deformationGradient, microDeformation, gradientMicroDeformation,
                                                                     A, B, C, D,
                                                                     PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                     dPK2StressdF, dPK2StressdChi, dPK2StressdGradChi,
                                                                     dReferenceMicroStressdF, dReferenceMicroStressdChi, dReferenceMicroStressdGradChi,
                                                                     dReferenceHigherOrderStressdF, dReferenceHigherOrderStressdGradChi ) );
    
            variableVector dCauchyStressdPK2Stress, dMicroStressdReferenceMicroStress, dHigherOrderStressdReferenceHigherOrderStress;
    
            TARDIGRADE_ERROR_TOOLS_CATCH( mapStressMeasuresToCurrent( deformationGradient, microDeformation, PK2Stress,
                                                                      referenceMicroStress, referenceHigherOrderStress,
                                                                      cauchyStress, microStress, higherOrderStress,
                                                                      dCauchyStressdF, dCauchyStressdPK2Stress,
                                                                      dMicroStressdF, dMicroStressdReferenceMicroStress,
                                                                      dHigherOrderStressdF, dHigherOrderStressdChi,
                                                                      dHigherOrderStressdReferenceHigherOrderStress ) );

            // Size the target vectors
            dCauchyStressdChi = floatVector( sot_dim * sot_dim, 0 );
            dCauchyStressdGradChi = floatVector( sot_dim * tot_dim, 0 );
            dMicroStressdChi = floatVector( sot_dim * sot_dim, 0 );
            dMicroStressdGradChi = floatVector( sot_dim * tot_dim, 0 );
            dHigherOrderStressdGradChi = floatVector( tot_dim * tot_dim, 0 );
 
            // Set up the Eigen Maps
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dCauchyStressdF_map(               dCauchyStressdF.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dCauchyStressdChi_map(             dCauchyStressdChi.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, tot_dim, Eigen::RowMajor > > dCauchyStressdGradChi_map(         dCauchyStressdGradChi.data( ), sot_dim, tot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dCauchyStressdPK2Stress_map(       dCauchyStressdPK2Stress.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dPK2StressdF_map(                  dPK2StressdF.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dPK2StressdChi_map(                dPK2StressdChi.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, tot_dim, Eigen::RowMajor > > dPK2StressdGradChi_map(            dPK2StressdGradChi.data( ), sot_dim, tot_dim );

            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dMicroStressdF_map(                    dMicroStressdF.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dMicroStressdChi_map(                  dMicroStressdChi.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, tot_dim, Eigen::RowMajor > > dMicroStressdGradChi_map(              dMicroStressdGradChi.data( ), sot_dim, tot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dMicroStressdReferenceMicroStress_map( dMicroStressdReferenceMicroStress.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dReferenceMicroStressdF_map(           dReferenceMicroStressdF.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dReferenceMicroStressdChi_map(         dReferenceMicroStressdChi.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, tot_dim, Eigen::RowMajor > > dReferenceMicroStressdGradChi_map(     dReferenceMicroStressdGradChi.data( ), sot_dim, tot_dim );

            Eigen::Map< Eigen::Matrix< variableType, tot_dim, sot_dim, Eigen::RowMajor > > dHigherOrderStressdF_map(                          dHigherOrderStressdF.data( ), tot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, tot_dim, tot_dim, Eigen::RowMajor > > dHigherOrderStressdGradChi_map(                    dHigherOrderStressdGradChi.data( ), tot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, tot_dim, tot_dim, Eigen::RowMajor > > dHigherOrderStressdReferenceHigherOrderStress_map( dHigherOrderStressdReferenceHigherOrderStress.data( ), tot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, tot_dim, sot_dim, Eigen::RowMajor > > dReferenceHigherOrderStressdF_map(                 dReferenceHigherOrderStressdF.data( ), tot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, tot_dim, tot_dim, Eigen::RowMajor > > dReferenceHigherOrderStressdGradChi_map(           dReferenceHigherOrderStressdGradChi.data( ), tot_dim, sot_dim );

            //Assemble the jacobians of the Cauchy stress
            dCauchyStressdF_map = ( dCauchyStressdF_map + dCauchyStressdPK2Stress_map * dPK2StressdF_map ).eval( );
            dCauchyStressdChi_map = ( dCauchyStressdPK2Stress_map * dPK2StressdChi_map ).eval( );
            dCauchyStressdGradChi_map = ( dCauchyStressdPK2Stress_map * dPK2StressdGradChi_map ).eval( );
    
            //Assemble the jacobians of the symmetric micro-stress
            dMicroStressdF_map = ( dMicroStressdF_map + dMicroStressdReferenceMicroStress_map * dReferenceMicroStressdF_map ).eval( );
            dMicroStressdChi_map = ( dMicroStressdReferenceMicroStress_map * dReferenceMicroStressdChi_map ).eval( );
            dMicroStressdGradChi_map = ( dMicroStressdReferenceMicroStress_map * dReferenceMicroStressdGradChi_map ).eval( );
    
            //Assemble the jacobians of the higher-order stress
            dHigherOrderStressdF_map       = ( dHigherOrderStressdF_map + dHigherOrderStressdReferenceHigherOrderStress_map * dReferenceHigherOrderStressdF_map ).eval( );
            dHigherOrderStressdGradChi_map = ( dHigherOrderStressdReferenceHigherOrderStress_map * dReferenceHigherOrderStressdGradChi_map ).eval( );
    
            return NULL;
        }

        errorOut linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                            const variableVector &gradientMicroDeformation,
                                            const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                            const parameterVector &D,
                                            variableVector &PK2Stress, variableVector &referenceMicroStress,
                                            variableVector &referenceHigherOrderStress ){
            /*!
             * Compute the stress measures in the reference configuration for micromorphic linear elasticity based off 
             * of a quadratic decomposition of the energy.
             *
             * \param &deformationGradient: The deformation gradient
             * \param &microDeformation: The micro-deformation
             * \param &gradientMicroDeformation: The spatial gradient of the micro-deformation
             * \param &A: The A stiffness matrix.
             * \param &B: The B stiffness matrix.
             * \param &C: The C stiffness matrix.
             * \param &D: The D stiffness matrix.
             * \param &PK2Stress: The second Piola-Kirchoff stress.
             * \param &referenceMicroStress: The symmetric micro-stress in the 
             * \   reference configuration.
             * \param &referenceHigherOrderStress: The higher-order stress in the 
             *     reference configuration.
             */
    
            //Compute the required deformation measures
            variableVector RCG, Psi, Gamma;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeDeformationMeasures( deformationGradient, microDeformation, gradientMicroDeformation,
                                                                      RCG, Psi, Gamma ) );
    
            TARDIGRADE_ERROR_TOOLS_CATCH( linearElasticityReferenceDerivedMeasures( RCG, Psi, Gamma, A, B, C, D,
                                                                                    PK2Stress, referenceMicroStress,
                                                                                    referenceHigherOrderStress ) );
    
            return NULL;
        }

        errorOut linearElasticityReference( const variableVector &deformationGradient, const variableVector &microDeformation,
                                            const variableVector &gradientMicroDeformation,
                                            const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                            const parameterVector &D,
                                            variableVector &PK2Stress, variableVector &referenceMicroStress,
                                            variableVector &referenceHigherOrderStress,
                                            variableVector &dPK2StressdF, variableVector &dPK2StressdChi, variableVector &dPK2StressdGradChi,
                                            variableVector &dReferenceMicroStressdF, variableVector &dReferenceMicroStressdChi,
                                            variableVector &dReferenceMicroStressdGradChi, variableVector &dMdF, variableVector &dMdGradChi ){
            /*!
             * Compute the stress measures in the reference configuration for micromorphic linear elasticity based off
             * of a quadratic decomposition of the energy.
             *
             * Also computes the Jacobians
             *
             * \param &deformationGradient: The deformation gradient
             * \param &microDeformation: The micro-deformation
             * \param &gradientMicroDeformation: The spatial gradient of the micro-deformation
             * \param &A: The A stiffness matrix.
             * \param &B: The B stiffness matrix.
             * \param &C: The C stiffness matrix.
             * \param &D: The D stiffness matrix.
             * \param &PK2Stress: The second Piola-Kirchoff stress.
             * \param &referenceMicroStress: The symmetric micro-stress in the 
             * \   reference configuration.
             * \param &referenceHigherOrderStress: The higher-order stress in the 
             * \   reference configuration.
             * \param &dPK2StressdF: The Jacobian of the PK2 stress w.r.t. the deformation gradient.
             * \param &dPK2StressdChi: The Jacobian of the PK2 stress w.r.t. the micro deformation.
             * \param &dPK2StressdGradChi: The Jacobian of the PK2 stress w.r.t. the gradient of the micro deformation.
             * \param &dReferenceMicroStressdF: The Jacobian of the Micro stress w.r.t. the deformation gradient.
             * \param &dReferenceMicroStressdChi: The Jacobian of the Micro stress w.r.t. the micro deformation.
             * \param &dReferenceMicroStressdGradChi: The Jacobian of the Micro stress w.r.t. the gradient of the micro deformation.
             * \param &dMdF: The Jacobian of the higher order stress w.r.t. the deformation gradient.
             * \param &dMdGradChi: The Jacobian of the higher order stress w.r.t. the gradient of the micro deformation.
             */
    
            //Assume 3d
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;
    
            //Compute the required deformation measures
            variableVector RCG, Psi, Gamma;
            variableVector dRCGdF, dPsidF, dPsidChi, dGammadF, dGammadGradChi;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeDeformationMeasures( deformationGradient, microDeformation, gradientMicroDeformation,
                                                                      RCG, Psi, Gamma, dRCGdF, dPsidF, dPsidChi, dGammadF, dGammadGradChi ) );
    
            variableVector dPK2StressdRCG, dPK2StressdPsi, dPK2StressdGamma;
            variableVector dReferenceMicroStressdRCG, dReferenceMicroStressdPsi, dReferenceMicroStressdGamma;
            variableVector dMdGamma;
    
            TARDIGRADE_ERROR_TOOLS_CATCH( linearElasticityReferenceDerivedMeasures( RCG, Psi, Gamma, A, B, C, D,
                                                                                    PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                                                    dPK2StressdRCG, dPK2StressdPsi, dPK2StressdGamma,
                                                                                    dReferenceMicroStressdRCG, dReferenceMicroStressdPsi,
                                                                                    dReferenceMicroStressdGamma, dMdGamma ) );
   
            // Size the arrays
            dPK2StressdF       = floatVector( sot_dim * sot_dim, 0 );
            dPK2StressdChi     = floatVector( sot_dim * sot_dim, 0 );
            dPK2StressdGradChi = floatVector( sot_dim * tot_dim, 0 );

            dReferenceMicroStressdF       = floatVector( sot_dim * sot_dim, 0 );
            dReferenceMicroStressdChi     = floatVector( sot_dim * sot_dim, 0 );
            dReferenceMicroStressdGradChi = floatVector( sot_dim * tot_dim, 0 );

            dMdF       = floatVector( tot_dim * sot_dim, 0 );
            dMdGradChi = floatVector( tot_dim * tot_dim, 0 );

            // Form the Eigen maps
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dPK2StressdF_map(       dPK2StressdF.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dPK2StressdChi_map(     dPK2StressdChi.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, tot_dim, Eigen::RowMajor > > dPK2StressdGradChi_map( dPK2StressdGradChi.data( ), sot_dim, tot_dim );

            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dPK2StressdRCG_map(   dPK2StressdRCG.data( ),   sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dPK2StressdPsi_map(   dPK2StressdPsi.data( ),   sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, tot_dim, Eigen::RowMajor > > dPK2StressdGamma_map( dPK2StressdGamma.data( ), sot_dim, tot_dim );

            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dReferenceMicroStressdF_map(       dReferenceMicroStressdF.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dReferenceMicroStressdChi_map(     dReferenceMicroStressdChi.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, tot_dim, Eigen::RowMajor > > dReferenceMicroStressdGradChi_map( dReferenceMicroStressdGradChi.data( ), sot_dim, tot_dim );

            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dReferenceMicroStressdRCG_map(   dReferenceMicroStressdRCG.data( ),   sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dReferenceMicroStressdPsi_map(   dReferenceMicroStressdPsi.data( ),   sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, tot_dim, Eigen::RowMajor > > dReferenceMicroStressdGamma_map( dReferenceMicroStressdGamma.data( ), sot_dim, tot_dim );

            Eigen::Map< Eigen::Matrix< variableType, tot_dim, sot_dim, Eigen::RowMajor > > dMdF_map(       dMdF.data( ), tot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, tot_dim, tot_dim, Eigen::RowMajor > > dMdGradChi_map( dMdGradChi.data( ), tot_dim, tot_dim );

            Eigen::Map< Eigen::Matrix< variableType, tot_dim, tot_dim, Eigen::RowMajor > > dMdGamma_map( dMdGamma.data( ), tot_dim, tot_dim );

            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dRCGdF_map(   dRCGdF.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dPsidF_map(   dPsidF.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, sot_dim, sot_dim, Eigen::RowMajor > > dPsidChi_map( dPsidChi.data( ), sot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, tot_dim, sot_dim, Eigen::RowMajor > > dGammadF_map( dGammadF.data( ), tot_dim, sot_dim );
            Eigen::Map< Eigen::Matrix< variableType, tot_dim, tot_dim, Eigen::RowMajor > > dGammadGradChi_map( dGammadGradChi.data( ), tot_dim, tot_dim );

            dPK2StressdF_map       = ( dPK2StressdRCG_map * dRCGdF_map + dPK2StressdPsi_map * dPsidF_map + dPK2StressdGamma_map * dGammadF_map ).eval( );
            dPK2StressdChi_map     = ( dPK2StressdPsi_map * dPsidChi_map ).eval( );
            dPK2StressdGradChi_map = ( dPK2StressdGamma_map * dGammadGradChi_map ).eval( );

            dReferenceMicroStressdF_map       = ( dReferenceMicroStressdRCG_map * dRCGdF_map + dReferenceMicroStressdPsi_map * dPsidF_map + dReferenceMicroStressdGamma_map * dGammadF_map ).eval( );
            dReferenceMicroStressdChi_map     = ( dReferenceMicroStressdPsi_map * dPsidChi_map ).eval( );
            dReferenceMicroStressdGradChi_map = ( dReferenceMicroStressdGamma_map * dGammadGradChi_map ).eval( );

            dMdF_map       = ( dMdGamma_map * dGammadF_map ).eval( );
            dMdGradChi_map = ( dMdGamma_map * dGammadGradChi_map ).eval( );
    
            return NULL;
        }

        errorOut linearElasticityReferenceDerivedMeasures( const variableVector &rightCauchyGreenDeformation,
                                                           const variableVector &Psi, const variableVector &Gamma,
                                                           const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                                           const parameterVector &D,
                                                           variableVector &PK2Stress, variableVector &referenceMicroStress,
                                                           variableVector &referenceHigherOrderStress ){
            /*!
             * Compute the stress measures in the reference configuration for micromorphic linear elasticity based off
             * of a quadratic decomposition of the energy.
             *
             * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation metric
             * \param &Psi: The micro-deformation measure Psi
             * \param &Gamma: The higher order deformation measure Gamma
             * \param &A: The A stiffness matrix.
             * \param &B: The B stiffness matrix.
             * \param &C: The C stiffness matrix.
             * \param &D: The D stiffness matrix.
             * \param &PK2Stress: The second Piola-Kirchoff stress.
             * \param &referenceMicroStress: The symmetric micro-stress in the 
             *     reference configuration.
             * \param &referenceHigherOrderStress: The higher-order stress in the 
             *     reference configuration.
             */
    
            //Assume 3D
            constexpr unsigned int dim = 3;
    
            variableVector invRCG = rightCauchyGreenDeformation;
            Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> > mat( invRCG.data(), 3, 3 );
            mat = mat.inverse( ).eval( );
    
            //Compute the strain measures
            variableVector greenLagrangeStrain = 0.5 * rightCauchyGreenDeformation;
            variableVector microStrain   = Psi;
            for ( unsigned int i = 0; i < dim; i++ ){
                greenLagrangeStrain[ dim * i + i ] -= 0.5;
                microStrain[dim * i + i ] -= 1;
            }
 
            //Compute the higher order stress
            TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStress( Gamma, C, referenceHigherOrderStress ) );
    
            //Compute the first common term for the PK2 and symmetric micro-stress
            variableVector term1;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm1( greenLagrangeStrain, microStrain, A, D, term1 ) );
    
            //Compute the second common term for the PK2 and symmetric micro-stress
            variableVector invRCGPsi;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeInvRCGPsi( invRCG, Psi, invRCGPsi ) );
    
            variableVector term2;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm2( greenLagrangeStrain, microStrain, invRCGPsi, B, D, term2 ) );
    
            //Compute the third common term for the PK2 and symmetric micro-stress
            variableVector invRCGGamma;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeInvRCGGamma( invRCG, Gamma, invRCGGamma ) );
    
            variableVector term3;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm3( invRCGGamma, referenceHigherOrderStress, term3 ) );
    
            //Construct the PK2 and reference symmetric stresses
            PK2Stress            = term1 + term2 + term3;
    
            variableVector symmTerm2Term3;
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::computeSymmetricPart( term2 + term3, symmTerm2Term3 ) );
            referenceMicroStress = term1 + 2 * symmTerm2Term3;
    
            return NULL;
        }

        errorOut linearElasticityReferenceDerivedMeasures( const variableVector &rightCauchyGreenDeformation, const variableVector &Psi,
                                                           const variableVector &Gamma,
                                                           const parameterVector &A, const parameterVector &B, const parameterVector &C,
                                                           const parameterVector &D,
                                                           variableVector &PK2Stress, variableVector &referenceMicroStress,
                                                           variableVector &referenceHigherOrderStress,
                                                           variableVector &dPK2StressdRCG, variableVector &dPK2StressdPsi,
                                                           variableVector &dPK2StressdGamma,
                                                           variableVector &dReferenceMicroStressdRCG,
                                                           variableVector &dReferenceMicroStressdPsi,
                                                           variableVector &dReferenceMicroStressdGamma,
                                                           variableVector &dMdGamma ){
            /*!
             * Compute the stress measures in the reference configuration for micromorphic linear elasticity based off
             * of a quadratic decomposition of the energy.
             *
             * \param &rightCauchyGreenDeformation: The right Cauchy-Green deformation metric
             * \param &Psi: The micro-deformation measure Psi
             * \param &Gamma: The higher order deformation measure Gamma
             * \param &A: The A stiffness matrix.
             * \param &B: The B stiffness matrix.
             * \param &C: The C stiffness matrix.
             * \param &D: The D stiffness matrix.
             * \param &PK2Stress: The second Piola-Kirchoff stress.
             * \param &referenceMicroStress: The symmetric micro-stress in the 
             *     reference configuration.
             * \param &referenceHigherOrderStress: The higher-order stress in the 
             *     reference configuration.
             * \param &dPK2StressdRCG: The Jacobian of the PK2 stress w.r.t. the right Cauchy-Green
             *     deformation metric.
             * \param &dPK2StressdPsi: The Jacobian of the PK2 stress w.r.t. the micro deformation 
             *     metric.
             * \param &dPK2StressdGamma: The Jacobian of the PK2 stress w.r.t. the higher order 
             *     deformation measure.
             * \param &dReferenceMicroStressdRCG: The Jacobian of the reference micro stress w.r.t. the 
             *     right Cacuhy-Green deformation metric.
             * \param &dReferenceMicroStressdPsi: The Jacobian of the reference micro stress w.r.t. the 
             *     micro deformation measure.
             * \param &dReferenceMicroStressdGamma: The Jacobian of the reference micro stress w.r.t. the 
             *     higher order deformation measure.
             * \param &dMdGamma: The Jacobian of the reference higher order stress w.r.t. 
             *     the higher order deformation measure.
             */

            //Assume 3d
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;   
 
            variableVector invRCG = rightCauchyGreenDeformation;
            Eigen::Map < Eigen::Matrix< floatType, 3, 3, Eigen::RowMajor> > mat( invRCG.data(), 3, 3 );
            mat = mat.inverse( ).eval( );
    
            //Compute the strain measures
            variableVector greenLagrangeStrain = 0.5 * rightCauchyGreenDeformation;
            variableVector microStrain   = Psi;
            for ( unsigned int i = 0; i < dim; i++ ){
                greenLagrangeStrain[ dim * i + i ] -= 0.5;
                microStrain[dim * i + i ] -= 1;
            }
    
            //Compute the higher order stress
            TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStress( Gamma, C, referenceHigherOrderStress, dMdGamma ) );
    
            //Compute the first common term for the PK2 and symmetric micro-stress
            variableVector term1;
    
            variableVector dTerm1dRCG, dTerm1dPsi;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm1( greenLagrangeStrain, microStrain, A, D, term1,
                                                                     dTerm1dRCG, dTerm1dPsi ) );
    
            //Assemble term1 jacobians w.r.t. F and Chi
            dTerm1dRCG *= 0.5;
    
            //Compute the second common term for the PK2 and symmetric micro-stress
            variableVector invRCGPsi;
            variableVector dInvRCGPsidRCG, dInvRCGPsidPsi;
    
            TARDIGRADE_ERROR_TOOLS_CATCH( computeInvRCGPsi( invRCG, Psi, invRCGPsi, dInvRCGPsidRCG, dInvRCGPsidPsi ) );
    
            variableVector term2;
            variableVector dTerm2dRCG, dTerm2dPsi, dTerm2dInvRCGPsi;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm2( greenLagrangeStrain, microStrain, invRCGPsi, B, D, term2,
                                                                     dTerm2dRCG, dTerm2dPsi, dTerm2dInvRCGPsi ) );

            dTerm2dRCG *= 0.5;
            dTerm2dRCG += tardigradeVectorTools::matrixMultiply( dTerm2dInvRCGPsi, dInvRCGPsidRCG, sot_dim, sot_dim, sot_dim, sot_dim );
    
            dTerm2dPsi += tardigradeVectorTools::matrixMultiply( dTerm2dInvRCGPsi, dInvRCGPsidPsi, sot_dim, sot_dim, sot_dim, sot_dim );
    
            //Compute the third common term for the PK2 and symmetric micro-stress
            variableVector invRCGGamma;
            variableVector dInvRCGGammadRCG, dInvRCGGammadGamma;
    
            TARDIGRADE_ERROR_TOOLS_CATCH( computeInvRCGGamma( invRCG, Gamma, invRCGGamma, dInvRCGGammadRCG, dInvRCGGammadGamma ) );
    
            variableVector term3;
            variableVector dTerm3dInvRCGGamma, dTerm3dM;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm3( invRCGGamma, referenceHigherOrderStress, term3, dTerm3dInvRCGGamma, dTerm3dM ) );
    
            variableVector dTerm3dRCG = tardigradeVectorTools::matrixMultiply( dTerm3dInvRCGGamma, dInvRCGGammadRCG, sot_dim, tot_dim, tot_dim, sot_dim );
            variableVector dTerm3dGamma = tardigradeVectorTools::matrixMultiply( dTerm3dInvRCGGamma, dInvRCGGammadGamma, sot_dim, tot_dim, tot_dim, tot_dim )
                                        + tardigradeVectorTools::matrixMultiply( dTerm3dM, dMdGamma, sot_dim, tot_dim, tot_dim, tot_dim );
    
            //Construct the PK2 and reference symmetric stresses
            PK2Stress            = term1 + term2 + term3;
    
            dPK2StressdRCG    = dTerm1dRCG + dTerm2dRCG + dTerm3dRCG;
            dPK2StressdPsi    = dTerm1dPsi + dTerm2dPsi;
            dPK2StressdGamma  = dTerm3dGamma;

            variableVector symmTerm2Term3;
            variableVector dSymmTerm2Term3dTerm2Term3;

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::computeSymmetricPart( term2 + term3, symmTerm2Term3, dSymmTerm2Term3dTerm2Term3 ) );
            referenceMicroStress = term1 + 2 * symmTerm2Term3;

            dReferenceMicroStressdRCG = dTerm1dRCG + 2 * ( tardigradeVectorTools::matrixMultiply( dSymmTerm2Term3dTerm2Term3, dTerm2dRCG, sot_dim, sot_dim, sot_dim, sot_dim )
                                                         + tardigradeVectorTools::matrixMultiply( dSymmTerm2Term3dTerm2Term3, dTerm3dRCG, sot_dim, sot_dim, sot_dim, sot_dim ) );
    
            dReferenceMicroStressdPsi = dTerm1dPsi + 2 * tardigradeVectorTools::matrixMultiply( dSymmTerm2Term3dTerm2Term3, dTerm2dPsi, sot_dim, sot_dim, sot_dim, sot_dim );
            dReferenceMicroStressdGamma = 2 * tardigradeVectorTools::matrixMultiply( dSymmTerm2Term3dTerm2Term3, dTerm3dGamma, sot_dim, sot_dim, sot_dim, tot_dim );

            return NULL;
        }

        errorOut mapStressMeasuresToCurrent( const variableVector &deformationGradient, const variableVector &microDeformation,
                                             const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                             const variableVector &referenceHigherOrderStress,
                                             variableVector &cauchyStress, variableVector &microStress,
                                             variableVector &higherOrderStress ){
            /*!
             * Map the stress measures in the reference configuration to the current configuration.
             *
             * \param &deformationGradient: The deformation gradient between the 
             *     reference configuration and the current configuration.
             * \param &microDeformation: The micro-deformation map between the 
             *     reference configuration and the current configuration.
             * \param &PK2Stress: The Second Piola-Kirchoff stress.
             * \param &referenceMicroStress: The symmetric micro-stress in the 
             *     reference configuration.
             * \param &referenceHigherOrderStress: The higher order stress in 
             *     the reference configuration.
             * \param &cauchyStress: The Cauchy stress (PK2 stress in the current configuration).
             * \param &microStress: The symmetric micro-stress in the current configuration.
             * \param &higherOrderStress: The higher order stress in the current configuration.
             */
    
            //Map the PK2 stress to the Cauchy stress
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient, cauchyStress ) );
    
            //Map the symmetric micro stress to the current configuration
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, microStress ) );
    
            //Map the higher order stress to the current configuration
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient,
                                                                                                     microDeformation, higherOrderStress ) );
    
            return NULL;
        }
    
        errorOut mapStressMeasuresToCurrent( const variableVector &deformationGradient, const variableVector &microDeformation,
                                             const variableVector &PK2Stress, const variableVector &referenceMicroStress,
                                             const variableVector &referenceHigherOrderStress,
                                             variableVector &cauchyStress, variableVector &microStress,
                                             variableVector &higherOrderStress,
                                             variableVector &dCauchyStressdF, variableVector &dCauchyStressdPK2Stress,
                                             variableVector &dMicroStressdF, variableVector &dMicroStressdReferenceMicroStress,
                                             variableVector &dHigherOrderStressdF, variableVector &dHigherOrderStressdChi,
                                             variableVector &dHigherOrderStressdReferenceHigherOrderStress ){
            /*!
             * Map the stress measures in the reference configuration to the current configuration.
             *
             * Also computes the Jacobians
             *
             * \param &deformationGradient: The deformation gradient between the 
             * \   reference configuration and the current configuration.
             * \param &microDeformation: The micro-deformation map between the 
             * \   reference configuration and the current configuration.
             * \param &PK2Stress: The Second Piola-Kirchoff stress.
             * \param &referenceMicroStress: The symmetric micro-stress in the 
             * \   reference configuration.
             * \param &referenceHigherOrderStress: The higher order stress in 
             * \   the reference configuration.
             * \param &cauchyStress: The Cauchy stress (PK2 stress in the current configuration).
             * \param &microStress: The symmetric micro-stress in the current configuration.
             * \param &higherOrderStress: The higher order stress in the current configuration.
             * \param &dCauchyStressdF: The Jacobian of the Cauchy stress w.r.t. the 
             * \   deformation gradient.
             * \param &dCauchyStressdPK2Stress: The Jacobian of the Cauchy stress w.r.t. the 
             * \   PK2 stress.
             * \param &dMicroStressdF: The Jacobian of the micro stress w.r.t. the 
             * \   deformation gradient.
             * \param &dMicroStressdReferenceMicroStress: The Jacobian of the micro-stress 
             * \   in the current configuration w.r.t. the micro-stress in the reference configuration.
             * \param &dHigherOrderStressdF: The Jacobian of the higher-order stress w.r.t.
             * \   the deformation gradient.
             * \param &dHigherOrderStressdChi: The Jacobian of the higher-order stress 
             * \   w.r.t. the micro-deformation.
             * \param &dHigherOrderStressdReferenceHigherOrderStress: The Jacobian of the 
             *     higher-order stress w.r.t. the higher order stress in the reference configuration.
             */
    
            //Map the PK2 stress to the Cauchy stress
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient, cauchyStress,
                                                                                             dCauchyStressdPK2Stress, dCauchyStressdF ) );
    
            //Map the symmetric micro stress to the current configuration
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, microStress,
                                                                                                        dMicroStressdReferenceMicroStress, dMicroStressdF ) );
    
            //Map the higher order stress to the current configuration
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient,
                                                                                                     microDeformation, higherOrderStress,
                                                                                                     dHigherOrderStressdReferenceHigherOrderStress,
                                                                                                     dHigherOrderStressdF,
                                                                                                     dHigherOrderStressdChi ) );
    
            return NULL;
        }

        errorOut computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                             const variableVector &gradientMicroDeformation,
                                             variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma ){
            /*!
             * Compute the deformation measures
             * \f$\begin{align}
             * C_{IJ} &= F_{iI} F_{iJ}\\
             * Psi_{IJ} &= F_{iI} \Chi_{iJ}\\
             * \Gamma_{IJK} &= F_{iI} \Chi_{iJ, K}
             * \end{align}\f$
             *
             * \param &deformationGradient: The deformation gradient
             * \param &microDeformation: The micro-deformation
             * \param &gradientMicroDeformation: The spatial gradient of the micro-deformation
             * \param &rightCauchyGreen: The Right Cauchy-Green deformation tensor
             * \param &Psi: The micro-deformation measure
             * \param &Gamma: The gradient micro-deformation measure
             */
    
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::computeRightCauchyGreen( deformationGradient, rightCauchyGreen ) );
    
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::computePsi( deformationGradient, microDeformation, Psi ) );
    
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::computeGamma( deformationGradient, gradientMicroDeformation, Gamma ) );
    
            return NULL;
    
        }
    
        errorOut computeDeformationMeasures( const variableVector &deformationGradient, const variableVector &microDeformation,
                                             const variableVector &gradientMicroDeformation,
                                             variableVector &rightCauchyGreen, variableVector &Psi, variableVector &Gamma,
                                             variableVector &dCdF, variableVector &dPsidF, variableVector &dPsidChi,
                                             variableVector &dGammadF, variableVector &dGammadGradChi ){
            /*!
             * Compute the deformation measures
             * \f$\begin{align}
             * C_{IJ} &= F_{iI} F_{iJ}\\
             * Psi_{IJ} &= F_{iI} \Chi_{iJ}\\
             * \Gamma_{IJK} &= F_{iI} \Chi_{iJ, K}
             * \end{align}\f$
             *
             * \param &deformationGradient: The deformation gradient
             * \param &microDeformation: The micro-deformation
             * \param &gradientMicroDeformation: The spatial gradient of the micro-deformation
             * \param &rightCauchyGreen: The Right Cauchy-Green deformation tensor
             * \param &Psi: The micro-deformation measure
             * \param &Gamma: The gradient micro-deformation measure
             * \param &dCdF: The gradient of the right Cauchy green deformation tensor w.r.t. 
             *    the deformation gradient.
             * \param &dPsidF: The gradient of Psi w.r.t. the deformation gradient.
             * \param &dPsidChi: The gradient of Psi w.r.t. the microDeformation.
             * \param &dGammadF: The gradient of Gamma w.r.t. the deformation gradient.
             * \param &dGammadGradChi: The gradient of Gamma w.r.t. the spatial gradient of Chi
             */
    
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeConstitutiveTools::computeRightCauchyGreen( deformationGradient, rightCauchyGreen, dCdF ) );
    
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::computePsi( deformationGradient, microDeformation, Psi, dPsidF, dPsidChi ) );
    
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::computeGamma( deformationGradient, gradientMicroDeformation, Gamma, dGammadF, dGammadGradChi ) );
    
            return NULL;
        }

        errorOut computeLinearElasticTerm1( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                            const parameterVector &A, const parameterVector &D, variableVector &term1 ){
            /*!
             * Compute the first term for the linear elastic model
             * \f$ term1_{IJ} = A_{IJKL} E_{KL} + D_{IJKL} * \mathcal{E}_{KL} \f$
             *
             * \param &greenLagrangeStrain: The Green-Lagrange strain.
             * \param &microStrain: The micro-strain
             * \param &A: The A stiffness matrix
             * \param &D: The D stiffness matrix
             * \param &term1: The first term.
             */
    
            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;
            constexpr unsigned int fot_dim = tot_dim * dim;
    
            TARDIGRADE_ERROR_TOOLS_CHECK( greenLagrangeStrain.size() == sot_dim, "The green lagrange strain must have a length of 9" );
    
            TARDIGRADE_ERROR_TOOLS_CHECK( microStrain.size() == sot_dim, "The micro-strain must have a length of 9" );
    
            TARDIGRADE_ERROR_TOOLS_CHECK( A.size() == fot_dim, "A must have a size of 3**4" );
    
            TARDIGRADE_ERROR_TOOLS_CHECK( D.size() == fot_dim, "D must have a size of 3**4" );
    
            //Compute the first common term for the PK2 and symmetric micro-stress
            term1 = variableVector( dim * dim, 0 );
            for ( unsigned int IJ = 0; IJ < sot_dim; IJ++ ){
                for ( unsigned int KL = 0; KL < sot_dim; KL++ ){
                    term1[ IJ ] += A[ sot_dim * IJ + KL ] * greenLagrangeStrain[ KL ]
                                 + D[ sot_dim * IJ + KL ] * microStrain[ KL ];
                }
            }
    
            return NULL;
        }
    
        errorOut computeLinearElasticTerm1( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                            const parameterVector &A, const parameterVector &D, variableVector &term1,
                                            variableVector &dTerm1dGreenLagrangeStrain, variableVector &dTerm1dMicroStrain ){
            /*!
             * Compute the first term for the linear elastic model
             * \f$ term1_{IJ} = A_{IJKL} E_{KL} + D_{IJKL} * \mathcal{E}_{KL} \f$
             *
             * Also return the Jacobian
             * \f$\begin{align}
             * \frac{\partial term^1_{IJ} }{ E_{MN} } &= A_{IJMN}\\
             * \frac{\partial term^1_{IJ} }{ \mathcal{E}_{MN} } &= D_{IJMN}
             * \end{align}\f$
             *
             * \param &greenLagrangeStrain: The Green-Lagrange strain.
             * \param &microStrain: The micro-strain
             * \param &A: The A stiffness matrix
             * \param &D: The D stiffness matrix
             * \param &term1: The first term.
             * \param &dTerm1dGreenLagrangeStrain: The derivative of the first term w.r.t. the Green-Lagrange strain
             * \param &dTerm1dMicroStrain: The derivative of the first term w.r.t. the micro strain
             */
    
            TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm1( greenLagrangeStrain, microStrain, A, D, term1 ) );
    
            //Compute the first common term for the PK2 and symmetric micro-stress
            dTerm1dGreenLagrangeStrain = A;
            dTerm1dMicroStrain = D;
    
            return NULL;
        }
    
        errorOut computeLinearElasticTerm2( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                            const variableVector &invCPsi, const parameterVector &B, const parameterVector &D,
                                            variableVector &term2 ){
            /*!
             * Compute the second term from the linear elastic constitutive model
             *
             * \f$term^2_{IJ} = \left( B_{IQKL} \mathcal{E}_{KL} + E_{KL} D_{KLIQ} \right) C_{JR}^{-1} \Psi_{RQ}\f$
             *
             * \param &greenLagrangeStrain: The Green-Lagrange strain \f$E_{IJ} = \frac{1}{2} \left( C_{IJ} - \delta_{IJ} \right)\f$
             * \param &microStrain: The micro-strain \f$\mathcal{E}_{IJ} = \Psi_{IJ} - \delta_{IJ}\f$
             * \param &invCPsi: The product \f$C_{JR}^{-1} \Psi_{RQ}\f$
             * \param &B: The B stiffness matrix
             * \param &D: The D stiffness matrix
             * \param &term2: The second term.
             */
    
            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
    
            TARDIGRADE_ERROR_TOOLS_CHECK( greenLagrangeStrain.size() == sot_dim, "The green lagrange strain must have a length of 9" );
    
            TARDIGRADE_ERROR_TOOLS_CHECK( microStrain.size() == sot_dim,  "The micro-strain must have a length of 9" );
    
            TARDIGRADE_ERROR_TOOLS_CHECK( invCPsi.size() == sot_dim, "invCPsi must have a size of 9" );
    
            TARDIGRADE_ERROR_TOOLS_CHECK( B.size() == sot_dim * sot_dim, "B must have a size of 3**4" );
    
            TARDIGRADE_ERROR_TOOLS_CHECK( D.size() == sot_dim * sot_dim, "D must have a size of 3**4" );
    
            term2 = variableVector( sot_dim, 0 );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int L = 0; L < dim; L++ ){
                            for ( unsigned int Q = 0; Q < dim; Q++ ){
                                term2[ dim * I + J] += ( B[ dim * dim * dim * I + dim * dim * Q + dim * K + L ] * microStrain[ dim * K + L ]
                                                     + greenLagrangeStrain[ dim * K + L ] * D[ dim * dim * dim * K + dim * dim * L + dim * I + Q ] )
                                                     * invCPsi[ dim * J + Q ];
                            }
                        }
                    }
                }
            }
    
            return NULL;
        }
    
       errorOut computeLinearElasticTerm2( const variableVector &greenLagrangeStrain, const variableVector &microStrain,
                                            const variableVector &invCPsi, const parameterVector &B, const parameterVector &D,
                                            variableVector &term2, variableVector &dTerm2dGreenLagrangeStrain,
                                            variableVector &dTerm2dMicroStrain, variableVector &dTerm2dInvCPsi ){
            /*!
             * Compute the second term from the linear elastic constitutive model
             *
             * \f$term^2_{IJ} = \left( B_{IQKL} \mathcal{E}_{KL} + E_{KL} D_{KLIQ} \right) C_{JR}^{-1} \Psi_{RQ}\f$
             *
             * Also return the Jacobians
             * \f$\begin{align}
             * \frac{ \partial term^2_{IJ} }{ \partial E_{MN} } &= D_{MNIK} C_{JR}^{-1} \Psi_{RK}\\
             * \frac{ \partial term^2_{IJ} }{ \partial \mathcal{E}_{MN} } &= B_{IKMN} C_{JR}^{-1} \Psi_{RK}\\
             * \frac{ \partial term^2_{IJ} }{ \partial C_{MO}^{-1} \Psi_{ON} } &= \left( B_{INKL} \mathcal{E}_{KL} + E_{KL} D_{KLIN} \right) \delta_{JM}
             * \end{align}\f$
             *
             * \param &greenLagrangeStrain: The Green-Lagrange strain \f$E_{IJ} = \frac{1}{2} \left( C_{IJ} - \delta_{IJ} \right)\f$
             * \param &microStrain: The micro-strain \f$\mathcal{E}_{IJ} = \Psi_{IJ} - \delta_{IJ}\f$
             * \param &invCPsi: The product \f$C_{JR}^{-1} \Psi_{RQ}\f$
             * \param &B: The B stiffness matrix
             * \param &D: The D stiffness matrix
             * \param &term2: The second term.
             * \param &dTerm2dGreenLagrangeStrain: The jacobian of term 2 w.r.t. the Green-Lagrange strain.
             * \param &dTerm2dMicroStrain: The jacobian of term 2 w.r.t. the microStrain.
             * \param &dTerm2dInvCPsi: The jacobian of term 2 w.r.t. \f$C_{JR}^{-1} \Psi_{RQ}\f$
             */
    
            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
    
            TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm2( greenLagrangeStrain, microStrain, invCPsi, B, D, term2 ) );
    
            //Compute the Jacobians
            dTerm2dGreenLagrangeStrain = variableVector( sot_dim * sot_dim, 0 );
            dTerm2dMicroStrain         = variableVector( sot_dim * sot_dim, 0 );
            dTerm2dInvCPsi             = variableVector( sot_dim * sot_dim, 0 );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int M = 0; M < dim; M++ ){
                        for ( unsigned int N = 0; N < dim; N++ ){
                            for ( unsigned int K = 0; K < dim; K++ ){
                                dTerm2dGreenLagrangeStrain[ dim * sot_dim * I + sot_dim * J + dim * M + N ] += D[ dim * dim * dim * M + dim * dim * N + dim * I + K] * invCPsi[ dim * J + K ];
                                dTerm2dMicroStrain[ dim * sot_dim * I + sot_dim * J + dim * M + N ] += B[ dim * dim * dim * I + dim * dim * K + dim * M + N] * invCPsi[ dim * J + K ];
                                dTerm2dInvCPsi[ dim * sot_dim * I + sot_dim * J + dim * J + M ] += ( B[ dim * dim * dim * I + dim * dim * M + dim * N + K ] * microStrain[ dim * N + K ] + greenLagrangeStrain[ dim * N + K ] * D[ dim * dim * dim * N + dim * dim * K + dim * I + M ] );
                            }
                        }
                    }
                }
            }
    
            return NULL;
        }

        errorOut computeReferenceHigherOrderStress( const variableVector &Gamma, const parameterVector &C,
                                                    variableVector &referenceHigherOrderStress ){
            /*!
             * Compute the higher order stress in the reference configuration.
             * \f$M_{IJK} = C_{JKILMN} Gamma_{LMN}\f$
             *
             * \param &Gamma: The micro-gradient deformation measure.
             * \param &C: The C stiffness tensor.
             * \param &referenceHigherOrderStress: The higher order stress in the reference 
             *     configuration.
             */
    
            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;
    
            TARDIGRADE_ERROR_TOOLS_CHECK( Gamma.size() == tot_dim, "Gamma must have a length of 27" );
    
            TARDIGRADE_ERROR_TOOLS_CHECK( C.size() == tot_dim * tot_dim, "The C stiffness tensor must have a length of 3**6.\nThe current size is " + std::to_string( C.size( ) ) );
    
            referenceHigherOrderStress = variableVector( tot_dim, 0 );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int LMN = 0; LMN < tot_dim; LMN++ ){
                                    referenceHigherOrderStress[ dim * dim * I + dim * J + K ] += C[ dim * dim * dim * dim * dim * J + dim * dim * dim * dim * K + dim * dim * dim * I + LMN ] * Gamma[ LMN ];
                        }
                    }
                }
            }
            return NULL;
        }
    
        errorOut computeReferenceHigherOrderStress( const variableVector &Gamma, const parameterVector &C,
                                                    variableVector &referenceHigherOrderStress,
                                                    variableVector &dReferenceHigherOrderStressdGamma ){
            /*!
             * Compute the higher order stress in the reference configuration.
             * \f$M_{IJK} = C_{JKILMN} \Gamma_{LMN}\f$
             *
             * Also compute the Jacobian
             * \f$\frac{ \partial M_{IJK} }{\partial \Gamma_{OPQ} } = C_{JKIOPQ}\f$
             *
             * \param &Gamma: The micro-gradient deformation measure.
             * \param &C: The C stiffness tensor.
             * \param &referenceHigherOrderStress: The higher order stress in the reference 
             *     configuration.
             * \param &dReferenceHigherOrderStressdGamma: The derivative of the higher order stress in the reference configuration w.r.t. \f$ \bf{\Gamma} \f$
             */
    
            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;
    
            TARDIGRADE_ERROR_TOOLS_CATCH( computeReferenceHigherOrderStress( Gamma, C, referenceHigherOrderStress ) );
    
            //Assemble the Jacobian
            dReferenceHigherOrderStressdGamma = variableVector( tot_dim * tot_dim, 0 );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int OPQ = 0; OPQ < tot_dim; OPQ++ ){
                            dReferenceHigherOrderStressdGamma[ dim * dim * tot_dim * I + dim * tot_dim * J + tot_dim * K + OPQ ] +=
                                C[ dim * dim * dim * dim * dim * J + dim * dim * dim * dim * K + dim * dim * dim * I + OPQ ];
                        }
                    }
                }
            }
    
            return NULL;
        }

        errorOut computeLinearElasticTerm3( const variableVector &invCGamma,
                                            const variableVector &referenceHigherOrderStress, variableVector &term3 ){
            /*!
             * Compute the value of the third term in the micromorphic linear elasticity formulation.
             * \f$ term3_{IJ} = M_{IQR} C_{JS}^{-1} \Gamma_{SQR} \f$
             *
             * \param &invCGamma: \f$ C_{JS}^{-1} \Gamma_{SQR} \f$
             * \param &referenceHigherOrderStress: The higher order stress in the reference configuration.
             * \param &term3: The third term in the linear elastic equation.
             */
    
            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;
    
            TARDIGRADE_ERROR_TOOLS_CHECK( invCGamma.size() == tot_dim, "invCGamma must have a size of 27" );
    
            TARDIGRADE_ERROR_TOOLS_CHECK( referenceHigherOrderStress.size() == tot_dim, "The referenceHigherOrder stress must have a size of 27" );
    
            term3 = variableVector( sot_dim, 0 );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int QR = 0; QR < sot_dim; QR++ ){
                        term3[ dim * I + J ] += referenceHigherOrderStress[ dim * dim * I + QR ] * invCGamma[ dim * dim * J + QR ];
                    }
                }
            }
    
            return NULL;
        }

        errorOut computeLinearElasticTerm3( const variableVector &invCGamma,
                                            const variableVector &referenceHigherOrderStress, variableVector &term3,
                                            variableVector &dTerm3dInvCGamma, variableVector &dTerm3dReferenceHigherOrderStress ){
            /*!
             * Compute the value of the third term in the micromorphic linear elasticity formulation.
             * \f$ term3_{IJ} = M_{IQR} C_{JS}^{-1} \Gamma_{SQR} \f$
             *
             * Also returns the Jacobians
             * \f$\begin{align}
             * \frac{ \partial term3_{IJ} }{ \partial M_{TUV} } &= \delta_{IT} C_{JS}^{-1} \Gamma_{SUV}\\
             * \frac{ \partial term3_{IJ} }{ \partial C_{TW}^{-1} \Gamma_{WUV} &= M_{IUV} \delta_{JT}
             * \end{align}\f$
             *
             * \param &invCGamma: \f$ C_{JS}^{-1} \Gamma_{SQR} \f$
             * \param &referenceHigherOrderStress: The higher order stress in the reference configuration.
             * \param &term3: The third term in the linear elastic equation.
             * \param &dTerm3dInvCGamma: The derivative of the third term w.r.t. \f$ C_{TW}^{-1} \Gamma_{WUV} \f$
             * \param &dTerm3dReferenceHigherOrderStress: The derivative of the third term w.r.t. the reference higher order stress
             */
    
            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;
    
            TARDIGRADE_ERROR_TOOLS_CATCH( computeLinearElasticTerm3( invCGamma, referenceHigherOrderStress, term3 ) );
    
            dTerm3dInvCGamma = variableVector( sot_dim * tot_dim, 0 );
            dTerm3dReferenceHigherOrderStress = variableVector( sot_dim * tot_dim, 0 );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int TU = 0; TU < sot_dim; TU++ ){
                            dTerm3dReferenceHigherOrderStress[ dim * tot_dim * I + tot_dim * J + dim * dim * I + TU ] += invCGamma[ dim * dim * J + TU ];
                            dTerm3dInvCGamma[ dim * tot_dim * I + tot_dim * J + dim * dim * J + TU ] += referenceHigherOrderStress[ dim * dim * I + TU ];
                    }
                }
            }
            return NULL;
        }

        errorOut computeInvRCGPsi( const variableVector &invRCG, const variableVector &Psi, variableVector &invRCGPsi ){
            /*!
             * Compute the product \f$ C_{IK}^{-1} \Psi_{KJ} \f$
             *
             * \param &invRCG: The inverse of the right cauchy green deformation tensor.
             * \param &Psi: The micro-deformation measure.
             * \param &invRCGPsi: the product.
             */
    
            //Assume 3d
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
    
            TARDIGRADE_ERROR_TOOLS_CHECK( invRCG.size() == sot_dim, "invRCG has an improper dimension" );
    
            TARDIGRADE_ERROR_TOOLS_CHECK( Psi.size() == sot_dim, "Psi has an improper dimension" );
    
            invRCGPsi = tardigradeVectorTools::matrixMultiply( invRCG, Psi, dim, dim, dim, dim );
    
            return NULL;
        }
    
        errorOut computeInvRCGPsi( const variableVector &invRCG, const variableVector &Psi, variableVector &invRCGPsi,
                                   variableVector &dInvRCGPsidRCG, variableVector &dInvRCGPsidPsi ){
            /*!
             * Compute the product \f$ C_{IK}^{-1} \Psi_{KJ} \f$
             *
             * Also compute the Jacobians
             * \f$\begin{align}
             * \frac{ \partial C_{IO}^{-1} \Psi_{OJ} } { \partial C_{KL} } &= -C_{IK}^{-1} C_{LO}^{-1} \Psi_{OJ}\\
             * \frac{ \partial C_{IO}^{-1} \Psi_{OJ} } { \partial \Psi_{KL} } &= C_{IK}^{-1} \delta_{JL}
             * \end{align}\f$
             *
             * \param &invRCG: The inverse of the right cauchy green deformation tensor.
             * \param &Psi: The micro-deformation measure.
             * \param &invRCGPsi: the product.
             * \param &dInvRCGPsidRCG: The Jacobian of the product w.r.t. the right cauchy green
             *     deformation tensor.
             * \param &dInvRCGPsidPsi: The Jacobian of the product w.r.t. the micro-deformation measure.
             */
    
            //Assume 3D
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
    
            TARDIGRADE_ERROR_TOOLS_CATCH( computeInvRCGPsi( invRCG, Psi, invRCGPsi ) );
    
            //Construct the jacobians
            dInvRCGPsidRCG = variableVector( sot_dim * sot_dim, 0 );
            dInvRCGPsidPsi = variableVector( sot_dim * sot_dim, 0 );
    
            for ( unsigned int I = 0; I < 3; I++ ){
                for ( unsigned int J = 0; J < 3; J++ ){
                    for ( unsigned int K = 0; K < 3; K++ ){
                        dInvRCGPsidPsi[ dim * sot_dim * I + sot_dim * J + dim * K + J ] = invRCG[ dim * I + K ];
                        for ( unsigned int L = 0; L < 3; L++ ){
                            dInvRCGPsidRCG[ dim * sot_dim * I + sot_dim * J + dim * K + L ] = -invRCG[ dim * I + K ] * invRCGPsi[ dim * L + J ];
                        }
                    }
                }
            }
    
            return NULL;
        }
    
        errorOut computeInvRCGGamma( const variableVector &invRCG, const variableVector &Gamma, variableVector &invRCGGamma ){
            /*!
             * Compute the product \f$ C_{IS}^{-1} \Gamma_{SQR} \f$
             *
             * \param &invRCG: The inverse of the right Cauchy Green deformation tensor.
             * \param &Gamma: The gradient of the micro-deformation deformation tensor.
             * \param &invRCGGamma: The product.
             */
    
            //Assume 3d
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;
    
            TARDIGRADE_ERROR_TOOLS_CHECK( invRCG.size() == sot_dim, "invRCG has an improper dimension" );
    
            TARDIGRADE_ERROR_TOOLS_CHECK( Gamma.size() == tot_dim, "Gamma has an improper dimension" );
    
            invRCGGamma = variableVector( tot_dim, 0 );
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int S = 0; S < dim; S++ ){
                    for ( unsigned int Q = 0; Q < dim; Q++ ){
                        for ( unsigned int R = 0; R < dim; R++ ){
                            invRCGGamma[ dim * dim * J + dim * Q + R ] += invRCG[ dim * J + S ] * Gamma[ dim * dim * S + dim * Q + R ];
                        }
                    }
                }
            }
    
            return NULL;
        }
    
        errorOut computeInvRCGGamma( const variableVector &invRCG, const variableVector &Gamma, variableVector &invRCGGamma,
                                     variableVector &dInvRCGGammadRCG, variableVector &dInvRCGGammadGamma ){
            /*!
             * Compute the product \f$ C_{IS}^{-1} \Gamma_{SQR} \f$
             *
             * Also compute the Jacobians
             * 
             * \f$\begin{align}
             * \frac{\partial C_{JS}^{-1} \Gamma_{SQR} }{ \partial C_{TU} } &= -C_{JT}^{-1} C_{US}^{-1} \Gamma_{SQR}\\
             * \frac{\partial C_{JS}^{-1} \Gamma_{SQR} }{ \partial \Gamma_{TUV} } &= C_{JT}^{-1} \delta_{QU} \delta_{RV}
             * \end{align}\f$
             *
             * \param &invRCG: The inverse of the right Cauchy-Green deformation tensor.
             * \param &Gamma: The gradient of the micro-deformation deformation tensor.
             * \param &invRCGGamma: The product.
             * \param &dInvRCGGammadRCG: The derivative of the product w.r.t. the right Cauchy-Green deformation tensor
             * \param &dInvRCGGammadGamma: The derivative of the product w.r.t. \f$ \bf{\Gamma} \f$
             */
    
            //Assume 3d
            constexpr unsigned int dim = 3;
            constexpr unsigned int sot_dim = dim * dim;
            constexpr unsigned int tot_dim = sot_dim * dim;
    
            TARDIGRADE_ERROR_TOOLS_CATCH( computeInvRCGGamma( invRCG, Gamma, invRCGGamma ) );
    
            //Assemble jacobians of invCGamma w.r.t. C and Gamma
            dInvRCGGammadRCG = variableVector( tot_dim * sot_dim, 0 );
            dInvRCGGammadGamma = variableVector( tot_dim * tot_dim, 0 );
    
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int Q = 0; Q < dim; Q++ ){
                    for ( unsigned int R = 0; R < dim; R++ ){
                        for ( unsigned int T = 0; T < dim; T++ ){
                            dInvRCGGammadGamma[ dim * dim * tot_dim * J + dim * tot_dim * Q + tot_dim * R + dim * dim * T + dim * Q + R]
                                = invRCG[ dim * J + T ];
                            for ( unsigned int U = 0; U < dim; U++ ){
                                dInvRCGGammadRCG[ dim * dim * sot_dim * J + dim * sot_dim * Q + sot_dim * R + dim * T + U ]
                                    = -invRCG[ dim * J + T] * invRCGGamma[ dim * dim * U + dim * Q + R ];
                            }
                        }
                    }
                }
            }
    
            return NULL;
        }

        errorOut formIsotropicA( const parameterType &lambda, const parameterType &mu, parameterVector &A ){
            /*!
             * Form the isotropic A stiffness tensor.
             * \f$\begin{align}
             * A_{KLMN} &= \lambda \delta_{KL} \delta_{MN} + \mu \left( \delta_{KM} \delta_{LN} + \delta_{KN} \delta_{LM} \right)
             * \end{align}\f$
             *
             * \param &lambda: The micromorphic lambda parameter.
             * \param &mu: The micromorphic mu parameter.
             * \param &A: The isotropic A stiffness tensor.
             */
    
            //Assume 3D
            constexpr unsigned int dim = 3;
    
            A = parameterVector( dim * dim * dim * dim, 0 );
    
            for ( unsigned int K = 0; K < dim; K++ ){
                for ( unsigned int L = 0; L < dim; L++ ){
                    A[ dim * dim * dim * K + dim * dim * K + dim * L + L ] += lambda;
                    A[ dim * dim * dim * K + dim * dim * L + dim * K + L ] += mu;
                    A[ dim * dim * dim * K + dim * dim * L + dim * L + K ] += mu;
                }
            }
    
            return NULL;
        }
    
        errorOut formIsotropicB( const parameterType &eta, const parameterType &tau,   const parameterType &kappa,
                                 const parameterType &nu,  const parameterType &sigma, parameterVector &B ){
            /*!
             * Form the isotropic B stiffness tensor.
             * \f$\begin{align}
             * B_{KLMN} &= ( eta - tau ) \delta_{KL} \delta_{MN} + \kappa \delta_{KM} \delta_{LN} + \nu \delta_{KN} \delta_{LM}\\
             *          &- \sigma \left( \delta_{KM} \delta_{LN} + \delta_{KN} \delta_{LM} \right)
             * \end{align}\f$
             *
             * \param &eta: The micromorphic eta parameter.
             * \param &tau: The micromorphic tau parameter.
             * \param &kappa: The micromorphic kappa parameter.
             * \param &nu: The micromorphic nu parameter.
             * \param &sigma: The micromorphic sigma parameter
             * \param &B: The isotropic B stiffnes tensor.
             */
    
            //Assume 3D
            constexpr unsigned int dim = 3;
    
            B = parameterVector( dim * dim * dim * dim, 0 );
    
            for ( unsigned int K = 0; K < dim; K++ ){
                for ( unsigned int L = 0; L < dim; L++ ){
                    B[ dim * dim * dim * K + dim * dim * K + dim * L + L ] += ( eta - tau );
                    B[ dim * dim * dim * K + dim * dim * L + dim * K + L ] += kappa;
                    B[ dim * dim * dim * K + dim * dim * L + dim * L + K ] += nu;
                    B[ dim * dim * dim * K + dim * dim * L + dim * K + L ] -= sigma;;
                    B[ dim * dim * dim * K + dim * dim * L + dim * L + K ] -= sigma;
                }
            }
    
            return NULL;
        }
    
        errorOut formIsotropicC( const parameterVector &taus, parameterVector &C ){
            /*!
             * Form the isotropic C stiffness tensor.
             * \f$\begin{align}
             * C_{KLMNPQ} &= \tau_1 \left( \delta_{KL} \delta_{MN} \delta_{PQ} + \delta_{KQ} \delta_{LM} \delta_{NP} \right)\\
             *            &+ \tau_2 \left( \delta_{KL} \delta_{MP} \delta_{NQ} + \delta_{KM} \delta_{LQ} \delta_{NP} \right)\\
             *            &+ \tau_3 \delta_{KL} \delta_{MQ} \delta_{NP}\\
             *            &+ \tau_4 \delta_{KN} \delta_{LM} \delta_{PQ}\\
             *            &+ \tau_5 \left( \delta_{KM} \delta_{LN} \delta_{PQ} + \delta_{KP} \delta_{LM} \delta_{NQ} )\\
             *            &+ \tau_6 \delta_{KM} \delta_{LP} \delta_{NQ}\\
             *            &+ \tau_7 \delta_{KN} \delta_{LP} \delta_{MQ}\\
             *            &+ \tau_8 \left( \delta_{KP} \delta_{LQ} \delta_{MN} + \delta_{KQ} \delta_{LN} \delta_{MP} )\\
             *            &+ \tau_9 \delta_{KN} \delta_{LQ} \delta_{MP}\\
             *            &+ \tau_{10} \delta_{KP} \delta_{LN} \delta_{MQ}\\
             *            &+ \tau_{11} \delta_{KQ} \delta_{LP} \delta_{MN}
             * \end{align}\f$
             *
             * \param &taus: The moduli (11 independent terms)
             * \param &C: The isotropic C stiffness tensor.
             */
    
            //Assume 3D
            constexpr unsigned int dim = 3;
    
            TARDIGRADE_ERROR_TOOLS_CHECK( taus.size() == 11, "11 moduli required to form C" );
    
            C = parameterVector( dim * dim * dim * dim * dim * dim, 0 );
    
            for ( unsigned int K = 0; K < dim; K++ ){
                for ( unsigned int L = 0; L < dim; L++ ){
                    for ( unsigned int M = 0; M < dim; M++ ){
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * K + dim * dim * dim * L + dim * dim * L + dim * M + M ]
                                                         += taus[0];
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * L + dim * dim * M + dim * M + K ]
                                                         += taus[0];
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * K + dim * dim * dim * L + dim * dim * M + dim * L + M ]
                                                         += taus[1];
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * K + dim * dim * M + dim * M + L ]
                                                         += taus[1];
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * K + dim * dim * dim * L + dim * dim * M + dim * M + L ]
                                                         += taus[2];
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * L + dim * dim * K + dim * M + M ]
                                                         += taus[3];
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * K + dim * dim * L + dim * M + M ]
                                                         += taus[4];
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * L + dim * dim * M + dim * K + M ]
                                                         += taus[4];
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * K + dim * dim * M + dim * L + M ]
                                                         += taus[5];
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * M + dim * dim * K + dim * L + M ]
                                                         += taus[6];
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * M + dim * dim * M + dim * K + L ]
                                                         += taus[7];
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * M + dim * dim * L + dim * M + K ]
                                                         += taus[7];
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * M + dim * dim * K + dim * M + L ]
                                                         += taus[8];
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * M + dim * dim * L + dim * K + M ]
                                                         += taus[9];
                        C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * M + dim * dim * M + dim * L + K ]
                                                         += taus[10];
                    }
                }
            }
    
            return NULL;
        }
    
        errorOut formIsotropicD( const parameterType &tau, const parameterType &sigma, parameterVector &D ) {
            /*!
             * Form the isotropic tensor D.
             * \f$ D_{KLMN} = \tau \delta_{KL} \delta_{MN} + \sigma \left( \delta_{KM} \delta_{LN} + \delta_{KN} \delta_{LM} \right) \f$
             *
             * \param &tau: The micromorphic tau parameter.
             * \param &sigma: The micromorphic sigma parameter.
             * \param &D: The D stiffness tensor.
             */
    
            //Assume 3D
            constexpr unsigned int dim = 3;
    
            D = parameterVector( dim * dim * dim * dim, 0 );
            for ( unsigned int K = 0; K < dim; K++ ){
                for ( unsigned int L = 0; L < dim; L++ ){
                    D[ dim * dim * dim * K + dim * dim * K + dim * L + L ] += tau;
                    D[ dim * dim * dim * K + dim * dim * L + dim * K + L ]
                        += sigma;
                    D[ dim * dim * dim * K + dim * dim * L + dim * L + K ]
                        += sigma;
                }
            }
    
            return NULL;
        }

        errorOut assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                         const double ( &grad_phi )[ 9 ][ 3 ],
                                                         variableVector &deformationGradient, variableVector &microDeformation,
                                                         variableVector &gradientMicroDeformation ){
            /*!
             * Assemble the fundamental deformation meaures from the degrees of freedom.
             *
             * \param &grad_u: The macro displacement gradient w.r.t. the reference configuration.
             * \param &phi: The micro displacement.
             * \param &grad_phi: The gradient of the micro displacement w.r.t. the reference configuration.
             * \param &deformationGradient: The deformation gradient
             * \param &microDeformation: The micro deformation
             * \param &gradientMicroDeformation: The gradient of the micro deformation.
             */
    
    
            //Extract the degrees of freedom
            variableVector displacementGradient = { grad_u[ 0 ][ 0 ], grad_u[ 0 ][ 1 ], grad_u[ 0 ][ 2 ],
                                                    grad_u[ 1 ][ 0 ], grad_u[ 1 ][ 1 ], grad_u[ 1 ][ 2 ],
                                                    grad_u[ 2 ][ 0 ], grad_u[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] };
    
            variableVector microDisplacement = { phi[ 0 ], phi[ 1 ], phi[ 2 ],
                                                 phi[ 3 ], phi[ 4 ], phi[ 5 ],
                                                 phi[ 6 ], phi[ 7 ], phi[ 8 ] };
    
            variableVector gradientMicroDisplacement = { grad_phi[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ],
                                                         grad_phi[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ],
                                                         grad_phi[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ],
                                                         grad_phi[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ],
                                                         grad_phi[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ],
                                                         grad_phi[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ],
                                                         grad_phi[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ],
                                                         grad_phi[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ],
                                                         grad_phi[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] };
    
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleDeformationGradient( displacementGradient, deformationGradient ) );
    
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleMicroDeformation( microDisplacement, microDeformation ) );
    
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, gradientMicroDeformation ) );
    
            return NULL;
        }

        errorOut assembleFundamentalDeformationMeasures( const double ( &grad_u )[ 3 ][ 3 ], const double ( &phi )[ 9 ],
                                                         const double ( &grad_phi )[ 9 ][ 3 ],
                                                         variableVector &deformationGradient, variableVector &microDeformation,
                                                         variableVector &gradientMicroDeformation, variableVector &dFdGradU,
                                                         variableVector &dChidPhi, variableVector &dGradChidGradPhi ){
            /*!
             * Assemble the fundamental deformation meaures from the degrees of freedom.
             *
             * \param &grad_u: The macro displacement gradient w.r.t. the reference configuration.
             * \param &phi: The micro displacement.
             * \param &grad_phi: The gradient of the micro displacement w.r.t. the reference configuration.
             * \param &deformationGradient: The deformation gradient
             * \param &microDeformation: The micro deformation
             * \param &gradientMicroDeformation: The gradient of the micro deformation.
             * \param &dFdGradU: The Jacobian of the deformation gradient w.r.t. the gradient of the displacement
             * \param &dChidPhi: The Jacobian of the micro deformation w.r.t. the micro displacement
             * \param &dGradChidGradPhi: The Jacobian of the gradient of the micro deformation w.r.t.
             *      the gradient of the micro displacement
             */
    
    
            //Extract the degrees of freedom
            variableVector displacementGradient = { grad_u[ 0 ][ 0 ], grad_u[ 0 ][ 1 ], grad_u[ 0 ][ 2 ],
                                                    grad_u[ 1 ][ 0 ], grad_u[ 1 ][ 1 ], grad_u[ 1 ][ 2 ],
                                                    grad_u[ 2 ][ 0 ], grad_u[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] };
    
            variableVector microDisplacement = { phi[ 0 ], phi[ 1 ], phi[ 2 ],
                                                 phi[ 3 ], phi[ 4 ], phi[ 5 ],
                                                 phi[ 6 ], phi[ 7 ], phi[ 8 ] };
    
            variableVector gradientMicroDisplacement = { grad_phi[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ],
                                                         grad_phi[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ],
                                                         grad_phi[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ],
                                                         grad_phi[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ],
                                                         grad_phi[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ],
                                                         grad_phi[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ],
                                                         grad_phi[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ],
                                                         grad_phi[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ],
                                                         grad_phi[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] };
    
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleDeformationGradient( displacementGradient, deformationGradient, dFdGradU ) );
    
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleMicroDeformation( microDisplacement, microDeformation, dChidPhi ) );
    
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeMicromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, gradientMicroDeformation,
                                                                                                         dGradChidGradPhi ) );
    
            return NULL;
        }

        errorOut extractMaterialParameters( const std::vector< double > &fparams,
                                            parameterVector &Amatrix, parameterVector &Bmatrix,
                                            parameterVector &Cmatrix, parameterVector &Dmatrix ){
            /*!
             * Extract the parameters from the parameter vector
             *
             * :param const std::vector< double > &fparams: The incoming parameter vector
             * :param parameterVector &Amatrix: The A stiffness matrix.
             * :param parameterVector &Bmatrix: The B stiffness matrix.
             * :param parameterVector &Cmatrix: The C stiffness matrix.
             * :param parameterVector &Dmatrix: The D stiffness matrix.
             */
    
            TARDIGRADE_ERROR_TOOLS_CHECK( fparams.size() != 0, "The material parameters vector has a length of 0" );
    
            unsigned int start = 0;
            unsigned int span;
    
            std::vector< parameterVector > outputs( 4 );
    
            //Extract the material parameters
            for ( unsigned int i = 0; i < outputs.size(); i++ ){
                span = ( unsigned int )std::floor( fparams[ start ]  + 0.5 ); //Extract the span of the parameter set

                TARDIGRADE_ERROR_TOOLS_CHECK( fparams.size( ) >= start + 1 + span, "fparams is not long enough to contain all of the required parameters:\n    filling variable " + std::to_string( i ) + "\n    size =          "  + std::to_string( fparams.size() ) + "\n    required size = "  + std::to_string( start + 1 + span ) + "\n" );
    
                outputs[ i ] = parameterVector( fparams.begin() + start + 1, fparams.begin() + start + 1 + span );
    
                start = start + 1 + span;
            }
    
            //Form the stiffness tensors

            TARDIGRADE_ERROR_TOOLS_CHECK( outputs[ 0 ].size( ) == 2, "Unrecognized number of parameters ( " + std::to_string( outputs[ 0 ].size() ) + " ) for the A stiffness tensor" );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeHydra::micromorphicLinearElasticity::formIsotropicA( outputs[ 0 ][ 0 ], outputs[ 0 ][ 1 ], Amatrix ) );

            TARDIGRADE_ERROR_TOOLS_CHECK( outputs[ 1 ].size() == 5, "Unrecognized number of parameters ( " + std::to_string( outputs[ 1 ].size() ) + " ) for the B stiffness tensor" );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeHydra::micromorphicLinearElasticity::formIsotropicB( outputs[ 1 ][ 0 ], outputs[ 1 ][ 1 ], outputs[ 1 ][ 2 ],
                                                                                                         outputs[ 1 ][ 3 ], outputs[ 1 ][ 4 ], Bmatrix ) );
    
            TARDIGRADE_ERROR_TOOLS_CHECK( outputs[ 2 ].size() == 11, "Unrecognized number of parameters ( " + std::to_string( outputs[ 2 ].size() ) + " ) for the C stiffness tensor" );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeHydra::micromorphicLinearElasticity::formIsotropicC( outputs[ 2 ], Cmatrix ) );
    
            TARDIGRADE_ERROR_TOOLS_CHECK( outputs[ 3 ].size() == 2, "Unrecognized number of parameters ( " + std::to_string( outputs[ 3 ].size() ) + " ) for the D stiffness tensor" );

            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeHydra::micromorphicLinearElasticity::formIsotropicD( outputs[ 3 ][ 0 ], outputs[ 3 ][ 1 ], Dmatrix ) );
    
            return NULL;

        }

        void residual::setRightCauchyGreen( ){
            /*!
             * Set the value of the right Cauchy-Green deformation tensor
             */

            setDeformation( false );

        }

        void residual::setPsi( ){
            /*!
             * Set the value of the micro deformation measure psi
             */

            setDeformation( false );

        }

        void residual::setGamma( ){
            /*!
             * Set the value of the micro deformation measure gamma
             */

            setDeformation( false );

        }

        void residual::setPreviousRightCauchyGreen( ){
            /*!
             * Set the value of the previous right Cauchy-Green deformation tensor
             */

            setDeformation( true );

        }

        void residual::setPreviousPsi( ){
            /*!
             * Set the value of the previous micro deformation measure psi
             */

            setDeformation( true );
        }

        void residual::setPreviousGamma( ){
            /*!
             * Set the value of the previous micro deformation measure gamma
             */

            setDeformation( true );

        }

        void residual::setDeformation( const bool isPrevious ){
            /*!
             * Evaluate the derived deformation measures
             * 
             * We assume that the first configuration in hydra.get_configurations is the elastic one
             *
             * \param isPrevious: Flag for whether the measures to be calculated are in the current or previous configuration
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int tot_dim = hydra->getTOTDimension( );

            floatVector deformationGradient1;

            floatVector microDeformation1;

            floatVector gradientMicroDeformation1;

            if ( isPrevious ){

                deformationGradient1 = floatVector( hydra->get_previousConfigurations( )->begin( ),
                                                    hydra->get_previousConfigurations( )->begin( ) + sot_dim );

                microDeformation1 = floatVector( hydra->get_previousMicroConfigurations( )->begin( ),
                                                 hydra->get_previousMicroConfigurations( )->begin( ) + sot_dim );

                gradientMicroDeformation1 = floatVector( hydra->get_previousGradientMicroConfigurations( )->begin( ),
                                                         hydra->get_previousGradientMicroConfigurations( )->begin( ) + tot_dim );

            }
            else{

                deformationGradient1 = floatVector( hydra->get_configurations( )->begin( ),
                                                    hydra->get_configurations( )->begin( ) + sot_dim );

                microDeformation1 = floatVector( hydra->get_microConfigurations( )->begin( ),
                                                 hydra->get_microConfigurations( )->begin( ) + sot_dim );

                gradientMicroDeformation1 = floatVector( hydra->get_gradientMicroConfigurations( )->begin( ),
                                                         hydra->get_gradientMicroConfigurations( )->begin( ) + tot_dim );

            }

            floatVector rightCauchyGreen;

            floatVector Psi;

            floatVector Gamma;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( computeDeformationMeasures( deformationGradient1, microDeformation1, gradientMicroDeformation1,
                                                                                   rightCauchyGreen, Psi, Gamma ) );

            if ( isPrevious ){

                set_previousRightCauchyGreen( rightCauchyGreen );

                set_previousPsi( Psi );

                set_previousGamma( Gamma );

            }
            else{

                set_rightCauchyGreen( rightCauchyGreen );

                set_psi( Psi );

                set_gamma( Gamma );

            }

        }

        void residual::setdRightCauchyGreendF( ){
            /*!
             * Set the jacobian of the right Cauchy-Green deformation tensor w.r.t. the total deformation gradient
             */

            setDeformationJacobians( false );

        }

        void residual::setdRightCauchyGreendFn( ){
            /*!
             * Set the jacobian of the right Cauchy-Green deformation tensor w.r.t. the remaining sub-deformation gradients
             */

            setDeformationJacobians( false );

        }

        void residual::setdPsidF( ){
            /*!
             * Set the jacobian of the micro deformation measure psi w.r.t. the total deformation gradient
             */

            setDeformationJacobians( false );

        }

        void residual::setdPsidFn( ){
            /*!
             * Set the jacobian of the micro deformation tensor psi w.r.t. the remaining sub-deformation gradients
             */

            setDeformationJacobians( false );

        }

        void residual::setdPsidChi( ){
            /*!
             * Set the jacobian of the micro deformation measure psi w.r.t. the total micro-deformation
             */

            setDeformationJacobians( false );

        }

        void residual::setdPsidChin( ){
            /*!
             * Set the jacobian of the micro deformation tensor psi w.r.t. the remaining sub-micro deformations
             */

            setDeformationJacobians( false );

        }

        void residual::setdGammadF( ){
            /*!
             * Set the jacobian of the micro deformation measure gamma w.r.t. the total deformation gradient
             */

            setDeformationJacobians( false );

        }

        void residual::setdGammadFn( ){
            /*!
             * Set the jacobian of the micro deformation tensor gamma w.r.t. the remaining sub-deformation gradients
             */

            setDeformationJacobians( false );

        }

        void residual::setdGammadChi( ){
            /*!
             * Set the jacobian of the micro deformation measure gamma w.r.t. the total micro-deformation
             */

            setDeformationJacobians( false );

        }

        void residual::setdGammadChin( ){
            /*!
             * Set the jacobian of the micro deformation tensor gamma w.r.t. the remaining sub-micro deformations
             */

            setDeformationJacobians( false );

        }

        void residual::setdGammadGradChi( ){
            /*!
             * Set the jacobian of the micro deformation measure gamma w.r.t. the reference spatial gradient of the total micro-deformation
             */

            setDeformationJacobians( false );

        }

        void residual::setdGammadGradChin( ){
            /*!
             * Set the jacobian of the micro deformation tensor gamma w.r.t. the local reference spatial gradient of the remaining sub-micro deformations
             */

            setDeformationJacobians( false );

        }

        void residual::setPreviousdRightCauchyGreendF( ){
            /*!
             * Set the jacobian of the previous right Cauchy-Green deformation tensor w.r.t. the total deformation gradient
             */

            setDeformationJacobians( true );

        }

        void residual::setPreviousdRightCauchyGreendFn( ){
            /*!
             * Set the jacobian of the previous right Cauchy-Green deformation tensor w.r.t. the remaining sub-deformation gradients
             */

            setDeformationJacobians( true );

        }

        void residual::setPreviousdPsidF( ){
            /*!
             * Set the jacobian of the previous micro deformation measure psi w.r.t. the total deformation gradient
             */

            setDeformationJacobians( true );

        }

        void residual::setPreviousdPsidFn( ){
            /*!
             * Set the jacobian of the previous micro deformation tensor psi w.r.t. the remaining sub-deformation gradients
             */

            setDeformationJacobians( true );

        }

        void residual::setPreviousdPsidChi( ){
            /*!
             * Set the jacobian of the previous micro deformation measure psi w.r.t. the total micro-deformation
             */

            setDeformationJacobians( true );

        }

        void residual::setPreviousdPsidChin( ){
            /*!
             * Set the jacobian of the previous micro deformation tensor psi w.r.t. the remaining sub-micro deformations
             */

            setDeformationJacobians( true );

        }

        void residual::setPreviousdGammadF( ){
            /*!
             * Set the jacobian of the previous micro deformation measure gamma w.r.t. the total deformation gradient
             */

            setDeformationJacobians( true );

        }

        void residual::setPreviousdGammadFn( ){
            /*!
             * Set the jacobian of the previous micro deformation tensor gamma w.r.t. the remaining sub-deformation gradients
             */

            setDeformationJacobians( true );

        }

        void residual::setPreviousdGammadChi( ){
            /*!
             * Set the jacobian of the previous micro deformation measure gamma w.r.t. the total micro-deformation
             */

            setDeformationJacobians( true );

        }

        void residual::setPreviousdGammadChin( ){
            /*!
             * Set the jacobian of the previous micro deformation tensor gamma w.r.t. the remaining sub-micro deformations
             */

            setDeformationJacobians( true );

        }

        void residual::setPreviousdGammadGradChi( ){
            /*!
             * Set the jacobian of the previous micro deformation measure gamma w.r.t. the reference spatial gradient of the total micro-deformation
             */

            setDeformationJacobians( true );

        }

        void residual::setPreviousdGammadGradChin( ){
            /*!
             * Set the jacobian of the previous micro deformation tensor gamma w.r.t. the local reference spatial gradient of the remaining sub-micro deformations
             */

            setDeformationJacobians( true );

        }

        void residual::setPK2Stress( ){
            /*!
             * Set the value of the second Piola-Kirchhoff stress
             */

            return setReferenceStresses( false );

        }

        void residual::setReferenceSymmetricMicroStress( ){
            /*!
             * Set the value of the reference symmetric micro stress
             */

            return setReferenceStresses( false );

        }

        void residual::setReferenceHigherOrderStress( ){
            /*!
             * Set the value of the reference higher order stress
             */

            return setReferenceStresses( false );

        }

        void residual::setPreviousPK2Stress( ){
            /*!
             * Set the value of the previous second Piola-Kirchhoff stress
             */

            return setReferenceStresses( true );

        }

        void residual::setPreviousReferenceSymmetricMicroStress( ){
            /*!
             * Set the value of the previous reference symmetric micro stress
             */

            return setReferenceStresses( true );

        }

        void residual::setPreviousReferenceHigherOrderStress( ){
            /*!
             * Set the value of the previous reference higher order stress
             */

            return setReferenceStresses( true );

        }

        void residual::setReferenceStresses( const bool isPrevious ){
            /*!
             * Set the values of the reference stresses
             * 
             * \param isPrevious: Flag for if the stresses to be calculated are the current (false) or previous (true)
             */

            const variableVector *C;

            const variableVector *Psi;

            const variableVector *Gamma;

            variableVector followingConfiguration;

            variableVector followingMicroConfiguration;

            if ( isPrevious ){

                C     = get_previousRightCauchyGreen( );

                Psi   = get_previousPsi( );

                Gamma = get_previousGamma( );

                followingConfiguration = hydra->getPreviousFollowingConfiguration( 0 );

                followingMicroConfiguration = hydra->getPreviousFollowingMicroConfiguration( 0 );

            }
            else{

                C     = get_rightCauchyGreen( );

                Psi   = get_psi( );

                Gamma = get_gamma( );

                followingConfiguration = hydra->getFollowingConfiguration( 0 );
    
                followingMicroConfiguration = hydra->getFollowingMicroConfiguration( 0 );

            }

            floatVector localPK2Stress;

            floatVector localReferenceSymmetricMicroStress;

            floatVector localReferenceHigherOrderStress;

            // Compute the stresses in the local configuration
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( linearElasticityReferenceDerivedMeasures( *C, *Psi, *Gamma, *getAMatrix( ), *getBMatrix( ), *getCMatrix( ), *getDMatrix( ),
                                                                                                 localPK2Stress, localReferenceSymmetricMicroStress, localReferenceHigherOrderStress ) );

            floatVector PK2Stress;

            floatVector referenceSymmetricMicroStress;

            floatVector referenceHigherOrderStress;

            // Pull the stresses back to the true reference configuration
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pullBackCauchyStress( localPK2Stress, followingConfiguration, PK2Stress ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pullBackMicroStress( localReferenceSymmetricMicroStress, followingConfiguration, referenceSymmetricMicroStress ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pullBackHigherOrderStress( localReferenceHigherOrderStress, followingConfiguration, followingMicroConfiguration,
                                                                                                               referenceHigherOrderStress ) );

            if ( isPrevious ){

                set_previousPK2Stress( PK2Stress );

                set_previousReferenceSymmetricMicroStress( referenceSymmetricMicroStress );

                set_previousReferenceHigherOrderStress( referenceHigherOrderStress );

            }
            else{

                set_PK2Stress( PK2Stress );

                set_referenceSymmetricMicroStress( referenceSymmetricMicroStress );

                set_referenceHigherOrderStress( referenceHigherOrderStress );

            }

        }

        void residual::setdPK2dF( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress tensor w.r.t. the deformation gradient
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdPK2dFn( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress tensor w.r.t. the sub-deformation gradients
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdPK2dChi( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress tensor w.r.t. the micro deformation
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdPK2dChin( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress tensor w.r.t. the sub-micro deformations
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdPK2dGradChi( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress tensor w.r.t. the spatial reference gradient of the micro deformations
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdPK2dGradChin( ){
            /*!
             * Set the derivative of the second Piola-Kirchhoff stress tensor w.r.t. the spatial local reference gradient of the sub-micro deformations
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdSIGMAdF( ){
            /*!
             * Set the derivative of the reference symmetric micro stress tensor w.r.t. the deformation gradient
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdSIGMAdFn( ){
            /*!
             * Set the derivative of the reference symmetric micro stress tensor w.r.t. the sub-deformation gradients
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdSIGMAdChi( ){
            /*!
             * Set the derivative of the reference symmetric micro stress tensor w.r.t. the micro deformation
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdSIGMAdChin( ){
            /*!
             * Set the derivative of the reference symmetric micro stress tensor w.r.t. the sub-micro deformations
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdSIGMAdGradChi( ){
            /*!
             * Set the derivative of the reference symmetric micro stress tensor w.r.t. the spatial reference gradient of the micro deformations
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdSIGMAdGradChin( ){
            /*!
             * Set the derivative of the reference symmetric micro stress tensor w.r.t. the spatial local reference gradient of the sub-micro deformations
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdMdF( ){
            /*!
             * Set the derivative of the reference higher order stress tensor w.r.t. the deformation gradient
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdMdFn( ){
            /*!
             * Set the derivative of the reference higher order micro stress tensor w.r.t. the sub-deformation gradients
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdMdChi( ){
            /*!
             * Set the derivative of the reference higher order stress tensor w.r.t. the micro deformation
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdMdChin( ){
            /*!
             * Set the derivative of the reference higher order stress tensor w.r.t. the sub-micro deformations
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdMdGradChi( ){
            /*!
             * Set the derivative of the reference higher order stress tensor w.r.t. the spatial reference gradient of the micro deformations
             */

            setReferenceStressJacobians( false );

        }

        void residual::setdMdGradChin( ){
            /*!
             * Set the derivative of the reference higher order stress tensor w.r.t. the spatial local reference gradient of the sub-micro deformations
             */

            setReferenceStressJacobians( false );

        }

        void residual::setPreviousdPK2dF( ){
            /*!
             * Set the previous derivative of the second Piola-Kirchhoff stress tensor w.r.t. the deformation gradient
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdPK2dFn( ){
            /*!
             * Set the previous derivative of the second Piola-Kirchhoff stress tensor w.r.t. the sub-deformation gradients
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdPK2dChi( ){
            /*!
             * Set the previous derivative of the second Piola-Kirchhoff stress tensor w.r.t. the micro deformation
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdPK2dChin( ){
            /*!
             * Set the previous derivative of the second Piola-Kirchhoff stress tensor w.r.t. the sub-micro deformations
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdPK2dGradChi( ){
            /*!
             * Set the previous derivative of the second Piola-Kirchhoff stress tensor w.r.t. the spatial reference gradient of the micro deformations
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdPK2dGradChin( ){
            /*!
             * Set the previous derivative of the second Piola-Kirchhoff stress tensor w.r.t. the spatial local reference gradient of the sub-micro deformations
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdSIGMAdF( ){
            /*!
             * Set the previous derivative of the reference symmetric micro stress tensor w.r.t. the deformation gradient
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdSIGMAdFn( ){
            /*!
             * Set the previous derivative of the reference symmetric micro stress tensor w.r.t. the sub-deformation gradients
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdSIGMAdChi( ){
            /*!
             * Set the derivative of the reference symmetric micro stress tensor w.r.t. the micro deformation
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdSIGMAdChin( ){
            /*!
             * Set the previous derivative of the reference symmetric micro stress tensor w.r.t. the sub-micro deformations
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdSIGMAdGradChi( ){
            /*!
             * Set the previous derivative of the reference symmetric micro stress tensor w.r.t. the spatial reference gradient of the micro deformations
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdSIGMAdGradChin( ){
            /*!
             * Set the previous derivative of the reference symmetric micro stress tensor w.r.t. the spatial local reference gradient of the sub-micro deformations
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdMdF( ){
            /*!
             * Set the previous derivative of the reference higher order stress tensor w.r.t. the deformation gradient
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdMdFn( ){
            /*!
             * Set the previous derivative of the reference higher order micro stress tensor w.r.t. the sub-deformation gradients
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdMdChi( ){
            /*!
             * Set the previous derivative of the reference higher order stress tensor w.r.t. the micro deformation
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdMdChin( ){
            /*!
             * Set the previous derivative of the reference higher order stress tensor w.r.t. the sub-micro deformations
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdMdGradChi( ){
            /*!
             * Set the previous derivative of the reference higher order stress tensor w.r.t. the spatial reference gradient of the micro deformations
             */

            setReferenceStressJacobians( true );

        }

        void residual::setPreviousdMdGradChin( ){
            /*!
             * Set the previous derivative of the reference higher order stress tensor w.r.t. the spatial local reference gradient of the sub-micro deformations
             */

            setReferenceStressJacobians( true );

        }

        void residual::setReferenceStressJacobians( const bool isPrevious ){
            /*!
             * Set the values and Jacobians of the reference stresses
             * 
             * \param isPrevious: Flag for if the stresses to be calculated are the current (false) or previous (true)
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            constexpr unsigned int tot_dim = sot_dim * dim;

            constexpr unsigned int fot_dim = tot_dim * dim;

            constexpr unsigned int fiot_dim = fot_dim * dim;

            constexpr unsigned int siot_dim = fiot_dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const variableVector *C;

            const variableVector *Psi;

            const variableVector *Gamma;

            const variableVector *dCdF;

            const variableVector *dCdFn;

            const variableVector *dPsidF;

            const variableVector *dPsidFn;

            const variableVector *dPsidChi;

            const variableVector *dPsidChin;

            const variableVector *dGammadF;

            const variableVector *dGammadFn;

            const variableVector *dGammadChi;

            const variableVector *dGammadChin;

            const variableVector *dGammadGradChi;

            const variableVector *dGammadGradChin;

            variableVector dFFollowdFs;

            variableVector dChiFollowdChis;

            variableVector followingConfiguration;

            variableVector followingMicroConfiguration;

            setDataStorageBase< secondOrderTensor > PK2Stress;

            setDataStorageBase< secondOrderTensor > referenceSymmetricMicroStress;

            setDataStorageBase< thirdOrderTensor  > referenceHigherOrderStress;

            setDataStorageBase< fourthOrderTensor > dPK2dF;

            setDataStorageBase< floatVector >       dPK2dFn;

            setDataStorageBase< fourthOrderTensor > dPK2dChi;

            setDataStorageBase< floatVector >       dPK2dChin;

            setDataStorageBase< fifthOrderTensor >  dPK2dGradChi;

            setDataStorageBase< floatVector >       dPK2dGradChin;

            setDataStorageBase< fourthOrderTensor > dSIGMAdF;

            setDataStorageBase< floatVector >       dSIGMAdFn;

            setDataStorageBase< fourthOrderTensor > dSIGMAdChi;

            setDataStorageBase< floatVector >       dSIGMAdChin;

            setDataStorageBase< fifthOrderTensor >  dSIGMAdGradChi;

            setDataStorageBase< floatVector >       dSIGMAdGradChin;

            setDataStorageBase< fifthOrderTensor >  dMdF;

            setDataStorageBase< floatVector >       dMdFn;

            setDataStorageBase< fifthOrderTensor >  dMdChi;

            setDataStorageBase< floatVector >       dMdChin;

            setDataStorageBase< sixthOrderTensor >  dMdGradChi;

            setDataStorageBase< floatVector >       dMdGradChin;

            if ( isPrevious ){

                dCdF            = get_previousdRightCauchyGreendF( );

                dCdFn           = get_previousdRightCauchyGreendFn( );

                C               = get_previousRightCauchyGreen( );

                dPsidF          = get_previousdPsidF( );

                dPsidFn         = get_previousdPsidFn( );

                dPsidChi        = get_previousdPsidChi( );

                dPsidChin       = get_previousdPsidChin( );

                Psi             = get_previousPsi( );

                dGammadF        = get_previousdGammadF( );

                dGammadFn       = get_previousdGammadFn( );

                dGammadChi      = get_previousdGammadChi( );

                dGammadChin     = get_previousdGammadChin( );

                dGammadGradChi  = get_previousdGammadGradChi( );

                dGammadGradChin = get_previousdGammadGradChin( );

                Gamma           = get_previousGamma( );

                dFFollowdFs     = hydra->getPreviousFollowingConfigurationJacobian( 0 );

                dChiFollowdChis = hydra->getPreviousFollowingMicroConfigurationJacobian( 0 );

                followingConfiguration = hydra->getPreviousFollowingConfiguration( 0 );

                followingMicroConfiguration = hydra->getPreviousFollowingMicroConfiguration( 0 );

                PK2Stress                     = get_setDataStorage_previousPK2Stress( );

                referenceSymmetricMicroStress = get_setDataStorage_previousReferenceSymmetricMicroStress( );

                referenceHigherOrderStress    = get_setDataStorage_previousReferenceHigherOrderStress( );

                dPK2dF                        = get_setDataStorage_previousdPK2dF( );

                dPK2dFn                       = get_setDataStorage_previousdPK2dFn( );

                dPK2dChi                      = get_setDataStorage_previousdPK2dChi( );

                dPK2dChin                     = get_setDataStorage_previousdPK2dChin( );

                dPK2dGradChi                  = get_setDataStorage_previousdPK2dGradChi( );

                dPK2dGradChin                 = get_setDataStorage_previousdPK2dGradChin( );

                dSIGMAdF                      = get_setDataStorage_previousdSIGMAdF( );

                dSIGMAdFn                     = get_setDataStorage_previousdSIGMAdFn( );

                dSIGMAdChi                    = get_setDataStorage_previousdSIGMAdChi( );

                dSIGMAdChin                   = get_setDataStorage_previousdSIGMAdChin( );

                dSIGMAdGradChi                = get_setDataStorage_previousdSIGMAdGradChi( );

                dSIGMAdGradChin               = get_setDataStorage_previousdSIGMAdGradChin( );

                dMdF                          = get_setDataStorage_previousdMdF( );

                dMdFn                         = get_setDataStorage_previousdMdFn( );

                dMdChi                        = get_setDataStorage_previousdMdChi( );

                dMdChin                       = get_setDataStorage_previousdMdChin( );

                dMdGradChi                    = get_setDataStorage_previousdMdGradChi( );

                dMdGradChin                   = get_setDataStorage_previousdMdGradChin( );

            }
            else{

                dCdF            = get_dRightCauchyGreendF( );

                dCdFn           = get_dRightCauchyGreendFn( );

                C               = get_rightCauchyGreen( );

                dPsidF          = get_dPsidF( );

                dPsidFn         = get_dPsidFn( );

                dPsidChi        = get_dPsidChi( );

                dPsidChin       = get_dPsidChin( );

                Psi             = get_psi( );

                dGammadF        = get_dGammadF( );

                dGammadFn       = get_dGammadFn( );

                dGammadChi      = get_dGammadChi( );

                dGammadChin     = get_dGammadChin( );

                dGammadGradChi  = get_dGammadGradChi( );

                dGammadGradChin = get_dGammadGradChin( );

                Gamma           = get_gamma( );

                dFFollowdFs     = hydra->getFollowingConfigurationJacobian( 0 );

                dChiFollowdChis = hydra->getFollowingMicroConfigurationJacobian( 0 );

                followingConfiguration = hydra->getFollowingConfiguration( 0 );
    
                followingMicroConfiguration = hydra->getFollowingMicroConfiguration( 0 );

                PK2Stress                     = get_setDataStorage_PK2Stress( );

                referenceSymmetricMicroStress = get_setDataStorage_referenceSymmetricMicroStress( );

                referenceHigherOrderStress    = get_setDataStorage_referenceHigherOrderStress( );

                dPK2dF                        = get_setDataStorage_dPK2dF( );

                dPK2dFn                       = get_setDataStorage_dPK2dFn( );

                dPK2dChi                      = get_setDataStorage_dPK2dChi( );

                dPK2dChin                     = get_setDataStorage_dPK2dChin( );

                dPK2dGradChi                  = get_setDataStorage_dPK2dGradChi( );

                dPK2dGradChin                 = get_setDataStorage_dPK2dGradChin( );

                dSIGMAdF                      = get_setDataStorage_dSIGMAdF( );

                dSIGMAdFn                     = get_setDataStorage_dSIGMAdFn( );

                dSIGMAdChi                    = get_setDataStorage_dSIGMAdChi( );

                dSIGMAdChin                   = get_setDataStorage_dSIGMAdChin( );

                dSIGMAdGradChi                = get_setDataStorage_dSIGMAdGradChi( );

                dSIGMAdGradChin               = get_setDataStorage_dSIGMAdGradChin( );

                dMdF                          = get_setDataStorage_dMdF( );

                dMdFn                         = get_setDataStorage_dMdFn( );

                dMdChi                        = get_setDataStorage_dMdChi( );

                dMdChin                       = get_setDataStorage_dMdChin( );

                dMdGradChi                    = get_setDataStorage_dMdGradChi( );

                dMdGradChin                   = get_setDataStorage_dMdGradChin( );

            }

            variableVector dFFollowdFn( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            variableVector dChiFollowdChin( sot_dim * ( num_configs - 1 ) * sot_dim, 0 );

            for ( unsigned int i = 0; i < sot_dim; i++ ){

                for ( unsigned int j = 0; j < ( num_configs - 1 ) * sot_dim; j++ ){

                    dFFollowdFn[ ( num_configs - 1 ) * sot_dim * i + j ]     = dFFollowdFs[ num_configs * sot_dim * i + j + sot_dim ];

                    dChiFollowdChin[ ( num_configs - 1 ) * sot_dim * i + j ] = dChiFollowdChis[ num_configs * sot_dim * i + j + sot_dim ];

                }

            }

            floatVector localPK2Stress;

            floatVector localReferenceSymmetricMicroStress;

            floatVector localReferenceHigherOrderStress;

            floatVector dlocalPK2dC;

            floatVector dlocalPK2dPsi;

            floatVector dlocalPK2dGamma;

            floatVector dlocalSIGMAdC;

            floatVector dlocalSIGMAdPsi;

            floatVector dlocalSIGMAdGamma;

            floatVector dlocalMdGamma;

            // Compute the stresses in the local configuration
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( linearElasticityReferenceDerivedMeasures( *C, *Psi, *Gamma, *getAMatrix( ), *getBMatrix( ), *getCMatrix( ), *getDMatrix( ),
                                                                                                 localPK2Stress, localReferenceSymmetricMicroStress, localReferenceHigherOrderStress,
                                                                                                 dlocalPK2dC, dlocalPK2dPsi, dlocalPK2dGamma, dlocalSIGMAdC, dlocalSIGMAdPsi, dlocalSIGMAdGamma,
                                                                                                 dlocalMdGamma ) );

            floatVector dLocalPK2dF = tardigradeVectorTools::matrixMultiply(     dlocalPK2dC,     *dCdF, sot_dim, sot_dim, sot_dim, sot_dim )
                                    + tardigradeVectorTools::matrixMultiply(   dlocalPK2dPsi,   *dPsidF, sot_dim, sot_dim, sot_dim, sot_dim )
                                    + tardigradeVectorTools::matrixMultiply( dlocalPK2dGamma, *dGammadF, sot_dim, tot_dim, tot_dim, sot_dim );

            floatVector dLocalPK2dFn = tardigradeVectorTools::matrixMultiply(     dlocalPK2dC,     *dCdFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                     + tardigradeVectorTools::matrixMultiply(   dlocalPK2dPsi,   *dPsidFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                     + tardigradeVectorTools::matrixMultiply( dlocalPK2dGamma, *dGammadFn, sot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dLocalPK2dChi = tardigradeVectorTools::matrixMultiply(   dlocalPK2dPsi,   *dPsidChi, sot_dim, sot_dim, sot_dim, sot_dim )
                                      + tardigradeVectorTools::matrixMultiply( dlocalPK2dGamma, *dGammadChi, sot_dim, tot_dim, tot_dim, sot_dim );

            floatVector dLocalPK2dChin = tardigradeVectorTools::matrixMultiply(   dlocalPK2dPsi,   *dPsidChin, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                       + tardigradeVectorTools::matrixMultiply( dlocalPK2dGamma, *dGammadChin, sot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dLocalPK2dGradChi = tardigradeVectorTools::matrixMultiply( dlocalPK2dGamma, *dGammadGradChi, sot_dim, tot_dim, tot_dim, tot_dim );

            floatVector dLocalPK2dGradChin = tardigradeVectorTools::matrixMultiply( dlocalPK2dGamma, *dGammadGradChin, sot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim );

            floatVector dLocalSIGMAdF = tardigradeVectorTools::matrixMultiply(     dlocalSIGMAdC,     *dCdF, sot_dim, sot_dim, sot_dim, sot_dim )
                                      + tardigradeVectorTools::matrixMultiply(   dlocalSIGMAdPsi,   *dPsidF, sot_dim, sot_dim, sot_dim, sot_dim )
                                      + tardigradeVectorTools::matrixMultiply( dlocalSIGMAdGamma, *dGammadF, sot_dim, tot_dim, tot_dim, sot_dim );

            floatVector dLocalSIGMAdFn = tardigradeVectorTools::matrixMultiply(     dlocalSIGMAdC,     *dCdFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                       + tardigradeVectorTools::matrixMultiply(   dlocalSIGMAdPsi,   *dPsidFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                       + tardigradeVectorTools::matrixMultiply( dlocalSIGMAdGamma, *dGammadFn, sot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dLocalSIGMAdChi = tardigradeVectorTools::matrixMultiply(   dlocalSIGMAdPsi,   *dPsidChi, sot_dim, sot_dim, sot_dim, sot_dim )
                                        + tardigradeVectorTools::matrixMultiply( dlocalSIGMAdGamma, *dGammadChi, sot_dim, tot_dim, tot_dim, sot_dim );

            floatVector dLocalSIGMAdChin = tardigradeVectorTools::matrixMultiply(   dlocalSIGMAdPsi,   *dPsidChin, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                         + tardigradeVectorTools::matrixMultiply( dlocalSIGMAdGamma, *dGammadChin, sot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dLocalSIGMAdGradChi = tardigradeVectorTools::matrixMultiply( dlocalSIGMAdGamma, *dGammadGradChi, sot_dim, tot_dim, tot_dim, tot_dim );

            floatVector dLocalSIGMAdGradChin = tardigradeVectorTools::matrixMultiply( dlocalSIGMAdGamma, *dGammadGradChin, sot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim );

            floatVector dLocalMdF = tardigradeVectorTools::matrixMultiply( dlocalMdGamma, *dGammadF, tot_dim, tot_dim, tot_dim, sot_dim );

            floatVector dLocalMdFn = tardigradeVectorTools::matrixMultiply( dlocalMdGamma, *dGammadFn, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dLocalMdChi = tardigradeVectorTools::matrixMultiply( dlocalMdGamma, *dGammadChi, tot_dim, tot_dim, tot_dim, sot_dim );

            floatVector dLocalMdChin = tardigradeVectorTools::matrixMultiply( dlocalMdGamma, *dGammadChin, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dLocalMdGradChi = tardigradeVectorTools::matrixMultiply( dlocalMdGamma, *dGammadGradChi, tot_dim, tot_dim, tot_dim, tot_dim );

            floatVector dLocalMdGradChin = tardigradeVectorTools::matrixMultiply( dlocalMdGamma, *dGammadGradChin, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim );

            floatVector dPK2dlocalPK2;

            floatVector dPK2dFFollow;

            floatVector dSIGMAdlocalSIGMA;

            floatVector dSIGMAdFFollow;

            floatVector dMdlocalM;

            floatVector dMdFFollow;

            floatVector dMdChiFollow;

            // Pull the stresses back to the true reference configuration
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pullBackCauchyStress( localPK2Stress, followingConfiguration, *PK2Stress.value, dPK2dlocalPK2, dPK2dFFollow ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pullBackMicroStress( localReferenceSymmetricMicroStress, followingConfiguration, *referenceSymmetricMicroStress.value,
                                                                                                         dSIGMAdlocalSIGMA, dSIGMAdFFollow ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pullBackHigherOrderStress( localReferenceHigherOrderStress, followingConfiguration, followingMicroConfiguration,
                                                                                                               *referenceHigherOrderStress.value, dMdlocalM, dMdFFollow, dMdChiFollow ) );

            *dPK2dF.value = tardigradeVectorTools::matrixMultiply( dPK2dlocalPK2, dLocalPK2dF, sot_dim, sot_dim, sot_dim, sot_dim );

            *dPK2dFn.value = tardigradeVectorTools::matrixMultiply( dPK2dlocalPK2, dLocalPK2dFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                           + tardigradeVectorTools::matrixMultiply( dPK2dFFollow, dFFollowdFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            *dPK2dChi.value = tardigradeVectorTools::matrixMultiply( dPK2dlocalPK2, dLocalPK2dChi, sot_dim, sot_dim, sot_dim, sot_dim );

            *dPK2dChin.value = tardigradeVectorTools::matrixMultiply( dPK2dlocalPK2, dLocalPK2dChin, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            *dPK2dGradChi.value = tardigradeVectorTools::matrixMultiply( dPK2dlocalPK2, dLocalPK2dGradChi, sot_dim, sot_dim, sot_dim, tot_dim );

            *dPK2dGradChin.value = tardigradeVectorTools::matrixMultiply( dPK2dlocalPK2, dLocalPK2dGradChin, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * tot_dim );

            *dSIGMAdF.value = tardigradeVectorTools::matrixMultiply( dSIGMAdlocalSIGMA, dLocalSIGMAdF, sot_dim, sot_dim, sot_dim, sot_dim );

            *dSIGMAdFn.value = tardigradeVectorTools::matrixMultiply( dSIGMAdlocalSIGMA, dLocalSIGMAdFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                             + tardigradeVectorTools::matrixMultiply( dSIGMAdFFollow, dFFollowdFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            *dSIGMAdChi.value = tardigradeVectorTools::matrixMultiply( dSIGMAdlocalSIGMA, dLocalSIGMAdChi, sot_dim, sot_dim, sot_dim, sot_dim );

            *dSIGMAdChin.value = tardigradeVectorTools::matrixMultiply( dSIGMAdlocalSIGMA, dLocalSIGMAdChin, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            *dSIGMAdGradChi.value = tardigradeVectorTools::matrixMultiply( dSIGMAdlocalSIGMA, dLocalSIGMAdGradChi, sot_dim, sot_dim, sot_dim, tot_dim );

            *dSIGMAdGradChin.value = tardigradeVectorTools::matrixMultiply( dSIGMAdlocalSIGMA, dLocalSIGMAdGradChin, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * tot_dim );

            *dMdF.value = tardigradeVectorTools::matrixMultiply( dMdlocalM, dLocalMdF, tot_dim, tot_dim, tot_dim, sot_dim );

            *dMdFn.value = tardigradeVectorTools::matrixMultiply( dMdlocalM,  dLocalMdFn,  tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim )
                         + tardigradeVectorTools::matrixMultiply( dMdFFollow, dFFollowdFn, tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            *dMdChi.value = tardigradeVectorTools::matrixMultiply( dMdlocalM, dLocalMdChi, tot_dim, tot_dim, tot_dim, sot_dim );

            *dMdChin.value = tardigradeVectorTools::matrixMultiply( dMdlocalM,    dLocalMdChin, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim )
                           + tardigradeVectorTools::matrixMultiply( dMdChiFollow, dChiFollowdChin, tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            *dMdGradChi.value = tardigradeVectorTools::matrixMultiply( dMdlocalM, dLocalMdGradChi, tot_dim, tot_dim, tot_dim, tot_dim );

            *dMdGradChin.value = tardigradeVectorTools::matrixMultiply( dMdlocalM, dLocalMdGradChin, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim );

        }

        void residual::setCauchyStress( ){
            /*!
             * Set the value of the Cauchy stress
             */

            setStresses( false );

        }

        void residual::setSymmetricMicroStress( ){
            /*!
             * Set the value of the symmetric micro stress
             */

            setStresses( false );

        }

        void residual::setHigherOrderStress( ){
            /*!
             * Set the value of the higher order stress
             */

            setStresses( false );

        }

        void residual::setPreviousCauchyStress( ){
            /*!
             * Set the value of the previous Cauchy stress
             */

            setStresses( true );

        }

        void residual::setPreviousSymmetricMicroStress( ){
            /*!
             * Set the value of the previous symmetric micro stress
             */

            setStresses( true );

        }

        void residual::setPreviousHigherOrderStress( ){
            /*!
             * Set the value of the previous higher order stress
             */

            setStresses( true );

        }

        void residual::setStresses( const bool isPrevious ){
            /*!
             * Set the values of the stresses in the current configuration
             * 
             * \param isPrevious: Flag for whether to compute the previous (true) or current (false) stresses
             */

            const floatVector *PK2Stress;

            const floatVector *referenceSymmetricMicroStress;

            const floatVector *referenceHigherOrderStress;

            const floatVector *deformationGradient;

            const floatVector *microDeformation;

            setDataStorageBase< secondOrderTensor > cauchyStress;

            setDataStorageBase< secondOrderTensor > symmetricMicroStress;

            setDataStorageBase< thirdOrderTensor >  higherOrderStress;

            if ( isPrevious ){

                PK2Stress = get_previousPK2Stress( );

                referenceSymmetricMicroStress = get_previousReferenceSymmetricMicroStress( );

                referenceHigherOrderStress = get_previousReferenceHigherOrderStress( );

                deformationGradient = hydra->getPreviousDeformationGradient( );

                microDeformation = hydra->getPreviousMicroDeformation( );

                cauchyStress         = get_setDataStorage_previousCauchyStress( );

                symmetricMicroStress = get_setDataStorage_previousSymmetricMicroStress( );

                higherOrderStress    = get_setDataStorage_previousHigherOrderStress( );

            }
            else{

                PK2Stress = get_PK2Stress( );

                referenceSymmetricMicroStress = get_referenceSymmetricMicroStress( );

                referenceHigherOrderStress = get_referenceHigherOrderStress( );

                deformationGradient = hydra->getDeformationGradient( );

                microDeformation = hydra->getMicroDeformation( );

                cauchyStress         = get_setDataStorage_cauchyStress( );

                symmetricMicroStress = get_setDataStorage_symmetricMicroStress( );

                higherOrderStress    = get_setDataStorage_higherOrderStress( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( mapStressMeasuresToCurrent( *deformationGradient, *microDeformation, *PK2Stress, *referenceSymmetricMicroStress,
                                                                                   *referenceHigherOrderStress, *cauchyStress.value, *symmetricMicroStress.value, *higherOrderStress.value )  );

        }

        void residual::setdCauchyStressdF( ){
            /*!
             * Set the jacobian of the cauchy stress w.r.t. the deformation gradient
             */

            setStressesJacobians( false );

        }

        void residual::setdCauchyStressdFn( ){
            /*!
             * Set the jacobian of the cauchy stress w.r.t. the sub-deformation gradients
             */

            setStressesJacobians( false );

        }

        void residual::setdCauchyStressdChi( ){
            /*!
             * Set the jacobian of the cauchy stress w.r.t. the micro deformation
             */

            setStressesJacobians( false );

        }

        void residual::setdCauchyStressdChin( ){
            /*!
             * Set the jacobian of the cauchy stress w.r.t. the sub-micro deformations
             */

            setStressesJacobians( false );

        }

        void residual::setdCauchyStressdGradChi( ){
            /*!
             * Set the jacobian of the cauchy stress w.r.t. the spatial gradient of the micro deformation
             */

            setStressesJacobians( false );

        }

        void residual::setdCauchyStressdGradChin( ){
            /*!
             * Set the jacobian of the cauchy stress w.r.t. the local reference spatial gradient of the sub-micro deformations
             */

            setStressesJacobians( false );

        }

        void residual::setdSymmetricMicroStressdF( ){
            /*!
             * Set the jacobian of the symmetric micro stress w.r.t. the deformation gradient
             */

            setStressesJacobians( false );

        }

        void residual::setdSymmetricMicroStressdFn( ){
            /*!
             * Set the jacobian of the symmetric micro stress w.r.t. the sub-deformation gradients
             */

            setStressesJacobians( false );

        }

        void residual::setdSymmetricMicroStressdChi( ){
            /*!
             * Set the jacobian of the symmetric micro stress w.r.t. the micro deformation
             */

            setStressesJacobians( false );

        }

        void residual::setdSymmetricMicroStressdChin( ){
            /*!
             * Set the jacobian of the symmetric micro stress w.r.t. the sub-micro deformations
             */

            setStressesJacobians( false );

        }

        void residual::setdSymmetricMicroStressdGradChi( ){
            /*!
             * Set the jacobian of the symmetric micro stress w.r.t. the spatial gradient of the micro deformation
             */

            setStressesJacobians( false );

        }

        void residual::setdSymmetricMicroStressdGradChin( ){
            /*!
             * Set the jacobian of the symmetric micro stress w.r.t. the local reference spatial gradient of the sub micro deformations
             */

            setStressesJacobians( false );

        }

        void residual::setdHigherOrderStressdF( ){
            /*!
             * Set the jacobian of the higher order stress w.r.t. the deformation gradient
             */

            setStressesJacobians( false );

        }

        void residual::setdHigherOrderStressdFn( ){
            /*!
             * Set the jacobian of the higher order stress w.r.t. the sub-deformation gradients
             */

            setStressesJacobians( false );

        }

        void residual::setdHigherOrderStressdChi( ){
            /*!
             * Set the jacobian of the higher order stress w.r.t. the micro deformation
             */

            setStressesJacobians( false );

        }

        void residual::setdHigherOrderStressdChin( ){
            /*!
             * Set the jacobian of the higher order stress w.r.t. the sub micro deformations
             */

            setStressesJacobians( false );

        }

        void residual::setdHigherOrderStressdGradChi( ){
            /*!
             * Set the jacobian of the higher order stress w.r.t. the reference spatial gradient of the micro deformation
             */

            setStressesJacobians( false );

        }

        void residual::setdHigherOrderStressdGradChin( ){
            /*!
             * Set the jacobian of the higher order stress w.r.t. the local reference spatial gradient of the sub micro deformations
             */

            setStressesJacobians( false );

        }

        void residual::setPreviousdCauchyStressdF( ){
            /*!
             * Set the jacobian of the previous cauchy stress w.r.t. the deformation gradient
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdCauchyStressdFn( ){
            /*!
             * Set the jacobian of the previous cauchy stress w.r.t. the sub-deformation gradients
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdCauchyStressdChi( ){
            /*!
             * Set the jacobian of the previous cauchy stress w.r.t. the micro deformation
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdCauchyStressdChin( ){
            /*!
             * Set the jacobian of the previous cauchy stress w.r.t. the sub-micro deformations
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdCauchyStressdGradChi( ){
            /*!
             * Set the jacobian of the previous cauchy stress w.r.t. the spatial gradient of the micro deformation
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdCauchyStressdGradChin( ){
            /*!
             * Set the jacobian of the previous cauchy stress w.r.t. the local reference spatial gradient of the sub-micro deformations
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdSymmetricMicroStressdF( ){
            /*!
             * Set the jacobian of the previous symmetric micro stress w.r.t. the deformation gradient
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdSymmetricMicroStressdFn( ){
            /*!
             * Set the jacobian of the previous symmetric micro stress w.r.t. the sub-deformation gradients
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdSymmetricMicroStressdChi( ){
            /*!
             * Set the jacobian of the previous symmetric micro stress w.r.t. the micro deformation
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdSymmetricMicroStressdChin( ){
            /*!
             * Set the jacobian of the previous symmetric micro stress w.r.t. the sub-micro deformations
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdSymmetricMicroStressdGradChi( ){
            /*!
             * Set the jacobian of the previous symmetric micro stress w.r.t. the spatial gradient of the micro deformation
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdSymmetricMicroStressdGradChin( ){
            /*!
             * Set the jacobian of the previous symmetric micro stress w.r.t. the local reference spatial gradient of the sub micro deformations
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderStressdF( ){
            /*!
             * Set the jacobian of the previous higher order stress w.r.t. the deformation gradient
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderStressdFn( ){
            /*!
             * Set the jacobian of the previous higher order stress w.r.t. the sub-deformation gradients
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderStressdChi( ){
            /*!
             * Set the jacobian of the previous higher order stress w.r.t. the micro deformation
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderStressdChin( ){
            /*!
             * Set the jacobian of the previous higher order stress w.r.t. the sub micro deformations
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderStressdGradChi( ){
            /*!
             * Set the jacobian of the previous higher order stress w.r.t. the reference spatial gradient of the micro deformation
             */

            setStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderStressdGradChin( ){
            /*!
             * Set the jacobian of the previous higher order stress w.r.t. the local reference spatial gradient of the sub micro deformations
             */

            setStressesJacobians( true );

        }

        void residual::setStressesJacobians( const bool isPrevious ){
            /*!
             * Set the stresses and their jacobians in the current configuration
             * 
             * \param isPrevious: A flag for whether to compute the previous (true) or current (false) stresses and their Jacobians
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            constexpr unsigned int tot_dim = sot_dim * dim;

            constexpr unsigned int fot_dim = tot_dim * dim;

            constexpr unsigned int fiot_dim = fot_dim * dim;

            constexpr unsigned int siot_dim = fiot_dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const secondOrderTensor *PK2Stress;

            const secondOrderTensor *referenceSymmetricMicroStress;

            const thirdOrderTensor  *referenceHigherOrderStress;

            const secondOrderTensor *deformationGradient;

            const secondOrderTensor *microDeformation;

            const fourthOrderTensor *dPK2dF;

            const floatVector *dPK2dFn;

            const fourthOrderTensor *dPK2dChi;

            const floatVector *dPK2dChin;

            const fifthOrderTensor *dPK2dGradChi;

            const floatVector *dPK2dGradChin;

            const fourthOrderTensor *dSIGMAdF;

            const floatVector *dSIGMAdFn;

            const fourthOrderTensor *dSIGMAdChi;

            const floatVector *dSIGMAdChin;

            const fifthOrderTensor *dSIGMAdGradChi;

            const floatVector *dSIGMAdGradChin;

            const fifthOrderTensor *dMdF;

            const floatVector *dMdFn;

            const fifthOrderTensor *dMdChi;

            const floatVector *dMdChin;

            const sixthOrderTensor *dMdGradChi;

            const floatVector *dMdGradChin;

            setDataStorageBase< secondOrderTensor > cauchyStress;

            setDataStorageBase< secondOrderTensor > symmetricMicroStress;

            setDataStorageBase< thirdOrderTensor >  higherOrderStress;

            setDataStorageBase< fourthOrderTensor > dCauchyStressdF;

            setDataStorageBase< floatVector >       dCauchyStressdFn;

            setDataStorageBase< fourthOrderTensor > dCauchyStressdChi;

            setDataStorageBase< floatVector >       dCauchyStressdChin;

            setDataStorageBase< fifthOrderTensor >  dCauchyStressdGradChi;

            setDataStorageBase< floatVector >       dCauchyStressdGradChin;

            setDataStorageBase< fourthOrderTensor > dMicroStressdF;

            setDataStorageBase< floatVector >       dMicroStressdFn;

            setDataStorageBase< fourthOrderTensor > dMicroStressdChi;

            setDataStorageBase< floatVector >       dMicroStressdChin;

            setDataStorageBase< fifthOrderTensor >  dMicroStressdGradChi;

            setDataStorageBase< floatVector >       dMicroStressdGradChin;

            setDataStorageBase< fifthOrderTensor >  dHigherOrderStressdF;

            setDataStorageBase< floatVector >       dHigherOrderStressdFn;

            setDataStorageBase< fifthOrderTensor >  dHigherOrderStressdChi;

            setDataStorageBase< floatVector >       dHigherOrderStressdChin;

            setDataStorageBase< sixthOrderTensor >  dHigherOrderStressdGradChi;

            setDataStorageBase< floatVector >       dHigherOrderStressdGradChin;

            if ( isPrevious ){

                dPK2dF          = get_previousdPK2dF( );

                dPK2dFn         = get_previousdPK2dFn( );

                dPK2dChi        = get_previousdPK2dChi( );

                dPK2dChin       = get_previousdPK2dChin( );

                dPK2dGradChi    = get_previousdPK2dGradChi( );

                dPK2dGradChin   = get_previousdPK2dGradChin( );

                dSIGMAdF        = get_previousdSIGMAdF( );

                dSIGMAdFn       = get_previousdSIGMAdFn( );

                dSIGMAdChi      = get_previousdSIGMAdChi( );

                dSIGMAdChin     = get_previousdSIGMAdChin( );

                dSIGMAdGradChi  = get_previousdSIGMAdGradChi( );

                dSIGMAdGradChin = get_previousdSIGMAdGradChin( );

                dMdF            = get_previousdMdF( );

                dMdFn           = get_previousdMdFn( );

                dMdChi          = get_previousdMdChi( );

                dMdChin         = get_previousdMdChin( );

                dMdGradChi      = get_previousdMdGradChi( );

                dMdGradChin     = get_previousdMdGradChin( );

                PK2Stress = get_previousPK2Stress( );

                referenceSymmetricMicroStress = get_previousReferenceSymmetricMicroStress( );

                referenceHigherOrderStress = get_previousReferenceHigherOrderStress( );

                deformationGradient = hydra->getPreviousDeformationGradient( );

                microDeformation = hydra->getPreviousMicroDeformation( );

                cauchyStress                = get_setDataStorage_previousCauchyStress( );

                symmetricMicroStress        = get_setDataStorage_previousSymmetricMicroStress( );

                higherOrderStress           = get_setDataStorage_previousHigherOrderStress( );

                dCauchyStressdF             = get_setDataStorage_previousdCauchyStressdF( );

                dCauchyStressdFn            = get_setDataStorage_previousdCauchyStressdFn( );

                dCauchyStressdChi           = get_setDataStorage_previousdCauchyStressdChi( );

                dCauchyStressdChin          = get_setDataStorage_previousdCauchyStressdChin( );

                dCauchyStressdGradChi       = get_setDataStorage_previousdCauchyStressdGradChi( );

                dCauchyStressdGradChin      = get_setDataStorage_previousdCauchyStressdGradChin( );

                dMicroStressdF              = get_setDataStorage_previousdSymmetricMicroStressdF( );

                dMicroStressdFn             = get_setDataStorage_previousdSymmetricMicroStressdFn( );

                dMicroStressdChi            = get_setDataStorage_previousdSymmetricMicroStressdChi( );

                dMicroStressdChin           = get_setDataStorage_previousdSymmetricMicroStressdChin( );

                dMicroStressdGradChi        = get_setDataStorage_previousdSymmetricMicroStressdGradChi( );

                dMicroStressdGradChin       = get_setDataStorage_previousdSymmetricMicroStressdGradChin( );

                dHigherOrderStressdF        = get_setDataStorage_previousdHigherOrderStressdF( );

                dHigherOrderStressdFn       = get_setDataStorage_previousdHigherOrderStressdFn( );

                dHigherOrderStressdChi      = get_setDataStorage_previousdHigherOrderStressdChi( );

                dHigherOrderStressdChin     = get_setDataStorage_previousdHigherOrderStressdChin( );

                dHigherOrderStressdGradChi  = get_setDataStorage_previousdHigherOrderStressdGradChi( );

                dHigherOrderStressdGradChin = get_setDataStorage_previousdHigherOrderStressdGradChin( );

            }
            else{

                dPK2dF          = get_dPK2dF( );

                dPK2dFn         = get_dPK2dFn( );

                dPK2dChi        = get_dPK2dChi( );

                dPK2dChin       = get_dPK2dChin( );

                dPK2dGradChi    = get_dPK2dGradChi( );

                dPK2dGradChin   = get_dPK2dGradChin( );

                dSIGMAdF        = get_dSIGMAdF( );

                dSIGMAdFn       = get_dSIGMAdFn( );

                dSIGMAdChi      = get_dSIGMAdChi( );

                dSIGMAdChin     = get_dSIGMAdChin( );

                dSIGMAdGradChi  = get_dSIGMAdGradChi( );

                dSIGMAdGradChin = get_dSIGMAdGradChin( );

                dMdF            = get_dMdF( );

                dMdFn           = get_dMdFn( );

                dMdChi          = get_dMdChi( );

                dMdChin         = get_dMdChin( );

                dMdGradChi      = get_dMdGradChi( );

                dMdGradChin     = get_dMdGradChin( );

                PK2Stress = get_PK2Stress( );

                referenceSymmetricMicroStress = get_referenceSymmetricMicroStress( );

                referenceHigherOrderStress = get_referenceHigherOrderStress( );

                deformationGradient = hydra->getDeformationGradient( );

                microDeformation = hydra->getMicroDeformation( );

                cauchyStress                = get_setDataStorage_cauchyStress( );

                symmetricMicroStress        = get_setDataStorage_symmetricMicroStress( );

                higherOrderStress           = get_setDataStorage_higherOrderStress( );

                dCauchyStressdF             = get_setDataStorage_dCauchyStressdF( );

                dCauchyStressdFn            = get_setDataStorage_dCauchyStressdFn( );

                dCauchyStressdChi           = get_setDataStorage_dCauchyStressdChi( );

                dCauchyStressdChin          = get_setDataStorage_dCauchyStressdChin( );

                dCauchyStressdGradChi       = get_setDataStorage_dCauchyStressdGradChi( );

                dCauchyStressdGradChin      = get_setDataStorage_dCauchyStressdGradChin( );

                dMicroStressdF              = get_setDataStorage_dSymmetricMicroStressdF( );

                dMicroStressdFn             = get_setDataStorage_dSymmetricMicroStressdFn( );

                dMicroStressdChi            = get_setDataStorage_dSymmetricMicroStressdChi( );

                dMicroStressdChin           = get_setDataStorage_dSymmetricMicroStressdChin( );

                dMicroStressdGradChi        = get_setDataStorage_dSymmetricMicroStressdGradChi( );

                dMicroStressdGradChin       = get_setDataStorage_dSymmetricMicroStressdGradChin( );

                dHigherOrderStressdF        = get_setDataStorage_dHigherOrderStressdF( );

                dHigherOrderStressdFn       = get_setDataStorage_dHigherOrderStressdFn( );

                dHigherOrderStressdChi      = get_setDataStorage_dHigherOrderStressdChi( );

                dHigherOrderStressdChin     = get_setDataStorage_dHigherOrderStressdChin( );

                dHigherOrderStressdGradChi  = get_setDataStorage_dHigherOrderStressdGradChi( );

                dHigherOrderStressdGradChin = get_setDataStorage_dHigherOrderStressdGradChin( );

            }

            fourthOrderTensor dCauchyStressdPK2Stress;

            fourthOrderTensor dMicroStressdSIGMA;

            sixthOrderTensor dHigherOrderStressdM;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( mapStressMeasuresToCurrent( *deformationGradient, *microDeformation, *PK2Stress, *referenceSymmetricMicroStress,
                                                                                   *referenceHigherOrderStress, *cauchyStress.value, *symmetricMicroStress.value, *higherOrderStress.value,
                                                                                   *dCauchyStressdF.value, dCauchyStressdPK2Stress,
                                                                                   *dMicroStressdF.value,  dMicroStressdSIGMA,
                                                                                   *dHigherOrderStressdF.value, *dHigherOrderStressdChi.value, dHigherOrderStressdM )  );

            Eigen::Map<       Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dCauchyStressdF(         dCauchyStressdF.value->data( ),  sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dCauchyStressdPK2Stress( dCauchyStressdPK2Stress.data( ), sot_dim, sot_dim );

            Eigen::Map<       Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dMicroStressdF(          dMicroStressdF.value->data( ),   sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dMicroStressdSIGMA(      dMicroStressdSIGMA.data( ),      sot_dim, sot_dim );

            Eigen::Map<       Eigen::Matrix< floatType, tot_dim, sot_dim, Eigen::RowMajor > > map_dHigherOrderStressdF(    dHigherOrderStressdF.value->data( ),    tot_dim, sot_dim );

            Eigen::Map<       Eigen::Matrix< floatType, tot_dim, sot_dim, Eigen::RowMajor > > map_dHigherOrderStressdChi(  dHigherOrderStressdChi.value->data( ),  tot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, tot_dim, tot_dim, Eigen::RowMajor > > map_dHigherOrderStressdM(    dHigherOrderStressdM.data( ),    tot_dim, tot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dPK2dF(          dPK2dF->data( ),          sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dPK2dFn(         dPK2dFn->data( ),         sot_dim, sot_dim * ( num_configs - 1 ) );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dPK2dChi(        dPK2dChi->data( ),        sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dPK2dChin(       dPK2dChin->data( ),       sot_dim, sot_dim * ( num_configs - 1 ) );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, tot_dim, Eigen::RowMajor > > map_dPK2dGradChi(    dPK2dGradChi->data( ),    sot_dim, tot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dPK2dGradChin(   dPK2dGradChin->data( ),   sot_dim, tot_dim * ( num_configs - 1 ) );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dSIGMAdF(        dSIGMAdF->data( ),        sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dSIGMAdFn(       dSIGMAdFn->data( ),       sot_dim, sot_dim * ( num_configs - 1 ) );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dSIGMAdChi(      dSIGMAdChi->data( ),      sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dSIGMAdChin(     dSIGMAdChin->data( ),     sot_dim, sot_dim * ( num_configs - 1 ) );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, tot_dim, Eigen::RowMajor > > map_dSIGMAdGradChi(  dSIGMAdGradChi->data( ),  sot_dim, tot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dSIGMAdGradChin( dSIGMAdGradChin->data( ), sot_dim, tot_dim * ( num_configs - 1 ) );

            Eigen::Map< const Eigen::Matrix< floatType, tot_dim, sot_dim, Eigen::RowMajor > > map_dMdF(            dMdF->data( ),            tot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, tot_dim,      -1, Eigen::RowMajor > > map_dMdFn(           dMdFn->data( ),           tot_dim, sot_dim * ( num_configs - 1 ) );

            Eigen::Map< const Eigen::Matrix< floatType, tot_dim, sot_dim, Eigen::RowMajor > > map_dMdChi(          dMdChi->data( ),          tot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, tot_dim,      -1, Eigen::RowMajor > > map_dMdChin(         dMdChin->data( ),         tot_dim, sot_dim * ( num_configs - 1 ) );

            Eigen::Map< const Eigen::Matrix< floatType, tot_dim, tot_dim, Eigen::RowMajor > > map_dMdGradChi(      dMdGradChi->data( ),      tot_dim, tot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, tot_dim,      -1, Eigen::RowMajor > > map_dMdGradChin(     dMdGradChin->data( ),     tot_dim, tot_dim * ( num_configs - 1 ) );

            dCauchyStressdFn.zero( fot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dCauchyStressdFn( dCauchyStressdFn.value->data( ), sot_dim, sot_dim * ( num_configs - 1 ) );

            dCauchyStressdChi.zero( fot_dim );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dCauchyStressdChi( dCauchyStressdChi.value->data( ), sot_dim, sot_dim );

            dCauchyStressdChin.zero( fot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dCauchyStressdChin( dCauchyStressdChin.value->data( ), sot_dim, sot_dim * ( num_configs - 1 ) );

            dCauchyStressdGradChi.zero( fiot_dim );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim, tot_dim, Eigen::RowMajor > > map_dCauchyStressdGradChi( dCauchyStressdGradChi.value->data( ),   sot_dim, tot_dim );

            dCauchyStressdGradChin.zero( fiot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dCauchyStressdGradChin( dCauchyStressdGradChin.value->data( ), sot_dim, tot_dim * ( num_configs - 1 ) );

            dMicroStressdFn.zero( fot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dMicroStressdFn( dMicroStressdFn.value->data( ), sot_dim, sot_dim * ( num_configs - 1 ) );

            dMicroStressdChi.zero( fot_dim );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dMicroStressdChi( dMicroStressdChi.value->data( ), sot_dim, sot_dim );

            dMicroStressdChin.zero( fot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dMicroStressdChin( dMicroStressdChin.value->data( ), sot_dim, sot_dim * ( num_configs - 1 ) );

            dMicroStressdGradChi.zero( fiot_dim );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim, tot_dim, Eigen::RowMajor > > map_dMicroStressdGradChi( dMicroStressdGradChi.value->data( ), sot_dim, tot_dim );

            dMicroStressdGradChin.zero( fiot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dMicroStressdGradChin( dMicroStressdGradChin.value->data( ), sot_dim, tot_dim * ( num_configs - 1 ) );

            dHigherOrderStressdFn.zero( fiot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, tot_dim,      -1, Eigen::RowMajor > > map_dHigherOrderStressdFn( dHigherOrderStressdFn.value->data( ), tot_dim, sot_dim * ( num_configs - 1 ) );

            dHigherOrderStressdChin.zero( fiot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, tot_dim,      -1, Eigen::RowMajor > > map_dHigherOrderStressdChin( dHigherOrderStressdChin.value->data( ), tot_dim, sot_dim * ( num_configs - 1 ) );

            dHigherOrderStressdGradChi.zero( siot_dim );
            Eigen::Map< Eigen::Matrix< floatType, tot_dim, tot_dim, Eigen::RowMajor > > map_dHigherOrderStressdGradChi( dHigherOrderStressdGradChi.value->data( ), tot_dim, tot_dim );

            dHigherOrderStressdGradChin.zero( siot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, tot_dim,      -1, Eigen::RowMajor > > map_dHigherOrderStressdGradChin( dHigherOrderStressdGradChin.value->data( ), tot_dim, tot_dim * ( num_configs - 1 ) );



            map_dCauchyStressdF += ( map_dCauchyStressdPK2Stress * map_dPK2dF ).eval( );

            map_dCauchyStressdFn = ( map_dCauchyStressdPK2Stress * map_dPK2dFn ).eval( );

            map_dCauchyStressdChi  = ( map_dCauchyStressdPK2Stress * map_dPK2dChi ).eval( );

            map_dCauchyStressdChin = ( map_dCauchyStressdPK2Stress * map_dPK2dChin ).eval( );

            map_dCauchyStressdGradChi  = ( map_dCauchyStressdPK2Stress * map_dPK2dGradChi ).eval( );

            map_dCauchyStressdGradChin = ( map_dCauchyStressdPK2Stress * map_dPK2dGradChin ).eval( );


            map_dMicroStressdF += ( map_dMicroStressdSIGMA * map_dSIGMAdF ).eval( );

            map_dMicroStressdFn = ( map_dMicroStressdSIGMA * map_dSIGMAdFn ).eval( );

            map_dMicroStressdChi  = ( map_dMicroStressdSIGMA * map_dSIGMAdChi ).eval( );

            map_dMicroStressdChin = ( map_dMicroStressdSIGMA * map_dSIGMAdChin ).eval( );

            map_dMicroStressdGradChi  = ( map_dMicroStressdSIGMA * map_dSIGMAdGradChi ).eval( );

            map_dMicroStressdGradChin = ( map_dMicroStressdSIGMA * map_dSIGMAdGradChin ).eval( );


            map_dHigherOrderStressdF += ( map_dHigherOrderStressdM * map_dMdF ).eval( );

            map_dHigherOrderStressdFn = ( map_dHigherOrderStressdM * map_dMdFn ).eval( );

            map_dHigherOrderStressdChi += ( map_dHigherOrderStressdM * map_dMdChi ).eval( );

            map_dHigherOrderStressdChin = ( map_dHigherOrderStressdM * map_dMdChin ).eval( );

            map_dHigherOrderStressdGradChi  = ( map_dHigherOrderStressdM * map_dMdGradChi ).eval( );

            map_dHigherOrderStressdGradChin = ( map_dHigherOrderStressdM * map_dMdGradChin ).eval( );

        }

        void residual::setDeformationJacobians( const bool isPrevious ){
            /*!
             * Evaluate the derived deformation Jacobians
             * 
             * We assume that the first configuration in hydra.get_configurations is the elastic one
             *
             * \param isPrevious: Flag for whether the measures to be calculated are in the current or previous configuration
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            constexpr unsigned int tot_dim = sot_dim * dim;

            constexpr unsigned int fot_dim = tot_dim * dim;

            constexpr unsigned int fiot_dim = fot_dim * dim;

            constexpr unsigned int siot_dim = fiot_dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            floatVector deformationGradient1;

            floatVector microDeformation1;

            floatVector gradientMicroDeformation1;

            const floatVector *dF1dF;

            const floatVector *dF1dFn;

            const floatVector *dChi1dChi;

            const floatVector *dChi1dChin;

            const floatVector *dGradChi1dFn;

            const floatVector *dGradChi1dChi;

            const floatVector *dGradChi1dChin;

            const floatVector *dGradChi1dGradChi;

            const floatVector *dGradChi1dGradChin;

            setDataStorageBase< secondOrderTensor > rightCauchyGreen;

            setDataStorageBase< secondOrderTensor > Psi;

            setDataStorageBase< thirdOrderTensor >  Gamma;

            setDataStorageBase< fourthOrderTensor > dRightCauchyGreendF;

            setDataStorageBase< floatVector > dRightCauchyGreendFn;

            setDataStorageBase< fourthOrderTensor > dPsidF;

            setDataStorageBase< floatVector > dPsidFn;

            setDataStorageBase< fourthOrderTensor > dPsidChi;

            setDataStorageBase< floatVector > dPsidChin;

            setDataStorageBase< fifthOrderTensor > dGammadF;

            setDataStorageBase< floatVector > dGammadFn;

            setDataStorageBase< fifthOrderTensor > dGammadChi;

            setDataStorageBase< floatVector > dGammadChin;

            setDataStorageBase< sixthOrderTensor > dGammadGradChi;

            setDataStorageBase< floatVector > dGammadGradChin;

            if ( isPrevious ){

                deformationGradient1 = floatVector( hydra->get_previousConfigurations( )->begin( ),
                                                    hydra->get_previousConfigurations( )->begin( ) + sot_dim );

                microDeformation1 = floatVector( hydra->get_previousMicroConfigurations( )->begin( ),
                                                 hydra->get_previousMicroConfigurations( )->begin( ) + sot_dim );

                gradientMicroDeformation1 = floatVector( hydra->get_previousGradientMicroConfigurations( )->begin( ),
                                                         hydra->get_previousGradientMicroConfigurations( )->begin( ) + tot_dim );

                dF1dF              = hydra->get_previousdF1dF( );

                dF1dFn             = hydra->get_previousdF1dFn( );

                dChi1dChi          = hydra->get_previousdChi1dChi( );

                dChi1dChin         = hydra->get_previousdChi1dChin( );

                dGradChi1dFn       = hydra->get_previousdGradChi1dFn( );

                dGradChi1dChi      = hydra->get_previousdGradChi1dChi( );

                dGradChi1dChin     = hydra->get_previousdGradChi1dChin( );

                dGradChi1dGradChi  = hydra->get_previousdGradChi1dGradChi( );

                dGradChi1dGradChin = hydra->get_previousdGradChi1dGradChin( );

                rightCauchyGreen   = get_setDataStorage_previousRightCauchyGreen( );

                Psi                = get_setDataStorage_previousPsi( );

                Gamma              = get_setDataStorage_previousGamma( );

                dRightCauchyGreendF = get_setDataStorage_previousdRightCauchyGreendF( );

                dRightCauchyGreendFn = get_setDataStorage_previousdRightCauchyGreendFn( );

                dPsidF               = get_setDataStorage_previousdPsidF( );

                dPsidFn              = get_setDataStorage_previousdPsidFn( );

                dPsidChi             = get_setDataStorage_previousdPsidChi( );

                dPsidChin            = get_setDataStorage_previousdPsidChin( );

                dGammadF             = get_setDataStorage_previousdGammadF( );

                dGammadFn            = get_setDataStorage_previousdGammadFn( );

                dGammadChi           = get_setDataStorage_previousdGammadChi( );

                dGammadChin          = get_setDataStorage_previousdGammadChin( );

                dGammadGradChi       = get_setDataStorage_previousdGammadGradChi( );

                dGammadGradChin      = get_setDataStorage_previousdGammadGradChin( );

            }
            else{

                deformationGradient1 = floatVector( hydra->get_configurations( )->begin( ),
                                                    hydra->get_configurations( )->begin( ) + sot_dim );

                microDeformation1 = floatVector( hydra->get_microConfigurations( )->begin( ),
                                                 hydra->get_microConfigurations( )->begin( ) + sot_dim );

                gradientMicroDeformation1 = floatVector( hydra->get_gradientMicroConfigurations( )->begin( ),
                                                         hydra->get_gradientMicroConfigurations( )->begin( ) + tot_dim );

                dF1dF              = hydra->get_dF1dF( );

                dF1dFn             = hydra->get_dF1dFn( );

                dChi1dChi          = hydra->get_dChi1dChi( );

                dChi1dChin         = hydra->get_dChi1dChin( );

                dGradChi1dFn       = hydra->get_dGradChi1dFn( );

                dGradChi1dChi      = hydra->get_dGradChi1dChi( );

                dGradChi1dChin     = hydra->get_dGradChi1dChin( );

                dGradChi1dGradChi  = hydra->get_dGradChi1dGradChi( );

                dGradChi1dGradChin = hydra->get_dGradChi1dGradChin( );

                rightCauchyGreen   = get_setDataStorage_rightCauchyGreen( );

                Psi                = get_setDataStorage_psi( );

                Gamma              = get_setDataStorage_gamma( );

                dRightCauchyGreendF = get_setDataStorage_dRightCauchyGreendF( );

                dRightCauchyGreendFn = get_setDataStorage_dRightCauchyGreendFn( );

                dPsidF               = get_setDataStorage_dPsidF( );

                dPsidFn              = get_setDataStorage_dPsidFn( );

                dPsidChi             = get_setDataStorage_dPsidChi( );

                dPsidChin            = get_setDataStorage_dPsidChin( );

                dGammadF             = get_setDataStorage_dGammadF( );

                dGammadFn            = get_setDataStorage_dGammadFn( );

                dGammadChi           = get_setDataStorage_dGammadChi( );

                dGammadChin          = get_setDataStorage_dGammadChin( );

                dGammadGradChi       = get_setDataStorage_dGammadGradChi( );

                dGammadGradChin      = get_setDataStorage_dGammadGradChin( );

            }

            floatVector dCdF1;

            floatVector dPsidF1;

            floatVector dPsidChi1;

            floatVector dGammadF1;

            floatVector dGammadGradChi1;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( computeDeformationMeasures( deformationGradient1, microDeformation1, gradientMicroDeformation1,
                                                                                   *rightCauchyGreen.value, *Psi.value, *Gamma.value,
                                                                                   dCdF1, dPsidF1, dPsidChi1, dGammadF1, dGammadGradChi1 ) );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dCdF1( dCdF1.data( ), sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dPsidF1( dPsidF1.data( ), sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dPsidChi1( dPsidChi1.data( ), sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, tot_dim, sot_dim, Eigen::RowMajor > > map_dGammadF1( dGammadF1.data( ), tot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, tot_dim, tot_dim, Eigen::RowMajor > > map_dGammadGradChi1( dGammadGradChi1.data( ), tot_dim, tot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dF1dF( dF1dF->data( ), sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dF1dFn( dF1dFn->data( ), sot_dim, ( num_configs - 1 ) * sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dChi1dChi( dChi1dChi->data( ), sot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dChi1dChin( dChi1dChin->data( ), sot_dim, ( num_configs - 1 ) * sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, tot_dim,      -1, Eigen::RowMajor > > map_dGradChi1dFn( dGradChi1dFn->data( ), tot_dim, ( num_configs - 1 ) * sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, tot_dim, sot_dim, Eigen::RowMajor > > map_dGradChi1dChi( dGradChi1dChi->data( ), tot_dim, sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, tot_dim,      -1, Eigen::RowMajor > > map_dGradChi1dChin( dGradChi1dChin->data( ), tot_dim, ( num_configs - 1 ) * sot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, tot_dim, tot_dim, Eigen::RowMajor > > map_dGradChi1dGradChi( dGradChi1dGradChi->data( ), tot_dim, tot_dim );

            Eigen::Map< const Eigen::Matrix< floatType, tot_dim,      -1, Eigen::RowMajor > > map_dGradChi1dGradChin( dGradChi1dGradChin->data( ), tot_dim, ( num_configs - 1 ) * tot_dim );

            dRightCauchyGreendF.zero( fot_dim );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dRightCauchyGreendF( dRightCauchyGreendF.value->data( ), sot_dim, sot_dim );

            map_dRightCauchyGreendF = ( map_dCdF1 * map_dF1dF ).eval( );

            dRightCauchyGreendFn.zero( fot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dRightCauchyGreendFn( dRightCauchyGreendFn.value->data( ), sot_dim, ( num_configs - 1 ) * sot_dim );

            map_dRightCauchyGreendFn = ( map_dCdF1 * map_dF1dFn ).eval( );

            dPsidF.zero( fot_dim );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dPsidF( dPsidF.value->data( ), sot_dim, sot_dim );

            map_dPsidF = ( map_dPsidF1 * map_dF1dF ).eval( );

            dPsidFn.zero( fot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dPsidFn( dPsidFn.value->data( ), sot_dim, ( num_configs - 1 ) * sot_dim );

            map_dPsidFn = ( map_dPsidF1 * map_dF1dFn ).eval( );

            dPsidChi.zero( fot_dim );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim, sot_dim, Eigen::RowMajor > > map_dPsidChi( dPsidChi.value->data( ), sot_dim, sot_dim );

            map_dPsidChi = ( map_dPsidChi1 * map_dChi1dChi ).eval( );

            dPsidChin.zero( fot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, sot_dim,      -1, Eigen::RowMajor > > map_dPsidChin( dPsidChin.value->data( ), sot_dim, ( num_configs - 1 ) * sot_dim );

            map_dPsidChin = ( map_dPsidChi1 * map_dChi1dChin ).eval( );

            dGammadF.zero( fiot_dim );
            Eigen::Map< Eigen::Matrix< floatType, tot_dim, sot_dim, Eigen::RowMajor > > map_dGammadF( dGammadF.value->data( ), tot_dim, sot_dim );

            map_dGammadF = ( map_dGammadF1 * map_dF1dF ).eval( );

            dGammadFn.zero( fiot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, tot_dim,      -1, Eigen::RowMajor > > map_dGammadFn( dGammadFn.value->data( ), tot_dim, ( num_configs - 1 ) * sot_dim );

            map_dGammadFn  = ( map_dGammadF1 * map_dF1dFn ).eval( );
            map_dGammadFn += ( map_dGammadGradChi1 * map_dGradChi1dFn ).eval( );

            dGammadChi.zero( fiot_dim );
            Eigen::Map< Eigen::Matrix< floatType, tot_dim, sot_dim, Eigen::RowMajor > > map_dGammadChi( dGammadChi.value->data( ), tot_dim, sot_dim );

            map_dGammadChi = ( map_dGammadGradChi1 * map_dGradChi1dChi ).eval( );

            dGammadChin.zero( fiot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, tot_dim,      -1, Eigen::RowMajor > > map_dGammadChin( dGammadChin.value->data( ), tot_dim, ( num_configs - 1 ) * sot_dim );

            map_dGammadChin  = ( map_dGammadGradChi1 * map_dGradChi1dChin ).eval( );

            dGammadGradChi.zero( siot_dim );
            Eigen::Map< Eigen::Matrix< floatType, tot_dim, tot_dim, Eigen::RowMajor > > map_dGammadGradChi( dGammadGradChi.value->data( ), tot_dim, tot_dim );

            map_dGammadGradChi = ( map_dGammadGradChi1 * map_dGradChi1dGradChi ).eval( );

            dGammadGradChin.zero( siot_dim * ( num_configs - 1 ) );
            Eigen::Map< Eigen::Matrix< floatType, tot_dim,      -1, Eigen::RowMajor > > map_dGammadGradChin( dGammadGradChin.value->data( ), tot_dim, ( num_configs - 1 ) * tot_dim );

            map_dGammadGradChin  = ( map_dGammadGradChi1 * map_dGradChi1dGradChin ).eval( );

        }

        void residual::setResidual( ){
            /*!
             * Set the residual w.r.t. the unknown vector
             */

            auto residual = get_setDataStorage_residual( );

            const floatVector *stress = hydra->getStress( );

            TARDIGRADE_ERROR_TOOLS_CATCH( *residual.value = *stress
                                                          - tardigradeVectorTools::appendVectors( { *get_PK2Stress( ), *get_referenceSymmetricMicroStress( ), *get_referenceHigherOrderStress( ) } ) );

        }

        void residual::setJacobian( ){
            /*!
             * Set the Jacobian w.r.t. the unknown vector
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int tot_dim = hydra->getTOTDimension( );

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const unsigned int num_unknowns = hydra->getNumUnknowns( );

            auto jacobian = get_setDataStorage_jacobian( );
            jacobian.zero( *getNumEquations( ) * num_unknowns );

            //Get references to the stress Jacobians. Doing it this way to allow changing the residual to the current configuration in the future.
            const floatVector *dS1dFn       = get_dPK2dFn( );

            const floatVector *dS1dChin     = get_dPK2dChin( );

            const floatVector *dS1dGradChin = get_dPK2dGradChin( );

            const floatVector *dS2dFn       = get_dSIGMAdFn( );

            const floatVector *dS2dChin     = get_dSIGMAdChin( );

            const floatVector *dS2dGradChin = get_dSIGMAdGradChin( );

            const floatVector *dS3dFn       = get_dMdFn( );

            const floatVector *dS3dChin     = get_dMdChin( );

            const floatVector *dS3dGradChin = get_dMdGradChin( );

            std::vector< std::vector< const floatVector * > > stressReferences = { { dS1dFn, dS1dChin, dS1dGradChin },
                                                                                   { dS2dFn, dS2dChin, dS2dGradChin },
                                                                                   { dS3dFn, dS3dChin, dS3dGradChin } };

            const std::array< unsigned int, 3 > dims = { sot_dim, sot_dim, tot_dim };

            unsigned int row = 0;
            for ( auto S = stressReferences.begin( ); S != stressReferences.end( ); S++ ){

                // Loop through the stress Jacobians

                for ( unsigned int i = 0; i < dims[ ( unsigned int )( S - stressReferences.begin( ) ) ]; i++ ){

                    unsigned int col = 0;

                    // Jacobians w.r.t. the stress
                    ( *jacobian.value )[ num_unknowns * row + row ] += 1;

                    // Jacobians w.r.t. the sub configurations
                    col = *hydra->getConfigurationUnknownCount( );

                    for ( auto Sn = S->begin( ); Sn != S->end( ); Sn++ ){

                        for ( unsigned int j = 0; j < ( num_configs - 1 ) * dims[ ( unsigned int )( Sn - S->begin( ) ) ]; j++ ){

                            ( *jacobian.value )[ num_unknowns * row + col ] -= ( **Sn )[ ( num_configs - 1 ) * dims[ ( unsigned int ) ( Sn - S->begin( ) ) ] * i + j ];

                            col++;

                        }

                    }

                    row++;

                }

            }

        }

        void residual::setdRdD( ){
            /*!
             * Set the Jacobian w.r.t. the deformation
             */

            const unsigned int sot_dim = hydra->getSOTDimension( );

            const unsigned int tot_dim = hydra->getTOTDimension( );

            const unsigned int num_equations = *getNumEquations( );

            const unsigned int num_configurationUnknowns = *hydra->getConfigurationUnknownCount( );

            auto dRdD = get_setDataStorage_dRdD( );

            dRdD.zero( num_equations * num_configurationUnknowns );

            //Get references to the stress Jacobians. Doing it this way to allow changing the residual to the current configuration in the future.
            const floatVector *dS1dF       = get_dPK2dF( );

            const floatVector *dS1dChi     = get_dPK2dChi( );

            const floatVector *dS1dGradChi = get_dPK2dGradChi( );

            const floatVector *dS2dF       = get_dSIGMAdF( );

            const floatVector *dS2dChi     = get_dSIGMAdChi( );

            const floatVector *dS2dGradChi = get_dSIGMAdGradChi( );

            const floatVector *dS3dF       = get_dMdF( );

            const floatVector *dS3dChi     = get_dMdChi( );

            const floatVector *dS3dGradChi = get_dMdGradChi( );

            std::vector< std::vector< const floatVector * > > stressReferences = { { dS1dF, dS1dChi, dS1dGradChi },
                                                                                   { dS2dF, dS2dChi, dS2dGradChi },
                                                                                   { dS3dF, dS3dChi, dS3dGradChi } };

            const std::array< unsigned int, 3 > dims = { sot_dim, sot_dim, tot_dim };

            unsigned int row = 0;
            for ( auto S = stressReferences.begin( ); S != stressReferences.end( ); S++ ){

                // Loop through the stress Jacobians

                for ( unsigned int i = 0; i < dims[ ( unsigned int )( S - stressReferences.begin( ) ) ]; i++ ){

                    unsigned int col = 0;

                    // Jacobians w.r.t. the deformation

                    for ( auto Sn = S->begin( ); Sn != S->end( ); Sn++ ){

                        for ( unsigned int j = 0; j < dims[ ( unsigned int )( Sn - S->begin( ) ) ]; j++ ){

                            ( *dRdD.value )[ num_configurationUnknowns * row + col ] -= ( **Sn )[ dims[ ( unsigned int )( Sn - S->begin( ) ) ] * i + j ];

                            col++;

                        }

                    }

                    row++;

                }

            }

        }

        void residual::setStress( ){
            /*!
             * Set the stresses
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            constexpr unsigned int tot_dim = sot_dim * dim;

            auto stress = get_setDataStorage_stress( );

            stress.zero( 2 * sot_dim + tot_dim );

            std::copy( std::begin( *get_PK2Stress( ) ), std::end( *get_PK2Stress( ) ), std::begin( *stress.value ) );

            std::copy( std::begin( *get_referenceSymmetricMicroStress( ) ), std::end( *get_referenceSymmetricMicroStress( ) ), std::begin( *stress.value ) + sot_dim );

            std::copy( std::begin( *get_referenceHigherOrderStress( ) ), std::end( *get_referenceHigherOrderStress( ) ), std::begin( *stress.value ) + 2 * sot_dim );

        }

        void residual::setPreviousStress( ){
            /*!
             * Set the previous stresses
             */

            constexpr unsigned int dim = 3;

            constexpr unsigned int sot_dim = dim * dim;

            constexpr unsigned int tot_dim = sot_dim * dim;

            auto previousStress = get_setDataStorage_previousStress( );

            previousStress.zero( 2 * sot_dim + tot_dim );

            std::copy( std::begin( *get_previousPK2Stress( ) ), std::end( *get_previousPK2Stress( ) ), std::begin( *previousStress.value ) );

            std::copy( std::begin( *get_previousReferenceSymmetricMicroStress( ) ), std::end( *get_previousReferenceSymmetricMicroStress( ) ), std::begin( *previousStress.value ) + sot_dim );

            std::copy( std::begin( *get_previousReferenceHigherOrderStress( ) ), std::end( *get_previousReferenceHigherOrderStress( ) ), std::begin( *previousStress.value ) + 2 * sot_dim );

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            auto dRdT = get_setDataStorage_dRdT( );
            dRdT.zero( *getNumEquations( ) );

        }

    }

}
