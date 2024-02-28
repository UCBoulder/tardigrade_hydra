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
             * :param variableVector &CauchyStress: The Cauchy stress.
             * :param variableVector &microStress: The symmetric micro-stress.
             * :param variableVector &higherOrderStress: The higher-order stress.
             */
    
            variableVector PK2Stress, referenceMicroStress, referenceHigherOrderStress;
            errorOut error = linearElasticityReference( deformationGradient, microDeformation, gradientMicroDeformation,
                                                        A, B, C, D,
                                                        PK2Stress, referenceMicroStress, referenceHigherOrderStress );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticity",
                                                 "Error in the computation of the stresses in the reference configuration" );
                result->addNext( error );
                return result;
            }
    
            error = mapStressMeasuresToCurrent( deformationGradient, microDeformation, PK2Stress,
                                                referenceMicroStress, referenceHigherOrderStress,
                                                cauchyStress, microStress, higherOrderStress );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticity",
                                                 "Error in mapping the reference stresses to the current configuration" );
                result->addNext( error );
                return result;
            }
    
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
             * \param &CauchyStress: The Cauchy stress.
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

            const unsigned int dim = 3;
            const unsigned int sot_dim = dim * dim;
            const unsigned int tot_dim = sot_dim * dim;
 
            variableVector PK2Stress, referenceMicroStress, referenceHigherOrderStress;
    
            variableVector dPK2StressdF, dPK2StressdChi, dPK2StressdGradChi;
            variableVector dReferenceMicroStressdF, dReferenceMicroStressdChi, dReferenceMicroStressdGradChi;
            variableVector dReferenceHigherOrderStressdF, dReferenceHigherOrderStressdGradChi;
    
            errorOut error = linearElasticityReference( deformationGradient, microDeformation, gradientMicroDeformation,
                                                        A, B, C, D,
                                                        PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                        dPK2StressdF, dPK2StressdChi, dPK2StressdGradChi,
                                                        dReferenceMicroStressdF, dReferenceMicroStressdChi, dReferenceMicroStressdGradChi,
                                                        dReferenceHigherOrderStressdF, dReferenceHigherOrderStressdGradChi );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticity",
                                                 "Error in the computation of the stresses in the reference configuration" );
                result->addNext( error );
                return result;
            }
    
            variableVector dCauchyStressdPK2Stress, dMicroStressdReferenceMicroStress, dHigherOrderStressdReferenceHigherOrderStress;
    
            error = mapStressMeasuresToCurrent( deformationGradient, microDeformation, PK2Stress,
                                                referenceMicroStress, referenceHigherOrderStress,
                                                cauchyStress, microStress, higherOrderStress,
                                                dCauchyStressdF, dCauchyStressdPK2Stress,
                                                dMicroStressdF, dMicroStressdReferenceMicroStress,
                                                dHigherOrderStressdF, dHigherOrderStressdChi,
                                                dHigherOrderStressdReferenceHigherOrderStress );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticity",
                                                 "Error in mapping the reference stresses to the current configuration" );
                result->addNext( error );
                return result;
            }
    
            //Assemble the jacobians of the Cauchy stress
            dCauchyStressdF += tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, dPK2StressdF, sot_dim, sot_dim, sot_dim, sot_dim );
            dCauchyStressdChi = tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, dPK2StressdChi, sot_dim, sot_dim, sot_dim, sot_dim );
            dCauchyStressdGradChi = tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, dPK2StressdGradChi, sot_dim, sot_dim, sot_dim, tot_dim );
    
            //Assemble the jacobians of the symmetric micro-stress
            dMicroStressdF += tardigradeVectorTools::matrixMultiply( dMicroStressdReferenceMicroStress, dReferenceMicroStressdF, sot_dim, sot_dim, sot_dim, sot_dim );
            dMicroStressdChi = tardigradeVectorTools::matrixMultiply( dMicroStressdReferenceMicroStress, dReferenceMicroStressdChi, sot_dim, sot_dim, sot_dim, sot_dim );
            dMicroStressdGradChi = tardigradeVectorTools::matrixMultiply( dMicroStressdReferenceMicroStress, dReferenceMicroStressdGradChi, sot_dim, sot_dim, sot_dim, tot_dim );
    
            //Assemble the jacobians of the higher-order stress
            dHigherOrderStressdF += tardigradeVectorTools::matrixMultiply( dHigherOrderStressdReferenceHigherOrderStress,
                                                      dReferenceHigherOrderStressdF, tot_dim, tot_dim, tot_dim, sot_dim );
            dHigherOrderStressdGradChi = tardigradeVectorTools::matrixMultiply( dHigherOrderStressdReferenceHigherOrderStress,
                                                          dReferenceHigherOrderStressdGradChi, tot_dim, tot_dim, tot_dim, tot_dim );
    
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
            errorOut error = computeDeformationMeasures( deformationGradient, microDeformation, gradientMicroDeformation,
                                                         RCG, Psi, Gamma );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReference",
                                                 "Error in the computation of the deformation measures" );
                result->addNext( error );
                return result;
            }
    
            error = linearElasticityReferenceDerivedMeasures( RCG, Psi, Gamma, A, B, C, D,
                                                              PK2Stress, referenceMicroStress,
                                                              referenceHigherOrderStress );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReference",
                                                 "Error in the computation of the reference stresses" );
                result->addNext( error );
                return result;
            }
    
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
             * \param variableVector &referenceHigherOrderStress: The higher-order stress in the 
             * \   reference configuration.
             * \param &dPK2StressdF: The Jacobian of the PK2 stress w.r.t. the deformation gradient.
             * \param &dPK2StressdChi: The Jacobian of the PK2 stress w.r.t. the micro deformation.
             * \param &dPK2StressdGradChi: The Jacobian of the PK2 stress w.r.t. the gradient of the micro deformation.
             * \param &dReferenceMicroStressdF: The Jacobian of the Micro stress w.r.t. the deformation gradient.
             * \param &dReferenceMicroStressdChi: The Jacobian of the Micro stress w.r.t. the micro deformation.
             * \param &dReferenceStressdGradChi: The Jacobian of the Micro stress w.r.t. the gradient of the micro deformation.
             * \param &dMdF: The Jacobian of the higher order stress w.r.t. the deformation gradient.
             * \param &dMdGradChi: The Jacobian of the higher order stress w.r.t. the gradient of the micro deformation.
             */
    
            //Assume 3d
            const unsigned int dim = 3;
            const unsigned int sot_dim = dim * dim;
            const unsigned int tot_dim = sot_dim * dim;
    
            constantVector eye( sot_dim );
            tardigradeVectorTools::eye( eye );
    
            //Compute the required deformation measures
            variableVector RCG, Psi, Gamma;
            variableVector dRCGdF, dPsidF, dPsidChi, dGammadF, dGammadGradChi;
            errorOut error = computeDeformationMeasures( deformationGradient, microDeformation, gradientMicroDeformation,
                                                         RCG, Psi, Gamma, dRCGdF, dPsidF, dPsidChi, dGammadF, dGammadGradChi );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReference (jacobian)",
                                                 "Error in the computation of the deformation measures" );
                result->addNext( error );
                return result;
            }
    
            variableVector dPK2StressdRCG, dPK2StressdPsi, dPK2StressdGamma;
            variableVector dReferenceMicroStressdRCG, dReferenceMicroStressdPsi, dReferenceMicroStressdGamma;
            variableVector dMdGamma;
    
            error = linearElasticityReferenceDerivedMeasures( RCG, Psi, Gamma, A, B, C, D,
                                                              PK2Stress, referenceMicroStress, referenceHigherOrderStress,
                                                              dPK2StressdRCG, dPK2StressdPsi, dPK2StressdGamma,
                                                              dReferenceMicroStressdRCG, dReferenceMicroStressdPsi,
                                                              dReferenceMicroStressdGamma, dMdGamma );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReference (jacobian)",
                                                 "Error in the computation of the deformation measures" );
                result->addNext( error );
                return result;
            }
    
            dPK2StressdF = tardigradeVectorTools::matrixMultiply( dPK2StressdRCG, dRCGdF, sot_dim, sot_dim, sot_dim, sot_dim )
                         + tardigradeVectorTools::matrixMultiply( dPK2StressdPsi, dPsidF, sot_dim, sot_dim, sot_dim, sot_dim )
                         + tardigradeVectorTools::matrixMultiply( dPK2StressdGamma, dGammadF, sot_dim, tot_dim, tot_dim, sot_dim );
    
            dPK2StressdChi = tardigradeVectorTools::matrixMultiply( dPK2StressdPsi, dPsidChi, sot_dim, sot_dim, sot_dim, sot_dim );
    
            dPK2StressdGradChi = tardigradeVectorTools::matrixMultiply( dPK2StressdGamma, dGammadGradChi, sot_dim, tot_dim, tot_dim, tot_dim );
    
            dReferenceMicroStressdF = tardigradeVectorTools::matrixMultiply( dReferenceMicroStressdRCG, dRCGdF, sot_dim, sot_dim, sot_dim, sot_dim )
                                    + tardigradeVectorTools::matrixMultiply( dReferenceMicroStressdPsi, dPsidF, sot_dim, sot_dim, sot_dim, sot_dim )
                                    + tardigradeVectorTools::matrixMultiply( dReferenceMicroStressdGamma, dGammadF, sot_dim, tot_dim, tot_dim, sot_dim );
    
            dReferenceMicroStressdChi = tardigradeVectorTools::matrixMultiply( dReferenceMicroStressdPsi, dPsidChi, sot_dim, sot_dim, sot_dim, sot_dim );
    
            dReferenceMicroStressdGradChi = tardigradeVectorTools::matrixMultiply( dReferenceMicroStressdGamma, dGammadGradChi, sot_dim, tot_dim, tot_dim, tot_dim );
    
            dMdF = tardigradeVectorTools::matrixMultiply( dMdGamma, dGammadF, tot_dim, tot_dim, tot_dim, sot_dim );
            dMdGradChi = tardigradeVectorTools::matrixMultiply( dMdGamma, dGammadGradChi, tot_dim, tot_dim, tot_dim, tot_dim );
    
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
            const unsigned int dim = 3;
    
            constantVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );
    
            variableVector invRCG = tardigradeVectorTools::inverse( rightCauchyGreenDeformation, dim, dim );
    
            //Compute the strain measures
            variableVector greenLagrangeStrain = 0.5 * ( rightCauchyGreenDeformation - eye );
            variableVector microStrain   = Psi - eye;
    
            //Compute the higher order stress
            errorOut error = computeReferenceHigherOrderStress( Gamma, C, referenceHigherOrderStress );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMeasures",
                                                 "Error in computation of higher-order stress" );
                result->addNext( error );
                return result;
            }
    
            //Compute the first common term for the PK2 and symmetric micro-stress
            variableVector term1;
            error = computeLinearElasticTerm1( greenLagrangeStrain, microStrain, A, D, term1 );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMeasures",
                                                 "Error in computation of term 1" );
                result->addNext( error );
                return result;
            }
    
            //Compute the second common term for the PK2 and symmetric micro-stress
            variableVector invRCGPsi;
            error = computeInvRCGPsi( invRCG, Psi, invRCGPsi );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMeasures",
                                                 "Error in computation of invRCG Psi product" );
                result->addNext( error );
                return result;
            }
    
            variableVector term2;
            error = computeLinearElasticTerm2( greenLagrangeStrain, microStrain, invRCGPsi, B, D, term2 );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMesures",
                                                 "Error in computation of term 2" );
                result->addNext( error );
                return result;
            }
    
            //Compute the third common term for the PK2 and symmetric micro-stress
            variableVector invRCGGamma;
            error = computeInvRCGGamma( invRCG, Gamma, invRCGGamma );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMeasures",
                                                 "Error in computation of invRCG Gamma product" );
                result->addNext(error);
                return result;
            }
    
            variableVector term3;
            error = computeLinearElasticTerm3( invRCGGamma, referenceHigherOrderStress, term3 );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMeasures",
                                                 "Error in computation of term 3" );
                result->addNext( error );
                return result;
            }
    
            //Construct the PK2 and reference symmetric stresses
            PK2Stress            = term1 + term2 + term3;
    
            variableVector symmTerm2Term3;
            error = tardigradeConstitutiveTools::computeSymmetricPart( term2 + term3, symmTerm2Term3 );
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
             * \param &dReferenceMicroStrssdGamma: The Jacobian of the reference micro stress w.r.t. the 
             *     higher order deformation measure.
             * \param &dMdGamma: The Jacobian of the reference higher order stress w.r.t. 
             *     the higher order deformation measure.
             */

            //Assume 3d
            const unsigned int dim = 3;
            const unsigned int sot_dim = dim * dim;
            const unsigned int tot_dim = sot_dim * dim;   
 
            constantVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );
    
            variableVector invRCG = tardigradeVectorTools::inverse( rightCauchyGreenDeformation, dim, dim );
    
            //Compute the strain measures
            variableVector greenLagrangeStrain = 0.5 * ( rightCauchyGreenDeformation - eye );
            variableVector microStrain   = Psi - eye;
    
            //Compute the higher order stress
            errorOut error = computeReferenceHigherOrderStress( Gamma, C, referenceHigherOrderStress, dMdGamma );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityRefereneDerivedMetrics (jacobian)",
                                                 "Error in computation of higher-order stress" );
                result->addNext( error );
                return result;
            }
    
            //Compute the first common term for the PK2 and symmetric micro-stress
            variableVector term1;
    
            variableVector dTerm1dRCG, dTerm1dPsi;
            error = computeLinearElasticTerm1( greenLagrangeStrain, microStrain, A, D, term1,
                                               dTerm1dRCG, dTerm1dPsi );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMetrics (jacobian)",
                                                 "Error in computation of term 1" );
                result->addNext( error );
                return result;
            }
    
            //Assemble term1 jacobians w.r.t. F and Chi
            dTerm1dRCG *= 0.5;
    
            //Compute the second common term for the PK2 and symmetric micro-stress
            variableVector invRCGPsi;
            variableVector dInvRCGPsidRCG, dInvRCGPsidPsi;
    
            error = computeInvRCGPsi( invRCG, Psi, invRCGPsi, dInvRCGPsidRCG, dInvRCGPsidPsi );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMetrics (jacobian)",
                                                 "Error in computation of invRCG Psi product" );
                result->addNext( error );
                return result;
            }
    
            variableVector term2;
            variableVector dTerm2dRCG, dTerm2dPsi, dTerm2dInvRCGPsi;
            error = computeLinearElasticTerm2( greenLagrangeStrain, microStrain, invRCGPsi, B, D, term2,
                                               dTerm2dRCG, dTerm2dPsi, dTerm2dInvRCGPsi );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMetrics (jacobian)",
                                                 "Error in computation of term 2" );
                result->addNext( error );
                return result;
            }

            dTerm2dRCG *= 0.5;
            dTerm2dRCG += tardigradeVectorTools::matrixMultiply( dTerm2dInvRCGPsi, dInvRCGPsidRCG, sot_dim, sot_dim, sot_dim, sot_dim );
    
            dTerm2dPsi += tardigradeVectorTools::matrixMultiply( dTerm2dInvRCGPsi, dInvRCGPsidPsi, sot_dim, sot_dim, sot_dim, sot_dim );
    
            //Compute the third common term for the PK2 and symmetric micro-stress
            variableVector invRCGGamma;
            variableVector dInvRCGGammadRCG, dInvRCGGammadGamma;
    
            error = computeInvRCGGamma( invRCG, Gamma, invRCGGamma, dInvRCGGammadRCG, dInvRCGGammadGamma );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMetrics (jacobian)",
                                                 "Error in computation of invRCG Gamma product" );
                result->addNext( error );
                return result;
            }
    
            variableVector term3;
            variableVector dTerm3dInvRCGGamma, dTerm3dM;
            error = computeLinearElasticTerm3( invRCGGamma, referenceHigherOrderStress, term3, dTerm3dInvRCGGamma, dTerm3dM );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMetrics (jacobian)",
                                                 "Error in computation of term 3" );
                result->addNext( error );
                return result;
            }
    
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
            error = tardigradeConstitutiveTools::computeSymmetricPart( term2 + term3, symmTerm2Term3, dSymmTerm2Term3dTerm2Term3 );
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
             * \param const variableVector &deformationGradient: The deformation gradient between the 
             *     reference configuration and the current configuration.
             * \param const variableVector &microDeformation: The micro-deformation map between the 
             *     reference configuration and the current configuration.
             * \param const variableVector &PK2Stress: The Second Piola-Kirchoff stress.
             * \param const variableVector &referenceMicroStress: The symmetric micro-stress in the 
             *     reference configuration.
             * \param const variableVector &referenceHigherOrderStress: The higher order stress in 
             *     the reference configuration.
             * \param variableVector &cauchyStress: The Cauchy stress (PK2 stress in the current configuration).
             * \param variableVector &microStress: The symmetric micro-stress in the current configuration.
             * \param variableVector &higherOrderStress: The higher order stress in the current configuration.
             */
    
            //Map the PK2 stress to the Cauchy stress
            errorOut error = tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient, cauchyStress );
    
            if ( error ){
                errorOut result = new errorNode( "mapStressMeasuresToCurrent",
                                                 "Error in the map of the PK2 stress to the Cauchy stress" );
                result->addNext( error );
                return result;
            }
    
            //Map the symmetric micro stress to the current configuration
            error = tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, microStress );
    
            if ( error ){
                errorOut result = new errorNode( "mapStressMeasuresToCurrent",
                                                 "Error in the map of the micro-stress to the current configuation" );
                result->addNext( error );
                return result;
            }
    
            //Map the higher order stress to the current configuration
            error = tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient,
                                                                     microDeformation, higherOrderStress );
    
            if ( error ){
                errorOut result = new errorNode( "mapStressMeasuresToCurrent",
                                                 "Error in the map of the higher-order stress to the current configuation" );
                result->addNext( error );
                return result;
            }
    
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
            errorOut error = tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient, cauchyStress,
                                                                      dCauchyStressdPK2Stress, dCauchyStressdF );
    
            if ( error ){
                errorOut result = new errorNode( "mapStressMeasuresToCurrent",
                                                 "Error in the map of the PK2 stress to the Cauchy stress" );
                result->addNext( error );
                return result;
            }
    
            //Map the symmetric micro stress to the current configuration
            error = tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, microStress,
                                                                        dMicroStressdReferenceMicroStress, dMicroStressdF );
    
            if ( error ){
                errorOut result = new errorNode( "mapStressMeasuresToCurrent",
                                                 "Error in the map of the micro-stress to the current configuation" );
                result->addNext( error );
                return result;
            }
    
            //Map the higher order stress to the current configuration
            error = tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient,
                                                                     microDeformation, higherOrderStress,
                                                                     dHigherOrderStressdReferenceHigherOrderStress,
                                                                     dHigherOrderStressdF,
                                                                     dHigherOrderStressdChi );
    
            if ( error ){
                errorOut result = new errorNode( "mapStressMeasuresToCurrent",
                                                 "Error in the map of the higher-order stress to the current configuation" );
                result->addNext( error );
                return result;
            }
    
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
    
            errorOut error = tardigradeConstitutiveTools::computeRightCauchyGreen( deformationGradient, rightCauchyGreen );
    
            if ( error ){
                errorOut result = new errorNode( "computeDeformationMeasures",
                                                 "Error in the computation of the right Cauchy-Green Deformation measure" );
                result->addNext( error );
                return result;
            }
    
            error = tardigradeMicromorphicTools::computePsi( deformationGradient, microDeformation, Psi );
    
            if ( error ){
                errorOut result = new errorNode( "computeDeformationMeasures",
                                                 "Error in the computation of Psi" );
                result->addNext( error );
                return result;
            }
    
            error = tardigradeMicromorphicTools::computeGamma( deformationGradient, gradientMicroDeformation, Gamma );
    
            if ( error ){
                errorOut result = new errorNode( "computeDeformationMeasures",
                                                 "Error in the computation of Gamma" );
                result->addNext( error );
                return result;
            }
    
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
    
            errorOut error = tardigradeConstitutiveTools::computeRightCauchyGreen( deformationGradient, rightCauchyGreen, dCdF );
    
            if ( error ){
                errorOut result = new errorNode( "computeDeformationMeasures (jacobian)",
                                                 "Error in the computation of the right Cauchy-Green Deformation measure" );
                result->addNext( error );
                return result;
            }
    
            error = tardigradeMicromorphicTools::computePsi( deformationGradient, microDeformation, Psi, dPsidF, dPsidChi );
    
            if ( error ){
                errorOut result = new errorNode( "computeDeformationMeasures (jacobian)",
                                                 "Error in the computation of Psi" );
                result->addNext( error );
                return result;
            }
    
            error = tardigradeMicromorphicTools::computeGamma( deformationGradient, gradientMicroDeformation, Gamma, dGammadF, dGammadGradChi );
    
            if ( error ){
                errorOut result = new errorNode( "computeDeformationMeasures (jacobian)",
                                                 "Error in the computation of Gamma" );
                result->addNext( error );
                return result;
            }
    
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
            const unsigned int dim = 3;
    
            if ( greenLagrangeStrain.size() != dim * dim ){
                return new errorNode( "computeLinearElasticTerm1",
                                      "The green lagrange strain must have a length of 9" );
            }
    
            if ( microStrain.size() != dim * dim ){
                return new errorNode( "computeLinearElasticTerm1",
                                      "The micro-strain must have a length of 9" );
            }
    
            if ( A.size() != dim * dim * dim * dim ){
                return new errorNode( "computeLinearElasticTerm1",
                                      "A must have a size of 3**4" );
            }
    
            if ( D.size() != dim * dim * dim * dim ){
                return new errorNode( "computeLinearElasticTerm1",
                                      "D must have a size of 3**4" );
            }
    
            //Compute the first common term for the PK2 and symmetric micro-stress
            term1 = variableVector( dim * dim, 0 );
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int L = 0; L < dim; L++ ){
                            term1[ dim * I + J ] += A[ dim * dim * dim * I + dim * dim * J + dim * K + L ] * greenLagrangeStrain[ dim * K + L ]
                                                  + D[ dim * dim * dim * I + dim * dim * J + dim * K + L ] * microStrain[ dim * K + L ];
                        }
                    }
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
             */
    
            //Assume 3D
            const unsigned int dim = 3;
            const unsigned int sot_dim = dim * dim;
    
            errorOut error = computeLinearElasticTerm1( greenLagrangeStrain, microStrain, A, D, term1 );
    
            if ( error ){
                errorOut result = new errorNode( "computeLinearElasticTerm1 (jacobian)",
                                                 "Error in computation of linear elastic term1" );
                result->addNext( error );
                return result;
            }
    
            //Compute the first common term for the PK2 and symmetric micro-stress
            dTerm1dGreenLagrangeStrain = variableVector( sot_dim * sot_dim, 0 );
            dTerm1dMicroStrain = variableVector( sot_dim * sot_dim, 0 );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int M = 0; M < dim; M++ ){
                        for ( unsigned int N = 0; N < dim; N++ ){
                            dTerm1dGreenLagrangeStrain[ dim * sot_dim * I + sot_dim * J + dim * M + N ] = A[ dim * dim * dim * I + dim * dim * J + dim * M + N ];
                            dTerm1dMicroStrain[ dim * sot_dim * I + sot_dim * J + dim * M + N ] = D[ dim * dim * dim * I + dim * dim * J + dim * M + N ];
                        }
                    }
                }
            }
    
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
            const unsigned int dim = 3;
            const unsigned int sot_dim = dim * dim;
    
            if ( greenLagrangeStrain.size() != sot_dim ){
                return new errorNode( "computeLinearElasticTerm2",
                                      "The green lagrange strain must have a length of 9" );
            }
    
            if ( microStrain.size() != sot_dim ){
                return new errorNode( "computeLinearElasticTerm2",
                                      "The micro-strain must have a length of 9" );
            }
    
            if ( invCPsi.size() != sot_dim ){
                return new errorNode( "computeLinearElasticTerm2",
                                      "invCPsi must have a size of 9" );
            }
    
            if ( B.size() != sot_dim * sot_dim ){
                return new errorNode( "computeLinearElasticTerm2",
                                      "B must have a size of 3**4" );
            }
    
            if ( D.size() != sot_dim * sot_dim ){
                return new errorNode( "computeLinearElasticTerm2",
                                      "D must have a size of 3**4" );
            }
    
            term2 = variableVector( greenLagrangeStrain.size(), 0 );
    
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
            const unsigned int dim = 3;
            const unsigned int sot_dim = dim * dim;
    
            errorOut error = computeLinearElasticTerm2( greenLagrangeStrain, microStrain, invCPsi, B, D, term2 );
    
            if ( error ){
                errorOut result = new errorNode("computeLinearElasticTerm2 (jacobian)",
                                                "Error in computation of term 2" );
                result->addNext( error );
                return result;
            }
    
            //Compute the Jacobians
            constantVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );
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
                                for ( unsigned int L = 0; L < dim; L++ ){
                                    dTerm2dInvCPsi[ dim * sot_dim * I + sot_dim * J + dim * M + N ] += ( B[ dim * dim * dim * I + dim * dim * N + dim * K + L ] * microStrain[ dim * K + L ] + greenLagrangeStrain[ dim * K + L ] * D[ dim * dim * dim * K + dim * dim * L + dim * I + N ] ) * eye[ dim * J + M ];
                                }
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
            const unsigned int dim = 3;
            const unsigned int sot_dim = dim * dim;
            const unsigned int tot_dim = sot_dim * dim;
    
            if ( Gamma.size() != tot_dim ){
                return new errorNode( "computeReferenceHigherOrderStress",
                                      "Gamma must have a length of 27" );
            }
    
            if ( C.size() != tot_dim * tot_dim ){
                return new errorNode( "computeReferenceHigherOrderStress",
                                      "The C stiffness tensor must have a length of 3**6.\nThe current size is " + std::to_string( C.size( ) ) );
            }
    
            referenceHigherOrderStress = variableVector( tot_dim, 0 );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int L = 0; L < dim; L++ ){
                            for ( unsigned int M = 0; M < dim; M++ ){
                                for ( unsigned int N = 0; N < dim; N++ ){
                                    referenceHigherOrderStress[ dim * dim * I + dim * J + K ] += C[ dim * dim * dim * dim * dim * J + dim * dim * dim * dim * K + dim * dim * dim * I + dim * dim * L + dim * M + N ] * Gamma[ dim * dim * L + dim * M + N ];
                                }
                            }
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
             * \f$M_{IJK} = C_{JKILMN} Gamma_{LMN}\f$
             *
             * Also compute the Jacobian
             * \f$\frac{ \partial M_{IJK} }{\partial \Gamma_{OPQ} } = C_{JKIOPQ}\f$
             *
             * \param &Gamma: The micro-gradient deformation measure.
             * \param &C: The C stiffness tensor.
             * \param &referenceHigherOrderStress: The higher order stress in the reference 
             *     configuration.
             */
    
            //Assume 3D
            const unsigned int dim = 3;
            const unsigned int sot_dim = dim * dim;
            const unsigned int tot_dim = sot_dim * dim;
    
            errorOut error = computeReferenceHigherOrderStress( Gamma, C, referenceHigherOrderStress );
    
            if ( error ){
                errorOut result = new errorNode( "computeReferenceHigherOrderStress (jacobian)",
                                                 "Error in computation of higher order stress" );
                result->addNext( error );
                return result;
            }
    
            //Assemble the Jacobian
            dReferenceHigherOrderStressdGamma = variableVector( tot_dim * tot_dim, 0 );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int O = 0; O < dim; O++ ){
                            for ( unsigned int P = 0; P < dim; P++ ){
                                for ( unsigned int Q = 0; Q < dim; Q++ ){
                                    dReferenceHigherOrderStressdGamma[ dim * dim * tot_dim * I + dim * tot_dim * J + tot_dim * K + dim * dim * O + dim * P + Q ] +=
                                        C[ dim * dim * dim * dim * dim * J + dim * dim * dim * dim * K + dim * dim * dim * I + dim * dim * O + dim * P + Q ];
                                }
                            }
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
            const unsigned int dim = 3;
            const unsigned int sot_dim = dim * dim;
    
            if ( invCGamma.size() != dim * dim * dim ){
                return new errorNode( "computeLinearElasticTerm3",
                                      "invCGamma must have a size of 27" );
            }
    
            if ( referenceHigherOrderStress.size() != dim * dim * dim ){
                return new errorNode( "computeLinearElasticTerm3",
                                      "The referenceHigherOrder stress must have a size of 27" );
            }
    
            term3 = variableVector( sot_dim, 0 );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int Q = 0; Q < dim; Q++ ){
                        for ( unsigned int R = 0; R < dim; R++ ){
                            term3[ dim * I + J ] += referenceHigherOrderStress[ dim * dim * I + dim * Q + R ] * invCGamma[ dim * dim * J + dim * Q + R ];
                        }
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
             */
    
            //Assume 3D
            const unsigned int dim = 3;
            const unsigned int sot_dim = dim * dim;
            const unsigned int tot_dim = sot_dim * dim;
    
            errorOut error = computeLinearElasticTerm3( invCGamma, referenceHigherOrderStress, term3 );
    
            if ( error ){
                errorOut result = new errorNode( "computeLinearElasticTerm3 (jacobian)",
                                                 "Error in computation of term 3" );
                result->addNext( error );
                return result;
            }
    
            constantVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );
    
            dTerm3dInvCGamma = variableVector( sot_dim * tot_dim, 0 );
            dTerm3dReferenceHigherOrderStress = variableVector( sot_dim * tot_dim, 0 );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int T = 0; T < dim; T++ ){
                        for ( unsigned int U = 0; U < dim; U++ ){
                            for ( unsigned int V = 0; V < dim; V++ ){
                                dTerm3dInvCGamma[ dim * tot_dim * I + tot_dim * J + dim * dim * T + dim * U + V ] = referenceHigherOrderStress[ dim * dim * I + dim * U + V ] * eye[ dim * J + T ];
                                dTerm3dReferenceHigherOrderStress[ dim * tot_dim * I + tot_dim * J + dim * dim * T + dim * U + V ] = eye[ dim * I + T ] * invCGamma[ dim * dim * J + dim * U + V ];
                            }
                        }
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
            const unsigned int dim = 3;
    
            if ( invRCG.size() != dim * dim ){
                return new errorNode( "computeInvRCGGamma", "invRCG has an improper dimension" );
            }
    
            if ( Psi.size() != dim * dim ){
                return new errorNode( "computeInvRCGGamma", "Psi has an improper dimension" );
            }
    
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
            const unsigned int dim = 3;
            const unsigned int sot_dim = dim * dim;
    
            errorOut error = computeInvRCGPsi( invRCG, Psi, invRCGPsi );
    
            if ( error ){
                errorOut result = new errorNode( "computeInvRCGPsi (jacobian)", "Error in computation of invRCG Psi product" );
                result->addNext( error );
                return result;
            }
    
            //Construct the jacobians
            variableVector eye( sot_dim );
            tardigradeVectorTools::eye( eye );
    
            dInvRCGPsidRCG = variableVector( sot_dim * sot_dim, 0 );
            dInvRCGPsidPsi = variableVector( sot_dim * sot_dim, 0 );
    
            for ( unsigned int I = 0; I < 3; I++ ){
                for ( unsigned int J = 0; J < 3; J++ ){
                    for ( unsigned int K = 0; K < 3; K++ ){
                        for ( unsigned int L = 0; L < 3; L++ ){
                            dInvRCGPsidRCG[ dim * sot_dim * I + sot_dim * J + dim * K + L ] = -invRCG[ dim * I + K ] * invRCGPsi[ dim * L + J ];
                            dInvRCGPsidPsi[ dim * sot_dim * I + sot_dim * J + dim * K + L ] = invRCG[ dim * I + K ] * eye[ dim * J + L ];
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
            const unsigned int dim = 3;
            const unsigned int sot_dim = dim * dim;
            const unsigned int tot_dim = sot_dim * dim;
    
            if ( invRCG.size() != sot_dim ){
                return new errorNode( "computeInvRCGGamma", "invRCG has an improper dimension" );
            }
    
            if ( Gamma.size() != tot_dim ){
                return new errorNode( "computeInvRCGGamma", "Gamma has an improper dimension" );
            }
    
            invRCGGamma = variableVector( tot_dim, 0 );
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int Q = 0; Q < dim; Q++ ){
                    for ( unsigned int R = 0; R < dim; R++ ){
                        for ( unsigned int S = 0; S < dim; S++ ){
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
             * \param &invRCG: The inverse of the right Cauchy Green deformation tensor.
             * \param &Gamma: The gradient of the micro-deformation deformation tensor.
             * \param &invRCGGamma: The product.
             */
    
            //Assume 3d
            const unsigned int dim = 3;
            const unsigned int sot_dim = dim * dim;
            const unsigned int tot_dim = sot_dim * dim;
    
            errorOut error = computeInvRCGGamma( invRCG, Gamma, invRCGGamma );
    
            if ( error ){
                errorOut result = new errorNode( "computeInvRCGGamma (jacobian)", "Error in computation of invRCG Gamma product" );
                result->addNext( error );
                return result;
            }
    
            //Assemble jacobians of invCGamma w.r.t. C and Gamma
            variableVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );
    
            dInvRCGGammadRCG = variableVector( tot_dim * sot_dim, 0 );
            dInvRCGGammadGamma = variableVector( tot_dim * tot_dim, 0 );
    
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int Q = 0; Q < dim; Q++ ){
                    for ( unsigned int R = 0; R < dim; R++ ){
                        for ( unsigned int T = 0; T < dim; T++ ){
                            for ( unsigned int U = 0; U < dim; U++ ){
                                dInvRCGGammadRCG[ dim * dim * sot_dim * J + dim * sot_dim * Q + sot_dim * R + dim * T + U ]
                                    = -invRCG[ dim * J + T] * invRCGGamma[ dim * dim * U + dim * Q + R ];
                                for ( unsigned int V = 0; V < dim; V++ ){
                                    dInvRCGGammadGamma[ dim * dim * tot_dim * J + dim * tot_dim * Q + tot_dim * R + dim * dim * T + dim * U + V]
                                        = invRCG[ dim * J + T ] * eye[ dim * Q + U ] * eye[ dim * R + V ];
                                }
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
            const unsigned int dim = 3;
    
            constantVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );
    
            A = parameterVector( dim * dim * dim * dim, 0 );
    
            for ( unsigned int K = 0; K < dim; K++ ){
                for ( unsigned int L = 0; L < dim; L++ ){
                    for ( unsigned int M = 0; M < dim; M++ ){
                        for ( unsigned int N = 0; N < dim; N++ ){
                            A[ dim * dim * dim * K + dim * dim * L + dim * M + N ] = lambda * eye[ dim * K + L ] * eye[ dim * M + N ]
                                                                                   + mu * ( eye[ dim * K + M ] * eye[ dim * L + N ]
                                                                                          + eye[ dim * K + N ] * eye[ dim * L + M ] );
                        }
                    }
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
            const unsigned int dim = 3;
    
            constantVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );
    
            B = parameterVector( dim * dim * dim * dim, 0 );
    
            for ( unsigned int K = 0; K < dim; K++ ){
                for ( unsigned int L = 0; L < dim; L++ ){
                    for ( unsigned int M = 0; M < dim; M++ ){
                        for ( unsigned int N = 0; N < dim; N++ ){
                            B[ dim * dim * dim * K + dim * dim * L + dim * M + N ] = ( eta - tau ) * eye[ dim * K + L ] * eye[ dim * M + N ]
                                                                                   + kappa * eye[ dim * K + M ] * eye[ dim * L + N ]
                                                                                   + nu * eye[ dim * K + N ] * eye[ dim * L + M ]
                                                                                   - sigma * ( eye[ dim * K + M ] * eye[ dim * L + N ]
                                                                                             + eye[ dim * K + N ] * eye[ dim * L + M ] );
                        }
                    }
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
            const unsigned int dim = 3;
    
            if ( taus.size() != 11 ){
                return new errorNode( "formIsotropicC", "11 moduli required to form C" );
            }
    
            constantVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );
    
            C = parameterVector( dim * dim * dim * dim * dim * dim, 0 );
    
            for ( unsigned int K = 0; K < dim; K++ ){
                for ( unsigned int L = 0; L < dim; L++ ){
                    for ( unsigned int M = 0; M < dim; M++ ){
                        for ( unsigned int N = 0; N < dim; N++ ){
                            for ( unsigned int P = 0; P < dim; P++ ){
                                for ( unsigned int Q = 0; Q < dim; Q++ ){
                                    C[ dim * dim * dim * dim * dim * K + dim * dim * dim * dim * L + dim * dim * dim * M
                                     + dim * dim * N + dim * P + Q ] = taus[0] * ( eye[ dim * K + L ] * eye[ dim * M + N ] * eye[ dim * P + Q ]
                                                                                 + eye[ dim * K + Q ] * eye[ dim * L + M ] * eye[ dim * N + P ] )
                                                                     + taus[1] * ( eye[ dim * K + L ] * eye[ dim * M + P ] * eye[ dim * N + Q ]
                                                                                 + eye[ dim * K + M ] * eye[ dim * L + Q ] * eye[ dim * N + P ] )
                                                                     + taus[2] * eye[ dim * K + L ] * eye[ dim * M + Q ] * eye[ dim * N + P]
                                                                     + taus[3] * eye[ dim * K + N ] * eye[ dim * L + M ] * eye[ dim * P + Q]
                                                                     + taus[4] * ( eye[ dim * K + M ] * eye[ dim * L + N ] * eye[ dim * P + Q ]
                                                                                 + eye[ dim * K + P ] * eye[ dim * L + M ] * eye[ dim * N + Q ] )
                                                                     + taus[5] * eye[ dim * K + M ] * eye[ dim * L + P ] * eye[ dim * N + Q ]
                                                                     + taus[6] * eye[ dim * K + N ] * eye[ dim * L + P ] * eye[ dim * M + Q ]
                                                                     + taus[7] * ( eye[ dim * K + P ] * eye[ dim * L + Q ] * eye[ dim * M + N ]
                                                                                 + eye[ dim * K + Q ] * eye[ dim * L + N ] * eye[ dim * M + P ] )
                                                                     + taus[8] * eye[ dim * K + N ] * eye[ dim * L + Q ] * eye[ dim * M + P ]
                                                                     + taus[9] * eye[ dim * K + P ] * eye[ dim * L + N ] * eye[ dim * M + Q ]
                                                                     + taus[10] * eye[ dim * K + Q ] * eye[ dim * L + P ] * eye[ dim * M + N ];
                                }
                            }
                        }
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
            const unsigned int dim = 3;
    
            constantVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );
    
            D = parameterVector( dim * dim * dim * dim, 0 );
            for ( unsigned int K = 0; K < dim; K++ ){
                for ( unsigned int L = 0; L < dim; L++ ){
                    for ( unsigned int M = 0; M < dim; M++ ){
                        for ( unsigned int N = 0; N < dim; N++ ){
                            D[ dim * dim * dim * K + dim * dim * L + dim * M + N ] = tau * eye[ dim * K + L ] * eye[ dim * M + N ]
                                + sigma * ( eye[ dim * K + M ] * eye[ dim * L + N ] + eye[ dim * K + N ] * eye[ dim * L + M ] );
                        }
                    }
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
    
            errorOut error = tardigradeMicromorphicTools::assembleDeformationGradient( displacementGradient, deformationGradient );
    
            if ( error ){
                errorOut result = new errorNode( "assembleFundamentalDeformationMeasures",
                                                 "Error in assembly of the deformation gradient" );
                result->addNext( error );
                return result;
            }
    
            error = tardigradeMicromorphicTools::assembleMicroDeformation( microDisplacement, microDeformation );
    
            if ( error ){
                errorOut result = new errorNode( "assembleFundamentalDeformationMeasures",
                                                 "Error in assembly of the micro deformation" );
                result->addNext( error );
                return result;
            }
    
            error = tardigradeMicromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, gradientMicroDeformation );
    
            if ( error ){
                errorOut result = new errorNode( "assembleFundamentalDeformationMeasures",
                                                 "Error in assembly of the gradient of the micro deformation" );
                result->addNext( error );
                return result;
            }
    
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
             * \param const &grad_u: The macro displacement gradient w.r.t. the reference configuration.
             * \param const &phi: The micro displacement.
             * \param const &grad_phi: The gradient of the micro displacement w.r.t. the reference configuration.
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
    
            errorOut error = tardigradeMicromorphicTools::assembleDeformationGradient( displacementGradient, deformationGradient, dFdGradU );
    
            if ( error ){
                errorOut result = new errorNode( "assembleFundamentalDeformationMeasures (jacobian)",
                                                 "Error in assembly of the deformation gradient" );
                result->addNext( error );
                return result;
            }
    
            error = tardigradeMicromorphicTools::assembleMicroDeformation( microDisplacement, microDeformation, dChidPhi );
    
            if ( error ){
                errorOut result = new errorNode( "assembleFundamentalDeformationMeasures (jacobian)",
                                                 "Error in assembly of the micro deformation" );
                result->addNext( error );
                return result;
            }
    
            error = tardigradeMicromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, gradientMicroDeformation,
                                                                         dGradChidGradPhi );
    
            if ( error ){
                errorOut result = new errorNode( "assembleFundamentalDeformationMeasures (jacobian)",
                                                 "Error in assembly of the gradient of the micro deformation" );
                result->addNext( error );
                return result;
            }
    
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
    
            if ( fparams.size() == 0 ){
                return new errorNode( "extractMaterialParameters",
                                      "The material parameters vector has a length of 0" );
            }
    
            unsigned int start = 0;
            unsigned int span;
    
            std::vector< parameterVector > outputs( 4 );
    
            //Extract the material parameters
            for ( unsigned int i = 0; i < outputs.size(); i++ ){
                span = ( unsigned int )std::floor( fparams[ start ]  + 0.5 ); //Extract the span of the parameter set
    
                if ( fparams.size() < start + 1 + span ){
                    std::string outstr = "fparams is not long enough to contain all of the required parameters:\n";
                    outstr +=            "    filling variable " + std::to_string( i ) + "\n";
                    outstr +=            "    size =          "  + std::to_string( fparams.size() ) + "\n";
                    outstr +=            "    required size = "  + std::to_string( start + 1 + span );
    
                    return new errorNode( "extractMaterialParameters",
                                          outstr.c_str() );
                }
    
                outputs[ i ] = parameterVector( fparams.begin() + start + 1, fparams.begin() + start + 1 + span );
    
                start = start + 1 + span;
            }
    
            //Form the stiffness tensors
            errorOut error;
            if ( outputs[ 0 ].size() == 2 ){
                error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicA( outputs[ 0 ][ 0 ], outputs[ 0 ][ 1 ], Amatrix );
            }
            else{
                std::string outstr = "Unrecognized number of parameters ( " + std::to_string( outputs[ 0 ].size() ) + " ) for the A stiffness tensor";
                return new errorNode( "extractMaterialParameters",
                                      outstr.c_str() );
            }
    
            if ( error ){
                errorOut result = new errorNode( "extractMaterialParameters", "Error in computation of the A stiffness tensor" );
                result->addNext( error );
                return result;
            }
    
            if ( outputs[ 1 ].size() == 5 ){
                error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicB( outputs[ 1 ][ 0 ], outputs[ 1 ][ 1 ], outputs[ 1 ][ 2 ],
                                                                                       outputs[ 1 ][ 3 ], outputs[ 1 ][ 4 ], Bmatrix );
            }
            else{
                std::string outstr = "Unrecognized number of parameters ( " + std::to_string( outputs[ 1 ].size() ) + " ) for the B stiffness tensor";
                return new errorNode( "extractMaterialParameters",
                                      outstr.c_str() );
            }
    
            if ( error ){
                errorOut result = new errorNode( "extractMaterialParameters", "Error in computation of the B stiffness tensor" );
                result->addNext( error );
                return result;
            }
    
            if ( outputs[ 2 ].size() == 11 ){
                error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicC( outputs[ 2 ], Cmatrix );
            }
            else{
                std::string outstr = "Unrecognized number of parameters ( " + std::to_string( outputs[ 2 ].size() ) + " ) for the C stiffness tensor";
                return new errorNode( "extractMaterialParameters",
                                      outstr.c_str() );
            }
    
            if ( error ){
                errorOut result = new errorNode( "extractMaterialParameters", "Error in computation of the C stiffness tensor" );
                result->addNext( error );
                return result;
            }
    
            if ( outputs[ 3 ].size() == 2 ){
                error = tardigradeHydra::micromorphicLinearElasticity::formIsotropicD( outputs[ 3 ][ 0 ], outputs[ 3 ][ 1 ], Dmatrix );
            }
            else{
                std::string outstr = "Unrecognized number of parameters ( " + std::to_string( outputs[ 3 ].size() ) + " ) for the D stiffness tensor";
                return new errorNode( "extractMaterialParameters",
                                      outstr.c_str() );
            }
    
            if ( error ){
                errorOut result = new errorNode( "extractMaterialParameters", "Error in computation of the D stiffness tensor" );
                result->addNext( error );
                return result;
            }
    
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

            floatVector deformationGradient1;

            floatVector microDeformation1;

            floatVector gradientMicroDeformation1;

            if ( isPrevious ){

                deformationGradient1 = ( *hydra->get_previousConfigurations( ) )[ 0 ];

                microDeformation1 = ( *hydra->get_previousMicroConfigurations( ) )[ 0 ];

                gradientMicroDeformation1 = ( *hydra->get_previousGradientMicroConfigurations( ) )[ 0 ];

            }
            else{

                deformationGradient1 = ( *hydra->get_configurations( ) )[ 0 ];

                microDeformation1 = ( *hydra->get_microConfigurations( ) )[ 0 ];

                gradientMicroDeformation1 = ( *hydra->get_gradientMicroConfigurations( ) )[ 0 ];

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

            const unsigned int dim = *hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int tot_dim = sot_dim * dim;

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

                dFFollowdFs     = tardigradeVectorTools::appendVectors( hydra->getPreviousFollowingConfigurationJacobian( 0 ) );

                dChiFollowdChis = tardigradeVectorTools::appendVectors( hydra->getPreviousFollowingMicroConfigurationJacobian( 0 ) );

                followingConfiguration = hydra->getPreviousFollowingConfiguration( 0 );

                followingMicroConfiguration = hydra->getPreviousFollowingMicroConfiguration( 0 );

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

                dFFollowdFs     = tardigradeVectorTools::appendVectors( hydra->getFollowingConfigurationJacobian( 0 ) );

                dChiFollowdChis = tardigradeVectorTools::appendVectors( hydra->getFollowingMicroConfigurationJacobian( 0 ) );

                followingConfiguration = hydra->getFollowingConfiguration( 0 );
    
                followingMicroConfiguration = hydra->getFollowingMicroConfiguration( 0 );

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

            floatVector PK2Stress;

            floatVector referenceSymmetricMicroStress;

            floatVector referenceHigherOrderStress;

            floatVector dPK2dlocalPK2;

            floatVector dPK2dFFollow;

            floatVector dSIGMAdlocalSIGMA;

            floatVector dSIGMAdFFollow;

            floatVector dMdlocalM;

            floatVector dMdFFollow;

            floatVector dMdChiFollow;

            // Pull the stresses back to the true reference configuration
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pullBackCauchyStress( localPK2Stress, followingConfiguration, PK2Stress, dPK2dlocalPK2, dPK2dFFollow ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pullBackMicroStress( localReferenceSymmetricMicroStress, followingConfiguration, referenceSymmetricMicroStress,
                                                                                                         dSIGMAdlocalSIGMA, dSIGMAdFFollow ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pullBackHigherOrderStress( localReferenceHigherOrderStress, followingConfiguration, followingMicroConfiguration,
                                                                                                               referenceHigherOrderStress, dMdlocalM, dMdFFollow, dMdChiFollow ) );

            floatVector dPK2dF = tardigradeVectorTools::matrixMultiply( dPK2dlocalPK2, tardigradeVectorTools::matrixMultiply(     dlocalPK2dC,     *dCdF, sot_dim, sot_dim, sot_dim, sot_dim )
                                                                                     + tardigradeVectorTools::matrixMultiply(   dlocalPK2dPsi,   *dPsidF, sot_dim, sot_dim, sot_dim, sot_dim )
                                                                                     + tardigradeVectorTools::matrixMultiply( dlocalPK2dGamma, *dGammadF, sot_dim, tot_dim, tot_dim, sot_dim ), sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dPK2dFn = tardigradeVectorTools::matrixMultiply( dPK2dlocalPK2, tardigradeVectorTools::matrixMultiply(     dlocalPK2dC,     *dCdFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                                                      + tardigradeVectorTools::matrixMultiply(   dlocalPK2dPsi,   *dPsidFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                                                      + tardigradeVectorTools::matrixMultiply( dlocalPK2dGamma, *dGammadFn, sot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim ), sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                + tardigradeVectorTools::matrixMultiply( dPK2dFFollow, dFFollowdFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dPK2dChi = tardigradeVectorTools::matrixMultiply( dPK2dlocalPK2, tardigradeVectorTools::matrixMultiply(   dlocalPK2dPsi,   *dPsidChi, sot_dim, sot_dim, sot_dim, sot_dim )
                                                                                       + tardigradeVectorTools::matrixMultiply( dlocalPK2dGamma, *dGammadChi, sot_dim, tot_dim, tot_dim, sot_dim ), sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dPK2dChin = tardigradeVectorTools::matrixMultiply( dPK2dlocalPK2, tardigradeVectorTools::matrixMultiply(   dlocalPK2dPsi,   *dPsidChin, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                                                        + tardigradeVectorTools::matrixMultiply( dlocalPK2dGamma, *dGammadChin, sot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim ), sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dPK2dGradChi = tardigradeVectorTools::matrixMultiply( dPK2dlocalPK2, tardigradeVectorTools::matrixMultiply( dlocalPK2dGamma, *dGammadGradChi, sot_dim, tot_dim, tot_dim, tot_dim ), sot_dim, sot_dim, sot_dim, tot_dim );

            floatVector dPK2dGradChin = tardigradeVectorTools::matrixMultiply( dPK2dlocalPK2, tardigradeVectorTools::matrixMultiply( dlocalPK2dGamma, *dGammadGradChin, sot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim ), sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * tot_dim );

            floatVector dSIGMAdF = tardigradeVectorTools::matrixMultiply( dSIGMAdlocalSIGMA, tardigradeVectorTools::matrixMultiply(     dlocalSIGMAdC,     *dCdF, sot_dim, sot_dim, sot_dim, sot_dim )
                                                                                           + tardigradeVectorTools::matrixMultiply(   dlocalSIGMAdPsi,   *dPsidF, sot_dim, sot_dim, sot_dim, sot_dim )
                                                                                           + tardigradeVectorTools::matrixMultiply( dlocalSIGMAdGamma, *dGammadF, sot_dim, tot_dim, tot_dim, sot_dim ), sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dSIGMAdFn = tardigradeVectorTools::matrixMultiply( dSIGMAdlocalSIGMA, tardigradeVectorTools::matrixMultiply(     dlocalSIGMAdC,     *dCdFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                                                            + tardigradeVectorTools::matrixMultiply(   dlocalSIGMAdPsi,   *dPsidFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                                                            + tardigradeVectorTools::matrixMultiply( dlocalSIGMAdGamma, *dGammadFn, sot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim ), sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                  + tardigradeVectorTools::matrixMultiply( dSIGMAdFFollow, dFFollowdFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dSIGMAdChi = tardigradeVectorTools::matrixMultiply( dSIGMAdlocalSIGMA, tardigradeVectorTools::matrixMultiply(   dlocalSIGMAdPsi,   *dPsidChi, sot_dim, sot_dim, sot_dim, sot_dim )
                                                                                             + tardigradeVectorTools::matrixMultiply( dlocalSIGMAdGamma, *dGammadChi, sot_dim, tot_dim, tot_dim, sot_dim ), sot_dim, sot_dim, sot_dim, sot_dim );

            floatVector dSIGMAdChin = tardigradeVectorTools::matrixMultiply( dSIGMAdlocalSIGMA, tardigradeVectorTools::matrixMultiply(   dlocalSIGMAdPsi,   *dPsidChin, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                                                                              + tardigradeVectorTools::matrixMultiply( dlocalSIGMAdGamma, *dGammadChin, sot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim ), sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dSIGMAdGradChi = tardigradeVectorTools::matrixMultiply( dSIGMAdlocalSIGMA, tardigradeVectorTools::matrixMultiply( dlocalSIGMAdGamma, *dGammadGradChi, sot_dim, tot_dim, tot_dim, tot_dim ), sot_dim, sot_dim, sot_dim, tot_dim );

            floatVector dSIGMAdGradChin = tardigradeVectorTools::matrixMultiply( dSIGMAdlocalSIGMA, tardigradeVectorTools::matrixMultiply( dlocalSIGMAdGamma, *dGammadGradChin, sot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim ), sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * tot_dim );

            floatVector dMdF = tardigradeVectorTools::matrixMultiply( dMdlocalM, tardigradeVectorTools::matrixMultiply( dlocalMdGamma, *dGammadF, tot_dim, tot_dim, tot_dim, sot_dim ), tot_dim, tot_dim, tot_dim, sot_dim );

            floatVector dMdFn = tardigradeVectorTools::matrixMultiply( dMdlocalM, tardigradeVectorTools::matrixMultiply( dlocalMdGamma, *dGammadFn, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim ), tot_dim, tot_dim, tot_dim, ( num_configs -1 ) * sot_dim )
                              + tardigradeVectorTools::matrixMultiply( dMdFFollow, dFFollowdFn, tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dMdChi = tardigradeVectorTools::matrixMultiply( dMdlocalM, tardigradeVectorTools::matrixMultiply( dlocalMdGamma, *dGammadChi, tot_dim, tot_dim, tot_dim, sot_dim ), tot_dim, tot_dim, tot_dim, sot_dim );

            floatVector dMdChin = tardigradeVectorTools::matrixMultiply( dMdlocalM, tardigradeVectorTools::matrixMultiply( dlocalMdGamma, *dGammadChin, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim ), tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim )
                                + tardigradeVectorTools::matrixMultiply( dMdChiFollow, dChiFollowdChin, tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim );

            floatVector dMdGradChi = tardigradeVectorTools::matrixMultiply( dMdlocalM, tardigradeVectorTools::matrixMultiply( dlocalMdGamma, *dGammadGradChi, tot_dim, tot_dim, tot_dim, tot_dim ), tot_dim, tot_dim, tot_dim, tot_dim );

            floatVector dMdGradChin = tardigradeVectorTools::matrixMultiply( dMdlocalM, tardigradeVectorTools::matrixMultiply( dlocalMdGamma, *dGammadGradChin, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim ), tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim );

            if ( isPrevious ){

                set_previousPK2Stress( PK2Stress );

                set_previousReferenceSymmetricMicroStress( referenceSymmetricMicroStress );

                set_previousReferenceHigherOrderStress( referenceHigherOrderStress );

                set_previousdPK2dF( dPK2dF );

                set_previousdPK2dFn( dPK2dFn );

                set_previousdPK2dChi( dPK2dChi );

                set_previousdPK2dChin( dPK2dChin );

                set_previousdPK2dGradChi( dPK2dGradChi );

                set_previousdPK2dGradChin( dPK2dGradChin );

                set_previousdSIGMAdF( dSIGMAdF );

                set_previousdSIGMAdFn( dSIGMAdFn );

                set_previousdSIGMAdChi( dSIGMAdChi );

                set_previousdSIGMAdChin( dSIGMAdChin );

                set_previousdSIGMAdGradChi( dSIGMAdGradChi );

                set_previousdSIGMAdGradChin( dSIGMAdGradChin );

                set_previousdMdF( dMdF );

                set_previousdMdFn( dMdFn );

                set_previousdMdChi( dMdChi );

                set_previousdMdChin( dMdChin );

                set_previousdMdGradChi( dMdGradChi );

                set_previousdMdGradChin( dMdGradChin );

            }
            else{

                set_PK2Stress( PK2Stress );

                set_referenceSymmetricMicroStress( referenceSymmetricMicroStress );

                set_referenceHigherOrderStress( referenceHigherOrderStress );

                set_dPK2dF( dPK2dF );

                set_dPK2dFn( dPK2dFn );

                set_dPK2dChi( dPK2dChi );

                set_dPK2dChin( dPK2dChin );

                set_dPK2dGradChi( dPK2dGradChi );

                set_dPK2dGradChin( dPK2dGradChin );

                set_dSIGMAdF( dSIGMAdF );

                set_dSIGMAdFn( dSIGMAdFn );

                set_dSIGMAdChi( dSIGMAdChi );

                set_dSIGMAdChin( dSIGMAdChin );

                set_dSIGMAdGradChi( dSIGMAdGradChi );

                set_dSIGMAdGradChin( dSIGMAdGradChin );

                set_dMdF( dMdF );

                set_dMdFn( dMdFn );

                set_dMdChi( dMdChi );

                set_dMdChin( dMdChin );

                set_dMdGradChi( dMdGradChi );

                set_dMdGradChin( dMdGradChin );

            }

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

            if ( isPrevious ){

                PK2Stress = get_previousPK2Stress( );

                referenceSymmetricMicroStress = get_previousReferenceSymmetricMicroStress( );

                referenceHigherOrderStress = get_previousReferenceHigherOrderStress( );

                deformationGradient = hydra->getPreviousDeformationGradient( );

                microDeformation = hydra->getPreviousMicroDeformation( );

            }
            else{

                PK2Stress = get_PK2Stress( );

                referenceSymmetricMicroStress = get_referenceSymmetricMicroStress( );

                referenceHigherOrderStress = get_referenceHigherOrderStress( );

                deformationGradient = hydra->getDeformationGradient( );

                microDeformation = hydra->getMicroDeformation( );

            }

            floatVector cauchyStress;

            floatVector symmetricMicroStress;

            floatVector higherOrderStress;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( mapStressMeasuresToCurrent( *deformationGradient, *microDeformation, *PK2Stress, *referenceSymmetricMicroStress,
                                                                                   *referenceHigherOrderStress, cauchyStress, symmetricMicroStress, higherOrderStress )  );

            if ( isPrevious ){

                set_previousCauchyStress( cauchyStress );

                set_previousSymmetricMicroStress( symmetricMicroStress );

                set_previousHigherOrderStress( higherOrderStress );

            }
            else{

                set_cauchyStress( cauchyStress );

                set_symmetricMicroStress( symmetricMicroStress );

                set_higherOrderStress( higherOrderStress );

            }

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

            const unsigned int dim = *hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int tot_dim = sot_dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            const floatVector *PK2Stress;

            const floatVector *referenceSymmetricMicroStress;

            const floatVector *referenceHigherOrderStress;

            const floatVector *deformationGradient;

            const floatVector *microDeformation;

            const floatVector *dPK2dF;

            const floatVector *dPK2dFn;

            const floatVector *dPK2dChi;

            const floatVector *dPK2dChin;

            const floatVector *dPK2dGradChi;

            const floatVector *dPK2dGradChin;

            const floatVector *dSIGMAdF;

            const floatVector *dSIGMAdFn;

            const floatVector *dSIGMAdChi;

            const floatVector *dSIGMAdChin;

            const floatVector *dSIGMAdGradChi;

            const floatVector *dSIGMAdGradChin;

            const floatVector *dMdF;

            const floatVector *dMdFn;

            const floatVector *dMdChi;

            const floatVector *dMdChin;

            const floatVector *dMdGradChi;

            const floatVector *dMdGradChin;

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

            }

            floatVector cauchyStress;

            floatVector symmetricMicroStress;

            floatVector higherOrderStress;

            floatVector dCauchyStressdF;

            floatVector dCauchyStressdPK2Stress;

            floatVector dMicroStressdF;

            floatVector dMicroStressdSIGMA;

            floatVector dHigherOrderStressdF;

            floatVector dHigherOrderStressdChi;

            floatVector dHigherOrderStressdM;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( mapStressMeasuresToCurrent( *deformationGradient, *microDeformation, *PK2Stress, *referenceSymmetricMicroStress,
                                                                                   *referenceHigherOrderStress, cauchyStress, symmetricMicroStress, higherOrderStress,
                                                                                    dCauchyStressdF, dCauchyStressdPK2Stress,
                                                                                    dMicroStressdF,  dMicroStressdSIGMA,
                                                                                    dHigherOrderStressdF, dHigherOrderStressdChi, dHigherOrderStressdM )  );

            if ( isPrevious ){

                set_previousdCauchyStressdF(                       tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, *dPK2dF         , sot_dim, sot_dim, sot_dim, sot_dim ) + dCauchyStressdF        );

                set_previousdCauchyStressdFn(                      tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, *dPK2dFn        , sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )    );

                set_previousdCauchyStressdChi(                     tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, *dPK2dChi       , sot_dim, sot_dim, sot_dim, sot_dim )                          );

                set_previousdCauchyStressdChin(                    tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, *dPK2dChin      , sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )    );

                set_previousdCauchyStressdGradChi(                 tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, *dPK2dGradChi   , sot_dim, sot_dim, sot_dim, tot_dim )                          );

                set_previousdCauchyStressdGradChin(                tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, *dPK2dGradChin  , sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * tot_dim )    );

                set_previousdSymmetricMicroStressdF(               tardigradeVectorTools::matrixMultiply( dMicroStressdSIGMA,      *dSIGMAdF       , sot_dim, sot_dim, sot_dim, sot_dim ) + dMicroStressdF         );

                set_previousdSymmetricMicroStressdFn(              tardigradeVectorTools::matrixMultiply( dMicroStressdSIGMA,      *dSIGMAdFn      , sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )    );

                set_previousdSymmetricMicroStressdChi(             tardigradeVectorTools::matrixMultiply( dMicroStressdSIGMA,      *dSIGMAdChi     , sot_dim, sot_dim, sot_dim, sot_dim )                          );

                set_previousdSymmetricMicroStressdChin(            tardigradeVectorTools::matrixMultiply( dMicroStressdSIGMA,      *dSIGMAdChin    , sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )    );

                set_previousdSymmetricMicroStressdGradChi(         tardigradeVectorTools::matrixMultiply( dMicroStressdSIGMA,      *dSIGMAdGradChi , sot_dim, sot_dim, sot_dim, tot_dim )                          );

                set_previousdSymmetricMicroStressdGradChin(        tardigradeVectorTools::matrixMultiply( dMicroStressdSIGMA,      *dSIGMAdGradChin, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * tot_dim )    );

                set_previousdHigherOrderStressdF(                  tardigradeVectorTools::matrixMultiply( dHigherOrderStressdM,    *dMdF           , tot_dim, tot_dim, tot_dim, sot_dim ) + dHigherOrderStressdF   );

                set_previousdHigherOrderStressdFn(                 tardigradeVectorTools::matrixMultiply( dHigherOrderStressdM,    *dMdFn          , tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim )    );

                set_previousdHigherOrderStressdChi(                tardigradeVectorTools::matrixMultiply( dHigherOrderStressdM,    *dMdChi         , tot_dim, tot_dim, tot_dim, sot_dim ) + dHigherOrderStressdChi );

                set_previousdHigherOrderStressdChin(               tardigradeVectorTools::matrixMultiply( dHigherOrderStressdM,    *dMdChin        , tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim )    );

                set_previousdHigherOrderStressdGradChi(            tardigradeVectorTools::matrixMultiply( dHigherOrderStressdM,    *dMdGradChi     , tot_dim, tot_dim, tot_dim, tot_dim )                          );

                set_previousdHigherOrderStressdGradChin(           tardigradeVectorTools::matrixMultiply( dHigherOrderStressdM,    *dMdGradChin    , tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim )    );

                set_previousCauchyStress( cauchyStress );

                set_previousSymmetricMicroStress( symmetricMicroStress );

                set_previousHigherOrderStress( higherOrderStress );

            }
            else{

                set_dCauchyStressdF(                       tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, *dPK2dF         , sot_dim, sot_dim, sot_dim, sot_dim ) + dCauchyStressdF        );

                set_dCauchyStressdFn(                      tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, *dPK2dFn        , sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )    );

                set_dCauchyStressdChi(                     tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, *dPK2dChi       , sot_dim, sot_dim, sot_dim, sot_dim )                          );

                set_dCauchyStressdChin(                    tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, *dPK2dChin      , sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )    );

                set_dCauchyStressdGradChi(                 tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, *dPK2dGradChi   , sot_dim, sot_dim, sot_dim, tot_dim )                          );

                set_dCauchyStressdGradChin(                tardigradeVectorTools::matrixMultiply( dCauchyStressdPK2Stress, *dPK2dGradChin  , sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * tot_dim )    );

                set_dSymmetricMicroStressdF(               tardigradeVectorTools::matrixMultiply( dMicroStressdSIGMA,      *dSIGMAdF       , sot_dim, sot_dim, sot_dim, sot_dim ) + dMicroStressdF         );

                set_dSymmetricMicroStressdFn(              tardigradeVectorTools::matrixMultiply( dMicroStressdSIGMA,      *dSIGMAdFn      , sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )    );

                set_dSymmetricMicroStressdChi(             tardigradeVectorTools::matrixMultiply( dMicroStressdSIGMA,      *dSIGMAdChi     , sot_dim, sot_dim, sot_dim, sot_dim )                          );

                set_dSymmetricMicroStressdChin(            tardigradeVectorTools::matrixMultiply( dMicroStressdSIGMA,      *dSIGMAdChin    , sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )    );

                set_dSymmetricMicroStressdGradChi(         tardigradeVectorTools::matrixMultiply( dMicroStressdSIGMA,      *dSIGMAdGradChi , sot_dim, sot_dim, sot_dim, tot_dim )                          );

                set_dSymmetricMicroStressdGradChin(        tardigradeVectorTools::matrixMultiply( dMicroStressdSIGMA,      *dSIGMAdGradChin, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * tot_dim )    );

                set_dHigherOrderStressdF(                  tardigradeVectorTools::matrixMultiply( dHigherOrderStressdM,    *dMdF           , tot_dim, tot_dim, tot_dim, sot_dim ) + dHigherOrderStressdF   );

                set_dHigherOrderStressdFn(                 tardigradeVectorTools::matrixMultiply( dHigherOrderStressdM,    *dMdFn          , tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim )    );

                set_dHigherOrderStressdChi(                tardigradeVectorTools::matrixMultiply( dHigherOrderStressdM,    *dMdChi         , tot_dim, tot_dim, tot_dim, sot_dim ) + dHigherOrderStressdChi );

                set_dHigherOrderStressdChin(               tardigradeVectorTools::matrixMultiply( dHigherOrderStressdM,    *dMdChin        , tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim )    );

                set_dHigherOrderStressdGradChi(            tardigradeVectorTools::matrixMultiply( dHigherOrderStressdM,    *dMdGradChi     , tot_dim, tot_dim, tot_dim, tot_dim )                          );

                set_dHigherOrderStressdGradChin(           tardigradeVectorTools::matrixMultiply( dHigherOrderStressdM,    *dMdGradChin    , tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim )    );

                set_cauchyStress( cauchyStress );

                set_symmetricMicroStress( symmetricMicroStress );

                set_higherOrderStress( higherOrderStress );

            }

        }

        void residual::setDeformationJacobians( const bool isPrevious ){
            /*!
             * Evaluate the derived deformation Jacobians
             * 
             * We assume that the first configuration in hydra.get_configurations is the elastic one
             *
             * \param isPrevious: Flag for whether the measures to be calculated are in the current or previous configuration
             */

            const unsigned int dim = *hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int tot_dim = sot_dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            floatVector deformationGradient1;

            floatVector microDeformation1;

            floatVector gradientMicroDeformation1;

            floatVector dF1dF;

            floatVector dF1dFn;

            floatVector dChi1dChi;

            floatVector dChi1dChin;

            floatVector dGradChi1dFn;

            floatVector dGradChi1dChi;

            floatVector dGradChi1dChin;

            floatVector dGradChi1dGradChi;

            floatVector dGradChi1dGradChin;

            if ( isPrevious ){

                deformationGradient1 = ( *hydra->get_previousConfigurations( ) )[ 0 ];

                microDeformation1 = ( *hydra->get_previousMicroConfigurations( ) )[ 0 ];

                gradientMicroDeformation1 = ( *hydra->get_previousGradientMicroConfigurations( ) )[ 0 ];

                dF1dF              = tardigradeVectorTools::appendVectors( *hydra->get_previousdF1dF( ) );

                dF1dFn             = tardigradeVectorTools::appendVectors( *hydra->get_previousdF1dFn( ) );

                dChi1dChi          = tardigradeVectorTools::appendVectors( *hydra->get_previousdChi1dChi( ) );

                dChi1dChin         = tardigradeVectorTools::appendVectors( *hydra->get_previousdChi1dChin( ) );

                dGradChi1dFn       = tardigradeVectorTools::appendVectors( *hydra->get_previousdGradChi1dFn( ) );

                dGradChi1dChi      = tardigradeVectorTools::appendVectors( *hydra->get_previousdGradChi1dChi( ) );

                dGradChi1dChin     = tardigradeVectorTools::appendVectors( *hydra->get_previousdGradChi1dChin( ) );

                dGradChi1dGradChi  = tardigradeVectorTools::appendVectors( *hydra->get_previousdGradChi1dGradChi( ) );

                dGradChi1dGradChin = tardigradeVectorTools::appendVectors( *hydra->get_previousdGradChi1dGradChin( ) );

            }
            else{

                deformationGradient1 = ( *hydra->get_configurations( ) )[ 0 ];

                microDeformation1 = ( *hydra->get_microConfigurations( ) )[ 0 ];

                gradientMicroDeformation1 = ( *hydra->get_gradientMicroConfigurations( ) )[ 0 ];

                dF1dF              = tardigradeVectorTools::appendVectors( *hydra->get_dF1dF( ) );

                dF1dFn             = tardigradeVectorTools::appendVectors( *hydra->get_dF1dFn( ) );

                dChi1dChi          = tardigradeVectorTools::appendVectors( *hydra->get_dChi1dChi( ) );

                dChi1dChin         = tardigradeVectorTools::appendVectors( *hydra->get_dChi1dChin( ) );

                dGradChi1dFn       = tardigradeVectorTools::appendVectors( *hydra->get_dGradChi1dFn( ) );

                dGradChi1dChi      = tardigradeVectorTools::appendVectors( *hydra->get_dGradChi1dChi( ) );

                dGradChi1dChin     = tardigradeVectorTools::appendVectors( *hydra->get_dGradChi1dChin( ) );

                dGradChi1dGradChi  = tardigradeVectorTools::appendVectors( *hydra->get_dGradChi1dGradChi( ) );

                dGradChi1dGradChin = tardigradeVectorTools::appendVectors( *hydra->get_dGradChi1dGradChin( ) );

            }

            floatVector rightCauchyGreen;

            floatVector Psi;

            floatVector Gamma;

            floatVector dCdF1;

            floatVector dPsidF1;

            floatVector dPsidChi1;

            floatVector dGammadF1;

            floatVector dGammadGradChi1;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( computeDeformationMeasures( deformationGradient1, microDeformation1, gradientMicroDeformation1,
                                                                                   rightCauchyGreen, Psi, Gamma,
                                                                                   dCdF1, dPsidF1, dPsidChi1, dGammadF1, dGammadGradChi1 ) );

            if ( isPrevious ){

                set_previousRightCauchyGreen( rightCauchyGreen );

                set_previousPsi( Psi );

                set_previousGamma( Gamma );

                set_previousdRightCauchyGreendF( tardigradeVectorTools::matrixMultiply( dCdF1, dF1dF, sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_previousdRightCauchyGreendFn( tardigradeVectorTools::matrixMultiply( dCdF1, dF1dFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_previousdPsidF( tardigradeVectorTools::matrixMultiply( dPsidF1, dF1dF, sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_previousdPsidFn( tardigradeVectorTools::matrixMultiply( dPsidF1, dF1dFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_previousdPsidChi( tardigradeVectorTools::matrixMultiply( dPsidChi1, dChi1dChi, sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_previousdPsidChin( tardigradeVectorTools::matrixMultiply( dPsidChi1, dChi1dChin, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_previousdGammadF( tardigradeVectorTools::matrixMultiply( dGammadF1, dF1dF, tot_dim, sot_dim, sot_dim, sot_dim ) );

                set_previousdGammadFn( tardigradeVectorTools::matrixMultiply( dGammadF1, dF1dFn, tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                                     + tardigradeVectorTools::matrixMultiply( dGammadGradChi1, dGradChi1dFn, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_previousdGammadChi( tardigradeVectorTools::matrixMultiply( dGammadGradChi1, dGradChi1dChi, tot_dim, tot_dim, tot_dim, sot_dim ) );

                set_previousdGammadChin( tardigradeVectorTools::matrixMultiply( dGammadGradChi1, dGradChi1dChin, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim  ) );

                set_previousdGammadGradChi( tardigradeVectorTools::matrixMultiply( dGammadGradChi1, dGradChi1dGradChi, tot_dim, tot_dim, tot_dim, tot_dim ) );

                set_previousdGammadGradChin( tardigradeVectorTools::matrixMultiply( dGammadGradChi1, dGradChi1dGradChin, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim ) );

            }
            else{

                set_rightCauchyGreen( rightCauchyGreen );

                set_psi( Psi );

                set_gamma( Gamma );

                set_dRightCauchyGreendF( tardigradeVectorTools::matrixMultiply( dCdF1, dF1dF, sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dRightCauchyGreendFn( tardigradeVectorTools::matrixMultiply( dCdF1, dF1dFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_dPsidF( tardigradeVectorTools::matrixMultiply( dPsidF1, dF1dF, sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dPsidFn( tardigradeVectorTools::matrixMultiply( dPsidF1, dF1dFn, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_dPsidChi( tardigradeVectorTools::matrixMultiply( dPsidChi1, dChi1dChi, sot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dPsidChin( tardigradeVectorTools::matrixMultiply( dPsidChi1, dChi1dChin, sot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_dGammadF( tardigradeVectorTools::matrixMultiply( dGammadF1, dF1dF, tot_dim, sot_dim, sot_dim, sot_dim ) );

                set_dGammadFn( tardigradeVectorTools::matrixMultiply( dGammadF1, dF1dFn, tot_dim, sot_dim, sot_dim, ( num_configs - 1 ) * sot_dim )
                             + tardigradeVectorTools::matrixMultiply( dGammadGradChi1, dGradChi1dFn, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim ) );

                set_dGammadChi( tardigradeVectorTools::matrixMultiply( dGammadGradChi1, dGradChi1dChi, tot_dim, tot_dim, tot_dim, sot_dim ) );

                set_dGammadChin( tardigradeVectorTools::matrixMultiply( dGammadGradChi1, dGradChi1dChin, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * sot_dim  ) );

                set_dGammadGradChi( tardigradeVectorTools::matrixMultiply( dGammadGradChi1, dGradChi1dGradChi, tot_dim, tot_dim, tot_dim, tot_dim ) );

                set_dGammadGradChin( tardigradeVectorTools::matrixMultiply( dGammadGradChi1, dGradChi1dGradChin, tot_dim, tot_dim, tot_dim, ( num_configs - 1 ) * tot_dim ) );

            }

        }

        void residual::setResidual( ){
            /*!
             * Set the residual w.r.t. the unknown vector
             */

            const floatVector *stress = hydra->getStress( );

            TARDIGRADE_ERROR_TOOLS_CATCH( setResidual( *stress 
                                                       - tardigradeVectorTools::appendVectors( { *get_PK2Stress( ), *get_referenceSymmetricMicroStress( ), *get_referenceHigherOrderStress( ) } ) ) );

        }

        void residual::setJacobian( ){
            /*!
             * Set the Jacobian w.r.t. the unknown vector
             */

            const unsigned int dim = *hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int tot_dim = dim * dim * dim;

            const unsigned int num_configs = *hydra->getNumConfigurations( );

            floatMatrix jacobian( *getNumEquations( ), floatVector( hydra->getUnknownVector( )->size( ), 0 ) );

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
                    jacobian[ row ][ row ] += 1;

                    // Jacobians w.r.t. the sub configurations
                    col = *hydra->getConfigurationUnknownCount( );

                    for ( auto Sn = S->begin( ); Sn != S->end( ); Sn++ ){

                        for ( unsigned int j = 0; j < ( num_configs - 1 ) * dims[ ( unsigned int )( Sn - S->begin( ) ) ]; j++ ){

                            jacobian[ row ][ col ] -= ( **Sn )[ ( num_configs - 1 ) * dims[ ( unsigned int ) ( Sn - S->begin( ) ) ] * i + j ];

                            col++;

                        }

                    }

                    row++;

                }

            }

            setJacobian( jacobian );

        }

        void residual::setdRdD( ){
            /*!
             * Set the Jacobian w.r.t. the deformation
             */

            const unsigned int dim = *hydra->getDimension( );

            const unsigned int sot_dim = dim * dim;

            const unsigned int tot_dim = dim * dim * dim;

            floatMatrix dRdD( *getNumEquations( ), floatVector( *hydra->getConfigurationUnknownCount( ), 0 ) );

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

                            dRdD[ row ][ col ] -= ( **Sn )[ dims[ ( unsigned int )( Sn - S->begin( ) ) ] * i + j ];

                            col++;

                        }

                    }

                    row++;

                }

            }

            setdRdD( dRdD );

        }

        void residual::setStress( ){
            /*!
             * Set the stresses
             */

            setStress( tardigradeVectorTools::appendVectors( { *get_PK2Stress( ), *get_referenceSymmetricMicroStress( ), *get_referenceHigherOrderStress( ) } ) );

        }

        void residual::setPreviousStress( ){
            /*!
             * Set the previous stresses
             */

            setPreviousStress( tardigradeVectorTools::appendVectors( { *get_previousPK2Stress( ), *get_previousReferenceSymmetricMicroStress( ), *get_previousReferenceHigherOrderStress( ) } ) );

        }

        void residual::setdRdT( ){
            /*!
             * Set the derivative of the residual w.r.t. the temperature
             */

            setdRdT( floatVector( *getNumEquations( ), 0 ) );

        }

    }

}
