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
                                   variableMatrix &dCauchyStressdF, variableMatrix &dCauchyStressdChi, variableMatrix &dCauchyStressdGradChi,
                                   variableMatrix &dMicroStressdF, variableMatrix &dMicroStressdChi, variableMatrix &dMicroStressdGradChi,
                                   variableMatrix &dHigherOrderStressdF, variableMatrix &dHigherOrderStressdChi,
                                   variableMatrix &dHigherOrderStressdGradChi ){
            /*!
             * Compute the stress measures in the current configuration for micromorphic linear elasticity based off 
             * of a quadratic decomposition of the energy.
             *
             * Also compute the Jacobians
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
             * :param variableMatrix &dCauchyStressdF: The Jacobian of the Cauchy stress w.r.t. the deformation gradient
             * :param variableMatrix &dCauchyStressdChi: The Jacobian of the Cauchy stress w.r.t. the micro deformation.
             * :param variableMatrix &dCauchyStressdGradChi: The Jacobian of the Cauchy stress w.r.t. the gradient of the 
             *     micro-deformation.
             * :param variableMatrix &dMicroStressdF: The Jacobian of the Micro stress w.r.t. the deformation gradient
             * :param variableMatrix &dMicroStressdChi: The Jacobian of the Micro stress w.r.t. the micro deformation.
             * :param variableMatrix &dMicroStressdGradChi: The Jacobian of the Micro stress w.r.t. the gradient of the 
             *     micro-deformation.
             * :param variableMatrix &dHigherOrderStressdF: The Jacobian of the Higher Order stress w.r.t. the deformation gradient
             * :param variableMatrix &dHigherOrderStressdChi: The Jacobian of the Higher Order stress w.r.t. the micro deformation.
             * :param variableMatrix &dHigherOrderStressdGradChi: The Jacobian of the Higher Order stress w.r.t. the gradient of the 
             *     micro-deformation.
             */
    
            variableVector PK2Stress, referenceMicroStress, referenceHigherOrderStress;
    
            variableMatrix dPK2StressdF, dPK2StressdChi, dPK2StressdGradChi;
            variableMatrix dReferenceMicroStressdF, dReferenceMicroStressdChi, dReferenceMicroStressdGradChi;
            variableMatrix dReferenceHigherOrderStressdF, dReferenceHigherOrderStressdGradChi;
    
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
    
            variableMatrix dCauchyStressdPK2Stress, dMicroStressdReferenceMicroStress, dHigherOrderStressdReferenceHigherOrderStress;
    
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
            dCauchyStressdF += tardigradeVectorTools::dot( dCauchyStressdPK2Stress, dPK2StressdF );
            dCauchyStressdChi = tardigradeVectorTools::dot( dCauchyStressdPK2Stress, dPK2StressdChi );
            dCauchyStressdGradChi = tardigradeVectorTools::dot( dCauchyStressdPK2Stress, dPK2StressdGradChi );
    
            //Assemble the jacobians of the symmetric micro-stress
            dMicroStressdF += tardigradeVectorTools::dot( dMicroStressdReferenceMicroStress, dReferenceMicroStressdF );
            dMicroStressdChi = tardigradeVectorTools::dot( dMicroStressdReferenceMicroStress, dReferenceMicroStressdChi );
            dMicroStressdGradChi = tardigradeVectorTools::dot( dMicroStressdReferenceMicroStress, dReferenceMicroStressdGradChi );
    
            //Assemble the jacobians of the higher-order stress
            dHigherOrderStressdF += tardigradeVectorTools::dot( dHigherOrderStressdReferenceHigherOrderStress,
                                                      dReferenceHigherOrderStressdF );
            dHigherOrderStressdGradChi = tardigradeVectorTools::dot( dHigherOrderStressdReferenceHigherOrderStress,
                                                          dReferenceHigherOrderStressdGradChi );
    
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
             * :param const variableVector &deformationGradient: The deformation gradient
             * :param const variableVector &microDeformation: The micro-deformation
             * :param const variableVector &gradientMicroDeformation: The spatial gradient of the micro-deformation
             * :param const parameterVector &A: The A stiffness matrix.
             * :param const parameterVector &B: The B stiffness matrix.
             * :param const parameterVector &C: The C stiffness matrix.
             * :param const parameterVector &D: The D stiffness matrix.
             * :param variableVector &PK2Stress: The second Piola-Kirchoff stress.
             * :param variableVector &referenceMicroStress: The symmetric micro-stress in the 
             *     reference configuration.
             * :param variableVector &referenceHigherOrderStress: The higher-order stress in the 
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
                                            variableMatrix &dPK2StressdF, variableMatrix &dPK2StressdChi, variableMatrix &dPK2StressdGradChi,
                                            variableMatrix &dReferenceMicroStressdF, variableMatrix &dReferenceMicroStressdChi,
                                            variableMatrix &dReferenceMicroStressdGradChi, variableMatrix &dMdF, variableMatrix &dMdGradChi ){
            /*!
             * Compute the stress measures in the reference configuration for micromorphic linear elasticity based off
             * of a quadratic decomposition of the energy.
             *
             * Also computes the Jacobians
             *
             * :param const variableVector &deformationGradient: The deformation gradient
             * :param const variableVector &microDeformation: The micro-deformation
             * :param const variableVector &gradientMicroDeformation: The spatial gradient of the micro-deformation
             * :param const parameterVector &A: The A stiffness matrix.
             * :param const parameterVector &B: The B stiffness matrix.
             * :param const parameterVector &C: The C stiffness matrix.
             * :param const parameterVector &D: The D stiffness matrix.
             * :param variableVector &PK2Stress: The second Piola-Kirchoff stress.
             * :param variableVector &referenceMicroStress: The symmetric micro-stress in the 
             *     reference configuration.
             * :param variableVector &referenceHigherOrderStress: The higher-order stress in the 
             *     reference configuration.
             * :param variableMatrix &dPK2StressdF: The Jacobian of the PK2 stress w.r.t. the deformation gradient.
             * :param variableMatrix &dPK2StressdChi: The Jacobian of the PK2 stress w.r.t. the micro deformation.
             * :param variableMatrix &dPK2StressdGradChi: The Jacobian of the PK2 stress w.r.t. the gradient of the micro deformation.
             * :param variableMatrix &dReferenceMicroStressdF: The Jacobian of the Micro stress w.r.t. the deformation gradient.
             * :param variableMatrix &dReferenceMicroStressdChi: The Jacobian of the Micro stress w.r.t. the micro deformation.
             * :param variableMatrix &dReferenceStressdGradChi: The Jacobian of the Micro stress w.r.t. the gradient of the micro deformation.
             * :param variableMatrix &dMdF: The Jacobian of the higher order stress w.r.t. the deformation gradient.
             * :param variableMatrix &dMdGradChi: The Jacobian of the higher order stress w.r.t. the gradient of the micro deformation.
             */
    
            //Assume 3d
            unsigned int dim = 3;
    
            constantVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );
    
            //Compute the required deformation measures
            variableVector RCG, Psi, Gamma;
            variableMatrix dRCGdF, dPsidF, dPsidChi, dGammadF, dGammadGradChi;
            errorOut error = computeDeformationMeasures( deformationGradient, microDeformation, gradientMicroDeformation,
                                                         RCG, Psi, Gamma, dRCGdF, dPsidF, dPsidChi, dGammadF, dGammadGradChi );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReference (jacobian)",
                                                 "Error in the computation of the deformation measures" );
                result->addNext( error );
                return result;
            }
    
            variableMatrix dPK2StressdRCG, dPK2StressdPsi, dPK2StressdGamma;
            variableMatrix dReferenceMicroStressdRCG, dReferenceMicroStressdPsi, dReferenceMicroStressdGamma;
            variableMatrix dMdGamma;
    
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
    
            dPK2StressdF = tardigradeVectorTools::dot( dPK2StressdRCG, dRCGdF )
                         + tardigradeVectorTools::dot( dPK2StressdPsi, dPsidF )
                         + tardigradeVectorTools::dot( dPK2StressdGamma, dGammadF );
    
            dPK2StressdChi = tardigradeVectorTools::dot( dPK2StressdPsi, dPsidChi );
    
            dPK2StressdGradChi = tardigradeVectorTools::dot( dPK2StressdGamma, dGammadGradChi );
    
            dReferenceMicroStressdF = tardigradeVectorTools::dot( dReferenceMicroStressdRCG, dRCGdF )
                                    + tardigradeVectorTools::dot( dReferenceMicroStressdPsi, dPsidF )
                                    + tardigradeVectorTools::dot( dReferenceMicroStressdGamma, dGammadF );
    
            dReferenceMicroStressdChi = tardigradeVectorTools::dot( dReferenceMicroStressdPsi, dPsidChi );
    
            dReferenceMicroStressdGradChi = tardigradeVectorTools::dot( dReferenceMicroStressdGamma, dGammadGradChi );
    
            dMdF = tardigradeVectorTools::dot( dMdGamma, dGammadF );
            dMdGradChi = tardigradeVectorTools::dot( dMdGamma, dGammadGradChi );
    
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
             * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation metric
             * :param const variableVector &Psi: The micro-deformation measure Psi
             * :param const variableVector &Gamma: The higher order deformation measure Gamma
             * :param const parameterVector &A: The A stiffness matrix.
             * :param const parameterVector &B: The B stiffness matrix.
             * :param const parameterVector &C: The C stiffness matrix.
             * :param const parameterVector &D: The D stiffness matrix.
             * :param variableVector &PK2Stress: The second Piola-Kirchoff stress.
             * :param variableVector &referenceMicroStress: The symmetric micro-stress in the 
             *     reference configuration.
             * :param variableVector &referenceHigherOrderStress: The higher-order stress in the 
             *     reference configuration.
             */
    
            //Assume 3D
            unsigned int dim = 3;
    
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
                                                           variableMatrix &dPK2StressdRCG, variableMatrix &dPK2StressdPsi,
                                                           variableMatrix &dPK2StressdGamma,
                                                           variableMatrix &dReferenceMicroStressdRCG,
                                                           variableMatrix &dReferenceMicroStressdPsi,
                                                           variableMatrix &dReferenceMicroStressdGamma,
                                                           variableMatrix &dMdGamma ){
            /*!
             * Compute the stress measures in the reference configuration for micromorphic linear elasticity based off
             * of a quadratic decomposition of the energy.
             *
             * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation metric
             * :param const variableVector &Psi: The micro-deformation measure Psi
             * :param const variableVector &Gamma: The higher order deformation measure Gamma
             * :param const parameterVector &A: The A stiffness matrix.
             * :param const parameterVector &B: The B stiffness matrix.
             * :param const parameterVector &C: The C stiffness matrix.
             * :param const parameterVector &D: The D stiffness matrix.
             * :param variableVector &PK2Stress: The second Piola-Kirchoff stress.
             * :param variableVector &referenceMicroStress: The symmetric micro-stress in the 
             *     reference configuration.
             * :param variableVector &referenceHigherOrderStress: The higher-order stress in the 
             *     reference configuration.
             * :param variableMatrix &dPK2StressdRCG: The Jacobian of the PK2 stress w.r.t. the right Cauchy-Green
             *     deformation metric.
             * :param variableMatrix &dPK2StressdPsi: The Jacobian of the PK2 stress w.r.t. the micro deformation 
             *     metric.
             * :param variableMatrix &dPK2StressdGamma: The Jacobian of the PK2 stress w.r.t. the higher order 
             *     deformation measure.
             * :param variableMatrix &dReferenceMicroStressdRCG: The Jacobian of the reference micro stress w.r.t. the 
             *     right Cacuhy-Green deformation metric.
             * :param variableMatrix &dReferenceMicroStressdPsi: The Jacobian of the reference micro stress w.r.t. the 
             *     micro deformation measure.
             * :param variableMatrix &dReferenceMicroStrssdGamma: The Jacobian of the reference micro stress w.r.t. the 
             *     higher order deformation measure.
             * :param variableMatrix &dMdGamma: The Jacobian of the reference higher order stress w.r.t. 
             *     the higher order deformation measure.
             */
    
            //Assume 3d
            unsigned int dim = 3;
    
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
    
            variableMatrix dTerm1dRCG, dTerm1dPsi;
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
            variableMatrix dInvRCGPsidRCG, dInvRCGPsidPsi;
    
            error = computeInvRCGPsi( invRCG, Psi, invRCGPsi, dInvRCGPsidRCG, dInvRCGPsidPsi );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMetrics (jacobian)",
                                                 "Error in computation of invRCG Psi product" );
                result->addNext( error );
                return result;
            }
    
            variableVector term2;
            variableMatrix dTerm2dRCG, dTerm2dPsi, dTerm2dInvRCGPsi;
            error = computeLinearElasticTerm2( greenLagrangeStrain, microStrain, invRCGPsi, B, D, term2,
                                               dTerm2dRCG, dTerm2dPsi, dTerm2dInvRCGPsi );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMetrics (jacobian)",
                                                 "Error in computation of term 2" );
                result->addNext( error );
                return result;
            }
    
            dTerm2dRCG *= 0.5;
            dTerm2dRCG += tardigradeVectorTools::dot( dTerm2dInvRCGPsi, dInvRCGPsidRCG );
    
            dTerm2dPsi += tardigradeVectorTools::dot( dTerm2dInvRCGPsi, dInvRCGPsidPsi );
    
            //Compute the third common term for the PK2 and symmetric micro-stress
            variableVector invRCGGamma;
            variableMatrix dInvRCGGammadRCG, dInvRCGGammadGamma;
    
            error = computeInvRCGGamma( invRCG, Gamma, invRCGGamma, dInvRCGGammadRCG, dInvRCGGammadGamma );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMetrics (jacobian)",
                                                 "Error in computation of invRCG Gamma product" );
                result->addNext( error );
                return result;
            }
    
            variableVector term3;
            variableMatrix dTerm3dInvRCGGamma, dTerm3dM;
            error = computeLinearElasticTerm3( invRCGGamma, referenceHigherOrderStress, term3, dTerm3dInvRCGGamma, dTerm3dM );
    
            if ( error ){
                errorOut result = new errorNode( "linearElasticityReferenceDerivedMetrics (jacobian)",
                                                 "Error in computation of term 3" );
                result->addNext( error );
                return result;
            }
    
            variableMatrix dTerm3dRCG = tardigradeVectorTools::dot( dTerm3dInvRCGGamma, dInvRCGGammadRCG );
            variableMatrix dTerm3dGamma = tardigradeVectorTools::dot( dTerm3dInvRCGGamma, dInvRCGGammadGamma )
                                        + tardigradeVectorTools::dot( dTerm3dM, dMdGamma );
    
            //Construct the PK2 and reference symmetric stresses
            PK2Stress            = term1 + term2 + term3;
    
            dPK2StressdRCG    = dTerm1dRCG + dTerm2dRCG + dTerm3dRCG;
            dPK2StressdPsi    = dTerm1dPsi + dTerm2dPsi;
            dPK2StressdGamma  = dTerm3dGamma;
    
            variableVector symmTerm2Term3;
            variableMatrix dSymmTerm2Term3dTerm2Term3;
            error = tardigradeConstitutiveTools::computeSymmetricPart( term2 + term3, symmTerm2Term3, dSymmTerm2Term3dTerm2Term3 );
            referenceMicroStress = term1 + 2 * symmTerm2Term3;
    
            dReferenceMicroStressdRCG = dTerm1dRCG + 2 * ( tardigradeVectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm2dRCG )
                                                         + tardigradeVectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm3dRCG ) );
    
            dReferenceMicroStressdPsi = dTerm1dPsi + 2 * tardigradeVectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm2dPsi );
            dReferenceMicroStressdGamma = 2 * tardigradeVectorTools::dot( dSymmTerm2Term3dTerm2Term3, dTerm3dGamma );
    
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
             * :param const variableVector &deformationGradient: The deformation gradient between the 
             *     reference configuration and the current configuration.
             * :param const variableVector &microDeformation: The micro-deformation map between the 
             *     reference configuration and the current configuration.
             * :param const variableVector &PK2Stress: The Second Piola-Kirchoff stress.
             * :param const variableVector &referenceMicroStress: The symmetric micro-stress in the 
             *     reference configuration.
             * :param const variableVector &referenceHigherOrderStress: The higher order stress in 
             *     the reference configuration.
             * :param variableVector &cauchyStress: The Cauchy stress (PK2 stress in the current configuration).
             * :param variableVector &microStress: The symmetric micro-stress in the current configuration.
             * :param variableVector &higherOrderStress: The higher order stress in the current configuration.
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
                                             variableMatrix &dCauchyStressdF, variableMatrix &dCauchyStressdPK2Stress,
                                             variableMatrix &dMicroStressdF, variableMatrix &dMicroStressdReferenceMicroStress,
                                             variableMatrix &dHigherOrderStressdF, variableMatrix &dHigherOrderStressdChi,
                                             variableMatrix &dHigherOrderStressdReferenceHigherOrderStress ){
            /*!
             * Map the stress measures in the reference configuration to the current configuration.
             *
             * Also computes the Jacobians
             *
             * :param const variableVector &deformationGradient: The deformation gradient between the 
             *     reference configuration and the current configuration.
             * :param const variableVector &microDeformation: The micro-deformation map between the 
             *     reference configuration and the current configuration.
             * :param const variableVector &PK2Stress: The Second Piola-Kirchoff stress.
             * :param const variableVector &referenceMicroStress: The symmetric micro-stress in the 
             *     reference configuration.
             * :param const variableVector &referenceHigherOrderStress: The higher order stress in 
             *     the reference configuration.
             * :param variableVector &cauchyStress: The Cauchy stress (PK2 stress in the current configuration).
             * :param variableVector &microStress: The symmetric micro-stress in the current configuration.
             * :param variableVector &higherOrderStress: The higher order stress in the current configuration.
             * :param variableMatrix &dCauchyStressdF: The Jacobian of the Cauchy stress w.r.t. the 
             *     deformation gradient.
             * :param variableMatrix &dCauchyStressdPK2Stress: The Jacobian of the Cauchy stress w.r.t. the 
             *     PK2 stress.
             * :param variableMatrix &dMicroStressdF: The Jacobian of the micro stress w.r.t. the 
             *     deformation gradient.
             * :param variableMatrix &dMicroStressdReferenceMicroStress: The Jacobian of the micro-stress 
             *     in the current configuration w.r.t. the micro-stress in the reference configuration.
             * :param variableMatrix &dHigherOrderStressdF: The Jacobian of the higher-order stress w.r.t.
             *     the deformation gradient.
             * :param variableMatrix &dHigherOrderStressdChi: The Jacobian of the higher-order stress 
             *     w.r.t. the micro-deformation.
             * :param variableMatrix &dHigherOrderStressdReferenceHigherOrderStress: The Jacobian of the 
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
             * :param const variableVector &deformationGradient: The deformation gradient
             * :param const variableVector &microDeformation: The micro-deformation
             * :param const variableVector &gradientMicroDeformation: The spatial gradient of the micro-deformation
             * :param variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor
             * :param variableVector &Psi: The micro-deformation measure
             * :param variableVector &Gamma: The gradient micro-deformation measure
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
                                             variableMatrix &dCdF, variableMatrix &dPsidF, variableMatrix &dPsidChi,
                                             variableMatrix &dGammadF, variableMatrix &dGammadGradChi ){
            /*!
             * Compute the deformation measures
             * \f$\begin{align}
             * C_{IJ} &= F_{iI} F_{iJ}\\
             * Psi_{IJ} &= F_{iI} \Chi_{iJ}\\
             * \Gamma_{IJK} &= F_{iI} \Chi_{iJ, K}
             * \end{align}\f$
             *
             * :param const variableVector &deformationGradient: The deformation gradient
             * :param const variableVector &microDeformation: The micro-deformation
             * :param const variableVector &gradientMicroDeformation: The spatial gradient of the micro-deformation
             * :param variableVector &rightCauchyGreen: The Right Cauchy-Green deformation tensor
             * :param variableVector &Psi: The micro-deformation measure
             * :param variableVector &Gamma: The gradient micro-deformation measure
             * :param variableMatrix &dCdF: The gradient of the right Cauchy green deformation tensor w.r.t. 
             *     the deformation gradient.
             * :param variableMatrix &dPsidF: The gradient of Psi w.r.t. the deformation gradient.
             * :param variableMatrix &dPsidChi: The gradient of Psi w.r.t. the microDeformation.
             * :param variableMatrix &dGammadF: The gradient of Gamma w.r.t. the deformation gradient.
             * :param variableMatrix &dGammadGradChi: The gradient of Gamma w.r.t. the spatial gradient of Chi
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
             * :param const variableVector &greenLagrangeStrain: The Green-Lagrange strain.
             * :param const variableVector &microStrain: The micro-strain
             * :param const parameterVector &A: The A stiffness matrix
             * :param const parameterVEctor &D: The D stiffness matrix
             * :param variableVector &term1: The first term.
             */
    
            //Assume 3D
            unsigned int dim = 3;
    
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
                                            variableMatrix &dTerm1dGreenLagrangeStrain, variableMatrix &dTerm1dMicroStrain ){
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
             * :param const variableVector &greenLagrangeStrain: The Green-Lagrange strain.
             * :param const variableVector &microStrain: The micro-strain
             * :param const parameterVector &A: The A stiffness matrix
             * :param const parameterVEctor &D: The D stiffness matrix
             * :param variableVector &term1: The first term.
             */
    
            //Assume 3D
            unsigned int dim = 3;
    
            errorOut error = computeLinearElasticTerm1( greenLagrangeStrain, microStrain, A, D, term1 );
    
            if ( error ){
                errorOut result = new errorNode( "computeLinearElasticTerm1 (jacobian)",
                                                 "Error in computation of linear elastic term1" );
                result->addNext( error );
                return result;
            }
    
            //Compute the first common term for the PK2 and symmetric micro-stress
            dTerm1dGreenLagrangeStrain = variableMatrix( term1.size(), variableVector( greenLagrangeStrain.size(), 0 ) );
            dTerm1dMicroStrain = variableMatrix( term1.size(), variableVector( microStrain.size(), 0 ) );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int M = 0; M < dim; M++ ){
                        for ( unsigned int N = 0; N < dim; N++ ){
                            dTerm1dGreenLagrangeStrain[ dim * I + J ][ dim * M + N ] = A[ dim * dim * dim * I + dim * dim * J + dim * M + N ];
                            dTerm1dMicroStrain[ dim * I + J ][ dim * M + N ] = D[ dim * dim * dim * I + dim * dim * J + dim * M + N ];
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
             * :param const variableVector &greenLagrangeStrain: The Green-Lagrange strain \f$E_{IJ} = \frac{1}{2} \left( C_{IJ} - \delta_{IJ} \right)\f$
             * :param const variableVector &microStrain: The micro-strain \f$\mathcal{E}_{IJ} = \Psi_{IJ} - \delta_{IJ}\f$
             * :param const variableVector &invCPsi: The product \f$C_{JR}^{-1} \Psi_{RQ}\f$
             * :param const variableVector &B: The B stiffness matrix
             * :param const variableVector &D: The D stiffness matrix
             * :param variableVector &term2: The second term.
             */
    
            //Assume 3D
            unsigned int dim = 3;
    
            if ( greenLagrangeStrain.size() != dim * dim ){
                return new errorNode( "computeLinearElasticTerm2",
                                      "The green lagrange strain must have a length of 9" );
            }
    
            if ( microStrain.size() != dim * dim ){
                return new errorNode( "computeLinearElasticTerm2",
                                      "The micro-strain must have a length of 9" );
            }
    
            if ( invCPsi.size() != dim * dim ){
                return new errorNode( "computeLinearElasticTerm2",
                                      "invCPsi must have a size of 9" );
            }
    
            if ( B.size() != dim * dim * dim * dim ){
                return new errorNode( "computeLinearElasticTerm2",
                                      "B must have a size of 3**4" );
            }
    
            if ( D.size() != dim * dim * dim * dim ){
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
                                            variableVector &term2, variableMatrix &dTerm2dGreenLagrangeStrain,
                                            variableMatrix &dTerm2dMicroStrain, variableMatrix &dTerm2dInvCPsi ){
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
             * :param const variableVector &greenLagrangeStrain: The Green-Lagrange strain \f$E_{IJ} = \frac{1}{2} \left( C_{IJ} - \delta_{IJ} \right)\f$
             * :param const variableVector &microStrain: The micro-strain \f$\mathcal{E}_{IJ} = \Psi_{IJ} - \delta_{IJ}\f$
             * :param const variableVector &invCPsi: The product \f$C_{JR}^{-1} \Psi_{RQ}\f$
             * :param const variableVector &B: The B stiffness matrix
             * :param const variableVector &D: The D stiffness matrix
             * :param variableVector &term2: The second term.
             * :param variableMatrix &dTerm2dGreenLagrangeStrain: The jacobian of term 2 w.r.t. the Green-Lagrange strain.
             * :param variableMatrix &dTerm2dMicroStrain: The jacobian of term 2 w.r.t. the microStrain.
             * :param variableMatrix &dTerm2dInvCPsi: The jacobian of term 2 w.r.t. \f$C_{JR}^{-1} \Psi_{RQ}\f$
             */
    
            //Assume 3D
            unsigned int dim = 3;
    
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
            dTerm2dGreenLagrangeStrain = variableMatrix( term2.size(), variableVector( greenLagrangeStrain.size(), 0 ) );
            dTerm2dMicroStrain         = variableMatrix( term2.size(), variableVector( microStrain.size(), 0 ) );
            dTerm2dInvCPsi             = variableMatrix( term2.size(), variableVector( invCPsi.size(), 0 ) );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int M = 0; M < dim; M++ ){
                        for ( unsigned int N = 0; N < dim; N++ ){
                            for ( unsigned int K = 0; K < dim; K++ ){
                                dTerm2dGreenLagrangeStrain[ dim * I + J ][ dim * M + N ] += D[ dim * dim * dim * M + dim * dim * N + dim * I + K] * invCPsi[ dim * J + K ];
                                dTerm2dMicroStrain[ dim * I + J ][ dim * M + N ] += B[ dim * dim * dim * I + dim * dim * K + dim * M + N] * invCPsi[ dim * J + K ];
                                for ( unsigned int L = 0; L < dim; L++ ){
                                    dTerm2dInvCPsi[ dim * I + J ][ dim * M + N ] += ( B[ dim * dim * dim * I + dim * dim * N + dim * K + L ] * microStrain[ dim * K + L ] + greenLagrangeStrain[ dim * K + L ] * D[ dim * dim * dim * K + dim * dim * L + dim * I + N ] ) * eye[ dim * J + M ];
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
             * :param const variableVector &Gamma: The micro-gradient deformation measure.
             * :param const parameterVector &C: The C stiffness tensor.
             * :param variableVector &referenceHigherOrderStress: The higher order stress in the reference 
             *     configuration.
             */
    
            //Assume 3D
            unsigned int dim = 3;
    
            if ( Gamma.size() != dim * dim * dim ){
                return new errorNode( "computeReferenceHigherOrderStress",
                                      "Gamma must have a length of 27" );
            }
    
            if ( C.size() != dim * dim * dim * dim * dim * dim ){
                return new errorNode( "computeReferenceHigherOrderStress",
                                      "The C stiffness tensor have a length of 3**6" );
            }
    
            referenceHigherOrderStress = variableVector( dim * dim * dim, 0 );
    
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
                                                    variableMatrix &dReferenceHigherOrderStressdGamma ){
            /*!
             * Compute the higher order stress in the reference configuration.
             * \f$M_{IJK} = C_{JKILMN} Gamma_{LMN}\f$
             *
             * Also compute the Jacobian
             * \f$\frac{ \partial M_{IJK} }{\partial \Gamma_{OPQ} } = C_{JKIOPQ}\f$
             *
             * :param const variableVector &Gamma: The micro-gradient deformation measure.
             * :param const parameterVector &C: The C stiffness tensor.
             * :param variableVector &referenceHigherOrderStress: The higher order stress in the reference 
             *     configuration.
             */
    
            //Assume 3D
            unsigned int dim = 3;
    
            errorOut error = computeReferenceHigherOrderStress( Gamma, C, referenceHigherOrderStress );
    
            if ( error ){
                errorOut result = new errorNode( "computeReferenceHigherOrderStress (jacobian)",
                                                 "Error in computation of higher order stress" );
                result->addNext( error );
                return result;
            }
    
            //Assemble the Jacobian
            dReferenceHigherOrderStressdGamma = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int O = 0; O < dim; O++ ){
                            for ( unsigned int P = 0; P < dim; P++ ){
                                for ( unsigned int Q = 0; Q < dim; Q++ ){
                                    dReferenceHigherOrderStressdGamma[ dim * dim * I + dim * J + K ][ dim * dim * O + dim * P + Q ] +=
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
             * :param const variableVector &invCGamma: \f$ C_{JS}^{-1} \Gamma_{SQR} \f$
             * :param const variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
             * :param variableVector &term3: The third term in the linear elastic equation.
             */
    
            //Assume 3D
            unsigned int dim = 3;
    
            if ( invCGamma.size() != dim * dim * dim ){
                return new errorNode( "computeLinearElasticTerm3",
                                      "invCGamma must have a size of 27" );
            }
    
            if ( referenceHigherOrderStress.size() != dim * dim * dim ){
                return new errorNode( "computeLinearElasticTerm3",
                                      "The referenceHigherOrder stress must have a size of 27" );
            }
    
            term3 = variableVector( dim * dim, 0 );
    
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
                                            variableMatrix &dTerm3dInvCGamma, variableMatrix &dTerm3dReferenceHigherOrderStress ){
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
             * :param const variableVector &invCGamma: \f$ C_{JS}^{-1} \Gamma_{SQR} \f$
             * :param const variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
             * :param variableVector &term3: The third term in the linear elastic equation.
             */
    
            //Assume 3D
            unsigned int dim = 3;
    
            errorOut error = computeLinearElasticTerm3( invCGamma, referenceHigherOrderStress, term3 );
    
            if ( error ){
                errorOut result = new errorNode( "computeLinearElasticTerm3 (jacobian)",
                                                 "Error in computation of term 3" );
                result->addNext( error );
                return result;
            }
    
            constantVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );
    
            dTerm3dInvCGamma = variableMatrix( dim * dim, variableVector( dim * dim * dim, 0 ) );
            dTerm3dReferenceHigherOrderStress = variableMatrix( dim * dim, variableVector( dim * dim * dim, 0 ) );
    
            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int T = 0; T < dim; T++ ){
                        for ( unsigned int U = 0; U < dim; U++ ){
                            for ( unsigned int V = 0; V < dim; V++ ){
                                dTerm3dInvCGamma[ dim * I + J ][ dim * dim * T + dim * U + V ] = referenceHigherOrderStress[ dim * dim * I + dim * U + V ] * eye[ dim * J + T ];
                                dTerm3dReferenceHigherOrderStress[ dim * I + J ][ dim * dim * T + dim * U + V ] = eye[ dim * I + T ] * invCGamma[ dim * dim * J + dim * U + V ];
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
             * :param const variableVector &invRCG: The inverse of the right cauchy green deformation tensor.
             * :param const variableVector &Psi: The micro-deformation measure.
             * :param variableVector &invRCGPsi: the product.
             */
    
            //Assume 3d
            unsigned int dim = 3;
    
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
                                   variableMatrix &dInvRCGPsidRCG, variableMatrix &dInvRCGPsidPsi ){
            /*!
             * Compute the product \f$ C_{IK}^{-1} \Psi_{KJ} \f$
             *
             * Also compute the Jacobians
             * \f$\begin{align}
             * \frac{ \partial C_{IO}^{-1} \Psi_{OJ} } { \partial C_{KL} } &= -C_{IK}^{-1} C_{LO}^{-1} \Psi_{OJ}\\
             * \frac{ \partial C_{IO}^{-1} \Psi_{OJ} } { \partial \Psi_{KL} } &= C_{IK}^{-1} \delta_{JL}
             * \end{align}\f$
             *
             * :param const variableVector &invRCG: The inverse of the right cauchy green deformation tensor.
             * :param const variableVector &Psi: The micro-deformation measure.
             * :param variableVector &invRCGPsi: the product.
             * :param variableMatrix &dInvRCGPsidRCG: The Jacobian of the product w.r.t. the right cauchy green
             *     deformation tensor.
             * :param variableMatrix &dInvRCGPsidPsi: The Jacobian of the product w.r.t. the micro-deformation measure.
             */
    
            //Assume 3D
            unsigned int dim = 3;
    
            errorOut error = computeInvRCGPsi( invRCG, Psi, invRCGPsi );
    
            if ( error ){
                errorOut result = new errorNode( "computeInvRCGPsi (jacobian)", "Error in computation of invRCG Psi product" );
                result->addNext( error );
                return result;
            }
    
            //Construct the jacobians
            variableVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );
    
            dInvRCGPsidRCG = variableMatrix( invRCGPsi.size(), variableVector( invRCG.size(), 0 ) );
            dInvRCGPsidPsi = variableMatrix( invRCGPsi.size(), variableVector( Psi.size(), 0 ) );
    
            for ( unsigned int I = 0; I < 3; I++ ){
                for ( unsigned int J = 0; J < 3; J++ ){
                    for ( unsigned int K = 0; K < 3; K++ ){
                        for ( unsigned int L = 0; L < 3; L++ ){
                            dInvRCGPsidRCG[ dim * I + J ][ dim * K + L ] = -invRCG[ dim * I + K ] * invRCGPsi[ dim * L + J ];
                            dInvRCGPsidPsi[ dim * I + J ][ dim * K + L ] = invRCG[ dim * I + K ] * eye[ dim * J + L ];
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
             * :param const variableVector &invRCG: The inverse of the right Cauchy Green deformation tensor.
             * :param const variableVector &Gamma: The gradient of the micro-deformation deformation tensor.
             * :param variableVector &invRCGGamma: The product.
             */
    
            //Assume 3d
            unsigned int dim = 3;
    
            if ( invRCG.size() != dim * dim ){
                return new errorNode( "computeInvRCGGamma", "invRCG has an improper dimension" );
            }
    
            if ( Gamma.size() != dim * dim * dim ){
                return new errorNode( "computeInvRCGGamma", "Gamma has an improper dimension" );
            }
    
            invRCGGamma = variableVector( dim * dim * dim, 0 );
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
                                     variableMatrix &dInvRCGGammadRCG, variableMatrix &dInvRCGGammadGamma ){
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
             * :param const variableVector &invRCG: The inverse of the right Cauchy Green deformation tensor.
             * :param const variableVector &Gamma: The gradient of the micro-deformation deformation tensor.
             * :param variableVector &invRCGGamma: The product.
             */
    
            //Assume 3d
            unsigned int dim = 3;
    
            errorOut error = computeInvRCGGamma( invRCG, Gamma, invRCGGamma );
    
            if ( error ){
                errorOut result = new errorNode( "computeInvRCGGamma (jacobian)", "Error in computation of invRCG Gamma product" );
                result->addNext( error );
                return result;
            }
    
            //Assemble jacobians of invCGamma w.r.t. C and Gamma
            variableVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );
    
            dInvRCGGammadRCG = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
            dInvRCGGammadGamma = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );
    
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int Q = 0; Q < dim; Q++ ){
                    for ( unsigned int R = 0; R < dim; R++ ){
                        for ( unsigned int T = 0; T < dim; T++ ){
                            for ( unsigned int U = 0; U < dim; U++ ){
                                dInvRCGGammadRCG[ dim * dim * J + dim * Q + R ][ dim * T + U ]
                                    = -invRCG[ dim * J + T] * invRCGGamma[ dim * dim * U + dim * Q + R ];
                                for ( unsigned int V = 0; V < dim; V++ ){
                                    dInvRCGGammadGamma[ dim * dim * J + dim * Q + R ][ dim * dim * T + dim * U + V]
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
             * :param const parameterType &lambda: The micromorphic lambda parameter.
             * :param const parameterType &mu: The micromorphic mu parameter.
             * :param parameterVector &A: The isotropic A stiffness tensor.
             */
    
            //Assume 3D
            unsigned int dim = 3;
    
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
             * :param const parameterType &eta: The micromorphic eta parameter.
             * :param const parameterType &tau: The micromorphic tau parameter.
             * :param const parameterType &kappa: The micromorphic kappa parameter.
             * :param const parameterType &nu: The micromorphic nu parameter.
             * :param const parameterType &sigma: The micromorphic sigma parameter
             * :param parameterVector &B: The isotropic B stiffnes tensor.
             */
    
            //Assume 3D
            unsigned int dim = 3;
    
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
             * :param const parameterVector &taus: The moduli (11 independent terms)
             * :param parameterVector &C: The isotropic C stiffness tensor.
             */
    
            //Assume 3D
            unsigned int dim = 3;
    
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
             * :param const parameterType &tau: The micromorphic tau parameter.
             * :param const parameterType &sigma: The micromorphic sigma parameter.
             * :param parameterVector &D: The D stiffness tensor.
             */
    
            //Assume 3D
            unsigned int dim = 3;
    
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
             * :param const double ( &grad_u )[ 3 ][ 3 ]: The macro displacement gradient w.r.t. the reference configuration.
             * :param const double ( &phi )[ 9 ]: The micro displacement.
             * :param const double ( &grad_phi )[ 9 ][ 3 ]: The gradient of the micro displacement w.r.t. the reference configuration.
             * :param variableVector &deformationGradient: The deformation gradient
             * :param variableVector &microDeformation: The micro deformation
             * :param variableVector &gradientMicroDeformation: The gradient of the micro deformation.
             */
    
    
            //Extract the degrees of freedom
            variableMatrix displacementGradient = { { grad_u[ 0 ][ 0 ], grad_u[ 0 ][ 1 ], grad_u[ 0 ][ 2 ] },
                                                    { grad_u[ 1 ][ 0 ], grad_u[ 1 ][ 1 ], grad_u[ 1 ][ 2 ] },
                                                    { grad_u[ 2 ][ 0 ], grad_u[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] } };
    
            variableVector microDisplacement = { phi[ 0 ], phi[ 1 ], phi[ 2 ],
                                                 phi[ 3 ], phi[ 4 ], phi[ 5 ],
                                                 phi[ 6 ], phi[ 7 ], phi[ 8 ] };
    
            variableMatrix gradientMicroDisplacement = { { grad_phi[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ] },
                                                         { grad_phi[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ] },
                                                         { grad_phi[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ] },
                                                         { grad_phi[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ] },
                                                         { grad_phi[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ] },
                                                         { grad_phi[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ] },
                                                         { grad_phi[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ] },
                                                         { grad_phi[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ] },
                                                         { grad_phi[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] } };
    
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
                                                         variableVector &gradientMicroDeformation, variableMatrix &dFdGradU,
                                                         variableMatrix &dChidPhi, variableMatrix &dGradChidGradPhi ){
            /*!
             * Assemble the fundamental deformation meaures from the degrees of freedom.
             *
             * :param const double ( &grad_u )[ 3 ][ 3 ]: The macro displacement gradient w.r.t. the reference configuration.
             * :param const double ( &phi )[ 9 ]: The micro displacement.
             * :param const double ( &grad_phi )[ 9 ][ 3 ]: The gradient of the micro displacement w.r.t. the reference configuration.
             * :param variableVector &deformationGradient: The deformation gradient
             * :param variableVector &microDeformation: The micro deformation
             * :param variableVector &gradientMicroDeformation: The gradient of the micro deformation.
             * :param variableMatrix &dFdGradU: The Jacobian of the deformation gradient w.r.t. the gradient of the displacement
             * :param variableMatrix &dChidPhi: The Jacobian of the micro deformation w.r.t. the micro displacement
             * :param variableMatrix &dGradChidGradPhi: The Jacobian of the gradient of the micro deformation w.r.t.
             *      the gradient of the micro displacement
             */
    
    
            //Extract the degrees of freedom
            variableMatrix displacementGradient = { { grad_u[ 0 ][ 0 ], grad_u[ 0 ][ 1 ], grad_u[ 0 ][ 2 ] },
                                                    { grad_u[ 1 ][ 0 ], grad_u[ 1 ][ 1 ], grad_u[ 1 ][ 2 ] },
                                                    { grad_u[ 2 ][ 0 ], grad_u[ 2 ][ 1 ], grad_u[ 2 ][ 2 ] } };
    
            variableVector microDisplacement = { phi[ 0 ], phi[ 1 ], phi[ 2 ],
                                                 phi[ 3 ], phi[ 4 ], phi[ 5 ],
                                                 phi[ 6 ], phi[ 7 ], phi[ 8 ] };
    
            variableMatrix gradientMicroDisplacement = { { grad_phi[ 0 ][ 0 ], grad_phi[ 0 ][ 1 ], grad_phi[ 0 ][ 2 ] },
                                                         { grad_phi[ 1 ][ 0 ], grad_phi[ 1 ][ 1 ], grad_phi[ 1 ][ 2 ] },
                                                         { grad_phi[ 2 ][ 0 ], grad_phi[ 2 ][ 1 ], grad_phi[ 2 ][ 2 ] },
                                                         { grad_phi[ 3 ][ 0 ], grad_phi[ 3 ][ 1 ], grad_phi[ 3 ][ 2 ] },
                                                         { grad_phi[ 4 ][ 0 ], grad_phi[ 4 ][ 1 ], grad_phi[ 4 ][ 2 ] },
                                                         { grad_phi[ 5 ][ 0 ], grad_phi[ 5 ][ 1 ], grad_phi[ 5 ][ 2 ] },
                                                         { grad_phi[ 6 ][ 0 ], grad_phi[ 6 ][ 1 ], grad_phi[ 6 ][ 2 ] },
                                                         { grad_phi[ 7 ][ 0 ], grad_phi[ 7 ][ 1 ], grad_phi[ 7 ][ 2 ] },
                                                         { grad_phi[ 8 ][ 0 ], grad_phi[ 8 ][ 1 ], grad_phi[ 8 ][ 2 ] } };
    
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

        void residual::setRightCauchyGreen( const variableVector &rightCauchyGreen ){
            /*!
             * Set the value of the right Cauchy-Green deformation tensor
             * 
             * \param &rightCauchyGreen: The value
             */

            _rightCauchyGreen.second = rightCauchyGreen;

            _rightCauchyGreen.first = true;

            addIterationData( &_rightCauchyGreen );

        }

        void residual::setPsi( const variableVector &psi ){
            /*!
             * Set the value of the micro deformation measure psi
             * 
             * \param &psi: The value
             */

            _psi.second = psi;

            _psi.first = true;

            addIterationData( &_psi );

        }

        void residual::setGamma( const variableVector &gamma ){
            /*!
             * Set the value of the micro deformation measure gamma
             * 
             * \param &gamma: The value
             */

            _gamma.second = gamma;

            _gamma.first = true;

            addIterationData( &_gamma );

        }

        void residual::setPreviousRightCauchyGreen( const variableVector &previousRightCauchyGreen ){
            /*!
             * Set the value of the previous right Cauchy-Green deformation tensor
             * 
             * \param &previousRightCauchyGreen: The value
             */

            _previousRightCauchyGreen.second = previousRightCauchyGreen;

            _previousRightCauchyGreen.first = true;

        }

        void residual::setPreviousPsi( const variableVector &previousPsi ){
            /*!
             * Set the value of the previous micro deformation measure psi
             * 
             * \param &previousPsi: The value
             */

            _previousPsi.second = previousPsi;

            _previousPsi.first = true;
        }

        void residual::setPreviousGamma( const variableVector &previousGamma ){
            /*!
             * Set the value of the previous micro deformation measure gamma
             * 
             * \param &previousGamma: The value
             */

            _previousGamma.second = previousGamma;

            _previousGamma.first = true;

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
             * We assume that the first configuration in hydra.getConfigurations is the elastic one
             *
             * \param isPrevious: Flag for whether the measures to be calculated are in the current or previous configuration
             */

            floatVector deformationGradient1;

            floatVector microDeformation1;

            floatVector gradientMicroDeformation1;

            if ( isPrevious ){

                deformationGradient1 = ( *hydra->getPreviousConfigurations( ) )[ 0 ];

                microDeformation1 = ( *hydra->getPreviousMicroConfigurations( ) )[ 0 ];

                gradientMicroDeformation1 = ( *hydra->getPreviousGradientMicroConfigurations( ) )[ 0 ];

            }
            else{

                deformationGradient1 = ( *hydra->getConfigurations( ) )[ 0 ];

                microDeformation1 = ( *hydra->getMicroConfigurations( ) )[ 0 ];

                gradientMicroDeformation1 = ( *hydra->getGradientMicroConfigurations( ) )[ 0 ];

            }

            floatVector rightCauchyGreen;

            floatVector Psi;

            floatVector Gamma;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( computeDeformationMeasures( deformationGradient1, microDeformation1, gradientMicroDeformation1,
                                                                                   rightCauchyGreen, Psi, Gamma ) );

            if ( isPrevious ){

                setPreviousRightCauchyGreen( rightCauchyGreen );

                setPreviousPsi( Psi );

                setPreviousGamma( Gamma );

            }
            else{

                setRightCauchyGreen( rightCauchyGreen );

                setPsi( Psi );

                setGamma( Gamma );

            }

        }

        void residual::setdRightCauchyGreendF( const floatMatrix &value ){
            /*!
             * Set the derivative of the right Cauchy-Green deformation measure w.r.t. the total deformation gradient
             * 
             * \param &value: The value of the Jacobian
             */

            _dRightCauchyGreendF.second = value;

            _dRightCauchyGreendF.first = true;

            addIterationData( &_dRightCauchyGreendF );

        }

        void residual::setdRightCauchyGreendFn( const floatMatrix &value ){
            /*!
             * Set the derivative of the right Cauchy-Green deformation measure w.r.t. the remaining sub-deformation gradients
             * 
             * \param &value: The value of the Jacobian
             */

            _dRightCauchyGreendFn.second = value;

            _dRightCauchyGreendFn.first = true;

            addIterationData( &_dRightCauchyGreendFn );

        }

        void residual::setdPsidF( const floatMatrix &value ){
            /*!
             * Set the derivative of the micro deformation measure w.r.t. the total deformation gradient
             * 
             * \param &value: The value of the Jacobian
             */

            _dPsidF.second = value;

            _dPsidF.first = true;

            addIterationData( &_dPsidF );

        }

        void residual::setdPsidFn( const floatMatrix &value ){
            /*!
             * Set the derivative of the micro deformation measure w.r.t. the remaining sub-deformation gradients
             * 
             * \param &value: The value of the Jacobian
             */

            _dPsidFn.second = value;

            _dPsidFn.first = true;

            addIterationData( &_dPsidFn );

        }

        void residual::setdPsidChi( const floatMatrix &value ){
            /*!
             * Set the derivative of the micro deformation measure w.r.t. the total micro deformation
             * 
             * \param &value: The value of the Jacobian
             */

            _dPsidChi.second = value;

            _dPsidChi.first = true;

            addIterationData( &_dPsidChi );

        }

        void residual::setdPsidChin( const floatMatrix &value ){
            /*!
             * Set the derivative of the micro deformation measure w.r.t. the remaining sub-micro deformation
             * 
             * \param &value: The value of the Jacobian
             */

            _dPsidChin.second = value;

            _dPsidChin.first = true;

            addIterationData( &_dPsidChin );

        }

        void residual::setdGammadF( const floatMatrix &value ){
            /*!
             * Set the derivative of the micro deformation measure gamma w.r.t. the total deformation gradient
             * 
             * \param &value: The value of the Jacobian
             */

            _dGammadF.second = value;

            _dGammadF.first = true;

            addIterationData( &_dGammadF );

        }

        void residual::setdGammadFn( const floatMatrix &value ){
            /*!
             * Set the derivative of the micro deformation measure gamma w.r.t. the remaining sub-deformation gradients
             * 
             * \param &value: The value of the Jacobian
             */

            _dGammadFn.second = value;

            _dGammadFn.first = true;

            addIterationData( &_dGammadFn );

        }

        void residual::setdGammadChi( const floatMatrix &value ){
            /*!
             * Set the derivative of the micro deformation measure gamma w.r.t. the total micro deformation
             * 
             * \param &value: The value of the Jacobian
             */

            _dGammadChi.second = value;

            _dGammadChi.first = true;

            addIterationData( &_dGammadChi );

        }

        void residual::setdGammadChin( const floatMatrix &value ){
            /*!
             * Set the derivative of the micro deformation measure gamma w.r.t. the remaining sub-micro deformation
             * 
             * \param &value: The value of the Jacobian
             */

            _dGammadChin.second = value;

            _dGammadChin.first = true;

            addIterationData( &_dGammadChin );

        }

        void residual::setdGammadGradChi( const floatMatrix &value ){
            /*!
             * Set the derivative of the micro deformation measure gamma w.r.t. the reference spatial gradient of the total micro deformation
             * 
             * \param &value: The value of the Jacobian
             */

            _dGammadGradChi.second = value;

            _dGammadGradChi.first = true;

            addIterationData( &_dGammadGradChi );

        }

        void residual::setdGammadGradChin( const floatMatrix &value ){
            /*!
             * Set the derivative of the micro deformation measure gamma w.r.t. the local reference spatial gradient of the remaining sub-micro deformation
             * 
             * \param &value: The value of the Jacobian
             */

            _dGammadGradChin.second = value;

            _dGammadGradChin.first = true;

            addIterationData( &_dGammadGradChin );

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

        const variableMatrix *residual::getdRightCauchyGreendF( ){
            /*!
             * Get the jacobian of the right Cauchy-Green deformation tensor w.r.t. the total deformation gradient
             */

            if ( !_dRightCauchyGreendF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdRightCauchyGreendF( ) );

            }

            return &_dRightCauchyGreendF.second;

        }

        const variableMatrix *residual::getdRightCauchyGreendFn( ){
            /*!
             * Get the jacobian of the right Cauchy-Green deformation tensor w.r.t. the remaining sub-deformation gradients
             */

            if ( !_dRightCauchyGreendFn.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdRightCauchyGreendFn( ) );

            }

            return &_dRightCauchyGreendFn.second;

        }

        const variableMatrix *residual::getdPsidF( ){
            /*!
             * Get the jacobian of the micro deformation measure psi w.r.t. the total deformation gradient
             */

            if ( !_dPsidF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPsidF( ) );

            }

            return &_dPsidF.second;

        }

        const variableMatrix *residual::getdPsidFn( ){
            /*!
             * Get the jacobian of the micro deformation tensor psi w.r.t. the remaining sub-deformation gradients
             */

            if ( !_dPsidFn.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPsidFn( ) );

            }

            return &_dPsidFn.second;

        }

        const variableMatrix *residual::getdPsidChi( ){
            /*!
             * Get the jacobian of the micro deformation measure psi w.r.t. the total micro-deformation
             */

            if ( !_dPsidChi.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPsidChi( ) );

            }

            return &_dPsidChi.second;

        }

        const variableMatrix *residual::getdPsidChin( ){
            /*!
             * Get the jacobian of the micro deformation tensor psi w.r.t. the remaining sub-micro deformations
             */

            if ( !_dPsidChin.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdPsidChin( ) );

            }

            return &_dPsidChin.second;

        }

        const variableMatrix *residual::getdGammadF( ){
            /*!
             * Get the jacobian of the micro deformation measure gamma w.r.t. the total deformation gradient
             */

            if ( !_dGammadF.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdGammadF( ) );

            }

            return &_dGammadF.second;

        }

        const variableMatrix *residual::getdGammadFn( ){
            /*!
             * Get the jacobian of the micro deformation tensor gamma w.r.t. the remaining sub-deformation gradients
             */

            if ( !_dGammadFn.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdGammadFn( ) );

            }

            return &_dGammadFn.second;

        }

        const variableMatrix *residual::getdGammadChi( ){
            /*!
             * Get the jacobian of the micro deformation measure gamma w.r.t. the total micro-deformation
             */

            if ( !_dGammadChi.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdGammadChi( ) );

            }

            return &_dGammadChi.second;

        }

        const variableMatrix *residual::getdGammadChin( ){
            /*!
             * Get the jacobian of the micro deformation tensor gamma w.r.t. the remaining sub-micro deformations
             */

            if ( !_dGammadChin.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdGammadChin( ) );

            }

            return &_dGammadChin.second;

        }

        const variableMatrix *residual::getdGammadGradChi( ){
            /*!
             * Get the jacobian of the micro deformation measure gamma w.r.t. the reference spatial gradient of the total micro-deformation
             */

            if ( !_dGammadGradChi.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdGammadGradChi( ) );

            }

            return &_dGammadGradChi.second;

        }

        const variableMatrix *residual::getdGammadGradChin( ){
            /*!
             * Get the jacobian of the micro deformation tensor gamma w.r.t. the local reference spatial gradient of the remaining sub-micro deformations
             */

            if ( !_dGammadGradChin.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setdGammadGradChin( ) );

            }

            return &_dGammadGradChin.second;

        }

        void residual::setDeformationJacobians( const bool isPrevious ){
            /*!
             * Evaluate the derived deformation Jacobians
             * 
             * We assume that the first configuration in hydra.getConfigurations is the elastic one
             *
             * \param isPrevious: Flag for whether the measures to be calculated are in the current or previous configuration
             */

            floatVector deformationGradient1;

            floatVector microDeformation1;

            floatVector gradientMicroDeformation1;

            if ( isPrevious ){

                deformationGradient1 = ( *hydra->getPreviousConfigurations( ) )[ 0 ];

                microDeformation1 = ( *hydra->getPreviousMicroConfigurations( ) )[ 0 ];

                gradientMicroDeformation1 = ( *hydra->getPreviousGradientMicroConfigurations( ) )[ 0 ];

            }
            else{

                deformationGradient1 = ( *hydra->getConfigurations( ) )[ 0 ];

                microDeformation1 = ( *hydra->getMicroConfigurations( ) )[ 0 ];

                gradientMicroDeformation1 = ( *hydra->getGradientMicroConfigurations( ) )[ 0 ];

            }

            floatVector rightCauchyGreen;

            floatVector Psi;

            floatVector Gamma;

            floatMatrix dCdF1;

            floatMatrix dPsidF1;

            floatMatrix dPsidChi1;

            floatMatrix dGammadF1;

            floatMatrix dGammadGradChi1;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( computeDeformationMeasures( deformationGradient1, microDeformation1, gradientMicroDeformation1,
                                                                                   rightCauchyGreen, Psi, Gamma,
                                                                                   dCdF1, dPsidF1, dPsidChi1, dGammadF1, dGammadGradChi1 ) );

            if ( isPrevious ){

                setPreviousRightCauchyGreen( rightCauchyGreen );

                setPreviousPsi( Psi );

                setPreviousGamma( Gamma );

            }
            else{

                setRightCauchyGreen( rightCauchyGreen );

                setPsi( Psi );

                setGamma( Gamma );

            }

        }

        const variableVector *residual::getRightCauchyGreen( ){
            /*!
             * Set the value of the right Cauchy-Green deformation tensor
             */

            if ( !_rightCauchyGreen.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setDeformation( false ) );

            }

            return &_rightCauchyGreen.second;

        }

        const variableVector *residual::getPsi( ){
            /*!
             * Set the value of the micro deformation measure Psi
             */

            if ( !_psi.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setDeformation( false ) );

            }

            return &_psi.second;

        }


        const variableVector *residual::getGamma( ){
            /*!
             * Set the value of the micro deformation measure Gamma
             */

            if ( !_gamma.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setDeformation( false ) );

            }

            return &_gamma.second;

        }

        const variableVector *residual::getPreviousRightCauchyGreen( ){
            /*!
             * Set the previous value of the right Cauchy-Green deformation tensor
             */

            if ( !_previousRightCauchyGreen.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setDeformation( true ) );

            }

            return &_previousRightCauchyGreen.second;

        }

        const variableVector *residual::getPreviousPsi( ){
            /*!
             * Set the value of the previous micro deformation measure Psi
             */

            if ( !_previousPsi.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setDeformation( true ) );

            }

            return &_previousPsi.second;

        }


        const variableVector *residual::getPreviousGamma( ){
            /*!
             * Set the value of the previous micro deformation measure Gamma
             */

            if ( !_previousGamma.first ){

                TARDIGRADE_ERROR_TOOLS_CATCH( setDeformation( true ) );

            }

            return &_previousGamma.second;

        }

    }

}
