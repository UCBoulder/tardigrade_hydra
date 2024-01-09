/**
  ******************************************************************************
  * \file tardigrade_hydraMicromorphicDruckerPragerPlasticity.cpp
  ******************************************************************************
  * An implementation of micromorphic linear elasticity using the hydra
  * framework. Used as an example and as the basis for more complex models.
  ******************************************************************************
  */

#include<tardigrade_hydraMicromorphicDruckerPragerPlasticity.h>
#include<tardigrade_constitutive_tools.h>
#include<tardigrade_micromorphic_tools.h>

namespace tardigradeHydra{

    namespace micromorphicDruckerPragerPlasticity{

        void computeDruckerPragerInternalParameters( const parameterType &frictionAngle, const parameterType &beta,
                                                     parameterType &A, parameterType &B ){
            /*!
             * Compute the Drucker-Prager internal parameters
             *
             * :param const parameterType &frictionAngle: The material friction angle ( 0 < frictionAngle < pi / 2 );
             * :param const parameterType &beta: The beta parameter.
             * :param parameterType &A: The A parameter.
             * :param parameterType &B: The B parameter.
             */
    
            //Make sure the parameters are within bounds
            TARDIGRADE_ERROR_TOOLS_CATCH( 
                if ( ( 0 > frictionAngle ) || ( frictionAngle > 1.570796 ) ){
                    throw std::runtime_error( "The friction angle must be betwen 0 and pi / 2 not " + std::to_string( frictionAngle ) );
                }
            )
    
            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( abs( beta ) > 1 ){
                    throw std::runtime_error( "Beta must be between -1 and 1 not " + std::to_string( beta ) );
                }
            )
    
            //Compute the parameters
            parameterType betaAngle = 2. * std::sqrt(6.) / ( 3. + beta * std::sin( frictionAngle ) );
    
            A = betaAngle * std::cos( frictionAngle );
    
            B = betaAngle * std::sin( frictionAngle );
    
        }

        void computeSecondOrderDruckerPragerYieldEquation( const variableVector &stressMeasure, const variableType &cohesion,
                                                               const variableVector &precedingDeformationGradient,
                                                               const parameterType &frictionAngle, const parameterType &beta,
                                                               variableType &yieldValue ){
            /*!
             * Compute the second-order Drucker Prager Yield equation
             *
             * \f$ F = ||dev ( stressMeasure ) || - \left( A^{\phi} \bar{c} - B^{\phi} \bar{p} \right) \leq 0
             * 
             * || dev ( stressMeasure ) || = \sqrt{ dev( referenceStressMeasure ) : dev( referenceStressMeasure ) }
             *  dev( referenceStressMeasure ) : dev( referenceStressMeasure ) = dev( referenceStressMeasure )_{IJ} dev( referenceStressMeasure )_{IJ}
             *  dev( referenceStressMeasure )_{IJ} = referenceStressMeasure_{IJ} - \bar{p} elasticRightCauchyGreen_{IJ}^{-1}
             *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} referenceStressMeasure_{IJ}
             *  A^{angle} = \beta^{angle} \cos( frictionAngle )
             *  B^{angle} = \beta^{angle} \sin( frictionAngle )
             *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }
             * \f$
             *
             * :param const variableVector &stressMeasure: The stress measure
             * :param const variableType &cohesion: The cohesion measure.
             * :param const variableVector &precedingDeformationGradient: The deformation gradients preceding the configuration of the stress measure
             * :param const parameterType &frictionAngle: The friction angle
             * :param const parameterType &beta: The beta parameter
             * :param variableType &yieldValue: The yield value.
             */

            parameterType AAngle, BAngle;
            TARDIGRADE_ERROR_TOOLS_CATCH( tardigradeHydra::micromorphicDruckerPragerPlasticity::computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle ) );

            //Compute the right Cauchy-Green deformation tensor
            floatVector rightCauchyGreen;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen ) );

            //Compute the decomposition of the stress
            variableType pressure;
            variableVector deviatoricReferenceStress;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( stressMeasure,
                                                                                                                                    rightCauchyGreen, deviatoricReferenceStress, pressure ) );

            //Compute the l2norm of the deviatoric stress
            variableType normDevStress = tardigradeVectorTools::l2norm( deviatoricReferenceStress );

            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        }

        void computeSecondOrderDruckerPragerYieldEquation( const variableVector &referenceStressMeasure, const variableType &cohesion,
                                                           const variableVector &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, variableVector &dFdStress, variableType &dFdc,
                                                           variableVector &dFdPrecedingF, double tol ){
            /*!
             * Compute the second-order Drucker Prager Yield equation
             *
             * \f$F = ||dev ( stressMeasure ) || - \left( A^{\phi} \bar{c} - B^{\phi} \bar{p} \right) \leq 0
             * 
             * || dev ( stressMeasure ) || = \sqrt{ dev( referenceStressMeasure ) : dev( referenceStressMeasure ) }
             *  dev( referenceStressMeasure ) : dev( referenceStressMeasure ) = dev( referenceStressMeasure )_{IJ} dev( referenceStressMeasure )_{IJ}
             *  dev( referenceStressMeasure )_{IJ} = referenceStressMeasure_{IJ} - \bar{p} elasticRightCauchyGreen_{IJ}^{-1}\f$
             *
             *  Also compute the Jacobians
             * \f$\frac{ \partial F }{ \partial stressMeasure_{IJ} } = \frac{ dev ( stressMeasure )_{AB} }{ || dev ( stressMeasure ) || } \frac{ \partial dev( stressMeasure ) \frac{ \partial dev( stressMeasure )_{AB} }{ \partial stressMeasure_{IJ} } + B^{\phi} \frac{ \partial \bar{p} }{ \partial stressMeasure_{IJ} }
             * \frac{ \partial F }{ \partial \bar{c} } = -A^{\phi}
             * \frac{ \partial F }{ \partial C_{IJ} } = \frac{ dev ( stressMeasure )_{AB} }{ || dev ( stressMeasure ) || } \frac{ \partial dev( stressMeasure ) \frac{ \partial dev( stressMeasure )_{AB} }{ \partial C_{IJ} } + B^{\phi} \frac{ \partial \bar{p} }{ \partial C_{IJ} }
             *
             *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} referenceStressMeasure_{IJ}
             *  A^{angle} = \beta^{angle} \cos( frictionAngle )
             *  B^{angle} = \beta^{angle} \sin( frictionAngle )
             *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }\f$
             *
             * :param const variableVector &referenceStressMeasure: The stress measure in the reference configuration
             * :param const variableType &cohesion: The cohesion measure.
             * :param const variableVector &precedingDeformationGradient: The preceding deformation gradient.
             * :param const parameterType &frictionAngle: The friction angle
             * :param const parameterType &beta: The beta parameter
             * :param variableType &yieldValue: The yield value.
             * :param variableVector &dFdStress: The Jacobian of the yield surface w.r.t. the stress measure.
             * :param variableType &dFdc: The Jacobian of the yield surface w.r.t. the cohesion.
             * :param variableVector &dFdPrecedingF: The Jacobian of the yield surface w.r.t. the preceding deformation gradient from the stress-measure's configuration
             * :param double tol: The tolerance used to prevent nans in the Jacobians
             */
    
            parameterType AAngle, BAngle;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle ) );
   
            //Compute the right Cauchy-Green deformation tensor
            floatVector rightCauchyGreen;
            floatMatrix dRCGdPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen, dRCGdPrecedingF ) );
 
            //Compute the decomposition of the stress
            variableType pressure;
            variableVector deviatoricReferenceStress;
    
            variableMatrix dDevStressdStress, dDevStressdRCG;
            variableVector dPressuredStress, dPressuredRCG;
    
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( referenceStressMeasure,
                                                       rightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                                                       dDevStressdRCG, dPressuredStress, dPressuredRCG ) );
   
            variableMatrix dDevStressdPrecedingF = tardigradeVectorTools::dot( dDevStressdRCG, dRCGdPrecedingF );
            variableVector dPressuredPrecedingF  = tardigradeVectorTools::Tdot( dRCGdPrecedingF, dPressuredRCG );
 
            //Compute the l2norm of the deviatoric stress
            variableType normDevStress = tardigradeVectorTools::l2norm( deviatoricReferenceStress );
    
            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );
    
            //Evaluate the jacobians
            variableVector devStressDirection = deviatoricReferenceStress / ( normDevStress + tol );
    
            dFdStress = tardigradeVectorTools::Tdot( dDevStressdStress, devStressDirection )
                      + BAngle * dPressuredStress;
    
            dFdc = - AAngle;
    
            dFdPrecedingF = tardigradeVectorTools::Tdot( dDevStressdPrecedingF, devStressDirection )
                          + BAngle * dPressuredPrecedingF;
        }

        void computeSecondOrderDruckerPragerYieldEquation( const variableVector &stressMeasure, const variableType &cohesion,
                                                           const variableVector &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableType &yieldValue, variableVector &dFdStress, variableType &dFdc,
                                                           variableVector &dFdPrecedingF, variableMatrix &d2FdStress2,
                                                           variableMatrix &d2FdStressdPrecedingF, double tol ){
            /*!
             * Compute the second-order Drucker Prager Yield equation
             *
             * \f$F = ||dev ( stressMeasure ) || - \left( A^{\phi} \bar{c} - B^{\phi} \bar{p} \right) \leq 0
             * 
             * || dev ( stressMeasure ) || = \sqrt{ dev( referenceStressMeasure ) : dev( referenceStressMeasure ) }
             *  dev( referenceStressMeasure ) : dev( referenceStressMeasure ) = dev( referenceStressMeasure )_{IJ} dev( referenceStressMeasure )_{IJ}
             *  dev( referenceStressMeasure )_{IJ} = referenceStressMeasure_{IJ} - \bar{p} elasticRightCauchyGreen_{IJ}^{-1}
             *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} referenceStressMeasure_{IJ}
             *  A^{angle} = \beta^{angle} \cos( frictionAngle )
             *  B^{angle} = \beta^{angle} \sin( frictionAngle )
             *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }\f$
             *
             *  Also compute the Jacobians
             * \f$\frac{ \partial F }{ \partial stressMeasure_{IJ} } = \frac{ dev ( stressMeasure )_{AB} }{ || dev ( stressMeasure ) || } \frac{ \partial dev( stressMeasure ) \frac{ \partial dev( stressMeasure )_{AB} }{ \partial stressMeasure_{IJ} } + B^{\phi} \frac{ \partial \bar{p} }{ \partial stressMeasure_{IJ} }
             * \frac{ \partial F }{ \partial \bar{c} } = -A^{\phi}
             * \frac{ \partial F }{ \partial C_{IJ} } = \frac{ dev ( stressMeasure )_{AB} }{ || dev ( stressMeasure ) || } \frac{ \partial dev( stressMeasure )_{AB} }{ \partial C_{IJ} } + B^{\phi} \frac{ \partial \bar{p} }{ \partial C_{IJ} }\f$
             *
             *  The second deriatives of \f$\frac{ \partial F }{ \partial \stressMeasure_{IJ}  }\f$ are
             *  \f$\frac{ \partial^2 F }{ \partial stressMeasure_{IJ} \partial stressMeasure_{KL} } = \frac{ \partial^2 || dev( stressMeasure ) || }{ \partial dev( stressMeasure )_{AB} \partial dev( stressMeasure )_{CD} } \frac{ \partial dev( stressMeasure )_{AB} } { \partial stressMeasure_{IJ} } \frac{ \partial dev( stressMeasure )_{CD} } { \partial stressMeasure_{KL} } + \frac{ dev ( stressMeasure )_{AB} }{ || dev( stressMeasure ) || } \frac{ \partial^2 dev( stressMeasure )_{AB} }{ \partial stressMeasure_{IJ} \partial stressMeasure_{KL} }
             *  \frac{ \partial^2 F }{ \partial stressMeasure_{IJ} \partial RCG_{KL} } = \frac{ \partial^2 || dev( stressMeasure ) || }{ \partial dev( stressMeasure )_{AB} \partial dev( stressMeasure )_{CD} } \frac{ \partial dev( stressMeasure )_{AB} } { \partial stressMeasure_{IJ} } \frac{ \partial dev( stressMeasure )_{CD} } { \partial C_{KL} } + \frac{ dev ( stressMeasure )_{AB} }{ || dev( stressMeasure ) || } \frac{ \partial^2 dev( stressMeasure )_{AB} }{ \partial stressMeasure_{IJ} \partial C_{KL} } + B^{\phi} \frac{ \partial^2 \bar{p} }{ \partial stressMeasure_{IJ} \partial C_{KL} } \f$
             *  
             * :param const variableVector &stressMeasure: The stress measure
             * :param const variableType &cohesion: The cohesion measure.
             * :param const variableVector &precedingDeformationGradient: The preceding deformation gradient
             * :param const parameterType &frictionAngle: The friction angle
             * :param const parameterType &beta: The beta parameter
             * :param variableType &yieldValue: The yield value.
             * :param variableVector &dFdStress: The Jacobian of the yield surface w.r.t. the stress measure.
             * :param variableType &dFdc: The Jacobian of the yield surface w.r.t. the cohesion.
             * :param variableVector &dFdElasticRCG: The Jacobian of the yield surface w.r.t. the elastic 
             *     right Cauchy-Green deformation tensor.
             * :param variableMatrix &d2FdStress2: The second derivative of the flow direction w.r.t. the stress. This 
             *     is useful if one is using this expression as the flow potential and wants the jacobian of the flow direction \f$\frac{ \partial G }{\partial Stress_{IJ} }\f$
             * :param variableMatrix &d2FdStressdPrecedingF: The second derivative of the flow direction w.r.t. the stress and the preceding deformation gradient
             * :param double tol: The tolerance used to prevent nans in the Jacobians
             */
            //Assume 3D
            unsigned int dim = 3;

            parameterType AAngle, BAngle;
            TARDIGRADE_ERROR_TOOLS_CATCH(  computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle ) );

            //Compute the right Cauchy-Green deformation tensor
            floatVector rightCauchyGreen;
            floatMatrix dRCGdPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen, dRCGdPrecedingF ) );

            //Compute the decomposition of the stress
            variableType pressure;
            variableVector deviatoricReferenceStress;

            variableMatrix dDevStressdStress, dDevStressdRCG;
            variableVector dPressuredStress, dPressuredRCG;

            variableMatrix d2DevStressdStressdRCG, d2PressuredStressdRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeSecondOrderReferenceStressDecomposition( stressMeasure,
                                                       rightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                                                       dDevStressdRCG, dPressuredStress, dPressuredRCG, d2DevStressdStressdRCG, d2PressuredStressdRCG ) )

            variableMatrix dDevStressdPrecedingF = tardigradeVectorTools::dot( dDevStressdRCG, dRCGdPrecedingF );
            variableVector dPressuredPrecedingF  = tardigradeVectorTools::Tdot( dRCGdPrecedingF, dPressuredRCG );

            variableMatrix d2DevStressdStressdPrecedingF( deviatoricReferenceStress.size( ), floatVector( deviatoricReferenceStress.size( ) * precedingDeformationGradient.size( ), 0 ) );

            variableMatrix d2PressuredStressdPrecedingF( deviatoricReferenceStress.size( ), floatVector( precedingDeformationGradient.size( ), 0 ) );

            for ( unsigned int I = 0; I < deviatoricReferenceStress.size( ); I++ ){

                for ( unsigned int J = 0; J < deviatoricReferenceStress.size( ); J++ ){

                    for ( unsigned int K = 0; K < precedingDeformationGradient.size( ); K++ ){

                        for ( unsigned int L = 0; L < rightCauchyGreen.size( ); L++ ){

                            d2DevStressdStressdPrecedingF[ I ][ precedingDeformationGradient.size( ) * J + K ]
                                += d2DevStressdStressdRCG[ I ][ rightCauchyGreen.size( ) * J + L ] * dRCGdPrecedingF[ L ][ K ];

                        }

                        d2PressuredStressdPrecedingF[ I ][ J ]
                            += d2PressuredStressdRCG[ I ][ K ] * dRCGdPrecedingF[ K ][ J ];

                    }

                }

            }

            //Compute the l2norm of the deviatoric stress
            variableType normDevStress = tardigradeVectorTools::l2norm( deviatoricReferenceStress );

            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

            //Evaluate the jacobians
            variableVector devStressDirection = deviatoricReferenceStress / ( normDevStress + tol );

            dFdStress = tardigradeVectorTools::Tdot( dDevStressdStress, devStressDirection )
                      + BAngle * dPressuredStress;

            dFdc = - AAngle;

            dFdPrecedingF = tardigradeVectorTools::Tdot( dDevStressdPrecedingF, devStressDirection )
                          + BAngle * dPressuredPrecedingF;

            //Evaluate the second-order jacobians
            constantMatrix EYE = tardigradeVectorTools::eye< constantType >( dim * dim );
            variableMatrix dDevStressDirectiondDevStress = ( EYE - tardigradeVectorTools::dyadic( devStressDirection, devStressDirection ) ) / ( normDevStress + tol );

            d2FdStress2 = tardigradeVectorTools::Tdot( dDevStressdStress, tardigradeVectorTools::dot( dDevStressDirectiondDevStress, dDevStressdStress ) );

            d2FdStressdPrecedingF = tardigradeVectorTools::Tdot( tardigradeVectorTools::dot( dDevStressDirectiondDevStress, dDevStressdStress ), dDevStressdPrecedingF )
                                  + BAngle * d2PressuredStressdPrecedingF;

            for ( unsigned int I = 0; I < dim; I++ ){
                for ( unsigned int J = 0; J < dim; J++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        for ( unsigned int L = 0; L < dim; L++ ){
                            for ( unsigned int A = 0; A < dim; A++ ){
                                for ( unsigned int B = 0; B < dim; B++ ){
                                    d2FdStressdPrecedingF[ dim * I + J ][ dim * K + L ] += devStressDirection[ dim * A + B ] * d2DevStressdStressdPrecedingF[ dim * A + B ][ dim * dim * dim * I + dim * dim * J + dim * K + L ];
                                }
                            }
                        }
                    }
                }
            }

        }

        void computeHigherOrderDruckerPragerYieldEquation( const variableVector &stressMeasure,
                                                               const variableVector &cohesion,
                                                               const variableVector &precedingDeformationGradient,
                                                               const parameterType &frictionAngle, const parameterType &beta,
                                                               variableVector &yieldValue ){
            /*!
             * Compute the higher-order Drucker Prager Yield equation
             *
             * \f$F_K = ||dev ( M ) ||_K - \left( A^{\phi} \bar{c}_K - B^{\phi} \bar{p}_K \right) \leq 0
             * 
             * || dev ( stressMeasure ) ||_K = \sqrt{ dev( M )_{IJK} : dev( M )_{IJK} }\f$
             * where the K's aren't summed.
             *  \f$dev( M )_{IJK} = M_{IJK} - \bar{p}_K elasticRightCauchyGreen_{IJ}^{-1}
             *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} M_{IJK}
             *  A^{angle} = \beta^{angle} \cos( frictionAngle )
             *  B^{angle} = \beta^{angle} \sin( frictionAngle )
             *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }\f$
             *
             * :param const variableVector &stressMeasure: The stress measure
             * :param const variableVector &cohesion: The cohesion measure.
             * :param const variableVector &precedingDeformationGradient: The preceding deformation gradient
             * :param const parameterType &frictionAngle: The friction angle
             * :param const parameterType &beta: The beta parameter
             * :param variableVector &yieldValue: The yield value.
             */

            parameterType AAngle, BAngle;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle ) )

            //Compute the right Cauchy-Green deformation tensor
            floatVector rightCauchyGreen;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen ) );

            //Compute the decomposition of the stress
            variableVector pressure;
            variableVector deviatoricReferenceStress;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( stressMeasure,
                                 rightCauchyGreen, deviatoricReferenceStress, pressure ) );

            //Compute the l2norm of the deviatoric stress
            variableVector normDevStress;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress ) );

            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

        }

        void computeHigherOrderDruckerPragerYieldEquation( const variableVector &stressMeasure,
                                                               const variableVector &cohesion,
                                                               const variableVector &precedingDeformationGradient,
                                                               const parameterType &frictionAngle, const parameterType &beta,
                                                               variableVector &yieldValue, variableMatrix &dFdStress, variableMatrix &dFdc,
                                                               variableMatrix &dFdPrecedingF ){
            /*!
             * Compute the higher-order Drucker Prager Yield equation
             *
             * \f$F_K = ||dev ( M ) ||_K - \left( A^{\phi} \bar{c}_K - B^{\phi} \bar{p}_K \right) \leq 0
             * 
             * || dev ( stressMeasure ) ||_K = \sqrt{ dev( M )_{IJK} : dev( M )_{IJK} }\f$
             * where the K's aren't summed.
             *  \f$dev( M )_{IJK} = M_{IJK} - \bar{p}_K elasticRightCauchyGreen_{IJ}^{-1}
             *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} M_{IJK}
             *  A^{angle} = \beta^{angle} \cos( frictionAngle )
             *  B^{angle} = \beta^{angle} \sin( frictionAngle )
             *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }\f$
             *
             *  Also computes the Jacobians
             *
             * :param const variableVector &stressMeasure: The stress measure
             * :param const variableVector &cohesion: The cohesion measure.
             * :param const variableVector &precedingDeformationGradient: The preceding deformation gradient
             * :param const parameterType &frictionAngle: The friction angle
             * :param const parameterType &beta: The beta parameter
             * :param variableVector &yieldValue: The yield value.
             * :param variableMatrix &dFdStress: The Jacobian of the yield function w.r.t. the reference higher order stress.
             * :param variableMatrix &dFdc: The Jacobian of the yield function w.r.t. the cohesion.
             * :param variableMatrix &dFdElasticRCG: The Jacobian of the yield function w.r.t. the preceding deformation gradient
             *     deformation tensor.
             */

            parameterType AAngle, BAngle;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle ) );

            //Compute the right Cauchy-Green deformation tensor
            floatVector rightCauchyGreen;
            floatMatrix dRCGdPrecedingF;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen, dRCGdPrecedingF ) );

            //Compute the decomposition of the stress
            variableVector pressure;
            variableVector deviatoricReferenceStress;

            variableMatrix dDevStressdStress, dDevStressdRCG;
            variableMatrix dPressuredStress, dPressuredRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( stressMeasure,
                                                       rightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                                                       dDevStressdRCG, dPressuredStress, dPressuredRCG ) );

            floatMatrix dDevStressdPrecedingF = tardigradeVectorTools::dot( dDevStressdRCG, dRCGdPrecedingF );
            floatMatrix dPressuredPrecedingF  = tardigradeVectorTools::dot( dPressuredRCG,  dRCGdPrecedingF );

            //Compute the l2norm of the deviatoric stress
            variableVector normDevStress;
            variableMatrix dNormDevStressdDevStress;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress, dNormDevStressdDevStress ) );

            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );

            //Construct the Jacobians
            dFdStress = tardigradeVectorTools::dot( dNormDevStressdDevStress, dDevStressdStress )
                      + BAngle * dPressuredStress;

            dFdc = -AAngle * tardigradeVectorTools::eye< constantType >( cohesion.size() );

            dFdPrecedingF = tardigradeVectorTools::dot( dNormDevStressdDevStress, dDevStressdPrecedingF )
                          + BAngle * dPressuredPrecedingF;

        }

        void computeHigherOrderDruckerPragerYieldEquation( const variableVector &stressMeasure,
                                                           const variableVector &cohesion,
                                                           const variableVector &precedingDeformationGradient,
                                                           const parameterType &frictionAngle, const parameterType &beta,
                                                           variableVector &yieldValue, variableMatrix &dFdStress, variableMatrix &dFdc,
                                                           variableMatrix &dFdPrecedingF, variableMatrix &d2FdStress2,
                                                           variableMatrix &d2FdStressdPrecedingF ){
            /*!
             * Compute the higher-order Drucker Prager Yield equation
             *
             * \f$F_K = ||dev ( M ) ||_K - \left( A^{\phi} \bar{c}_K - B^{\phi} \bar{p}_K \right) \leq 0
             * 
             * || dev ( stressMeasure ) ||_K = \sqrt{ dev( M )_{IJK} : dev( M )_{IJK} }\f$
             * where the K's aren't summed.
             * \f$dev( M )_{IJK} = M_{IJK} - \bar{p}_K elasticRightCauchyGreen_{IJ}^{-1}
             *  \bar{p} = \frac{1}{3} elasticRightCauchyGreen_{IJ} M_{IJK}
             *  A^{angle} = \beta^{angle} \cos( frictionAngle )
             *  B^{angle} = \beta^{angle} \sin( frictionAngle )
             *  \beta^{angle} = \frac{2 \sqrt{6} }{3 + \beta \sin( frictionAngle ) }\f$
             *
             *  Also computes the Jacobians
             *
             * :param const variableVector &stressMeasure: The stress measure
             * :param const variableVector &cohesion: The cohesion measure.
             * :param const variableVector &precedingDeformationGradient: The preceding deformation gradient
             * :param const parameterType &frictionAngle: The friction angle
             * :param const parameterType &beta: The beta parameter
             * :param variableVector &yieldValue: The yield value.
             * :param variableMatrix &dFdStress: The Jacobian of the yield function w.r.t. the reference higher order stress.
             * :param variableMatrix &dFdc: The Jacobian of the yield function w.r.t. the cohesion.
             * :param variableMatrix &dFdPrecedingF: The Jacobian of the yield function w.r.t. the preceding deformation gradient
             * :param variableMatrix &d2FdStress2: The second order Jacobian of the yield function w.r.t. the reference 
             *     higher order stress.
             * :param variableMatrix &d2FdStressdPrecedingF: The second order Jacobian of the yield function w.r.t. the 
             *     reference higher order stress and the preceding deformation gradient
             */

            //Assume 3D
            unsigned int dim = 3;

            parameterType AAngle, BAngle;
            TARDIGRADE_ERROR_TOOLS_CATCH( computeDruckerPragerInternalParameters( frictionAngle, beta, AAngle, BAngle ) );

            //Compute the right Cauchy-Green deformation tensor
            floatVector rightCauchyGreen;
            floatMatrix dRCGdPrecedingF;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( precedingDeformationGradient, rightCauchyGreen, dRCGdPrecedingF ) );

            //Compute the decomposition of the stress
            variableVector pressure;
            variableVector deviatoricReferenceStress;

            variableMatrix dDevStressdStress, dDevStressdRCG;
            variableMatrix dPressuredStress, dPressuredRCG;

            variableMatrix d2DevStressdStressdRCG, d2PressuredStressdRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderReferenceStressDecomposition( stressMeasure,
                                                       rightCauchyGreen, deviatoricReferenceStress, pressure, dDevStressdStress,
                                                       dDevStressdRCG, dPressuredStress, dPressuredRCG, d2DevStressdStressdRCG, d2PressuredStressdRCG ) )

            floatMatrix dDevStressdPrecedingF = tardigradeVectorTools::dot( dDevStressdRCG, dRCGdPrecedingF );
            floatMatrix dPressuredPrecedingF  = tardigradeVectorTools::dot( dPressuredRCG,  dRCGdPrecedingF );

            variableMatrix d2DevStressdStressdPrecedingF( deviatoricReferenceStress.size( ), floatVector( deviatoricReferenceStress.size( ) * precedingDeformationGradient.size( ), 0 ) );

            variableMatrix d2PressuredStressdPrecedingF( pressure.size( ), floatVector( deviatoricReferenceStress.size( ) * precedingDeformationGradient.size( ), 0 ) );

            for ( unsigned int I = 0; I < deviatoricReferenceStress.size( ); I++ ){

                for ( unsigned int J = 0; J < stressMeasure.size( ); J++ ){

                    for ( unsigned int K = 0; K < precedingDeformationGradient.size( ); K++ ){

                        for ( unsigned int L = 0; L < rightCauchyGreen.size( ); L++ ){

                            d2DevStressdStressdPrecedingF[ I ][ precedingDeformationGradient.size( ) * J + K ]
                                += d2DevStressdStressdRCG[ I ][ rightCauchyGreen.size( ) * J + L ] * dRCGdPrecedingF[ L ][ K ];

                        }

                    }

                }

            }

            for ( unsigned int I = 0; I < pressure.size( ); I++ ){

                for ( unsigned int J = 0; J < stressMeasure.size( ); J++ ){

                    for ( unsigned int K = 0; K < precedingDeformationGradient.size( ); K++ ){

                        for ( unsigned int L = 0; L < rightCauchyGreen.size( ); L++ ){

                            d2PressuredStressdPrecedingF[ I ][ precedingDeformationGradient.size( ) * J + K ]
                                += d2PressuredStressdRCG[ I ][ rightCauchyGreen.size( ) * J + L ] * dRCGdPrecedingF[ L ][ K ];

                        }

                    }

                }

            }

            //Compute the l2norm of the deviatoric stress
            variableVector normDevStress;
            variableMatrix dNormDevStressdDevStress;
            variableMatrix d2NormDevStressdDevStress2;
            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::computeHigherOrderStressNorm( deviatoricReferenceStress, normDevStress,
                                                                                                                  dNormDevStressdDevStress,
                                                                                                                  d2NormDevStressdDevStress2 ) )

            //Evaluate the yield equation
            yieldValue = normDevStress - ( AAngle * cohesion - BAngle * pressure );
    
            //Construct the Jacobians
            dFdStress = tardigradeVectorTools::dot( dNormDevStressdDevStress, dDevStressdStress )
                      + BAngle * dPressuredStress;
    
            dFdc = -AAngle * tardigradeVectorTools::eye< constantType >( cohesion.size() );
    
            dFdPrecedingF = tardigradeVectorTools::dot( dNormDevStressdDevStress, dDevStressdPrecedingF )
                          + BAngle * dPressuredPrecedingF;
    
            //Construct the second-order jacobians
            d2FdStress2 = variableMatrix( dim, variableVector( dim * dim * dim * dim * dim * dim, 0 ) );
            d2FdStressdPrecedingF = tardigradeVectorTools::dot( dNormDevStressdDevStress, d2DevStressdStressdPrecedingF )
                                  + BAngle * d2PressuredStressdPrecedingF;

            for ( unsigned int K = 0; K < 3; K++ ){
                for ( unsigned int L = 0; L < 3; L++ ){
                    for ( unsigned int M = 0; M < 3; M++ ){
                        for ( unsigned int N = 0; N < 3; N++ ){
                            for ( unsigned int O = 0; O < 3; O++ ){
                                for ( unsigned int P = 0; P < 3; P++ ){
                                    for ( unsigned int Q = 0; Q < 3; Q++ ){
                                        for ( unsigned int A = 0; A < 3; A++ ){
                                            for ( unsigned int B = 0; B < 3; B++ ){
                                                for ( unsigned int C = 0; C < 3; C++ ){
                                                    for ( unsigned int D = 0; D < 3; D++ ){
                                                        for ( unsigned int E = 0; E < 3; E++ ){
                                                            d2FdStressdPrecedingF[ K ][ dim * dim * dim * dim * L + dim * dim * dim * M + dim * dim * N + dim * O + P ]
                                                                += d2NormDevStressdDevStress2[ K ][ dim * dim * dim * dim * dim * Q + dim * dim * dim * dim * A + dim * dim * dim * B + dim * dim * C + dim * D + E ]
                                                                 * dDevStressdStress[ dim * dim * Q + dim * A + B ][ dim * dim * L + dim * M + N ]
                                                                 * dDevStressdPrecedingF[ dim * dim * C + dim * D + E ][ dim * O + P ];
                                                            for ( unsigned int F = 0; F < 3; F++ ){
                                                                d2FdStress2[ K ][ dim * dim * dim * dim * dim * L + dim * dim * dim * dim * M + dim * dim * dim * N + dim * dim * O + dim * P + Q ]
                                                                    += d2NormDevStressdDevStress2[ K ][ dim * dim * dim * dim * dim * A + dim * dim * dim * dim * B + dim * dim * dim * C + dim * dim * D + dim * E + F ]
                                                                     * dDevStressdStress[ dim * dim * A + dim * B + C ][ dim * dim * L + dim * M + N ]
                                                                     * dDevStressdStress[ dim * dim * D + dim * E + F ][ dim * dim * O + dim * P + Q ];
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

        }

        void residual::setMacroDrivingStress( ){
            /*!
             * Set the macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration)
             */

            setDrivingStresses( false );

        }

        void residual::setSymmetricMicroDrivingStress( ){
            /*!
             * Set the symmetric micro driving stress (stress in current configuration of plastic configuration)
             */

            setDrivingStresses( false );

        }

        void residual::setHigherOrderDrivingStress( ){
            /*!
             * Set the higher-order driving stress (stress in current configuration of plastic configuration)
             */

            setDrivingStresses( false );

        }

        void residual::setPreviousMacroDrivingStress( ){
            /*!
             * Set the previous macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration)
             */

            setDrivingStresses( true );

        }

        void residual::setPreviousSymmetricMicroDrivingStress( ){
            /*!
             * Set the previous symmetric micro driving stress (stress in current configuration of plastic configuration)
             */

            setDrivingStresses( true );

        }

        void residual::setPreviousHigherOrderDrivingStress( ){
            /*!
             * Set the previous higher-order driving stress (stress in current configuration of plastic configuration)
             */

            setDrivingStresses( true );

        }

        void residual::setDrivingStresses( const bool isPrevious ){
            /*!
             * Set the driving stresses for the plasticity
             * 
             * We here assume that the driving stresses are in the current configuration of
             * this residual's configuration.
             *
             * We also assume that the stress from hydra is in the reference configuration
             */

            const unsigned int *dim = hydra->getDimension( );

            const floatVector *stress;

            floatVector Fp;

            floatVector chip;

            if ( isPrevious ){

                stress = hydra->getPreviousStress( );

                Fp     = hydra->getPreviousFollowingConfiguration(      ( *getPlasticConfigurationIndex( ) ) - 1 );

                chip   = hydra->getPreviousFollowingMicroConfiguration( ( *getPlasticConfigurationIndex( ) ) - 1 );

            }
            else{

                stress = hydra->getStress( );

                Fp     = hydra->getFollowingConfiguration(      ( *getPlasticConfigurationIndex( ) ) - 1 );

                chip   = hydra->getFollowingMicroConfiguration( ( *getPlasticConfigurationIndex( ) ) - 1 );

            }

            // Extract the stresses from the stress vector
            floatVector PK2Stress(                     stress->begin( ),                           stress->begin( ) + 1 * ( *dim ) * ( *dim ) );

            floatVector referenceSymmetricMicroStress( stress->begin( ) + 1 * ( *dim ) * ( *dim ), stress->begin( ) + 2 * ( *dim ) * ( *dim ) );;

            floatVector referenceHigherOrderStress(    stress->begin( ) + 2 * ( *dim ) * ( *dim ), stress->begin( ) + 2 * ( *dim ) * ( *dim ) + ( *dim ) * ( *dim ) * ( *dim ) );

            // Push the stresses forward to the current configuration of the plastic configuration
            floatVector macroDrivingStress;

            floatVector symmetricMicroDrivingStress;

            floatVector higherOrderDrivingStress;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress, Fp, macroDrivingStress ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceSymmetricMicroStress, Fp, symmetricMicroDrivingStress ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, Fp, chip, higherOrderDrivingStress ) );

            if ( isPrevious ){

                set_previousMacroDrivingStress(          macroDrivingStress );

                set_previousSymmetricMicroDrivingStress( symmetricMicroDrivingStress );

                set_previousHigherOrderDrivingStress(    higherOrderDrivingStress );

            }
            else{

                set_macroDrivingStress(          macroDrivingStress );

                set_symmetricMicroDrivingStress( symmetricMicroDrivingStress );

                set_higherOrderDrivingStress(    higherOrderDrivingStress );

            }

        }

        void residual::setdMacroDrivingStressdMacroStress( ){
            /*!
             * Set the jacobian of the macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration) w.r.t. the macro stress
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdSymmetricMicroDrivingStressdMicroStress( ){
            /*!
             * Set the jacobian of the symmetric micro driving stress (stress in current configuration of plastic configuration) w.r.t. the micro stress
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdHigherOrderDrivingStressdHigherOrderStress( ){
            /*!
             * Set the jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdMacroDrivingStressdF( ){
            /*!
             * Set the jacobian of the macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration) w.r.t. the deformation gradient
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdSymmetricMicroDrivingStressdF( ){
            /*!
             * Set the jacobian of the symmetric micro driving stress (stress in current configuration of plastic configuration) w.r.t. the deformation gradient
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdHigherOrderDrivingStressdF( ){
            /*!
             * Set the jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the deformation gradient
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdHigherOrderDrivingStressdChi( ){
            /*!
             * Set the jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the micro deformation
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdMacroDrivingStressdFn( ){
            /*!
             * Set the jacobian of the macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration) w.r.t. the sub deformation gradients
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdSymmetricMicroDrivingStressdFn( ){
            /*!
             * Set the jacobian of the symmetric micro driving stress (stress in current configuration of plastic configuration) w.r.t. the sub deformation gradients
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdHigherOrderDrivingStressdFn( ){
            /*!
             * Set the jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the sub deformation gradients
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setdHigherOrderDrivingStressdChin( ){
            /*!
             * Set the jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the sub micro deformations
             */

            setDrivingStressesJacobians( false );

        }

        void residual::setPreviousdMacroDrivingStressdMacroStress( ){
            /*!
             * Set the previous jacobian of the macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration) w.r.t. the macro stress
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdSymmetricMicroDrivingStressdMicroStress( ){
            /*!
             * Set the previous jacobian of the symmetric micro driving stress (stress in current configuration of plastic configuration) w.r.t. the micro stress
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderDrivingStressdHigherOrderStress( ){
            /*!
             * Set the previous jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdMacroDrivingStressdF( ){
            /*!
             * Set the previous jacobian of the macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration) w.r.t. the deformation gradient
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdSymmetricMicroDrivingStressdF( ){
            /*!
             * Set the previous jacobian of the symmetric micro driving stress (stress in current configuration of plastic configuration) w.r.t. the deformation gradient
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderDrivingStressdF( ){
            /*!
             * Set the previous jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the deformation gradient
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderDrivingStressdChi( ){
            /*!
             * Set the previous jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the micro deformation
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdMacroDrivingStressdFn( ){
            /*!
             * Set the previous jacobian of the macro (i.e. the stress associated with the Cauchy stress) driving stress (stress in current configuration of plastic configuration) w.r.t. the sub deformation gradients
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdSymmetricMicroDrivingStressdFn( ){
            /*!
             * Set the previous jacobian of the symmetric micro driving stress (stress in current configuration of plastic configuration) w.r.t. the sub deformation gradients
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderDrivingStressdFn( ){
            /*!
             * Set the previous jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the sub deformation gradients
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setPreviousdHigherOrderDrivingStressdChin( ){
            /*!
             * Set the previous jacobian of the higher-order driving stress (stress in current configuration of plastic configuration) w.r.t. the higher order stress w.r.t. the sub micro deformations
             */

            setDrivingStressesJacobians( true );

        }

        void residual::setDrivingStressesJacobians( const bool isPrevious ){
            /*!
             * Set the driving stresses for the plasticity along with the Jacobians
             * 
             * We here assume that the driving stresses are in the current configuration of
             * this residual's configuration.
             *
             * We also assume that the stress from hydra is in the reference configuration
             */

            const unsigned int *dim = hydra->getDimension( );

            const floatVector *stress;

            floatMatrix dFpdSubFs;

            const floatMatrix *dF1dF;

            const floatMatrix *dF1dFn;

            floatMatrix dChipdSubChis;

            const floatMatrix *dChi1dChi;

            const floatMatrix *dChi1dChin;

            floatVector Fp;

            floatVector chip;

            if ( isPrevious ){

                stress = hydra->getPreviousStress( );

                dF1dF         = hydra->get_previousdF1dF( );

                dF1dFn        = hydra->get_previousdF1dFn( );

                dFpdSubFs     = hydra->getPreviousFollowingConfigurationJacobian( ( *getPlasticConfigurationIndex( ) ) - 1 );

                dChi1dChi     = hydra->get_previousdChi1dChi( );

                dChi1dChin    = hydra->get_previousdChi1dChin( );

                dChipdSubChis = hydra->getPreviousFollowingMicroConfigurationJacobian( ( *getPlasticConfigurationIndex( ) ) - 1 );

                Fp            = hydra->getPreviousFollowingConfiguration(         ( *getPlasticConfigurationIndex( ) ) - 1 );

                chip          = hydra->getPreviousFollowingMicroConfiguration(    ( *getPlasticConfigurationIndex( ) ) - 1 );

            }
            else{

                stress = hydra->getStress( );

                dF1dF         = hydra->get_dF1dF( );

                dF1dFn        = hydra->get_dF1dFn( );

                dFpdSubFs     = hydra->getFollowingConfigurationJacobian( ( *getPlasticConfigurationIndex( ) ) - 1 );

                dChi1dChi     = hydra->get_dChi1dChi( );

                dChi1dChin    = hydra->get_dChi1dChin( );

                dChipdSubChis = hydra->getFollowingMicroConfigurationJacobian( ( *getPlasticConfigurationIndex( ) ) - 1 );

                Fp            = hydra->getFollowingConfiguration(         ( *getPlasticConfigurationIndex( ) ) - 1 );

                chip          = hydra->getFollowingMicroConfiguration(    ( *getPlasticConfigurationIndex( ) ) - 1 );

            }

            // Assemble the derivatives of the deformation gradient map
            floatMatrix dFpdF(  Fp.size( ), floatVector( hydra->getDeformationGradient( )->size( ), 0 ) );

            floatMatrix dFpdFn( Fp.size( ), floatVector( ( ( *hydra->getNumConfigurations( ) ) - 1 ) * hydra->getDeformationGradient( )->size( ), 0 ) );

            floatMatrix dChipdChi(  chip.size( ), floatVector( hydra->getDeformationGradient( )->size( ), 0 ) );

            floatMatrix dChipdChin( chip.size( ), floatVector( ( ( *hydra->getNumConfigurations( ) ) - 1 ) * hydra->getMicroDeformation( )->size( ), 0 ) );

            for ( unsigned int i = 0; i < ( *dim ) * ( *dim ); i++ ){

                for ( unsigned int j = 0; j < ( *dim ) * ( *dim ); j++ ){

                    for ( unsigned int k = 0; k < ( *dim ) * ( *dim ); k++ ){

                        dFpdF[ i ][ j ] += dFpdSubFs[ i ][ k ] * ( *dF1dF )[ k ][ j ];

                        dChipdChi[ i ][ j ] += dChipdSubChis[ i ][ k ] * ( *dChi1dChi )[ k ][ j ];

                    }

                }

                for ( unsigned int j = 0; j < ( ( *hydra->getNumConfigurations( ) ) - 1 ) * ( *dim ) * ( *dim ); j++ ){

                    dFpdFn[ i ][ j ] += dFpdSubFs[ i ][ j + ( *dim ) * ( *dim ) ];

                    dChipdChin[ i ][ j ] += dChipdSubChis[ i ][ j + ( *dim ) * ( *dim ) ];

                    for ( unsigned int k = 0; k < ( *dim ) * ( *dim ); k++ ){

                        dFpdFn[ i ][ j ] += dFpdSubFs[ i ][ k ] * ( *dF1dFn )[ k ][ j ];

                        dChipdChin[ i ][ j ] += dChipdSubChis[ i ][ k ] * ( *dChi1dChin )[ k ][ j ];

                    }

                }

            }

            // Extract the stresses from the stress vector
            floatVector PK2Stress(                     stress->begin( ),                           stress->begin( ) + 1 * ( *dim ) * ( *dim ) );

            floatVector referenceSymmetricMicroStress( stress->begin( ) + 1 * ( *dim ) * ( *dim ), stress->begin( ) + 2 * ( *dim ) * ( *dim ) );;

            floatVector referenceHigherOrderStress(    stress->begin( ) + 2 * ( *dim ) * ( *dim ), stress->begin( ) + 2 * ( *dim ) * ( *dim ) + ( *dim ) * ( *dim ) * ( *dim ) );

            // Push the stresses forward to the current configuration of the plastic configuration
            floatVector macroDrivingStress;

            floatVector symmetricMicroDrivingStress;

            floatVector higherOrderDrivingStress;

            floatMatrix dMacrodFp;

            floatMatrix dMacrodPK2;

            floatMatrix dMicrodFp;

            floatMatrix dMicrodSigma;

            floatMatrix dHigherdFp;

            floatMatrix dHigherdChip;

            floatMatrix dHigherdM;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pushForwardPK2Stress( PK2Stress, Fp, macroDrivingStress,
                                                                                                          dMacrodPK2, dMacrodFp ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pushForwardReferenceMicroStress( referenceSymmetricMicroStress, Fp, symmetricMicroDrivingStress,
                                                                                                                     dMicrodSigma, dMicrodFp ) );

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeMicromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, Fp, chip, higherOrderDrivingStress,
                                                                                                                  dHigherdM, dHigherdFp, dHigherdChip ) );

            if ( isPrevious ){

                set_previousMacroDrivingStress(          macroDrivingStress );

                set_previousSymmetricMicroDrivingStress( symmetricMicroDrivingStress );

                set_previousHigherOrderDrivingStress(    higherOrderDrivingStress );

                set_previousdMacroDrivingStressdMacroStress( dMacrodPK2 );

                set_previousdMacroDrivingStressdF( tardigradeVectorTools::dot( dMacrodFp, dFpdF ) );

                set_previousdMacroDrivingStressdFn( tardigradeVectorTools::dot( dMacrodFp, dFpdFn ) );

                set_previousdSymmetricMicroDrivingStressdMicroStress( dMicrodSigma );

                set_previousdSymmetricMicroDrivingStressdF( tardigradeVectorTools::dot( dMicrodFp, dFpdF ) );

                set_previousdSymmetricMicroDrivingStressdFn( tardigradeVectorTools::dot( dMicrodFp, dFpdFn ) );

                set_previousdHigherOrderDrivingStressdHigherOrderStress( dHigherdM );

                set_previousdHigherOrderDrivingStressdF( tardigradeVectorTools::dot( dHigherdFp, dFpdF ) );

                set_previousdHigherOrderDrivingStressdFn( tardigradeVectorTools::dot( dHigherdFp, dFpdFn ) );

                set_previousdHigherOrderDrivingStressdChi( tardigradeVectorTools::dot( dHigherdChip, dChipdChi ) );

                set_previousdHigherOrderDrivingStressdChin( tardigradeVectorTools::dot( dHigherdChip, dChipdChin ) );

            }
            else{

                set_macroDrivingStress(          macroDrivingStress );

                set_symmetricMicroDrivingStress( symmetricMicroDrivingStress );

                set_higherOrderDrivingStress(    higherOrderDrivingStress );

                set_dMacroDrivingStressdMacroStress( dMacrodPK2 );

                set_dMacroDrivingStressdF( tardigradeVectorTools::dot( dMacrodFp, dFpdF ) );

                set_dMacroDrivingStressdFn( tardigradeVectorTools::dot( dMacrodFp, dFpdFn ) );

                set_dSymmetricMicroDrivingStressdMicroStress( dMicrodSigma );

                set_dSymmetricMicroDrivingStressdF( tardigradeVectorTools::dot( dMicrodFp, dFpdF ) );

                set_dSymmetricMicroDrivingStressdFn( tardigradeVectorTools::dot( dMicrodFp, dFpdFn ) );

                set_dHigherOrderDrivingStressdHigherOrderStress( dHigherdM );

                set_dHigherOrderDrivingStressdF( tardigradeVectorTools::dot( dHigherdFp, dFpdF ) );

                set_dHigherOrderDrivingStressdFn( tardigradeVectorTools::dot( dHigherdFp, dFpdFn ) );

                set_dHigherOrderDrivingStressdChi( tardigradeVectorTools::dot( dHigherdChip, dChipdChi ) );

                set_dHigherOrderDrivingStressdChin( tardigradeVectorTools::dot( dHigherdChip, dChipdChin ) );

            }

        }

        void residual::extractMaterialParameters( const parameterVector &parameters ){
            /*!
             * Extract the parameters from the parameter vector
             *
             * :param const std::vector< double > &parameters: The incoming parameter vector
             * :param parameterVector &macroHardeningParameters: The parameters used in the hardening of the macro Strain ISV
             *     (initial cohesion, hardening modulus)
             * :param parameterVector &microHardeningParameters: The parameters used in the hardening of the micro Strain ISV
             *     (initial cohesion, hardening modulus)
             * :param parameterVector &microGradientHardeningParameters: The parameters used in the hardening of the micro Gradient Strain ISV
             *     (initial cohesion, hardening modulus)
             * :param parameterVector &macroFlowParameters: The parameters used in the macro flow direction computation.
             *     (friction angle, beta )
             * :param parameterVector &microFlowParameters: The parameters used in the micro flow direction computation
             *     (friction angle, beta )
             * :param parameterVector &microGradientFlowParameters: The parameters used in the micro Gradient flow direction computation.
             *     (friction angle, beta )
             * :param parameterVector &macroYieldParameters: The parameters used in the macro yielding computation.
             *     (friction angle, beta )
             * :param parameterVector &microYieldParameters: The parameters used in the micro yielding computation
             *     (friction angle, beta )
             * :param parameterVector &microGradientYieldParameters: The parameters used in the micro Gradient yielding computation.
             *     (friction angle, beta )
             */
        
            if ( parameters.size() == 0 ){
        
                TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( "The parameter vector has a length of zero" ) );
        
            }
        
            unsigned int start = 0;
            unsigned int span;
        
            std::vector< parameterVector > outputs( 9 );
        
            //Extract the material parameters
            for ( unsigned int i = 0; i < outputs.size(); i++ ){
                span = ( unsigned int )std::floor( parameters[ start ]  + 0.5 ); //Extract the span of the parameter set
        
                if ( parameters.size() < start + 1 + span ){
                    std::string outstr = "fparams is not long enough to contain all of the required parameters:\n";
                    outstr +=            "    filling variable " + std::to_string( i ) + "\n";
                    outstr +=            "    size =          "  + std::to_string( parameters.size() ) + "\n";
                    outstr +=            "    required size = "  + std::to_string( start + 1 + span );
        
                    TARDIGRADE_ERROR_TOOLS_CATCH( throw std::runtime_error( outstr ); )
        
                }
        
                outputs[ i ] = parameterVector( parameters.begin() + start + 1, parameters.begin() + start + 1 + span );
        
                start = start + 1 + span;
            }

            //Set the output values
            set_macroHardeningParameters(         outputs[ 0 ] );
        
            set_microHardeningParameters(         outputs[ 1 ] );
        
            set_microGradientHardeningParameters( outputs[ 2 ] );
        
            set_macroFlowParameters(              outputs[ 3 ] );
        
            set_microFlowParameters(              outputs[ 4 ] );
        
            set_microGradientFlowParameters(      outputs[ 5 ] );
        
            set_macroYieldParameters(             outputs[ 6 ] );
        
            set_microYieldParameters(             outputs[ 7 ] );
        
            set_microGradientYieldParameters(     outputs[ 8 ] );
        
        }

        const unsigned int* residual::getPlasticConfigurationIndex( ){
            /*!
             * Get plastic configuration index
             */

            return &_plasticConfigurationIndex;

        }

        const std::vector< unsigned int >* residual::getStateVariableIndices( ){
            /*!
             * Get state variable indices
             */

            return &_stateVariableIndices;

        }

        const floatType* residual::getIntegrationParameter( ){
            /*!
             * Get the integration parameter
             */

            return &_integrationParameter;

        }

        void residual::setPlasticStateVariables( const bool isPrevious ){
            /*!
             * Set the plastic state variables
             * 
             * \param isPrevious: Flag for whether to set the current (false) or previous (true) values of the plastic state variables
             */

            floatVector plasticStateVariables( getStateVariableIndices( )->size( ), 0 );

            const floatVector *nonlinearISVs;

            if ( isPrevious ){

                nonlinearISVs = hydra->get_previousNonLinearSolveStateVariables( );

            }
            else{

                nonlinearISVs = hydra->get_nonLinearSolveStateVariables( );

            }

            for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->end( ); ind++ ){

                plasticStateVariables[ ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ] = ( *nonlinearISVs )[ *ind ];

            }

            if ( isPrevious ){

                set_previousPlasticStateVariables( plasticStateVariables );

            }
            else{

                set_plasticStateVariables( plasticStateVariables );

            }

        }

        void residual::setPlasticStateVariables( ){
            /*!
             * Set the plastic state variables
             */

            setPlasticStateVariables( false );

        }

        void residual::setPreviousPlasticStateVariables( ){
            /*!
             * Set the previous plastic state variables
             */

            setPlasticStateVariables( true );

        }

        void residual::setPlasticMultipliers( const bool isPrevious ){
            /*!
             * Set the plastic multipliers
             * 
             * \param isPrevious: Flag for whether to set the current (false) or previous (true) values of the plastic multipliers
             */

            floatVector plasticMultipliers( *getNumPlasticMultipliers( ), 0 );

            const floatVector *nonlinearISVs;

            if ( isPrevious ){

                nonlinearISVs = hydra->get_previousNonLinearSolveStateVariables( );

            }
            else{

                nonlinearISVs = hydra->get_nonLinearSolveStateVariables( );

            }

            for ( auto ind = getStateVariableIndices( )->begin( ); ind != getStateVariableIndices( )->begin( ) + *getNumPlasticMultipliers( ); ind++ ){

                plasticMultipliers[ ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) ] = ( *nonlinearISVs )[ *ind ];

            }

            if ( isPrevious ){

                set_previousPlasticMultipliers( plasticMultipliers );

            }
            else{

                set_plasticMultipliers( plasticMultipliers );

            }

        }

        void residual::setPlasticMultipliers( ){
            /*!
             * Set the plastic multipliers
             */

            setPlasticMultipliers( false );

        }

        void residual::setPreviousPlasticMultipliers( ){
            /*!
             * Set the previous plastic multipliers
             */

            setPlasticMultipliers( true );

        }

        void residual::setPlasticStrainLikeISVs( const bool isPrevious ){
            /*!
             * Set the plastic strain-like internal state variables
             * 
             * \param isPrevious: Flag for whether to set the current (false) or previous (true) values of the plastic multipliers
             */

            floatVector plasticStrainLikeISVs( getStateVariableIndices( )->size( ) - *getNumPlasticMultipliers( ), 0 );

            const floatVector *nonlinearISVs;

            if ( isPrevious ){

                nonlinearISVs = hydra->get_previousNonLinearSolveStateVariables( );

            }
            else{

                nonlinearISVs = hydra->get_nonLinearSolveStateVariables( );

            }

            for ( auto ind = getStateVariableIndices( )->begin( ) + *getNumPlasticMultipliers( ); ind != getStateVariableIndices( )->end( ); ind++ ){

                plasticStrainLikeISVs[ ( unsigned int )( ind - getStateVariableIndices( )->begin( ) ) - *getNumPlasticMultipliers( ) ] = ( *nonlinearISVs )[ *ind ];

            }

            if ( isPrevious ){

                set_previousPlasticStrainLikeISVs( plasticStrainLikeISVs );

            }
            else{

                set_plasticStrainLikeISVs( plasticStrainLikeISVs );

            }

        }

        void residual::setPlasticStrainLikeISVs( ){
            /*!
             * Set the plastic strain-like isvs
             */

            setPlasticStrainLikeISVs( false );

        }

        void residual::setPreviousPlasticStrainLikeISVs( ){
            /*!
             * Set the previous plastic strain-like isvs
             */

            setPlasticStrainLikeISVs( true );

        }

        void residual::setdMacroFlowdc( ){
            /*!
             * Set the derivative of the macro flow potential w.r.t. the cohesion
             */

            setFlowPotentialGradients( false );

        }

        void residual::setdMicroFlowdc( ){
            /*!
             * Set the derivative of the micro flow potential w.r.t. the cohesion
             */

            setFlowPotentialGradients( false );

        }

        void residual::setdMicroGradientFlowdc( ){
            /*!
             * Set the derivative of the micro gradient flow potential w.r.t. the cohesion
             */

            setFlowPotentialGradients( false );

        }

        void residual::setdMacroFlowdDrivingStress( ){
            /*!
             * Set the derivative of the macro flow potential w.r.t. the macro driving stress
             */

            setFlowPotentialGradients( false );

        }

        void residual::setdMicroFlowdDrivingStress( ){
            /*!
             * Set the derivative of the micro flow potential w.r.t. the micro driving stress
             */

            setFlowPotentialGradients( false );

        }

        void residual::setdMicroGradientFlowdDrivingStress( ){
            /*!
             * Set the derivative of the micro gradient flow potential w.r.t. the micro-gradient driving stress
             */

            setFlowPotentialGradients( false );

        }

        void residual::setPreviousdMacroFlowdc( ){
            /*!
             * Set the previous derivative of the macro flow potential w.r.t. the cohesion
             */

            setFlowPotentialGradients( true );

        }

        void residual::setPreviousdMicroFlowdc( ){
            /*!
             * Set the previous derivative of the micro flow potential w.r.t. the cohesion
             */

            setFlowPotentialGradients( true );

        }

        void residual::setPreviousdMicroGradientFlowdc( ){
            /*!
             * Set the previous derivative of the micro gradient flow potential w.r.t. the cohesion
             */

            setFlowPotentialGradients( true );

        }

        void residual::setPreviousdMacroFlowdDrivingStress( ){
            /*!
             * Set the previous derivative of the macro flow potential w.r.t. the macro driving stress
             */

            setFlowPotentialGradients( true );

        }

        void residual::setPreviousdMicroFlowdDrivingStress( ){
            /*!
             * Set the previous derivative of the micro flow potential w.r.t. the micro driving stress
             */

            setFlowPotentialGradients( true );

        }

        void residual::setPreviousdMicroGradientFlowdDrivingStress( ){
            /*!
             * Set the previous derivative of the micro gradient flow potential w.r.t. the micro-gradient driving stress
             */

            setFlowPotentialGradients( true );

        }

        void residual::setFlowPotentialGradients( const bool isPrevious ){
            /*!
             * Set the gradients of the flow potential
             *
             * \param &isPrevious: Flag for whether to set the current (false) or previous (true) gradients
             */

            const floatType   *macroCohesion;

            const floatType   *microCohesion;

            const floatVector *microGradientCohesion;

            const floatVector *macroDrivingStress;

            const floatVector *microDrivingStress;

            const floatVector *microGradientDrivingStress;

            floatVector precedingDeformationGradient;

            const floatVector *macroFlowParameters         = get_macroFlowParameters( );

            const floatVector *microFlowParameters         = get_microFlowParameters( );

            const floatVector *microGradientFlowParameters = get_microGradientFlowParameters( );

            if ( isPrevious ){

                precedingDeformationGradient = hydra->getPreviousPrecedingConfiguration( *getPlasticConfigurationIndex( ) );

                macroCohesion                = get_previousMacroCohesion( );

                microCohesion                = get_previousMicroCohesion( );

                microGradientCohesion        = get_previousMicroGradientCohesion( );

                macroDrivingStress           = get_previousMacroDrivingStress( );

                microDrivingStress           = get_previousSymmetricMicroDrivingStress( );

                microGradientDrivingStress   = get_previousHigherOrderDrivingStress( );

            }
            else{

                precedingDeformationGradient = hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) );

                macroCohesion                = get_macroCohesion( );

                microCohesion                = get_microCohesion( );

                microGradientCohesion        = get_microGradientCohesion( );

                macroDrivingStress           = get_macroDrivingStress( );

                microDrivingStress           = get_symmetricMicroDrivingStress( );

                microGradientDrivingStress   = get_higherOrderDrivingStress( );

            }

            floatType tempYield;

            floatVector tempVectorYield;

            floatType   dMacroFlowdCohesion, dMicroFlowdCohesion;

            floatVector dMacroFlowdDrivingStress, dMicroFlowdDrivingStress,
                        dMacroFlowdPrecedingF,    dMicroFlowdPrecedingF;

            floatMatrix dMicroGradientFlowdDrivingStress, dMicroGradientFlowdCohesion, dMicroGradientFlowdPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderDruckerPragerYieldEquation( *macroDrivingStress, *macroCohesion, precedingDeformationGradient,
                                                                                        ( *macroFlowParameters )[ 0 ], ( *macroFlowParameters )[ 1 ],
                                                                                        tempYield, dMacroFlowdDrivingStress, dMacroFlowdCohesion, dMacroFlowdPrecedingF ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderDruckerPragerYieldEquation( *microDrivingStress, *microCohesion, precedingDeformationGradient,
                                                                                        ( *microFlowParameters )[ 0 ], ( *microFlowParameters )[ 1 ],
                                                                                        tempYield, dMicroFlowdDrivingStress, dMicroFlowdCohesion, dMicroFlowdPrecedingF ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( computeHigherOrderDruckerPragerYieldEquation( *microGradientDrivingStress, *microGradientCohesion, precedingDeformationGradient,
                                                                                       ( *microGradientFlowParameters )[ 0 ], ( *microGradientFlowParameters )[ 1 ],
                                                                                       tempVectorYield, dMicroGradientFlowdDrivingStress, dMicroGradientFlowdCohesion, dMicroGradientFlowdPrecedingF ) );

            if ( isPrevious ){

                set_previousdMacroFlowdc( dMacroFlowdCohesion );

                set_previousdMicroFlowdc( dMicroFlowdCohesion );

                set_previousdMicroGradientFlowdc( dMicroGradientFlowdCohesion );

                set_previousdMacroFlowdDrivingStress( dMacroFlowdDrivingStress );

                set_previousdMicroFlowdDrivingStress( dMicroFlowdDrivingStress );

                set_previousdMicroGradientFlowdDrivingStress( dMicroGradientFlowdDrivingStress );

            }
            else{

                set_dMacroFlowdc( dMacroFlowdCohesion );

                set_dMicroFlowdc( dMicroFlowdCohesion );

                set_dMicroGradientFlowdc( dMicroGradientFlowdCohesion );

                set_dMacroFlowdDrivingStress( dMacroFlowdDrivingStress );

                set_dMicroFlowdDrivingStress( dMicroFlowdDrivingStress );

                set_dMicroGradientFlowdDrivingStress( dMicroGradientFlowdDrivingStress );

            }

        }

        void residual::setd2MacroFlowdDrivingStressdStress( ){
            /*!
             * Set the Jacobian of the derivative of the macro flow potential w.r.t. the macro driving stress w.r.t. the macro stress
             */

            setFlowPotentialGradientsJacobians( false );

        }

        void residual::setd2MacroFlowdDrivingStressdF( ){
            /*!
             * Set the Jacobian of the derivative of the macro flow potential w.r.t. the macro driving stress w.r.t. the deformation gradient
             */

            setFlowPotentialGradientsJacobians( false );

        }

        void residual::setd2MacroFlowdDrivingStressdFn( ){
            /*!
             * Set the Jacobian of the derivative of the macro flow potential w.r.t. the macro driving stress w.r.t. the sub deformation gradients
             */

            setFlowPotentialGradientsJacobians( false );

        }

        void residual::setd2MicroFlowdDrivingStressdStress( ){
            /*!
             * Set the Jacobian of the derivative of the micro flow potential w.r.t. the micro driving stress w.r.t. the micro stress
             */

            setFlowPotentialGradientsJacobians( false );

        }

        void residual::setd2MicroFlowdDrivingStressdF( ){
            /*!
             * Set the Jacobian of the derivative of the micro flow potential w.r.t. the micro driving stress w.r.t. the deformation gradient
             */

            setFlowPotentialGradientsJacobians( false );

        }

        void residual::setd2MicroFlowdDrivingStressdFn( ){
            /*!
             * Set the Jacobian of the derivative of the micro flow potential w.r.t. the micro driving stress w.r.t. the sub deformation gradients
             */

            setFlowPotentialGradientsJacobians( false );

        }

        void residual::setd2MicroGradientFlowdDrivingStressdStress( ){
            /*!
             * Set the Jacobian of the derivative of the micro gradient flow potential w.r.t. the micro-gradient driving stress w.r.t. the micro gradient stress
             */

            setFlowPotentialGradientsJacobians( false );

        }

        void residual::setd2MicroGradientFlowdDrivingStressdF( ){
            /*!
             * Set the Jacobian of the derivative of the micro gradient flow potential w.r.t. the micro-gradient driving stress w.r.t. the deformation gradient
             */

            setFlowPotentialGradientsJacobians( false );

        }

        void residual::setd2MicroGradientFlowdDrivingStressdFn( ){
            /*!
             * Set the Jacobian of the derivative of the micro gradient flow potential w.r.t. the micro-gradient driving stress w.r.t. the sub deformation gradients
             */

            setFlowPotentialGradientsJacobians( false );

        }

        void residual::setd2MicroGradientFlowdDrivingStressdChi( ){
            /*!
             * Set the Jacobian of the derivative of the micro gradient flow potential w.r.t. the micro-gradient driving stress w.r.t. the micro deformation
             */

            setFlowPotentialGradientsJacobians( false );

        }

        void residual::setd2MicroGradientFlowdDrivingStressdChin( ){
            /*!
             * Set the Jacobian of the derivative of the micro gradient flow potential w.r.t. the micro-gradient driving stress w.r.t. the sub micro deformation
             */

            setFlowPotentialGradientsJacobians( false );

        }

        void residual::setPreviousd2MacroFlowdDrivingStressdStress( ){
            /*!
             * Set the Jacobian of the previous derivative of the macro flow potential w.r.t. the macro driving stress w.r.t. the macro stress
             */

            setFlowPotentialGradientsJacobians( true );

        }

        void residual::setPreviousd2MacroFlowdDrivingStressdF( ){
            /*!
             * Set the Jacobian of the previous derivative of the macro flow potential w.r.t. the macro driving stress w.r.t. the deformation gradient
             */

            setFlowPotentialGradientsJacobians( true );

        }

        void residual::setPreviousd2MacroFlowdDrivingStressdFn( ){
            /*!
             * Set the Jacobian of the previous derivative of the macro flow potential w.r.t. the macro driving stress w.r.t. the sub deformation gradients
             */

            setFlowPotentialGradientsJacobians( true );

        }

        void residual::setPreviousd2MicroFlowdDrivingStressdStress( ){
            /*!
             * Set the Jacobian of the previous derivative of the micro flow potential w.r.t. the micro driving stress w.r.t. the micro stress
             */

            setFlowPotentialGradientsJacobians( true );

        }

        void residual::setPreviousd2MicroFlowdDrivingStressdF( ){
            /*!
             * Set the Jacobian of the previous derivative of the micro flow potential w.r.t. the micro driving stress w.r.t. the deformation gradient
             */

            setFlowPotentialGradientsJacobians( true );

        }

        void residual::setPreviousd2MicroFlowdDrivingStressdFn( ){
            /*!
             * Set the Jacobian of the previous derivative of the micro flow potential w.r.t. the micro driving stress w.r.t. the sub deformation gradients
             */

            setFlowPotentialGradientsJacobians( true );

        }

        void residual::setPreviousd2MicroGradientFlowdDrivingStressdStress( ){
            /*!
             * Set the Jacobian of the previous derivative of the micro gradient flow potential w.r.t. the micro-gradient driving stress w.r.t. the micro gradient stress
             */

            setFlowPotentialGradientsJacobians( true );

        }

        void residual::setPreviousd2MicroGradientFlowdDrivingStressdF( ){
            /*!
             * Set the Jacobian of the previous derivative of the micro gradient flow potential w.r.t. the micro-gradient driving stress w.r.t. the deformation gradient
             */

            setFlowPotentialGradientsJacobians( true );

        }

        void residual::setPreviousd2MicroGradientFlowdDrivingStressdFn( ){
            /*!
             * Set the Jacobian of the previous derivative of the micro gradient flow potential w.r.t. the micro-gradient driving stress w.r.t. the sub deformation gradients
             */

            setFlowPotentialGradientsJacobians( true );

        }

        void residual::setPreviousd2MicroGradientFlowdDrivingStressdChi( ){
            /*!
             * Set the Jacobian of the previous derivative of the micro gradient flow potential w.r.t. the micro-gradient driving stress w.r.t. the micro deformation
             */

            setFlowPotentialGradientsJacobians( true );

        }

        void residual::setPreviousd2MicroGradientFlowdDrivingStressdChin( ){
            /*!
             * Set the Jacobian of the previous derivative of the micro gradient flow potential w.r.t. the micro-gradient driving stress w.r.t. the sub micro deformation
             */

            setFlowPotentialGradientsJacobians( true );

        }

        void residual::setFlowPotentialGradientsJacobians( const bool isPrevious ){
            /*!
             * Set the Jacobians of the flow potential gradients
             *
             * \param isPrevious: A flag for whether to set the current (false) or previous (true) derivatives
             */

            const floatType   *macroCohesion;

            const floatType   *microCohesion;

            const floatVector *microGradientCohesion;

            const floatVector *macroDrivingStress;

            const floatVector *microDrivingStress;

            const floatVector *microGradientDrivingStress;

            const floatMatrix *dMacroDrivingStressdStress;

            const floatMatrix *dMacroDrivingStressdF;

            const floatMatrix *dMacroDrivingStressdFn;

            const floatMatrix *dMicroDrivingStressdStress;

            const floatMatrix *dMicroDrivingStressdF;

            const floatMatrix *dMicroDrivingStressdFn;

            const floatMatrix *dMicroGradientDrivingStressdStress;

            const floatMatrix *dMicroGradientDrivingStressdF;

            const floatMatrix *dMicroGradientDrivingStressdFn;

            const floatMatrix *dMicroGradientDrivingStressdChi;

            const floatMatrix *dMicroGradientDrivingStressdChin;

            floatVector precedingDeformationGradient;

            floatMatrix dPrecedingFdSubFs;

            const floatMatrix *dF1dF;

            const floatMatrix *dF1dFn;

            const floatVector *macroFlowParameters         = get_macroFlowParameters( );

            const floatVector *microFlowParameters         = get_microFlowParameters( );

            const floatVector *microGradientFlowParameters = get_microGradientFlowParameters( );

            if ( isPrevious ){

                precedingDeformationGradient = hydra->getPreviousPrecedingConfiguration( *getPlasticConfigurationIndex( ) );

                dPrecedingFdSubFs = hydra->getPreviousPrecedingConfigurationJacobian( *getPlasticConfigurationIndex( ) );

                dF1dF = hydra->get_previousdF1dF( );

                dF1dFn = hydra->get_previousdF1dFn( );

                macroCohesion                      = get_previousMacroCohesion( );

                microCohesion                      = get_previousMicroCohesion( );

                microGradientCohesion              = get_previousMicroGradientCohesion( );

                dMacroDrivingStressdStress         = get_previousdMacroDrivingStressdMacroStress( );

                dMicroDrivingStressdStress         = get_previousdSymmetricMicroDrivingStressdMicroStress( );

                dMicroGradientDrivingStressdStress = get_previousdHigherOrderDrivingStressdHigherOrderStress( );

                dMacroDrivingStressdF              = get_previousdMacroDrivingStressdF( );

                dMicroDrivingStressdF              = get_previousdSymmetricMicroDrivingStressdF( );

                dMicroGradientDrivingStressdF      = get_previousdHigherOrderDrivingStressdF( );

                dMicroGradientDrivingStressdChi    = get_previousdHigherOrderDrivingStressdChi( );

                dMacroDrivingStressdFn             = get_previousdMacroDrivingStressdFn( );

                dMicroDrivingStressdFn             = get_previousdSymmetricMicroDrivingStressdFn( );

                dMicroGradientDrivingStressdFn     = get_previousdHigherOrderDrivingStressdFn( );

                dMicroGradientDrivingStressdChin   = get_previousdHigherOrderDrivingStressdChin( );

                macroDrivingStress                 = get_previousMacroDrivingStress( );

                microDrivingStress                 = get_previousSymmetricMicroDrivingStress( );

                microGradientDrivingStress         = get_previousHigherOrderDrivingStress( );

            }
            else{

                precedingDeformationGradient = hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) );

                dPrecedingFdSubFs = hydra->getPrecedingConfigurationJacobian( *getPlasticConfigurationIndex( ) );

                dF1dF = hydra->get_dF1dF( );

                dF1dFn = hydra->get_dF1dFn( );

                macroCohesion                      = get_macroCohesion( );

                microCohesion                      = get_microCohesion( );

                microGradientCohesion              = get_microGradientCohesion( );

                dMacroDrivingStressdStress         = get_dMacroDrivingStressdMacroStress( );

                dMicroDrivingStressdStress         = get_dSymmetricMicroDrivingStressdMicroStress( );

                dMicroGradientDrivingStressdStress = get_dHigherOrderDrivingStressdHigherOrderStress( );

                dMacroDrivingStressdF              = get_dMacroDrivingStressdF( );

                dMicroDrivingStressdF              = get_dSymmetricMicroDrivingStressdF( );

                dMicroGradientDrivingStressdF      = get_dHigherOrderDrivingStressdF( );

                dMicroGradientDrivingStressdChi    = get_dHigherOrderDrivingStressdChi( );

                dMacroDrivingStressdFn             = get_dMacroDrivingStressdFn( );

                dMicroDrivingStressdFn             = get_dSymmetricMicroDrivingStressdFn( );

                dMicroGradientDrivingStressdFn     = get_dHigherOrderDrivingStressdFn( );

                dMicroGradientDrivingStressdChin   = get_dHigherOrderDrivingStressdChin( );

                macroDrivingStress                 = get_macroDrivingStress( );

                microDrivingStress                 = get_symmetricMicroDrivingStress( );

                microGradientDrivingStress         = get_higherOrderDrivingStress( );

            }

            // Construct the derivatives of the preceding F

            floatMatrix dPrecedingFdF( precedingDeformationGradient.size( ), floatVector( hydra->getDeformationGradient( )->size( ), 0 ) );

            floatMatrix dPrecedingFdFn( precedingDeformationGradient.size( ), floatVector( ( ( *hydra->getNumConfigurations( ) ) - 1 ) * hydra->getDeformationGradient( )->size( ), 0 ) );

            for ( unsigned int i = 0; i < hydra->getDeformationGradient( )->size( ); i++ ){

                for ( unsigned int j = 0; j < hydra->getDeformationGradient( )->size( ); j++ ){

                    for ( unsigned int k = 0; k < hydra->getDeformationGradient( )->size( ); k++ ){

                        dPrecedingFdF[ i ][ j ] += dPrecedingFdSubFs[ i ][ k ] * ( *dF1dF )[ k ][ j ];

                    }

                }

                for ( unsigned int j = 0; j < ( ( *hydra->getNumConfigurations( ) ) - 1 ) * hydra->getDeformationGradient( )->size( ); j++ ){

                    dPrecedingFdFn[ i ][ j ] = dPrecedingFdSubFs[ i ][ hydra->getDeformationGradient( )->size( ) + j ];

                    for ( unsigned int k = 0; k < hydra->getDeformationGradient( )->size( ); k++ ){

                        dPrecedingFdFn[ i ][ j ] += dPrecedingFdSubFs[ i ][ k ] * ( *dF1dFn )[ k ][ j ];

                    }

                }

            }

            floatType tempYield;

            floatVector tempVectorYield;

            floatType   dMacroFlowdCohesion, dMicroFlowdCohesion;

            floatVector dMacroFlowdDrivingStress, dMicroFlowdDrivingStress,
                        dMacroFlowdPrecedingF,    dMicroFlowdPrecedingF;

            floatMatrix dMicroGradientFlowdDrivingStress, dMicroGradientFlowdCohesion, dMicroGradientFlowdPrecedingF;

            floatMatrix d2MacroFlowdDrivingStress2,         d2MacroFlowdDrivingStressdPrecedingF,
                        d2MicroFlowdDrivingStress2,         d2MicroFlowdDrivingStressdPrecedingF,
                        d2MicroGradientFlowdDrivingStress2, d2MicroGradientFlowdDrivingStressdPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderDruckerPragerYieldEquation( *macroDrivingStress, *macroCohesion, precedingDeformationGradient,
                                                                                        ( *macroFlowParameters )[ 0 ], ( *macroFlowParameters )[ 1 ],
                                                                                        tempYield, dMacroFlowdDrivingStress, dMacroFlowdCohesion, dMacroFlowdPrecedingF,
                                                                                        d2MacroFlowdDrivingStress2, d2MacroFlowdDrivingStressdPrecedingF ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderDruckerPragerYieldEquation( *microDrivingStress, *microCohesion, precedingDeformationGradient,
                                                                                        ( *microFlowParameters )[ 0 ], ( *microFlowParameters )[ 1 ],
                                                                                        tempYield, dMicroFlowdDrivingStress, dMicroFlowdCohesion, dMicroFlowdPrecedingF,
                                                                                        d2MicroFlowdDrivingStress2, d2MicroFlowdDrivingStressdPrecedingF ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( computeHigherOrderDruckerPragerYieldEquation( *microGradientDrivingStress, *microGradientCohesion, precedingDeformationGradient,
                                                                                       ( *microGradientFlowParameters )[ 0 ], ( *microGradientFlowParameters )[ 1 ],
                                                                                       tempVectorYield, dMicroGradientFlowdDrivingStress, dMicroGradientFlowdCohesion, dMicroGradientFlowdPrecedingF,
                                                                                       d2MicroGradientFlowdDrivingStress2, d2MicroGradientFlowdDrivingStressdPrecedingF ) );

            floatMatrix d2MacroFlowdDrivingStressdMacroStress( dMacroFlowdDrivingStress.size( ), floatVector( macroDrivingStress->size( ) * macroDrivingStress->size( ), 0 ) );

            floatMatrix d2MicroFlowdDrivingStressdMicroStress( dMicroFlowdDrivingStress.size( ), floatVector( microDrivingStress->size( ) * microDrivingStress->size( ), 0 ) );

            floatMatrix d2MicroGradientFlowdDrivingStressdMicroGradientStress( dMicroGradientFlowdDrivingStress.size( ), floatVector( microGradientDrivingStress->size( ) * microGradientDrivingStress->size( ), 0 ) );

            floatMatrix d2MacroFlowdDrivingStressdFn( dMacroFlowdDrivingStress.size( ), floatVector( macroDrivingStress->size( ) * ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient.size( ), 0 ) );

            floatMatrix d2MicroFlowdDrivingStressdFn( dMicroFlowdDrivingStress.size( ), floatVector( microDrivingStress->size( ) * ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient.size( ), 0 ) );

            floatMatrix d2MicroGradientFlowdDrivingStressdFn( dMicroGradientFlowdDrivingStress.size( ), floatVector( microGradientDrivingStress->size( ) * ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient.size( ), 0 ) );

            floatMatrix d2MicroGradientFlowdDrivingStressdChin( dMicroGradientFlowdDrivingStress.size( ), floatVector( microGradientDrivingStress->size( ) * ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient.size( ), 0 ) );

            for ( unsigned int I = 0; I < d2MicroFlowdDrivingStress2.size( ); I++ ){

                for ( unsigned int J = 0; J < macroDrivingStress->size( ); J++ ){

                    for ( unsigned int K = 0; K < macroDrivingStress->size( ); K++ ){

                        for ( unsigned int L = 0; L < macroDrivingStress->size( ); L++ ){

                            d2MacroFlowdDrivingStressdMacroStress[ I ][ macroDrivingStress->size( ) * J + K ]
                                += d2MacroFlowdDrivingStress2[ I ][ macroDrivingStress->size( ) * J + L ] * ( *dMacroDrivingStressdStress )[ L ][ K ];

                            d2MicroFlowdDrivingStressdMicroStress[ I ][ microDrivingStress->size( ) * J + K ]
                                += d2MicroFlowdDrivingStress2[ I ][ microDrivingStress->size( ) * J + L ] * ( *dMicroDrivingStressdStress )[ L ][ K ];

                        }

                        for ( unsigned int L = 0; L < ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient.size( ); L++ ){

                            d2MacroFlowdDrivingStressdFn[ I ][ ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient.size( ) * J + L ]
                                += d2MacroFlowdDrivingStressdPrecedingF[ I ][ macroDrivingStress->size( ) * J + K ] * dPrecedingFdFn[ K ][ L ]
                                 + d2MacroFlowdDrivingStress2[ I ][ macroDrivingStress->size( ) * J + K ] * ( *dMacroDrivingStressdFn ) [ K ][ L ];

                            d2MicroFlowdDrivingStressdFn[ I ][ ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient.size( ) * J + L ]
                                += d2MicroFlowdDrivingStressdPrecedingF[ I ][ precedingDeformationGradient.size( ) * J + K ] * dPrecedingFdFn[ K ][ L ]
                                 + d2MicroFlowdDrivingStress2[ I ][ microDrivingStress->size( ) * J + K ] * ( *dMicroDrivingStressdFn ) [ K ][ L ];

                        }

                    }

                }

            }

            for ( unsigned int I = 0; I < d2MicroGradientFlowdDrivingStress2.size( ); I++ ){

                for ( unsigned int J = 0; J < microGradientDrivingStress->size( ); J++ ){

                    for ( unsigned int K = 0; K < microGradientDrivingStress->size( ); K++ ){

                        for ( unsigned int L = 0; L < microGradientDrivingStress->size( ); L++ ){

                            d2MicroGradientFlowdDrivingStressdMicroGradientStress[ I ][ microGradientDrivingStress->size( ) * J + K ]
                                += d2MicroGradientFlowdDrivingStress2[ I ][ microGradientDrivingStress->size( ) * J + L ] * ( *dMicroGradientDrivingStressdStress )[ L ][ K ];

                        }

                        for ( unsigned int L = 0; L < ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient.size( ); L++ ){

                            d2MicroGradientFlowdDrivingStressdFn[ I ][ ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient.size( ) * J + L ]
                                += d2MicroGradientFlowdDrivingStress2[ I ][ microGradientDrivingStress->size( ) * J + K ] * ( *dMicroGradientDrivingStressdFn ) [ K ][ L ];

                            d2MicroGradientFlowdDrivingStressdChin[ I ][ ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient.size( ) * J + L ]
                                += d2MicroGradientFlowdDrivingStress2[ I ][ microGradientDrivingStress->size( ) * J + K ] * ( *dMicroGradientDrivingStressdChin ) [ K ][ L ];

                        }

                    }

                    for ( unsigned int K = 0; K < precedingDeformationGradient.size( ); K++ ){

                        for ( unsigned int L = 0; L < ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient.size( ); L++ ){

                            d2MicroGradientFlowdDrivingStressdFn[ I ][ ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient.size( ) * J + L ]
                                += d2MicroGradientFlowdDrivingStressdPrecedingF[ I ][ precedingDeformationGradient.size( ) * J + K ] * dPrecedingFdFn[ K ][ L ];

                        }

                    }

                }

            }

            if ( isPrevious ){

                set_previousdMacroFlowdc( dMacroFlowdCohesion );

                set_previousdMicroFlowdc( dMicroFlowdCohesion );

                set_previousdMicroGradientFlowdc( dMicroGradientFlowdCohesion );

                set_previousdMacroFlowdDrivingStress( dMacroFlowdDrivingStress );

                set_previousdMicroFlowdDrivingStress( dMicroFlowdDrivingStress );

                set_previousdMicroGradientFlowdDrivingStress( dMicroGradientFlowdDrivingStress );

                set_previousd2MacroFlowdDrivingStressdStress( d2MacroFlowdDrivingStressdMacroStress );

                set_previousd2MicroFlowdDrivingStressdStress( d2MicroFlowdDrivingStressdMicroStress );

                set_previousd2MicroGradientFlowdDrivingStressdStress( d2MicroGradientFlowdDrivingStressdMicroGradientStress );

                set_previousd2MacroFlowdDrivingStressdFn( d2MacroFlowdDrivingStressdFn );

                set_previousd2MicroFlowdDrivingStressdFn( d2MicroFlowdDrivingStressdFn );

                set_previousd2MicroGradientFlowdDrivingStressdFn( d2MicroGradientFlowdDrivingStressdFn );

                set_previousd2MicroGradientFlowdDrivingStressdChin( d2MicroGradientFlowdDrivingStressdChin );

            }
            else{

                set_dMacroFlowdc( dMacroFlowdCohesion );

                set_dMicroFlowdc( dMicroFlowdCohesion );

                set_dMicroGradientFlowdc( dMicroGradientFlowdCohesion );

                set_dMacroFlowdDrivingStress( dMacroFlowdDrivingStress );

                set_dMicroFlowdDrivingStress( dMicroFlowdDrivingStress );

                set_dMicroGradientFlowdDrivingStress( dMicroGradientFlowdDrivingStress );

                set_d2MacroFlowdDrivingStressdStress( d2MacroFlowdDrivingStressdMacroStress );

                set_d2MicroFlowdDrivingStressdStress( d2MicroFlowdDrivingStressdMicroStress );

                set_d2MicroGradientFlowdDrivingStressdStress( d2MicroGradientFlowdDrivingStressdMicroGradientStress );

                set_d2MacroFlowdDrivingStressdFn( d2MacroFlowdDrivingStressdFn );

                set_d2MicroFlowdDrivingStressdFn( d2MicroFlowdDrivingStressdFn );

                set_d2MicroGradientFlowdDrivingStressdFn( d2MicroGradientFlowdDrivingStressdFn );

                set_d2MicroGradientFlowdDrivingStressdChin( d2MicroGradientFlowdDrivingStressdChin );

            }

        }

        void residual::setMacroCohesion( ){
            /*!
             * Set the macro cohesion
             */

            setCohesions( false );

        }

        void residual::setMicroCohesion( ){
            /*!
             * Set the micro cohesion
             */

            setCohesions( false );

        }

        void residual::setMicroGradientCohesion( ){
            /*!
             * Set the micro gradient cohesion
             */

            setCohesions( false );

        }

        void residual::setPreviousMacroCohesion( ){
            /*!
             * Set the previous macro cohesion
             */

            setCohesions( true );

        }

        void residual::setPreviousMicroCohesion( ){
            /*!
             * Set the previous macro cohesion
             */

            setCohesions( true );

        }

        void residual::setPreviousMicroGradientCohesion( ){
            /*!
             * Set the micro gradient cohesion
             */

            setCohesions( true );

        }

        void residual::setCohesions( const bool isPrevious ){
            /*!
             * Set the values of the cohesion
             * 
             * \param isPrevious: Flag for whether to compute the current (false) or previous (true) cohesions
             */

            const floatVector *plasticStrainLikeISVs;

            if ( isPrevious ){

                plasticStrainLikeISVs = get_previousPlasticStrainLikeISVs( );

            }
            else{

                plasticStrainLikeISVs = get_plasticStrainLikeISVs( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( get_macroHardeningParameters( )->size( ) != 2 ){
    
                    throw std::runtime_error( "The micro hardening parameters must have a length of 2 rather than " + std::to_string( get_macroHardeningParameters( )->size( ) ) );
    
                }
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( get_microHardeningParameters( )->size( ) != 2 ){
    
                    throw std::runtime_error( "The micro hardening parameters must have a length of 2 rather than " + std::to_string( get_microHardeningParameters( )->size( ) ) );
    
                }
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( get_microGradientHardeningParameters( )->size( ) != 2 ){
    
                    throw std::runtime_error( "The micro hardening parameters must have a length of 2 rather than " + std::to_string( get_microGradientHardeningParameters( )->size( ) ) );
    
                }
            )

            floatType macroCohesion           = ( *get_macroHardeningParameters( ) )[ 0 ] + ( *get_macroHardeningParameters( ) )[ 1 ] * ( *plasticStrainLikeISVs )[ 0 ];

            floatType microCohesion           = ( *get_microHardeningParameters( ) )[ 0 ] + ( *get_microHardeningParameters( ) )[ 1 ] * ( *plasticStrainLikeISVs )[ 1 ];

            floatVector microGradientCohesion = ( *get_microGradientHardeningParameters( ) )[ 0 ] + ( *get_microGradientHardeningParameters( ) )[ 1 ] * floatVector( plasticStrainLikeISVs->begin( ) + 2,
                                                                                                                                                                     plasticStrainLikeISVs->end( ) );

            if ( isPrevious ){

                set_previousMacroCohesion( macroCohesion );

                set_previousMicroCohesion( microCohesion );

                set_previousMicroGradientCohesion( microGradientCohesion );

            }
            else{

                set_macroCohesion( macroCohesion );

                set_microCohesion( microCohesion );

                set_microGradientCohesion( microGradientCohesion );

            }

        }

        void residual::setdMacroCohesiondStateVariables( ){
            /*!
             * Set the jacobian of the macro cohesion w.r.t. the nonlinear state variables
             */

            setCohesionsJacobians( false );

        }

        void residual::setdMicroCohesiondStateVariables( ){
            /*!
             * Set the jacobian of the micro cohesion w.r.t. the nonlinear state variables
             */

            setCohesionsJacobians( false );

        }

        void residual::setdMicroGradientCohesiondStateVariables( ){
            /*!
             * Set the jacobian of the micro gradient cohesion w.r.t. the nonlinear state variables
             */

            setCohesionsJacobians( false );

        }

        void residual::setPreviousdMacroCohesiondStateVariables( ){
            /*!
             * Set the jacobians of the previous macro cohesion w.r.t. the nonlinear state variables
             */

            setCohesionsJacobians( true );

        }

        void residual::setPreviousdMicroCohesiondStateVariables( ){
            /*!
             * Set the jacobians of the previous macro cohesion w.r.t. the nonlinear state variables
             */

            setCohesionsJacobians( true );

        }

        void residual::setPreviousdMicroGradientCohesiondStateVariables( ){
            /*!
             * Set the jacobians of the micro gradient cohesion w.r.t. the nonlinear state variables
             */

            setCohesionsJacobians( true );

        }

        void residual::setCohesionsJacobians( const bool isPrevious ){
            /*!
             * Set the values of the cohesion and their Jacobians
             * 
             * \param isPrevious: Flag for whether to compute the current (false) or previous (true) cohesions
             */

            const floatVector *plasticStrainLikeISVs;

            if ( isPrevious ){

                plasticStrainLikeISVs = get_previousPlasticStrainLikeISVs( );

            }
            else{

                plasticStrainLikeISVs = get_plasticStrainLikeISVs( );

            }

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( get_macroHardeningParameters( )->size( ) != 2 ){
    
                    throw std::runtime_error( "The micro hardening parameters must have a length of 2 rather than " + std::to_string( get_macroHardeningParameters( )->size( ) ) );
    
                }
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( get_microHardeningParameters( )->size( ) != 2 ){
    
                    throw std::runtime_error( "The micro hardening parameters must have a length of 2 rather than " + std::to_string( get_microHardeningParameters( )->size( ) ) );
    
                }
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( get_microGradientHardeningParameters( )->size( ) != 2 ){
    
                    throw std::runtime_error( "The micro hardening parameters must have a length of 2 rather than " + std::to_string( get_microGradientHardeningParameters( )->size( ) ) );
    
                }
            )

            floatVector dMacroCohesiondISVs( get_plasticStateVariables( )->size( ), 0 );

            floatVector dMicroCohesiondISVs( get_plasticStateVariables( )->size( ), 0 );

            floatMatrix dMicroGradientCohesiondISVs( get_plasticStrainLikeISVs( )->size( ) - 2, floatVector( get_plasticStateVariables( )->size( ), 0 ) );

            floatType macroCohesion           = ( *get_macroHardeningParameters( ) )[ 0 ] + ( *get_macroHardeningParameters( ) )[ 1 ] * ( *plasticStrainLikeISVs )[ 0 ];

            dMacroCohesiondISVs[ get_plasticMultipliers( )->size( ) + 0 ] = ( *get_macroHardeningParameters( ) )[ 1 ];

            floatType microCohesion           = ( *get_microHardeningParameters( ) )[ 0 ] + ( *get_microHardeningParameters( ) )[ 1 ] * ( *plasticStrainLikeISVs )[ 1 ];

            dMicroCohesiondISVs[ get_plasticMultipliers( )->size( ) + 1 ] = ( *get_microHardeningParameters( ) )[ 1 ];

            floatVector microGradientCohesion = ( *get_microGradientHardeningParameters( ) )[ 0 ] + ( *get_microGradientHardeningParameters( ) )[ 1 ] * floatVector( plasticStrainLikeISVs->begin( ) + 2,
                                                                                                                                                                     plasticStrainLikeISVs->end( ) );

            for ( unsigned int i = 2; i < plasticStrainLikeISVs->size( ); i++ ){

                dMicroGradientCohesiondISVs[ i - 2 ][ get_plasticMultipliers( )->size( ) + i ] = ( *get_microGradientHardeningParameters( ) )[ 1 ];

            }

            if ( isPrevious ){

                set_previousMacroCohesion( macroCohesion );

                set_previousMicroCohesion( microCohesion );

                set_previousMicroGradientCohesion( microGradientCohesion );

                set_previousdMacroCohesiondStateVariables( dMacroCohesiondISVs );

                set_previousdMicroCohesiondStateVariables( dMicroCohesiondISVs );

                set_previousdMicroGradientCohesiondStateVariables( dMicroGradientCohesiondISVs );

            }
            else{

                set_macroCohesion( macroCohesion );

                set_microCohesion( microCohesion );

                set_microGradientCohesion( microGradientCohesion );

                set_dMacroCohesiondStateVariables( dMacroCohesiondISVs );

                set_dMicroCohesiondStateVariables( dMicroCohesiondISVs );

                set_dMicroGradientCohesiondStateVariables( dMicroGradientCohesiondISVs );

            }

        }

    }

}
