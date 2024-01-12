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

        void computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMacroVelocityGradient ){
            /*!
             * Compute the plastic macro velocity gradient in the intermediate configuration.
             *
             * \f$\bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]\f$
             *
             * \param &macroGamma: The macro plastic multiplier.
             * \param &microGamma: The micro plastic multiplier.
             * \param &inverseElasticRightCauchyGreen: The inverse of the elastic right Cauchy-Green deformation tensor.
             * \param &macroFlowDirection: The flow direction of the macro plasticity.
             * \param &microFlowDirection: The flow direction of the micro plasticity.
             * \param &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
             */

            //Assume 3D
            unsigned int dim = 3;

            TARDIGRADE_ERROR_TOOLS_CATCH( 
                if ( inverseElasticRightCauchyGreen.size() != dim * dim ){
                    throw std::runtime_error( "The inverse elastic right Cauchy-Green deformation tensor must be 3D" );
                }
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( macroFlowDirection.size() != dim * dim ){
                    throw std::runtime_error( "The macro flow direction tensor must be 3D" );
                }
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( microFlowDirection.size() != dim * dim ){
                    throw std::runtime_error( "The micro flow direction tensor must be 3D" );
                }
            )

            //Compute the macro-scale velocity gradient
            plasticMacroVelocityGradient = variableVector( dim * dim, 0 );

            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                    for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                        plasticMacroVelocityGradient[ dim * Bb + Kb ]
                            += inverseElasticRightCauchyGreen[ dim * Bb + Lb ]
                             * ( macroGamma * macroFlowDirection[ dim * Kb + Lb ]
                             +   microGamma * microFlowDirection[ dim * Kb + Lb ] );

                    }

                }

            }

        }

        void computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMacroVelocityGradient,
                                                  variableVector &dPlasticMacroLdMacroGamma,
                                                  variableVector &dPlasticMacroLdMicroGamma ){
            /*!
             * Compute the plastic macro velocity gradient in the intermediate configuration.
             *
             * \f$ \bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right] \f$
             *
             * \param &macroGamma: The macro plastic multiplier.
             * \param &microGamma: The micro plastic multiplier.
             * \param &inverseElasticRightCauchyGreen: The inverse of the elastic right Cauchy-Green deformation tensor.
             * \param &macroFlowDirection: The flow direction of the macro plasticity.
             * \param &microFlowDirection: The flow direction of the micro plasticity.
             * \param &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
             * \param &dPlasticMacroLdMacroGamma: The Jacobian of the plastic velocity gradient w.r.t. the 
             *     macro plastic multiplier.
             * \param &dPlasticMacroLdMicroGamma: The Jacobian of the plastic velocity gradient w.r.t. the 
             *     micro plastic multiplier.
             */

            //Assume 3D
            unsigned int dim = 3;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseElasticRightCauchyGreen,
                                                     macroFlowDirection, microFlowDirection, plasticMacroVelocityGradient )
            )

            dPlasticMacroLdMacroGamma = variableVector( dim * dim, 0 );
            dPlasticMacroLdMicroGamma = variableVector( dim * dim, 0 );

            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                    for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                        dPlasticMacroLdMacroGamma[ dim * Bb + Kb ] += inverseElasticRightCauchyGreen[ dim * Bb + Lb ]
                                                                    * macroFlowDirection[ dim * Kb + Lb ];

                        dPlasticMacroLdMicroGamma[ dim * Bb + Kb ] += inverseElasticRightCauchyGreen[ dim * Bb + Lb ]
                                                                    * microFlowDirection[ dim * Kb + Lb ];

                    }

                }

            }

        }

        void computePlasticMacroVelocityGradient( const variableType &macroGamma, const variableType &microGamma,
                                                  const variableVector &inverseElasticRightCauchyGreen,
                                                  const variableVector &macroFlowDirection,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMacroVelocityGradient,
                                                  variableVector &dPlasticMacroLdMacroGamma,
                                                  variableVector &dPlasticMacroLdMicroGamma,
                                                  variableMatrix &dPlasticMacroLdElasticRCG,
                                                  variableMatrix &dPlasticMacroLdMacroFlowDirection,
                                                  variableMatrix &dPlasticMacroLdMicroFlowDirection ){
            /*!
             * Compute the plastic macro velocity gradient in the intermediate configuration.
             *
             * \f$\bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]\f$
             *
             * :param const variableType &macroGamma: The macro plastic multiplier.
             * :param const variableType &microGamma: The micro plastic multiplier.
             * :param const variableVector &inverseElasticRightCauchyGreen: The inverse of the elastic right Cauchy-Green deformation tensor.
             * :param const variableVector &macroFlowDirection: The flow direction of the macro plasticity.
             * :param const variableVector &microFlowDirection: The flow direction of the micro plasticity.
             * :param variableVector &plasticMacroVelocityGradient: The plastic velocity gradient for the macro plastic deformation.
             * :param variableVector &dPlasticMacroLdMacroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
             *     macro plastic multiplier.
             * :param variableVector &dPlasticMacroLdMicroGamma: The Jacobian of the plastic macro velocity gradient w.r.t. the 
             *     micro plastic multiplier.
             * :param variableMatrix &dPlasticMacroLdElasticRCG: The Jacobian of the plastic macro velocity gradient w.r.t. the 
             *     elastic right Cauchy-Green deformation tensor.
             * :param variableMatrix &dPlasticMacroLdMacroFlowDirection: The Jacobian of the plastic macro velocity gradient w.r.t. the 
             *     macro flow direction tensor.
             * :param variableMatrix &dPlasticMacroLdMicroFlowDirection: The Jacobian of the plastic macro velocity gradient w.r.t. the 
             *     micro flow direction tensor.
             */

            //Assume 3D
            unsigned int dim = 3;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMacroVelocityGradient( macroGamma, microGamma, inverseElasticRightCauchyGreen,
                                                     macroFlowDirection, microFlowDirection, plasticMacroVelocityGradient,
                                                     dPlasticMacroLdMacroGamma, dPlasticMacroLdMicroGamma );
            )

            constantVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );

            dPlasticMacroLdElasticRCG = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
            dPlasticMacroLdMacroFlowDirection = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );
            dPlasticMacroLdMicroFlowDirection = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );

            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                    for ( unsigned int Ob = 0; Ob < dim; Ob++ ){

                        for ( unsigned int Pb = 0; Pb < dim; Pb++ ){

                            dPlasticMacroLdElasticRCG[ dim * Bb + Kb ][ dim * Ob + Pb ]
                                -= inverseElasticRightCauchyGreen[ dim * Bb + Ob ]
                                 * plasticMacroVelocityGradient[ dim * Pb + Kb ];

                            dPlasticMacroLdMacroFlowDirection[ dim * Bb + Kb ][ dim * Ob + Pb ]
                                += macroGamma * inverseElasticRightCauchyGreen[ dim * Bb + Pb ] * eye[ dim * Kb + Ob ];

                            dPlasticMacroLdMicroFlowDirection[ dim * Bb + Kb ][ dim * Ob + Pb ]
                                += microGamma * inverseElasticRightCauchyGreen[ dim * Bb + Pb ] * eye[ dim * Kb + Ob ];

                        }

                    }

                }

            }

        }

        void computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient ){
            /*!
             * Compute the plastic micro velocity gradient
             *
             *  \f$\bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]\f$
             *
             *  Note: This function is used in conjunction with other functions. If it is used by itself, the user must guarantee 
             *        that elasticPsi and inverseElasticPsi are actually inverses of each-other. This is not checked in code.
             *
             * :param const variableType &microGamma: The micro plastic multiplier.
             * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
             * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
             * :param const variableVector &inverseElasticPsi: The inverse of the elastic micro deformation measure Psi.
             * :param const variableVector &microFlowDirection: The micro plastic flow direction.
             * :param variableVector &plasticMicroVelocityGradient: The plastic micro velocity gradient.
             */

            //Assume 3D
            unsigned int dim = 3;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( elasticMicroRightCauchyGreen.size() != dim * dim ){
                    throw std::runtime_error( "The elastic micro right Cauchy-Green deformation tensor is not 3D" );
                }
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( elasticPsi.size() != dim * dim ){
                    throw std::runtime_error( "The elastic micro deformation tensor Psi is not 3D" );
                }
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( inverseElasticPsi.size() != dim * dim ){
                    throw std::runtime_error( "The inverse of the elastic micro deformation tensor Psi is not 3D" );
                }
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( microFlowDirection.size() != dim * dim ){
                    throw std::runtime_error( "The micro flow direction of the elastic micro plastic flow direction is not 3D" );
                }
            )

            plasticMicroVelocityGradient = variableVector( dim * dim, 0 );

            //NOTE: I'm making the second inverse elastic Psi be the transpose of what was done previously.
            //      I think the way it was is a bug since it isn't consistent with the form in my dissertation.
            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                    for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                        for ( unsigned int Nb = 0; Nb < dim; Nb++ ){

                            for ( unsigned int Eb = 0; Eb < dim; Eb++ ){

                                plasticMicroVelocityGradient[ dim * Bb + Kb ]
                                    += microGamma
                                     * inverseElasticPsi[ dim * Bb + Lb ]
                                     * microFlowDirection[ dim * Eb + Lb ]
                                     * inverseElasticPsi[ dim * Nb + Eb ]
                                     * elasticMicroRightCauchyGreen[ dim * Nb + Kb ];

                            }

                        }

                    }

                }

            }

        }

        void computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient,
                                                  variableVector &dPlasticMicroLdMicroGamma ){
            /*!
             * Compute the plastic micro velocity gradient
             *
             *  \f$\bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]\f$
             *
             *  Note: This function is used in conjunction with other functions. If it is used by itself, the user must guarantee 
             *        that elasticPsi and inverseElasticPsi are actually inverses of each-other. This is not checked in code.
             *
             * :param const variableType &microGamma: The micro plastic multiplier.
             * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
             * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
             * :param const variableVector &inverseElasticPsi: The inverse of the elastic micro deformation measure Psi.
             * :param const variableVector &microFlowDirection: The micro plastic flow direction.
             * :param variableVector &plasticMicroVelocityGradient: The plastic micro velocity gradient.
             * :param variableVector &dPlasticMicroLdMicroGamma: The Jacobian of the plastic micro velocity gradient
             *     w.r.t. the micro plastic multiplier.
             */

            //Assume 3D
            unsigned int dim = 3;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( elasticMicroRightCauchyGreen.size() != dim * dim ){
                    throw std::runtime_error( "The elastic micro right Cauchy-Green deformation tensor is not 3D" );
                }
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( elasticPsi.size() != dim * dim ){
                    throw std::runtime_error( "The elastic micro deformation tensor Psi is not 3D" );
                }
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( inverseElasticPsi.size() != dim * dim ){
                    throw std::runtime_error( "The inverse of the elastic micro deformation tensor Psi is not 3D" );
                }
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                if ( microFlowDirection.size() != dim * dim ){
                    throw std::runtime_error( "The micro flow direction of the elastic micro plastic flow direction is not 3D" );
                }
            )

            plasticMicroVelocityGradient = variableVector( dim * dim, 0 );
            dPlasticMicroLdMicroGamma = variableVector( dim * dim, 0 );

            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                    for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                        for ( unsigned int Nb = 0; Nb < dim; Nb++ ){

                            for ( unsigned int Eb = 0; Eb < dim; Eb++ ){

                                dPlasticMicroLdMicroGamma[ dim * Bb + Kb ]
                                    += inverseElasticPsi[ dim * Bb + Lb ]
                                     * microFlowDirection[ dim * Eb + Lb ]
                                     * inverseElasticPsi[ dim * Nb + Eb ]
                                     * elasticMicroRightCauchyGreen[ dim * Nb + Kb ];

                            }

                        }

                    }

                    plasticMicroVelocityGradient[ dim * Bb + Kb ] = microGamma * dPlasticMicroLdMicroGamma[ dim * Bb + Kb ];

                }

            }

        }

        void computePlasticMicroVelocityGradient( const variableType &microGamma, const variableVector &elasticMicroRightCauchyGreen,
                                                  const variableVector &elasticPsi, const variableVector &inverseElasticPsi,
                                                  const variableVector &microFlowDirection,
                                                  variableVector &plasticMicroVelocityGradient,
                                                  variableVector &dPlasticMicroLdMicroGamma,
                                                  variableMatrix &dPlasticMicroLdElasticMicroRCG,
                                                  variableMatrix &dPlasticMicroLdElasticPsi,
                                                  variableMatrix &dPlasticMicroLdMicroFlowDirection ){
            /*!
             * Compute the plastic micro velocity gradient
             *
             *  \f$\bar{ L }_{ \bar{B} \bar{K} }^p = \bar{ C }_{ \bar{B} \bar{L} }^{e, -1} \left[ \dot{ \bar{ \gamma } } \frac{ \partial \bar{G}^{\text{MACRO}} }{ \partial \bar{S}_{ \bar{K} \bar{L} } + \dot{ \bar{ \gamma } }^{\chi} \frac{ \partial \bar{G}^{\chi} }{ \partial \bar{ \Sigma }_{ \bar{K} \bar{L} } \right]\f$
             *
             *  Note: This function is used in conjunction with other functions. If it is used by itself, the user must guarantee 
             *        that elasticPsi and inverseElasticPsi are actually inverses of each-other. This is not checked in code.
             *
             * :param const variableType &microGamma: The micro plastic multiplier.
             * :param const variableVector &elasticMicroRightCauchyGreen: The elastic micro right Cauchy-Green deformation tensor.
             * :param const variableVector &elasticPsi: The elastic micro deformation measure Psi.
             * :param const variableVector &inverseElasticPsi: The inverse of the elastic micro deformation measure Psi.
             * :param const variableVector &microFlowDirection: The micro plastic flow direction.
             * :param variableVector &plasticMicroVelocityGradient: The plastic micro velocity gradient.
             * :param variableVector &dPlasticMicroLdMicroGamma: The Jacobian of the plastic micro velocity gradient
             *     w.r.t. the micro plastic multiplier.
             * :param variableMatrix &dPlasticMicroLdElasticMicroRCG: The Jacobian of the plastic micro velocity gradient
             *     w.r.t. the micro right Cauchy-Green deformation tensor.
             * :param variableMatrix &dPlasticMicroLdElasticPsi: The Jacobian of the plastic micro velocity gradient
             *     w.r.t. the micro deformation measure Psi.
             * :param variableMatrix &dPlasticMicroLdMicroFlowDirection: The Jacobian of the plastic micro velocity gradient
             *     w.r.t. the micro flow direction.
             */

            //Assume 3D
            unsigned int dim = 3;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMicroVelocityGradient( microGamma, elasticMicroRightCauchyGreen,
                                                     elasticPsi, inverseElasticPsi, microFlowDirection,
                                                     plasticMicroVelocityGradient, dPlasticMicroLdMicroGamma );
            )

            constantVector eye( dim * dim );
            tardigradeVectorTools::eye( eye );

            //Assemble the Jacobians
            dPlasticMicroLdElasticMicroRCG = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );

            dPlasticMicroLdElasticPsi = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );

            dPlasticMicroLdMicroFlowDirection = variableMatrix( dim * dim, variableVector( dim * dim, 0 ) );

            for ( unsigned int Bb = 0; Bb < dim; Bb++ ){

                for ( unsigned int Kb = 0; Kb < dim; Kb++ ){

                    for ( unsigned int Ob = 0; Ob < dim; Ob++ ){

                        for ( unsigned int Pb = 0; Pb < dim; Pb++ ){

                            dPlasticMicroLdElasticPsi[ dim * Bb + Kb ][ dim * Ob + Pb ]
                                -= inverseElasticPsi[ dim * Bb + Ob ] * plasticMicroVelocityGradient[ dim * Pb + Kb ];

                            for ( unsigned int Lb = 0; Lb < dim; Lb++ ){

                                dPlasticMicroLdMicroFlowDirection[ dim * Bb + Kb ][ dim * Ob + Pb ]
                                    += microGamma * inverseElasticPsi[ dim * Bb + Pb ]
                                     * inverseElasticPsi[ dim * Lb + Ob ]
                                     * elasticMicroRightCauchyGreen[ dim * Lb + Kb ];

                                for ( unsigned int Eb = 0; Eb < dim; Eb++ ){

                                    dPlasticMicroLdElasticMicroRCG[ dim * Bb + Kb ][ dim * Ob + Pb ]
                                        += microGamma * inverseElasticPsi[ dim * Bb + Lb ] * microFlowDirection[ dim * Eb + Lb ]
                                         * inverseElasticPsi[ dim * Ob + Eb ] * eye[ dim * Kb + Pb ];
                                    for ( unsigned int Nb = 0; Nb < dim; Nb++ ){

                                        dPlasticMicroLdElasticPsi[ dim * Bb + Kb ][ dim * Ob + Pb ]
                                            -= microGamma * inverseElasticPsi[ dim * Bb + Lb ] * microFlowDirection[ dim * Eb + Lb ]
                                             * inverseElasticPsi[ dim * Nb + Ob ] * inverseElasticPsi[ dim * Pb + Eb ]
                                             * elasticMicroRightCauchyGreen[ dim * Nb + Kb ];
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

            const floatVector *precedingDeformationGradient;

            const floatVector *macroFlowParameters         = get_macroFlowParameters( );

            const floatVector *microFlowParameters         = get_microFlowParameters( );

            const floatVector *microGradientFlowParameters = get_microGradientFlowParameters( );

            if ( isPrevious ){

                precedingDeformationGradient = get_previousPrecedingDeformationGradient( );

                macroCohesion                = get_previousMacroCohesion( );

                microCohesion                = get_previousMicroCohesion( );

                microGradientCohesion        = get_previousMicroGradientCohesion( );

                macroDrivingStress           = get_previousMacroDrivingStress( );

                microDrivingStress           = get_previousSymmetricMicroDrivingStress( );

                microGradientDrivingStress   = get_previousHigherOrderDrivingStress( );

            }
            else{

                precedingDeformationGradient = get_precedingDeformationGradient( );

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

            TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderDruckerPragerYieldEquation( *macroDrivingStress, *macroCohesion, *precedingDeformationGradient,
                                                                                        ( *macroFlowParameters )[ 0 ], ( *macroFlowParameters )[ 1 ],
                                                                                        tempYield, dMacroFlowdDrivingStress, dMacroFlowdCohesion, dMacroFlowdPrecedingF ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderDruckerPragerYieldEquation( *microDrivingStress, *microCohesion, *precedingDeformationGradient,
                                                                                        ( *microFlowParameters )[ 0 ], ( *microFlowParameters )[ 1 ],
                                                                                        tempYield, dMicroFlowdDrivingStress, dMicroFlowdCohesion, dMicroFlowdPrecedingF ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( computeHigherOrderDruckerPragerYieldEquation( *microGradientDrivingStress, *microGradientCohesion, *precedingDeformationGradient,
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

            const floatVector *precedingDeformationGradient;

            const floatMatrix *dPrecedingFdF;

            const floatMatrix *dPrecedingFdFn;

            const floatVector *macroFlowParameters         = get_macroFlowParameters( );

            const floatVector *microFlowParameters         = get_microFlowParameters( );

            const floatVector *microGradientFlowParameters = get_microGradientFlowParameters( );

            if ( isPrevious ){

                precedingDeformationGradient = get_previousPrecedingDeformationGradient( );

                dPrecedingFdF                      = get_previousdPrecedingDeformationGradientdF( );

                dPrecedingFdFn                     = get_previousdPrecedingDeformationGradientdFn( );

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

                precedingDeformationGradient = get_precedingDeformationGradient( );

                dPrecedingFdF                      = get_dPrecedingDeformationGradientdF( );

                dPrecedingFdFn                     = get_dPrecedingDeformationGradientdFn( );

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

            floatType tempYield;

            floatVector tempVectorYield;

            floatType   dMacroFlowdCohesion, dMicroFlowdCohesion;

            floatVector dMacroFlowdDrivingStress, dMicroFlowdDrivingStress,
                        dMacroFlowdPrecedingF,    dMicroFlowdPrecedingF;

            floatMatrix dMicroGradientFlowdDrivingStress, dMicroGradientFlowdCohesion, dMicroGradientFlowdPrecedingF;

            floatMatrix d2MacroFlowdDrivingStress2,         d2MacroFlowdDrivingStressdPrecedingF,
                        d2MicroFlowdDrivingStress2,         d2MicroFlowdDrivingStressdPrecedingF,
                        d2MicroGradientFlowdDrivingStress2, d2MicroGradientFlowdDrivingStressdPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderDruckerPragerYieldEquation( *macroDrivingStress, *macroCohesion, *precedingDeformationGradient,
                                                                                        ( *macroFlowParameters )[ 0 ], ( *macroFlowParameters )[ 1 ],
                                                                                        tempYield, dMacroFlowdDrivingStress, dMacroFlowdCohesion, dMacroFlowdPrecedingF,
                                                                                        d2MacroFlowdDrivingStress2, d2MacroFlowdDrivingStressdPrecedingF ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderDruckerPragerYieldEquation( *microDrivingStress, *microCohesion, *precedingDeformationGradient,
                                                                                        ( *microFlowParameters )[ 0 ], ( *microFlowParameters )[ 1 ],
                                                                                        tempYield, dMicroFlowdDrivingStress, dMicroFlowdCohesion, dMicroFlowdPrecedingF,
                                                                                        d2MicroFlowdDrivingStress2, d2MicroFlowdDrivingStressdPrecedingF ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( computeHigherOrderDruckerPragerYieldEquation( *microGradientDrivingStress, *microGradientCohesion, *precedingDeformationGradient,
                                                                                       ( *microGradientFlowParameters )[ 0 ], ( *microGradientFlowParameters )[ 1 ],
                                                                                       tempVectorYield, dMicroGradientFlowdDrivingStress, dMicroGradientFlowdCohesion, dMicroGradientFlowdPrecedingF,
                                                                                       d2MicroGradientFlowdDrivingStress2, d2MicroGradientFlowdDrivingStressdPrecedingF ) );

            floatMatrix d2MacroFlowdDrivingStressdMacroStress( dMacroFlowdDrivingStress.size( ), floatVector( macroDrivingStress->size( ), 0 ) );

            floatMatrix d2MicroFlowdDrivingStressdMicroStress( dMicroFlowdDrivingStress.size( ), floatVector( microDrivingStress->size( ), 0 ) );

            floatMatrix d2MicroGradientFlowdDrivingStressdMicroGradientStress( dMicroGradientFlowdDrivingStress.size( ), floatVector( microGradientDrivingStress->size( ) * microGradientDrivingStress->size( ), 0 ) );

            floatMatrix d2MacroFlowdDrivingStressdF( dMacroFlowdDrivingStress.size( ), floatVector( precedingDeformationGradient->size( ), 0 ) );

            floatMatrix d2MicroFlowdDrivingStressdF( dMicroFlowdDrivingStress.size( ), floatVector( precedingDeformationGradient->size( ), 0 ) );

            floatMatrix d2MicroGradientFlowdDrivingStressdF( dMicroGradientFlowdDrivingStress.size( ), floatVector( microGradientDrivingStress->size( ) * precedingDeformationGradient->size( ), 0 ) );

            floatMatrix d2MicroGradientFlowdDrivingStressdChi( dMicroGradientFlowdDrivingStress.size( ), floatVector( microGradientDrivingStress->size( ) * precedingDeformationGradient->size( ), 0 ) );

            floatMatrix d2MacroFlowdDrivingStressdFn( dMacroFlowdDrivingStress.size( ), floatVector( ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient->size( ), 0 ) );

            floatMatrix d2MicroFlowdDrivingStressdFn( dMicroFlowdDrivingStress.size( ), floatVector( ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient->size( ), 0 ) );

            floatMatrix d2MicroGradientFlowdDrivingStressdFn( dMicroGradientFlowdDrivingStress.size( ), floatVector( microGradientDrivingStress->size( ) * ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient->size( ), 0 ) );

            floatMatrix d2MicroGradientFlowdDrivingStressdChin( dMicroGradientFlowdDrivingStress.size( ), floatVector( microGradientDrivingStress->size( ) * ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient->size( ), 0 ) );

            for ( unsigned int I = 0; I < d2MicroFlowdDrivingStress2.size( ); I++ ){

                for ( unsigned int J = 0; J < macroDrivingStress->size( ); J++ ){

                    for ( unsigned int K = 0; K < macroDrivingStress->size( ); K++ ){

                        d2MacroFlowdDrivingStressdMacroStress[ I ][ J ]
                            += d2MacroFlowdDrivingStress2[ I ][ K ] * ( *dMacroDrivingStressdStress )[ K ][ J ];

                        d2MicroFlowdDrivingStressdMicroStress[ I ][ J ]
                            += d2MicroFlowdDrivingStress2[ I ][ K ] * ( *dMicroDrivingStressdStress )[ K ][ J ];

                        d2MacroFlowdDrivingStressdF[ I ][ J ]
                            += d2MacroFlowdDrivingStress2[ I ][ K ] * ( *dMacroDrivingStressdF )[ K ][ J ]
                             + d2MacroFlowdDrivingStressdPrecedingF[ I ][ K ] * ( *dPrecedingFdF )[ K ][ J ];

                        d2MicroFlowdDrivingStressdF[ I ][ J ]
                            += d2MicroFlowdDrivingStress2[ I ][ K ] * ( *dMicroDrivingStressdF )[ K ][ J ]
                             + d2MicroFlowdDrivingStressdPrecedingF[ I ][ K ] * ( *dPrecedingFdF )[ K ][ J ];

                    }

                    for ( unsigned int K = 0; K < ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient->size( ); K++ ){

                        d2MacroFlowdDrivingStressdFn[ I ][ K ]
                            += d2MacroFlowdDrivingStressdPrecedingF[ I ][ J ] * ( *dPrecedingFdFn )[ J ][ K ]
                             + d2MacroFlowdDrivingStress2[ I ][ J ] * ( *dMacroDrivingStressdFn )[ J ][ K ];

                        d2MicroFlowdDrivingStressdFn[ I ][ K ]
                            += d2MicroFlowdDrivingStressdPrecedingF[ I ][ J ] * ( *dPrecedingFdFn )[ J ][ K ]
                             + d2MicroFlowdDrivingStress2[ I ][ J ] * ( *dMicroDrivingStressdFn ) [ J ][ K ];

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

                        for ( unsigned int L = 0; L < precedingDeformationGradient->size( ); L++ ){

                            d2MicroGradientFlowdDrivingStressdF[ I ][ precedingDeformationGradient->size( ) * J + L ]
                                += d2MicroGradientFlowdDrivingStress2[ I ][ microGradientDrivingStress->size( ) * J + K ] * ( *dMicroGradientDrivingStressdF )[ K ][ L ];

                            d2MicroGradientFlowdDrivingStressdChi[ I ][ precedingDeformationGradient->size( ) * J + L ]
                                += d2MicroGradientFlowdDrivingStress2[ I ][ microGradientDrivingStress->size( ) * J + K ] * ( *dMicroGradientDrivingStressdChi )[ K ][ L ];

                        }

                        for ( unsigned int L = 0; L < ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient->size( ); L++ ){

                            d2MicroGradientFlowdDrivingStressdFn[ I ][ ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient->size( ) * J + L ]
                                += d2MicroGradientFlowdDrivingStress2[ I ][ microGradientDrivingStress->size( ) * J + K ] * ( *dMicroGradientDrivingStressdFn ) [ K ][ L ];

                            d2MicroGradientFlowdDrivingStressdChin[ I ][ ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient->size( ) * J + L ]
                                += d2MicroGradientFlowdDrivingStress2[ I ][ microGradientDrivingStress->size( ) * J + K ] * ( *dMicroGradientDrivingStressdChin ) [ K ][ L ];

                        }

                    }

                    for ( unsigned int K = 0; K < precedingDeformationGradient->size( ); K++ ){

                        for ( unsigned int L = 0; L < precedingDeformationGradient->size( ); L++ ){

                            d2MicroGradientFlowdDrivingStressdF[ I ][ precedingDeformationGradient->size( ) * J + L ]
                                += d2MicroGradientFlowdDrivingStressdPrecedingF[ I ][ precedingDeformationGradient->size( ) * J + K ] * ( *dPrecedingFdF )[ K ][ L ];

                        }

                        for ( unsigned int L = 0; L < ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient->size( ); L++ ){

                            d2MicroGradientFlowdDrivingStressdFn[ I ][ ( ( *hydra->getNumConfigurations( ) ) - 1 ) * precedingDeformationGradient->size( ) * J + L ]
                                += d2MicroGradientFlowdDrivingStressdPrecedingF[ I ][ precedingDeformationGradient->size( ) * J + K ] * ( *dPrecedingFdFn )[ K ][ L ];

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

                set_previousd2MacroFlowdDrivingStressdF( d2MacroFlowdDrivingStressdF );

                set_previousd2MicroFlowdDrivingStressdF( d2MicroFlowdDrivingStressdF );

                set_previousd2MicroGradientFlowdDrivingStressdF( d2MicroGradientFlowdDrivingStressdF );

                set_previousd2MicroGradientFlowdDrivingStressdChi( d2MicroGradientFlowdDrivingStressdChi );

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

                set_d2MacroFlowdDrivingStressdF( d2MacroFlowdDrivingStressdF );

                set_d2MicroFlowdDrivingStressdF( d2MicroFlowdDrivingStressdF );

                set_d2MicroGradientFlowdDrivingStressdF( d2MicroGradientFlowdDrivingStressdF );

                set_d2MicroGradientFlowdDrivingStressdChi( d2MicroGradientFlowdDrivingStressdChi );

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

        void residual::setPlasticStrainLikeISVEvolutionRates( ){
            /*!
             * Set the evolution rates of the strain-like ISVs
             */

            setPlasticStrainLikeISVEvolutionRates( false );

        }

        void residual::setPreviousPlasticStrainLikeISVEvolutionRates( ){
            /*!
             * Set the previous evolution rates of the strain-like ISVs
             */

            setPlasticStrainLikeISVEvolutionRates( true );

        }

        void residual::setPlasticStrainLikeISVEvolutionRates( const bool isPrevious ){
            /*!
             * Set the evolution rates of the strain-like ISVs
             * 
             * \param isPrevious: A flag for if we should set the current (false) or previous (true) ISV evolution rates
             */

            const floatType   *dMacroFlowdc;

            const floatType   *dMicroFlowdc;

            const floatMatrix *dMicroGradientFlowdc;

            const floatVector *plasticMultipliers;

            if ( isPrevious ){

                plasticMultipliers    = get_previousPlasticMultipliers( );

                dMacroFlowdc          = get_previousdMacroFlowdc( );

                dMicroFlowdc          = get_previousdMicroFlowdc( );

                dMicroGradientFlowdc = get_previousdMicroGradientFlowdc( );

            }
            else{

                plasticMultipliers    = get_plasticMultipliers( );

                dMacroFlowdc          = get_dMacroFlowdc( );

                dMicroFlowdc          = get_dMicroFlowdc( );

                dMicroGradientFlowdc = get_dMicroGradientFlowdc( );

            }

            floatVector evolutionRates( plasticMultipliers->size( ), 0 );

            evolutionRates[ 0 ] = -( *dMacroFlowdc ) * ( *plasticMultipliers )[ 0 ];

            evolutionRates[ 1 ] = -( *dMicroFlowdc ) * ( *plasticMultipliers )[ 1 ];

            for ( unsigned int i = 2; i < plasticMultipliers->size( ); i++ ){

                for ( unsigned int j = 2; j < plasticMultipliers->size( ); j++ ){

                    evolutionRates[ i ] -= ( *dMicroGradientFlowdc )[ i - 2 ][ j - 2 ] * ( *plasticMultipliers )[ j ];

                }

            }

            if ( isPrevious ){

                set_previousPlasticStrainLikeISVEvolutionRates( evolutionRates );

            }
            else{

                set_plasticStrainLikeISVEvolutionRates( evolutionRates );

            }

        }

        void residual::setdPlasticStrainLikeISVEvolutionRatesdStateVariables( ){
            /*!
             * Set the Jacobian of the evolution rates of the strain-like ISVs w.r.t. the nonlinear state variables
             */

            setPlasticStrainLikeISVEvolutionRatesJacobians( false );

        }

        void residual::setPreviousdPlasticStrainLikeISVEvolutionRatesdStateVariables( ){
            /*!
             * Set the Jacobian of the previous evolution rates of the strain-like ISVs w.r.t. the nonlinear state variables
             */

            setPlasticStrainLikeISVEvolutionRatesJacobians( true );

        }

        void residual::setPlasticStrainLikeISVEvolutionRatesJacobians( const bool isPrevious ){
            /*!
             * Set the evolution rates and Jacobians of the strain-like ISVs
             * 
             * \param isPrevious: A flag for if we should set the current (false) or previous (true) ISV evolution rates
             */

            const floatType   *dMacroFlowdc;

            const floatType   *dMicroFlowdc;

            const floatMatrix *dMicroGradientFlowdc;

            const floatVector *plasticMultipliers;

            if ( isPrevious ){

                plasticMultipliers    = get_previousPlasticMultipliers( );

                dMacroFlowdc          = get_previousdMacroFlowdc( );

                dMicroFlowdc          = get_previousdMicroFlowdc( );

                dMicroGradientFlowdc = get_previousdMicroGradientFlowdc( );

            }
            else{

                plasticMultipliers    = get_plasticMultipliers( );

                dMacroFlowdc          = get_dMacroFlowdc( );

                dMicroFlowdc          = get_dMicroFlowdc( );

                dMicroGradientFlowdc = get_dMicroGradientFlowdc( );

            }

            floatVector evolutionRates( plasticMultipliers->size( ), 0 );

            floatMatrix dEvolutionRatesdStateVariables( plasticMultipliers->size( ), floatVector( get_plasticStateVariables( )->size( ), 0 ) );

            evolutionRates[ 0 ] = -( *dMacroFlowdc ) * ( *plasticMultipliers )[ 0 ];

            evolutionRates[ 1 ] = -( *dMicroFlowdc ) * ( *plasticMultipliers )[ 1 ];

            dEvolutionRatesdStateVariables[ 0 ][ 0 ] = -( *dMacroFlowdc );

            dEvolutionRatesdStateVariables[ 1 ][ 1 ] = -( *dMicroFlowdc );

            for ( unsigned int i = 2; i < plasticMultipliers->size( ); i++ ){

                for ( unsigned int j = 2; j < plasticMultipliers->size( ); j++ ){

                    evolutionRates[ i ] -= ( *dMicroGradientFlowdc )[ i - 2 ][ j - 2 ] * ( *plasticMultipliers )[ j ];

                    dEvolutionRatesdStateVariables[ i ][ j ] = -( *dMicroGradientFlowdc )[ i - 2 ][ j - 2 ];

                }

            }

            if ( isPrevious ){

                set_previousPlasticStrainLikeISVEvolutionRates( evolutionRates );

                set_previousdPlasticStrainLikeISVEvolutionRatesdStateVariables( dEvolutionRatesdStateVariables );

            }
            else{

                set_plasticStrainLikeISVEvolutionRates( evolutionRates );

                set_dPlasticStrainLikeISVEvolutionRatesdStateVariables( dEvolutionRatesdStateVariables );

            }

        }

        void residual::setUpdatedPlasticStrainLikeISVs( ){
            /*!
             * Set the updated strain like ISVs
             */

            const floatVector *previousPlasticStrainLikeISVs = get_previousPlasticStrainLikeISVs( );

            const floatVector *evolutionRates                = get_plasticStrainLikeISVEvolutionRates( );

            const floatVector *previousEvolutionRates        = get_previousPlasticStrainLikeISVEvolutionRates( );

            floatVector dISVs, updatedISVs;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolution( *hydra->getDeltaTime( ), *previousPlasticStrainLikeISVs, *previousEvolutionRates, *evolutionRates, dISVs, updatedISVs, *getIntegrationParameter( ) ) );

            set_updatedPlasticStrainLikeISVs( updatedISVs );

        }

        void residual::setdUpdatedPlasticStrainLikeISVsdStateVariables( ){
            /*!
             * Set the updated plastic strain-like ISVs jacobians
             */

            setUpdatedPlasticStrainLikeISVsJacobians( false );

        }

        void residual::setdUpdatedPlasticStrainLikeISVsdPreviousStateVariables( ){
            /*!
             * Set the previous updated plastic strain-like ISVs jacobians
             */

            setUpdatedPlasticStrainLikeISVsJacobians( true );

        }

        void residual::setUpdatedPlasticStrainLikeISVsJacobians( const bool addPrevious ){
            /*!
             * Set the updated plastic strain like ISVs Jacobians
             * 
             * \param addPrevious: Flag for whether to add the Jacobians w.r.t. the previous state variables (true) or not (false)
             */

            const floatVector *previousPlasticStrainLikeISVs = get_previousPlasticStrainLikeISVs( );

            const floatVector *evolutionRates                = get_plasticStrainLikeISVEvolutionRates( );

            const floatVector *previousEvolutionRates        = get_previousPlasticStrainLikeISVEvolutionRates( );

            floatVector dISVs, updatedISVs;

            floatMatrix dISVsdEvolutionRates;

            if ( addPrevious ){

                floatMatrix dISVsdPreviousEvolutionRates;

                TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolution( *hydra->getDeltaTime( ), *previousPlasticStrainLikeISVs, *previousEvolutionRates, *evolutionRates, dISVs, updatedISVs, dISVsdEvolutionRates, dISVsdPreviousEvolutionRates, *getIntegrationParameter( ) ) );

                floatMatrix dISVsdStateVariables = tardigradeVectorTools::dot( dISVsdPreviousEvolutionRates, *get_previousdPlasticStrainLikeISVEvolutionRatesdStateVariables( ) );

                for ( unsigned int i = 0; i < updatedISVs.size( ); i++ ){

                    dISVsdStateVariables[ i ][ i + ( *getNumPlasticMultipliers( ) ) ] += 1;

                }

                set_dUpdatedPlasticStrainLikeISVsdPreviousStateVariables( dISVsdStateVariables );

            }
            else{

                TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::midpointEvolution( *hydra->getDeltaTime( ), *previousPlasticStrainLikeISVs, *previousEvolutionRates, *evolutionRates, dISVs, updatedISVs, dISVsdEvolutionRates, *getIntegrationParameter( ) ) );

            }

            set_dUpdatedPlasticStrainLikeISVsdStateVariables( tardigradeVectorTools::dot( dISVsdEvolutionRates, *get_dPlasticStrainLikeISVEvolutionRatesdStateVariables( ) ) );

            set_updatedPlasticStrainLikeISVs( updatedISVs );

        }

        void residual::setMacroYield( ){
            /*!
             * Set the value of the macro-yield equation
             */

            setYield( false );

        }

        void residual::setMicroYield( ){
            /*!
             * Set the value of the micro-yield equation
             */

            setYield( false );

        }

        void residual::setMicroGradientYield( ){
            /*!
             * Set the value of the micro gradient-yield equation
             */

            setYield( false );

        }

        void residual::setPreviousMacroYield( ){
            /*!
             * Set the previous value of the macro-yield equation
             */

            setYield( true );

        }

        void residual::setPreviousMicroYield( ){
            /*!
             * Set the previous value of the micro-yield equation
             */

            setYield( true );

        }

        void residual::setPreviousMicroGradientYield( ){
            /*!
             * Set the previous value of the micro gradient-yield equation
             */

            setYield( true );

        }

        void residual::setYield( const bool isPrevious ){
            /*!
             * Set the values of the yield equations
             * 
             * \param isPrevious: A flag for if the current (false) or previous (true) values should be calculated
             */

            const floatVector *macroDrivingStress;

            const floatVector *microDrivingStress;

            const floatVector *microGradientDrivingStress;

            const floatType   *macroCohesion;

            const floatType   *microCohesion;

            const floatVector *microGradientCohesion;

            const floatVector *precedingDeformationGradient;

            const floatVector *macroYieldParameters = get_macroYieldParameters( );

            const floatVector *microYieldParameters = get_microYieldParameters( );

            const floatVector *microGradientYieldParameters = get_microGradientYieldParameters( );

            if ( isPrevious ){

                precedingDeformationGradient = get_previousPrecedingDeformationGradient( );

                macroDrivingStress         = get_previousMacroDrivingStress( );

                microDrivingStress         = get_previousSymmetricMicroDrivingStress( );

                microGradientDrivingStress = get_previousHigherOrderDrivingStress( );

                macroCohesion              = get_previousMacroCohesion( );

                microCohesion              = get_previousMicroCohesion( );

                microGradientCohesion      = get_previousMicroGradientCohesion( );

            }
            else{

                precedingDeformationGradient = get_precedingDeformationGradient( );

                macroDrivingStress         = get_macroDrivingStress( );

                microDrivingStress         = get_symmetricMicroDrivingStress( );

                microGradientDrivingStress = get_higherOrderDrivingStress( );

                macroCohesion              = get_macroCohesion( );

                microCohesion              = get_microCohesion( );

                microGradientCohesion      = get_microGradientCohesion( );

            }

            floatType macroYield;

            floatType microYield;

            floatVector microGradientYield;

            TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderDruckerPragerYieldEquation( *macroDrivingStress, *macroCohesion, *precedingDeformationGradient,
                                                                                        ( *macroYieldParameters )[ 0 ], ( *macroYieldParameters )[ 1 ],
                                                                                        macroYield ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderDruckerPragerYieldEquation( *microDrivingStress, *microCohesion, *precedingDeformationGradient,
                                                                                        ( *microYieldParameters )[ 0 ], ( *microYieldParameters )[ 1 ],
                                                                                        microYield ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( computeHigherOrderDruckerPragerYieldEquation( *microGradientDrivingStress, *microGradientCohesion, *precedingDeformationGradient,
                                                                                        ( *microGradientYieldParameters )[ 0 ], ( *microGradientYieldParameters )[ 1 ],
                                                                                        microGradientYield ) );

            if ( isPrevious ){

                set_previousMacroYield( macroYield );

                set_previousMicroYield( microYield );

                set_previousMicroGradientYield( microGradientYield );

            }
            else{

                set_macroYield( macroYield );

                set_microYield( microYield );

                set_microGradientYield( microGradientYield );

            }

        }

        void residual::setdMacroYielddStress( ){
            /*!
             * Set the Jacobian of the macro yield stress w.r.t. the macro stress measure
             */

            setYieldJacobians( false );

        }

        void residual::setdMacroYielddStateVariables( ){
            /*!
             * Set the Jacobian of the macro yield stress w.r.t. the state variables
             */

            setYieldJacobians( false );

        }

        void residual::setdMacroYielddF( ){
            /*!
             * Set the Jacobian of the macro yield stress w.r.t. the deformation gradient
             */

            setYieldJacobians( false );

        }

        void residual::setdMacroYielddFn( ){
            /*!
             * Set the Jacobian of the macro yield stress w.r.t. the sub deformation gradients
             */

            setYieldJacobians( false );

        }

        void residual::setdMicroYielddStress( ){
            /*!
             * Set the Jacobian of the micro yield stress w.r.t. the micro stress measure
             */

            setYieldJacobians( false );

        }

        void residual::setdMicroYielddStateVariables( ){
            /*!
             * Set the Jacobian of the micro yield stress w.r.t. the state variables
             */

            setYieldJacobians( false );

        }

        void residual::setdMicroYielddF( ){
            /*!
             * Set the Jacobian of the micro yield stress w.r.t. the deformation gradient
             */

            setYieldJacobians( false );

        }

        void residual::setdMicroYielddFn( ){
            /*!
             * Set the Jacobian of the micro yield stress w.r.t. the sub deformation gradients
             */

            setYieldJacobians( false );

        }

        void residual::setdMicroGradientYielddStress( ){
            /*!
             * Set the Jacobian of the micro gradient yield stress w.r.t. the micro gradient stress measure
             */

            setYieldJacobians( false );

        }

        void residual::setdMicroGradientYielddStateVariables( ){
            /*!
             * Set the Jacobian of the micro gradient yield stress w.r.t. the state variables
             */

            setYieldJacobians( false );

        }

        void residual::setdMicroGradientYielddF( ){
            /*!
             * Set the Jacobian of the micro gradient yield stress w.r.t. the deformation gradient
             */

            setYieldJacobians( false );

        }

        void residual::setdMicroGradientYielddFn( ){
            /*!
             * Set the Jacobian of the micro gradient yield stress w.r.t. the sub deformation gradient
             */

            setYieldJacobians( false );

        }

        void residual::setdMicroGradientYielddChi( ){
            /*!
             * Set the Jacobian of the micro gradient yield stress w.r.t. the micro deformation
             */

            setYieldJacobians( false );

        }

        void residual::setdMicroGradientYielddChin( ){
            /*!
             * Set the Jacobian of the micro gradient yield stress w.r.t. the sub micro deformations
             */

            setYieldJacobians( false );

        }

        void residual::setPreviousdMacroYielddStress( ){
            /*!
             * Set the previous Jacobian of the macro yield stress w.r.t. the macro stress measure
             */

            setYieldJacobians( true );

        }

        void residual::setPreviousdMacroYielddStateVariables( ){
            /*!
             * Set the previous Jacobian of the macro yield stress w.r.t. the state variables
             */

            setYieldJacobians( true );

        }

        void residual::setPreviousdMacroYielddF( ){
            /*!
             * Set the previous Jacobian of the macro yield stress w.r.t. the deformation gradient
             */

            setYieldJacobians( true );

        }

        void residual::setPreviousdMacroYielddFn( ){
            /*!
             * Set the previous Jacobian of the macro yield stress w.r.t. the sub deformation gradients
             */

            setYieldJacobians( true );

        }

        void residual::setPreviousdMicroYielddStress( ){
            /*!
             * Set the previous Jacobian of the micro yield stress w.r.t. the micro stress measure
             */

            setYieldJacobians( true );

        }

        void residual::setPreviousdMicroYielddStateVariables( ){
            /*!
             * Set the previous Jacobian of the micro yield stress w.r.t. the state variables
             */

            setYieldJacobians( true );

        }

        void residual::setPreviousdMicroYielddF( ){
            /*!
             * Set the previous Jacobian of the micro yield stress w.r.t. the deformation gradient
             */

            setYieldJacobians( true );

        }

        void residual::setPreviousdMicroYielddFn( ){
            /*!
             * Set the previous Jacobian of the micro yield stress w.r.t. the sub deformation gradients
             */

            setYieldJacobians( true );

        }

        void residual::setPreviousdMicroGradientYielddStress( ){
            /*!
             * Set the previous Jacobian of the micro gradient yield stress w.r.t. the micro gradient stress measure
             */

            setYieldJacobians( true );

        }

        void residual::setPreviousdMicroGradientYielddStateVariables( ){
            /*!
             * Set the previous Jacobian of the micro gradient yield stress w.r.t. the state variables
             */

            setYieldJacobians( true );

        }

        void residual::setPreviousdMicroGradientYielddF( ){
            /*!
             * Set the previous Jacobian of the micro gradient yield stress w.r.t. the deformation gradient
             */

            setYieldJacobians( true );

        }

        void residual::setPreviousdMicroGradientYielddFn( ){
            /*!
             * Set the previous Jacobian of the micro gradient yield stress w.r.t. the sub deformation gradient
             */

            setYieldJacobians( true );

        }

        void residual::setPreviousdMicroGradientYielddChi( ){
            /*!
             * Set the previous Jacobian of the micro gradient yield stress w.r.t. the micro deformation
             */

            setYieldJacobians( true );

        }

        void residual::setPreviousdMicroGradientYielddChin( ){
            /*!
             * Set the previous Jacobian of the micro gradient yield stress w.r.t. the sub micro deformations
             */

            setYieldJacobians( true );

        }

        void residual::setYieldJacobians( const bool isPrevious ){
            /*!
             * Set the value and jacobians of the yield functions
             * 
             * \param isPrevious: Flag for whether to get the jacobians of the current (false) or previous (true) yield functions
             */

            const floatVector *macroDrivingStress;

            const floatVector *microDrivingStress;

            const floatVector *microGradientDrivingStress;

            const floatType   *macroCohesion;

            const floatType   *microCohesion;

            const floatVector *microGradientCohesion;

            const floatVector *dMacroCohesiondStateVariables;

            const floatMatrix *dMacroDrivingStressdStress;

            const floatMatrix *dMacroDrivingStressdF;

            const floatMatrix *dMacroDrivingStressdFn;

            const floatVector *dMicroCohesiondStateVariables;

            const floatMatrix *dMicroDrivingStressdStress;

            const floatMatrix *dMicroDrivingStressdF;

            const floatMatrix *dMicroDrivingStressdFn;

            const floatMatrix *dMicroGradientCohesiondStateVariables;

            const floatMatrix *dMicroGradientDrivingStressdStress;

            const floatMatrix *dMicroGradientDrivingStressdF;

            const floatMatrix *dMicroGradientDrivingStressdFn;

            const floatMatrix *dMicroGradientDrivingStressdChi;

            const floatMatrix *dMicroGradientDrivingStressdChin;

            const floatVector *precedingDeformationGradient;

            const floatMatrix *dPrecedingFdF;

            const floatMatrix *dPrecedingFdFn;

            const floatVector *macroYieldParameters = get_macroYieldParameters( );

            const floatVector *microYieldParameters = get_microYieldParameters( );

            const floatVector *microGradientYieldParameters = get_microGradientYieldParameters( );

            if ( isPrevious ){

                precedingDeformationGradient = get_previousPrecedingDeformationGradient( );

                dPrecedingFdF  = get_previousdPrecedingDeformationGradientdF( );

                dPrecedingFdFn = get_previousdPrecedingDeformationGradientdFn( );

                dMacroCohesiondStateVariables         = get_previousdMacroCohesiondStateVariables( );

                dMicroCohesiondStateVariables         = get_previousdMicroCohesiondStateVariables( );

                dMicroGradientCohesiondStateVariables = get_previousdMicroGradientCohesiondStateVariables( );

                dMacroDrivingStressdStress            = get_previousdMacroDrivingStressdMacroStress( );

                dMicroDrivingStressdStress            = get_previousdSymmetricMicroDrivingStressdMicroStress( );

                dMicroGradientDrivingStressdStress    = get_previousdHigherOrderDrivingStressdHigherOrderStress( );

                dMacroDrivingStressdF                 = get_previousdMacroDrivingStressdF( );

                dMicroDrivingStressdF                 = get_previousdSymmetricMicroDrivingStressdF( );

                dMicroGradientDrivingStressdF         = get_previousdHigherOrderDrivingStressdF( );

                dMicroGradientDrivingStressdChi       = get_previousdHigherOrderDrivingStressdChi( );

                dMacroDrivingStressdFn                = get_previousdMacroDrivingStressdFn( );

                dMicroDrivingStressdFn                = get_previousdSymmetricMicroDrivingStressdFn( );

                dMicroGradientDrivingStressdFn        = get_previousdHigherOrderDrivingStressdFn( );

                dMicroGradientDrivingStressdChin      = get_previousdHigherOrderDrivingStressdChin( );

                macroDrivingStress         = get_previousMacroDrivingStress( );

                microDrivingStress         = get_previousSymmetricMicroDrivingStress( );

                microGradientDrivingStress = get_previousHigherOrderDrivingStress( );

                macroCohesion              = get_previousMacroCohesion( );

                microCohesion              = get_previousMicroCohesion( );

                microGradientCohesion      = get_previousMicroGradientCohesion( );

            }
            else{

                precedingDeformationGradient = get_precedingDeformationGradient( );

                dPrecedingFdF  = get_dPrecedingDeformationGradientdF( );

                dPrecedingFdFn = get_dPrecedingDeformationGradientdFn( );

                dMacroCohesiondStateVariables         = get_dMacroCohesiondStateVariables( );

                dMicroCohesiondStateVariables         = get_dMicroCohesiondStateVariables( );

                dMicroGradientCohesiondStateVariables = get_dMicroGradientCohesiondStateVariables( );

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

                macroDrivingStress         = get_macroDrivingStress( );

                microDrivingStress         = get_symmetricMicroDrivingStress( );

                microGradientDrivingStress = get_higherOrderDrivingStress( );

                macroCohesion              = get_macroCohesion( );

                microCohesion              = get_microCohesion( );

                microGradientCohesion      = get_microGradientCohesion( );

            }

            floatType   macroYield;

            floatType   microYield;

            floatVector microGradientYield;

            floatVector dMacroYielddDrivingStress;

            floatType   dMacroYielddCohesion;

            floatVector dMacroYielddPrecedingF;

            floatVector dMicroYielddDrivingStress;

            floatType   dMicroYielddCohesion;

            floatVector dMicroYielddPrecedingF;

            floatMatrix dMicroGradientYielddDrivingStress;

            floatMatrix dMicroGradientYielddCohesion;

            floatMatrix dMicroGradientYielddPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderDruckerPragerYieldEquation( *macroDrivingStress, *macroCohesion, *precedingDeformationGradient,
                                                                                        ( *macroYieldParameters )[ 0 ], ( *macroYieldParameters )[ 1 ],
                                                                                        macroYield, dMacroYielddDrivingStress, dMacroYielddCohesion, dMacroYielddPrecedingF ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( computeSecondOrderDruckerPragerYieldEquation( *microDrivingStress, *microCohesion, *precedingDeformationGradient,
                                                                                        ( *microYieldParameters )[ 0 ], ( *microYieldParameters )[ 1 ],
                                                                                        microYield, dMicroYielddDrivingStress, dMicroYielddCohesion, dMicroYielddPrecedingF ) );

            TARDIGRADE_ERROR_TOOLS_CATCH( computeHigherOrderDruckerPragerYieldEquation( *microGradientDrivingStress, *microGradientCohesion, *precedingDeformationGradient,
                                                                                        ( *microGradientYieldParameters )[ 0 ], ( *microGradientYieldParameters )[ 1 ],
                                                                                        microGradientYield, dMicroGradientYielddDrivingStress, dMicroGradientYielddCohesion,
                                                                                        dMicroGradientYielddPrecedingF ) );

            floatVector dMacroYielddStress = tardigradeVectorTools::Tdot( *dMacroDrivingStressdStress, dMacroYielddDrivingStress );

            floatVector dMacroYielddStateVariables = dMacroYielddCohesion * ( *dMacroCohesiondStateVariables );

            floatVector dMacroYielddF = tardigradeVectorTools::Tdot( *dMacroDrivingStressdF, dMacroYielddDrivingStress )
                                      + tardigradeVectorTools::Tdot( *dPrecedingFdF, dMacroYielddPrecedingF );

            floatVector dMacroYielddFn = tardigradeVectorTools::Tdot( *dMacroDrivingStressdFn, dMacroYielddDrivingStress )
                                       + tardigradeVectorTools::Tdot( *dPrecedingFdFn, dMacroYielddPrecedingF );

            floatVector dMicroYielddStress = tardigradeVectorTools::Tdot( *dMicroDrivingStressdStress, dMicroYielddDrivingStress );

            floatVector dMicroYielddStateVariables = dMicroYielddCohesion * ( *dMicroCohesiondStateVariables );

            floatVector dMicroYielddF = tardigradeVectorTools::Tdot( *dMicroDrivingStressdF, dMicroYielddDrivingStress )
                                      + tardigradeVectorTools::Tdot( *dPrecedingFdF, dMicroYielddPrecedingF );

            floatVector dMicroYielddFn = tardigradeVectorTools::Tdot( *dMicroDrivingStressdFn, dMicroYielddDrivingStress )
                                       + tardigradeVectorTools::Tdot( *dPrecedingFdFn, dMicroYielddPrecedingF );

            floatMatrix dMicroGradientYielddStress = tardigradeVectorTools::dot( dMicroGradientYielddDrivingStress, *dMicroGradientDrivingStressdStress );

            floatMatrix dMicroGradientYielddStateVariables = tardigradeVectorTools::dot( dMicroGradientYielddCohesion, *dMicroGradientCohesiondStateVariables );

            floatMatrix dMicroGradientYielddF = tardigradeVectorTools::dot( dMicroGradientYielddDrivingStress, *dMicroGradientDrivingStressdF )
                                              + tardigradeVectorTools::dot( dMicroGradientYielddPrecedingF, *dPrecedingFdF );

            floatMatrix dMicroGradientYielddFn = tardigradeVectorTools::dot( dMicroGradientYielddDrivingStress, *dMicroGradientDrivingStressdFn )
                                               + tardigradeVectorTools::dot( dMicroGradientYielddPrecedingF, *dPrecedingFdFn );

            floatMatrix dMicroGradientYielddChi = tardigradeVectorTools::dot( dMicroGradientYielddDrivingStress, *dMicroGradientDrivingStressdChi );

            floatMatrix dMicroGradientYielddChin = tardigradeVectorTools::dot( dMicroGradientYielddDrivingStress, *dMicroGradientDrivingStressdChin );

            if ( isPrevious ){

                set_previousMacroYield( macroYield );

                set_previousMicroYield( microYield );

                set_previousMicroGradientYield( microGradientYield );

                set_previousdMacroYielddStress( dMacroYielddStress );

                set_previousdMacroYielddStateVariables( dMacroYielddStateVariables );

                set_previousdMacroYielddF( dMacroYielddF );

                set_previousdMacroYielddFn( dMacroYielddFn );

                set_previousdMicroYielddStress( dMicroYielddStress );

                set_previousdMicroYielddStateVariables( dMicroYielddStateVariables );

                set_previousdMicroYielddF( dMicroYielddF );

                set_previousdMicroYielddFn( dMicroYielddFn );

                set_previousdMicroGradientYielddStress( dMicroGradientYielddStress );

                set_previousdMicroGradientYielddStateVariables( dMicroGradientYielddStateVariables );

                set_previousdMicroGradientYielddF( dMicroGradientYielddF );

                set_previousdMicroGradientYielddFn( dMicroGradientYielddFn );

                set_previousdMicroGradientYielddChi( dMicroGradientYielddChi );

                set_previousdMicroGradientYielddChin( dMicroGradientYielddChin );

            }
            else{

                set_macroYield( macroYield );

                set_microYield( microYield );

                set_microGradientYield( microGradientYield );

                set_dMacroYielddStress( dMacroYielddStress );

                set_dMacroYielddStateVariables( dMacroYielddStateVariables );

                set_dMacroYielddF( dMacroYielddF );

                set_dMacroYielddFn( dMacroYielddFn );

                set_dMicroYielddStress( dMicroYielddStress );

                set_dMicroYielddStateVariables( dMicroYielddStateVariables );

                set_dMicroYielddF( dMicroYielddF );

                set_dMicroYielddFn( dMicroYielddFn );

                set_dMicroGradientYielddStress( dMicroGradientYielddStress );

                set_dMicroGradientYielddStateVariables( dMicroGradientYielddStateVariables );

                set_dMicroGradientYielddF( dMicroGradientYielddF );

                set_dMicroGradientYielddFn( dMicroGradientYielddFn );

                set_dMicroGradientYielddChi( dMicroGradientYielddChi );

                set_dMicroGradientYielddChin( dMicroGradientYielddChin );

            }

        }

        void residual::setPrecedingDeformationGradient( ){
            /*!
             * Set the preceding deformation gradient
             */

            setPrecedingDeformationGradient( false );

        }

        void residual::setPreviousPrecedingDeformationGradient( ){
            /*!
             * Set the previous preceding deformation gradient
             */

            setPrecedingDeformationGradient( true );

        }

        void residual::setPrecedingDeformationGradient( const bool isPrevious ){
            /*!
             * Set the preceding deformation gradient to the plastic configuration
             * 
             * \param isPrevious: Whether to set the current (false) or previous (true) preceding deformation gradient
             */

            if ( isPrevious ){

                set_previousPrecedingDeformationGradient( hydra->getPreviousPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

            }
            else{

                set_precedingDeformationGradient( hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

            }

        }

        void residual::setdPrecedingDeformationGradientdF( ){
            /*!
             * Set the jacobian of the preceding deformation gradient w.r.t. the deformation gradient
             */

            setPrecedingDeformationGradientJacobians( false );

        }

        void residual::setdPrecedingDeformationGradientdFn( ){
            /*!
             * Set the jacobian of the preceding deformation gradient w.r.t. the sub-deformation gradients
             */

            setPrecedingDeformationGradientJacobians( false );

        }

        void residual::setPreviousdPrecedingDeformationGradientdF( ){
            /*!
             * Set the previous jacobian of the preceding deformation gradient w.r.t. the deformation gradient
             */

            setPrecedingDeformationGradientJacobians( true );

        }

        void residual::setPreviousdPrecedingDeformationGradientdFn( ){
            /*!
             * Set the previous jacobian of the preceding deformation gradient w.r.t. the sub-deformation gradients
             */

            setPrecedingDeformationGradientJacobians( true );

        }

        void residual::setPrecedingDeformationGradientJacobians( const bool isPrevious ){
            /*!
             * Set the preceding deformation gradient to the plastic configuration and its jacobians
             * 
             * \param isPrevious: Whether to set the current (false) or previous (true) preceding deformation gradient
             */

            floatMatrix dPrecedingFdSubFs;

            const floatMatrix *dF1dF;

            const floatMatrix *dF1dFn;

            if ( isPrevious ){

                set_previousPrecedingDeformationGradient( hydra->getPreviousPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                dPrecedingFdSubFs = hydra->getPreviousPrecedingConfigurationJacobian( *getPlasticConfigurationIndex( ) );

                dF1dF = hydra->get_previousdF1dF( );

                dF1dFn = hydra->get_previousdF1dFn( );

            }
            else{

                set_precedingDeformationGradient( hydra->getPrecedingConfiguration( *getPlasticConfigurationIndex( ) ) );

                dPrecedingFdSubFs = hydra->getPrecedingConfigurationJacobian( *getPlasticConfigurationIndex( ) );

                dF1dF = hydra->get_dF1dF( );

                dF1dFn = hydra->get_dF1dFn( );

            }

            // Construct the derivatives of the preceding F

            floatMatrix dPrecedingFdF( get_precedingDeformationGradient( )->size( ), floatVector( hydra->getDeformationGradient( )->size( ), 0 ) );

            floatMatrix dPrecedingFdFn( get_precedingDeformationGradient( )->size( ), floatVector( ( ( *hydra->getNumConfigurations( ) ) - 1 ) * hydra->getDeformationGradient( )->size( ), 0 ) );

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

            if ( isPrevious ){

                set_previousdPrecedingDeformationGradientdF( dPrecedingFdF );

                set_previousdPrecedingDeformationGradientdFn( dPrecedingFdFn );

            }
            else{

                set_dPrecedingDeformationGradientdF( dPrecedingFdF );

                set_dPrecedingDeformationGradientdFn( dPrecedingFdFn );

            }

        }

        void residual::setPrecedingMicroDeformation( ){
            /*!
             * Set the preceding micro deformation
             */

            setPrecedingMicroDeformation( false );

        }

        void residual::setPreviousPrecedingMicroDeformation( ){
            /*!
             * Set the previous preceding micro deformation
             */

            setPrecedingMicroDeformation( true );

        }

        void residual::setPrecedingMicroDeformation( const bool isPrevious ){
            /*!
             * Set the preceding micro deformation to the plastic configuration
             * 
             * \param isPrevious: Whether to set the current (false) or previous (true) preceding micro deformation
             */

            if ( isPrevious ){

                set_previousPrecedingMicroDeformation( hydra->getPreviousPrecedingMicroConfiguration( *getPlasticConfigurationIndex( ) ) );

            }
            else{

                set_precedingMicroDeformation( hydra->getPrecedingMicroConfiguration( *getPlasticConfigurationIndex( ) ) );

            }

        }

        void residual::setdPrecedingMicroDeformationdChi( ){
            /*!
             * Set the jacobian of the preceding micro deformation w.r.t. the micro deformation
             */

            setPrecedingMicroDeformationJacobians( false );

        }

        void residual::setdPrecedingMicroDeformationdChin( ){
            /*!
             * Set the jacobian of the preceding micro deformation w.r.t. the sub micro deformation
             */

            setPrecedingMicroDeformationJacobians( false );

        }

        void residual::setPreviousdPrecedingMicroDeformationdChi( ){
            /*!
             * Set the previous jacobian of the preceding micro deformation w.r.t. the micro deformation
             */

            setPrecedingMicroDeformationJacobians( true );

        }

        void residual::setPreviousdPrecedingMicroDeformationdChin( ){
            /*!
             * Set the previous jacobian of the preceding micro deformation w.r.t. the sub-micro deformation
             */

            setPrecedingMicroDeformationJacobians( true );

        }

        void residual::setPrecedingMicroDeformationJacobians( const bool isPrevious ){
            /*!
             * Set the preceding micro deformation to the plastic configuration and its jacobians
             * 
             * \param isPrevious: Whether to set the current (false) or previous (true) preceding micro deformation
             */

            floatMatrix dPrecedingChidSubChis;

            const floatMatrix *dChi1dChi;

            const floatMatrix *dChi1dChin;

            if ( isPrevious ){

                set_previousPrecedingMicroDeformation( hydra->getPreviousPrecedingMicroConfiguration( *getPlasticConfigurationIndex( ) ) );

                dPrecedingChidSubChis = hydra->getPreviousPrecedingMicroConfigurationJacobian( *getPlasticConfigurationIndex( ) );

                dChi1dChi  = hydra->get_previousdChi1dChi( );

                dChi1dChin = hydra->get_previousdChi1dChin( );

            }
            else{

                set_precedingMicroDeformation( hydra->getPrecedingMicroConfiguration( *getPlasticConfigurationIndex( ) ) );

                dPrecedingChidSubChis = hydra->getPrecedingMicroConfigurationJacobian( *getPlasticConfigurationIndex( ) );

                dChi1dChi  = hydra->get_dChi1dChi( );

                dChi1dChin = hydra->get_dChi1dChin( );

            }

            // Construct the derivatives of the preceding F

            floatMatrix dPrecedingChidChi( get_precedingMicroDeformation( )->size( ), floatVector( hydra->getMicroDeformation( )->size( ), 0 ) );

            floatMatrix dPrecedingChidChin( get_precedingMicroDeformation( )->size( ), floatVector( ( ( *hydra->getNumConfigurations( ) ) - 1 ) * hydra->getMicroDeformation( )->size( ), 0 ) );

            for ( unsigned int i = 0; i < hydra->getMicroDeformation( )->size( ); i++ ){

                for ( unsigned int j = 0; j < hydra->getMicroDeformation( )->size( ); j++ ){

                    for ( unsigned int k = 0; k < hydra->getMicroDeformation( )->size( ); k++ ){

                        dPrecedingChidChi[ i ][ j ] += dPrecedingChidSubChis[ i ][ k ] * ( *dChi1dChi )[ k ][ j ];

                    }

                }

                for ( unsigned int j = 0; j < ( ( *hydra->getNumConfigurations( ) ) - 1 ) * hydra->getMicroDeformation( )->size( ); j++ ){

                    dPrecedingChidChin[ i ][ j ] = dPrecedingChidSubChis[ i ][ hydra->getMicroDeformation( )->size( ) + j ];

                    for ( unsigned int k = 0; k < hydra->getMicroDeformation( )->size( ); k++ ){

                        dPrecedingChidChin[ i ][ j ] += dPrecedingChidSubChis[ i ][ k ] * ( *dChi1dChin )[ k ][ j ];

                    }

                }

            }

            if ( isPrevious ){

                set_previousdPrecedingMicroDeformationdChi( dPrecedingChidChi );

                set_previousdPrecedingMicroDeformationdChin( dPrecedingChidChin );

            }
            else{

                set_dPrecedingMicroDeformationdChi( dPrecedingChidChi );

                set_dPrecedingMicroDeformationdChin( dPrecedingChidChin );

            }

        }

        void residual::setPlasticMacroVelocityGradient( ){
            /*!
             * Set the plastic macro velocity gradient
             */

            setPlasticVelocityGradients( false );

        }

        void residual::setPreviousPlasticMacroVelocityGradient( ){
            /*!
             * Set the previous plastic macro velocity gradient
             */

            setPlasticVelocityGradients( true );

        }

        void residual::setPlasticMicroVelocityGradient( ){
            /*!
             * Set the plastic micro velocity gradient
             */

            setPlasticVelocityGradients( false );

        }

        void residual::setPreviousPlasticMicroVelocityGradient( ){
            /*!
             * Set the previous plastic macro velocity gradient
             */

            setPlasticVelocityGradients( true );

        }

        void residual::setPlasticVelocityGradients( const bool isPrevious ){
            /*!
             * Set the plastic macro velocity gradient
             * 
             * \param isPrevious: Flag for if the previous (true) or current (false) velocity gradient should be calculated
             */

            const unsigned int *dim = hydra->getDimension( );

            const floatVector *precedingDeformationGradient;

            const floatVector *precedingMicroDeformation;

            const floatVector *precedingGradientMicroDeformation;

            const floatVector *plasticMultipliers;

            const floatVector *dMacroFlowdDrivingStress;

            const floatVector *dMicroFlowdDrivingStress;

            if ( isPrevious ){

                precedingDeformationGradient = get_previousPrecedingDeformationGradient( );

                precedingMicroDeformation = get_previousPrecedingMicroDeformation( );

//                precedingGradientMicroDeformation = get_previousPrecedingGradientMicroDeformation( );

                plasticMultipliers = get_previousPlasticMultipliers( );

                dMacroFlowdDrivingStress = get_previousdMacroFlowdDrivingStress( );

                dMicroFlowdDrivingStress = get_previousdMicroFlowdDrivingStress( );

            }
            else{

                precedingDeformationGradient = get_precedingDeformationGradient( );

                precedingMicroDeformation = get_precedingMicroDeformation( );

//                precedingGradientMicroDeformation = get_precedingGradientMicroDeformation( );

                plasticMultipliers = get_plasticMultipliers( );

                dMacroFlowdDrivingStress = get_dMacroFlowdDrivingStress( );

                dMicroFlowdDrivingStress = get_dMicroFlowdDrivingStress( );

            }

            // Form the preceding RCG and its inverse
            floatVector precedingRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( *precedingDeformationGradient, precedingRCG ) );

            floatVector inversePrecedingRCG = tardigradeVectorTools::inverse( precedingRCG, *dim, *dim );

            // Form the precedingPsi and its inverse
            floatVector precedingPsi;

            TARDIGRADE_ERROR_TOOLS_CATCH( precedingPsi = tardigradeVectorTools::matrixMultiply( *precedingDeformationGradient, *precedingMicroDeformation, *dim, *dim, *dim, *dim, true, false ) );

            floatVector inversePrecedingPsi = tardigradeVectorTools::inverse( precedingPsi, *dim, *dim );

            // Form the preceding micro RCG and its inverse
            floatVector precedingMicroRCG;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( *precedingMicroDeformation, precedingMicroRCG ) );

//            // Form Gamma
//            floatVector precedingGamma( ( *dim ) * ( *dim ) * ( *dim ), 0 );
//            for ( unsigned int I = 0; I < *dim; I++ ){
//
//                for ( unsigned int J = 0; J < *dim; J++ ){
//
//                    for ( unsigned int K = 0; K < *dim; K++ ){
//
//                        for ( unsigned int i = 0; i < *dim; i++ ){
//
//                            precedingGamma[ ( *dim ) * ( *dim ) * I + ( *dim ) * J + K ]
//                                += ( *precedingDeformationGradient )[ ( *dim ) * i + I ] * ( *precedingGradientMicroDeformation )[ ( *dim ) * ( *dim ) * i + ( *dim ) * J + K ];
//
//                        }
//
//                    }
//
//                }
//
//            }

            floatVector macroVelocityGradient;

            floatVector microVelocityGradient;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMacroVelocityGradient( ( *plasticMultipliers )[ 0 ], ( *plasticMultipliers )[ 1 ],
                                                     inversePrecedingRCG, *dMacroFlowdDrivingStress, *dMicroFlowdDrivingStress,
                                                     macroVelocityGradient )
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMicroVelocityGradient( ( *plasticMultipliers )[ 1 ], precedingMicroRCG, precedingPsi, inversePrecedingPsi,
                                                     *dMicroFlowdDrivingStress, microVelocityGradient );
            )

            if ( isPrevious ){

                set_previousPlasticMacroVelocityGradient( macroVelocityGradient );

                set_previousPlasticMicroVelocityGradient( microVelocityGradient );

            }
            else{

                set_plasticMacroVelocityGradient( macroVelocityGradient );

                set_plasticMicroVelocityGradient( microVelocityGradient );

            }

        }

        void residual::setdPlasticMacroVelocityGradientdMacroStress( ){
            /*!
             * Set the Jacobian of the plastic macro velocity gradient w.r.t. the macro stress
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticMacroVelocityGradientdMacroStress( ){
            /*!
             * Set the Jacobian of the previous plastic macro velocity gradient w.r.t. the macro stress
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticMacroVelocityGradientdMicroStress( ){
            /*!
             * Set the Jacobian of the plastic macro velocity gradient w.r.t. the micro stress
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticMacroVelocityGradientdMicroStress( ){
            /*!
             * Set the Jacobian of the previous plastic macro velocity gradient w.r.t. the micro stress
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticMacroVelocityGradientdF( ){
            /*!
             * Set the Jacobian of the plastic macro velocity gradient w.r.t. the deformation gradient
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticMacroVelocityGradientdF( ){
            /*!
             * Set the Jacobian of the previous plastic macro velocity gradient w.r.t. the deformation gradient
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticMacroVelocityGradientdFn( ){
            /*!
             * Set the Jacobian of the plastic macro velocity gradient w.r.t. the sub deformation gradients
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticMacroVelocityGradientdFn( ){
            /*!
             * Set the Jacobian of the previous plastic macro velocity gradient w.r.t. the sub deformation gradients
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticMacroVelocityGradientdStateVariables( ){
            /*!
             * Set the Jacobian of the plastic macro velocity gradient w.r.t. the state variables
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticMacroVelocityGradientdStateVariables( ){
            /*!
             * Set the Jacobian of the previous plastic macro velocity gradient w.r.t. the state variables
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticMicroVelocityGradientdMicroStress( ){
            /*!
             * Set the Jacobian of the plastic micro velocity gradient w.r.t. the micro stress
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticMicroVelocityGradientdMicroStress( ){
            /*!
             * Set the Jacobian of the previous plastic micro velocity gradient w.r.t. the micro stress
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticMicroVelocityGradientdFn( ){
            /*!
             * Set the Jacobian of the plastic micro velocity gradient w.r.t. the sub deformation gradients
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticMicroVelocityGradientdFn( ){
            /*!
             * Set the Jacobian of the previous plastic micro velocity gradient w.r.t. the sub deformation gradients
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setdPlasticMicroVelocityGradientdChin( ){
            /*!
             * Set the Jacobian of the plastic micro velocity gradient w.r.t. the sub micro deformations
             */

            setPlasticVelocityGradientsJacobians( false );

        }

        void residual::setPreviousdPlasticMicroVelocityGradientdChin( ){
            /*!
             * Set the Jacobian of the previous plastic micro velocity gradient w.r.t. the sub micro deformations
             */

            setPlasticVelocityGradientsJacobians( true );

        }

        void residual::setPlasticVelocityGradientsJacobians( const bool isPrevious ){
            /*!
             * Set the plastic macro velocity gradient and their Jacobians
             * 
             * \param isPrevious: Flag for if the previous (true) or current (false) velocity gradient should be calculated
             */

            const unsigned int *dim = hydra->getDimension( );

            const floatVector *precedingDeformationGradient;

            const floatMatrix *dPrecedingFdF;

            const floatMatrix *dPrecedingFdFn;

            const floatMatrix *dPrecedingChidChi;

            const floatMatrix *dPrecedingChidChin;

            const floatVector *precedingMicroDeformation;

            const floatVector *precedingGradientMicroDeformation;

            const floatVector *plasticMultipliers;

            const floatVector *dMacroFlowdDrivingStress;

            const floatVector *dMicroFlowdDrivingStress;

            const floatMatrix *d2MacroFlowdDrivingStressdStress;

            const floatMatrix *d2MacroFlowdDrivingStressdF;

            const floatMatrix *d2MacroFlowdDrivingStressdFn;

            const floatMatrix *d2MicroFlowdDrivingStressdF;

            const floatMatrix *d2MicroFlowdDrivingStressdFn;

            const floatMatrix *d2MicroFlowdDrivingStressdStress;

            if ( isPrevious ){

                dPrecedingFdF = get_previousdPrecedingDeformationGradientdF( );

                dPrecedingFdFn = get_previousdPrecedingDeformationGradientdFn( );

                dPrecedingChidChi = get_previousdPrecedingMicroDeformationdChi( );

                dPrecedingChidChin = get_previousdPrecedingMicroDeformationdChin( );

                precedingDeformationGradient = get_previousPrecedingDeformationGradient( );

                precedingMicroDeformation = get_previousPrecedingMicroDeformation( );

//                precedingGradientMicroDeformation = get_previousPrecedingGradientMicroDeformation( );

                plasticMultipliers = get_previousPlasticMultipliers( );

                d2MacroFlowdDrivingStressdStress = get_previousd2MacroFlowdDrivingStressdStress( );

                d2MacroFlowdDrivingStressdF      = get_previousd2MacroFlowdDrivingStressdF( );

                d2MacroFlowdDrivingStressdFn     = get_previousd2MacroFlowdDrivingStressdFn( );

                d2MicroFlowdDrivingStressdStress = get_previousd2MicroFlowdDrivingStressdStress( );

                d2MicroFlowdDrivingStressdF      = get_previousd2MicroFlowdDrivingStressdF( );

                d2MicroFlowdDrivingStressdFn     = get_previousd2MicroFlowdDrivingStressdFn( );

                dMacroFlowdDrivingStress         = get_previousdMacroFlowdDrivingStress( );

                dMicroFlowdDrivingStress         = get_previousdMicroFlowdDrivingStress( );

            }
            else{

                dPrecedingFdF = get_dPrecedingDeformationGradientdF( );

                dPrecedingFdFn = get_dPrecedingDeformationGradientdFn( );

                dPrecedingChidChi = get_dPrecedingMicroDeformationdChi( );

                dPrecedingChidChin = get_dPrecedingMicroDeformationdChin( );

                precedingDeformationGradient = get_precedingDeformationGradient( );

                precedingMicroDeformation = get_precedingMicroDeformation( );

//                precedingGradientMicroDeformation = get_precedingGradientMicroDeformation( );

                plasticMultipliers = get_plasticMultipliers( );

                d2MacroFlowdDrivingStressdStress = get_d2MacroFlowdDrivingStressdStress( );

                d2MacroFlowdDrivingStressdF      = get_d2MacroFlowdDrivingStressdF( );

                d2MacroFlowdDrivingStressdFn     = get_d2MacroFlowdDrivingStressdFn( );

                d2MicroFlowdDrivingStressdStress = get_d2MicroFlowdDrivingStressdStress( );

                d2MicroFlowdDrivingStressdF      = get_d2MicroFlowdDrivingStressdF( );

                d2MicroFlowdDrivingStressdFn     = get_d2MicroFlowdDrivingStressdFn( );

                dMacroFlowdDrivingStress         = get_dMacroFlowdDrivingStress( );

                dMicroFlowdDrivingStress         = get_dMicroFlowdDrivingStress( );

            }

            // Form the preceding RCG and its inverse
            floatVector precedingRCG;

            floatMatrix dRCGdPrecedingF;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( *precedingDeformationGradient, precedingRCG, dRCGdPrecedingF ) );

            floatVector inversePrecedingRCG = tardigradeVectorTools::inverse( precedingRCG, *dim, *dim );

            floatMatrix dRCGdF = tardigradeVectorTools::dot( dRCGdPrecedingF,  *dPrecedingFdF );

            floatMatrix dRCGdFn = tardigradeVectorTools::dot( dRCGdPrecedingF, *dPrecedingFdFn );

            // Form the precedingPsi and its inverse
            floatVector precedingPsi;

            TARDIGRADE_ERROR_TOOLS_CATCH( precedingPsi = tardigradeVectorTools::matrixMultiply( *precedingDeformationGradient, *precedingMicroDeformation, *dim, *dim, *dim, *dim, true, false ) );

            floatMatrix dPsidPrecedingF( precedingPsi.size( ), floatVector( precedingDeformationGradient->size( ), 0 ) );

            floatMatrix dPsidPrecedingChi( precedingPsi.size( ), floatVector( precedingDeformationGradient->size( ), 0 ) );

            floatVector eye( precedingDeformationGradient->size( ), 0 );
            tardigradeVectorTools::eye( eye );

            for ( unsigned int I = 0; I < *dim; I++ ){

                for ( unsigned int J = 0; J < *dim; J++ ){

                    for ( unsigned int A = 0; A < *dim; A++ ){

                        for ( unsigned int B = 0; B < *dim; B++ ){

                            dPsidPrecedingF[   ( *dim ) * I + J ][ ( *dim ) * A + B ] += eye[ ( *dim ) * I + B ] * ( *precedingMicroDeformation )[ ( *dim ) * A + J ];

                            dPsidPrecedingChi[ ( *dim ) * I + J ][ ( *dim ) * A + B ] += ( *precedingDeformationGradient )[ ( *dim ) * A + I ] * eye[ ( *dim ) * J + B ];

                        }

                    }

                }

            }

            floatMatrix dPsidF    = tardigradeVectorTools::dot( dPsidPrecedingF,   *dPrecedingFdF  );

            floatMatrix dPsidFn   = tardigradeVectorTools::dot( dPsidPrecedingF,   *dPrecedingFdFn );

            floatMatrix dPsidChi  = tardigradeVectorTools::dot( dPsidPrecedingChi, *dPrecedingChidChi  );

            floatMatrix dPsidChin = tardigradeVectorTools::dot( dPsidPrecedingChi, *dPrecedingChidChin );

            floatVector inversePrecedingPsi = tardigradeVectorTools::inverse( precedingPsi, *dim, *dim );

            // Form the preceding micro RCG and its inverse
            floatVector precedingMicroRCG;

            floatMatrix dMicroRCGdPrecedingChi;

            TARDIGRADE_ERROR_TOOLS_CATCH_NODE_POINTER( tardigradeConstitutiveTools::computeRightCauchyGreen( *precedingMicroDeformation, precedingMicroRCG, dMicroRCGdPrecedingChi ) );

            floatMatrix dMicroRCGdChi  = tardigradeVectorTools::dot( dMicroRCGdPrecedingChi, *dPrecedingChidChi );

            floatMatrix dMicroRCGdChin = tardigradeVectorTools::dot( dMicroRCGdPrecedingChi, *dPrecedingChidChin );

//            // Form Gamma
//            floatVector precedingGamma( ( *dim ) * ( *dim ) * ( *dim ), 0 );
//            for ( unsigned int I = 0; I < *dim; I++ ){
//
//                for ( unsigned int J = 0; J < *dim; J++ ){
//
//                    for ( unsigned int K = 0; K < *dim; K++ ){
//
//                        for ( unsigned int i = 0; i < *dim; i++ ){
//
//                            precedingGamma[ ( *dim ) * ( *dim ) * I + ( *dim ) * J + K ]
//                                += ( *precedingDeformationGradient )[ ( *dim ) * i + I ] * ( *precedingGradientMicroDeformation )[ ( *dim ) * ( *dim ) * i + ( *dim ) * J + K ];
//
//                        }
//
//                    }
//
//                }
//
//            }

            floatVector macroVelocityGradient;

            floatVector microVelocityGradient;

            floatVector dPlasticMacroLdMacroGamma;

            floatVector dPlasticMacroLdMicroGamma;

            floatMatrix dPlasticMacroLdPrecedingRCG;

            floatMatrix dPlasticMacroLdMacroFlowDirection;

            floatMatrix dPlasticMacroLdMicroFlowDirection;

            floatVector dPlasticMicroLdMicroGamma;

            floatMatrix dPlasticMicroLdPrecedingMicroRCG;

            floatMatrix dPlasticMicroLdPrecedingPsi;

            floatMatrix dPlasticMicroLdMicroFlowDirection;

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMacroVelocityGradient( ( *plasticMultipliers )[ 0 ], ( *plasticMultipliers )[ 1 ],
                                                     inversePrecedingRCG, *dMacroFlowdDrivingStress, *dMicroFlowdDrivingStress,
                                                     macroVelocityGradient, dPlasticMacroLdMacroGamma, dPlasticMacroLdMicroGamma,
                                                     dPlasticMacroLdPrecedingRCG, dPlasticMacroLdMacroFlowDirection, dPlasticMacroLdMicroFlowDirection )
            )

            TARDIGRADE_ERROR_TOOLS_CATCH(
                computePlasticMicroVelocityGradient( ( *plasticMultipliers )[ 1 ], precedingMicroRCG, precedingPsi, inversePrecedingPsi,
                                                     *dMicroFlowdDrivingStress, microVelocityGradient, dPlasticMicroLdMicroGamma,
                                                     dPlasticMicroLdPrecedingMicroRCG, dPlasticMicroLdPrecedingPsi, dPlasticMicroLdMicroFlowDirection );
            )

            // Assemble the Jacobians
            floatMatrix dPlasticMacroLdMacroStress = tardigradeVectorTools::dot( dPlasticMacroLdMacroFlowDirection, *d2MacroFlowdDrivingStressdStress );

            floatMatrix dPlasticMacroLdMicroStress = tardigradeVectorTools::dot( dPlasticMacroLdMicroFlowDirection, *d2MicroFlowdDrivingStressdStress );

            floatMatrix dPlasticMicroLdMicroStress = tardigradeVectorTools::dot( dPlasticMicroLdMicroFlowDirection, *d2MicroFlowdDrivingStressdStress );

            floatMatrix dPlasticMacroLdF           = tardigradeVectorTools::dot( dPlasticMacroLdMacroFlowDirection, *d2MacroFlowdDrivingStressdF )
                                                   + tardigradeVectorTools::dot( dPlasticMacroLdMicroFlowDirection, *d2MicroFlowdDrivingStressdF )
                                                   + tardigradeVectorTools::dot( dPlasticMacroLdPrecedingRCG, dRCGdF );

            floatMatrix dPlasticMacroLdFn          = tardigradeVectorTools::dot( dPlasticMacroLdMacroFlowDirection, *d2MacroFlowdDrivingStressdFn )
                                                   + tardigradeVectorTools::dot( dPlasticMacroLdMicroFlowDirection, *d2MicroFlowdDrivingStressdFn )
                                                   + tardigradeVectorTools::dot( dPlasticMacroLdPrecedingRCG, dRCGdFn );

            floatMatrix dPlasticMicroLdFn          = tardigradeVectorTools::dot( dPlasticMicroLdMicroFlowDirection, *d2MicroFlowdDrivingStressdFn )
                                                   + tardigradeVectorTools::dot( dPlasticMicroLdPrecedingPsi, dPsidFn );

            floatMatrix dPlasticMicroLdChin        = tardigradeVectorTools::dot( dPlasticMicroLdPrecedingPsi, dPsidChin );
                                                     tardigradeVectorTools::dot( dPlasticMicroLdPrecedingMicroRCG, dMicroRCGdChin );

            floatMatrix dPlasticMacroLdISVs( macroVelocityGradient.size( ), floatVector( get_plasticStateVariables( )->size( ), 0 ) );

            for ( unsigned int i = 0; i < ( *dim ) * ( *dim ); i++ ){

                dPlasticMacroLdISVs[ i ][ 0 ] = dPlasticMacroLdMacroGamma[ i ];

                dPlasticMacroLdISVs[ i ][ 1 ] = dPlasticMacroLdMicroGamma[ i ];

            }

            if ( isPrevious ){

                set_previousPlasticMacroVelocityGradient( macroVelocityGradient );

                set_previousPlasticMicroVelocityGradient( microVelocityGradient );

                set_previousdPlasticMacroVelocityGradientdMacroStress( dPlasticMacroLdMacroStress );

                set_previousdPlasticMacroVelocityGradientdMicroStress( dPlasticMacroLdMicroStress );

                set_previousdPlasticMicroVelocityGradientdMicroStress( dPlasticMicroLdMicroStress );

                set_previousdPlasticMacroVelocityGradientdF( dPlasticMacroLdF );

                set_previousdPlasticMacroVelocityGradientdFn( dPlasticMacroLdFn );

                set_previousdPlasticMicroVelocityGradientdFn( dPlasticMicroLdFn );

                set_previousdPlasticMicroVelocityGradientdChin( dPlasticMicroLdChin );

                set_previousdPlasticMacroVelocityGradientdStateVariables( dPlasticMacroLdISVs );

            }
            else{

                set_plasticMacroVelocityGradient( macroVelocityGradient );

                set_plasticMicroVelocityGradient( microVelocityGradient );

                set_dPlasticMacroVelocityGradientdMacroStress( dPlasticMacroLdMacroStress );

                set_dPlasticMacroVelocityGradientdMicroStress( dPlasticMacroLdMicroStress );

                set_dPlasticMicroVelocityGradientdMicroStress( dPlasticMicroLdMicroStress );

                set_dPlasticMacroVelocityGradientdF( dPlasticMacroLdF );

                set_dPlasticMacroVelocityGradientdFn( dPlasticMacroLdFn );

                set_dPlasticMicroVelocityGradientdFn( dPlasticMicroLdFn );

                set_dPlasticMicroVelocityGradientdChin( dPlasticMicroLdChin );

                set_dPlasticMacroVelocityGradientdStateVariables( dPlasticMacroLdISVs );

            }

        }

    }

}
